/* arcbatch.c -- fits and bootstraps that take MANY runs in one call and share one parallel pool.
 *
 *   arcbt_k4eq       kappa curve with drift, shape held on the equivalence locus
 *   arcbt_host       beta-companion curve with drift, free and on its equivalence curve
 *   arcbt_k4_boot    residual bootstrap of a free (or h = 0) kappa fit
 *   arcbt_k4eq_boot  residual bootstrap of the constrained kappa fit
 *   arcbt_locus_k    the equivalence locus k(h) itself
 *
 * The mean curve is g0 + m x + g1 F((x - location)/scale; shape). It is linear in (g0, m, g1), which
 * are solved exactly at every trial shape; the search is Levenberg-Marquardt on the projected
 * residual (the Jacobian is -g1 P dF/dq, P the projector off the span of 1, x and F), from a grid of
 * starts, the best of which are carried to convergence and polished by a simplex.
 *
 * Every parallel loop writes its own slots, the selections are serial, and every bootstrap replicate
 * has its own splitmix64 stream seeded by (seed, replicate), so nothing depends on the number of
 * threads.
 */
#if defined(__clang__)
#pragma clang fp contract(off)
#elif defined(__GNUC__)
#pragma GCC optimize ("fp-contract=off")
#else
#pragma STDC FP_CONTRACT OFF
#endif

#include <math.h>
#include <stdlib.h>
#include <string.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "arck4.h"
#include "arceq2.h"
#include "arcef.h"
#include "arcbatch.h"
#include "arcthreads.h"

#define BT_BAD 1e12
#define BT_DMAX 4

/* ================= the kappa distribution function and its equivalence locus ================= */

/* F and t dF/dt-free derivative dF/dt at t = (x - mu)/sg; log1p throughout, so that the pass
 * through h = 0 and k = 0 is smooth */
static double bt_k4F(double t, double k, double h, double *dFdt) {
  double z = k * t, ly;
  if (dFdt) *dFdt = 0.0;
  if (k == 0.0) ly = -t;
  else {
    if (z >= 1.0) return (k > 0.0) ? 1.0 : 0.0;
    ly = log1p(-z) / k;
  }
  double y = exp(ly);
  if (!(y < HUGE_VAL)) return 0.0;
  double hy = h * y, F;
  if (hy >= 1.0) return 0.0;
  F = (h == 0.0) ? exp(-y) : exp(log1p(-hy) / h);
  if (dFdt && F > 0.0) {
    double d = F * y / ((1.0 - hy) * (1.0 - z));
    if (isfinite(d)) *dFdt = d;
  }
  return F;
}

/* tangent minus third-derivative reading of the standard kappa, as a function of one log
 * coordinate s that keeps full precision where the locus runs into a boundary:
 *   h < -1:  s = log(1 - 2hk), the distance to the frontier hk = 1/2 (the locus approaches it like
 *            exp(c h), far below the spacing of doubles in k for h < -15);
 *   h >= -1: s = log(1 - k), the distance to k = 1, which the locus reaches as h -> 1/2.
 * The readings are those of arck4_readings. */
static double bt_locus_disc(double s, double h, double *kout) {
  double d = 1.0 - h, c, k, e, hk;
  if (h < -1.0) { e = exp(s); hk = 0.5 * (1.0 - e); k = hk / h; c = 1.0 - k; }
  else         { c = exp(s); k = 1.0 - c; hk = h * k; e = 1.0 - 2.0 * hk; }
  *kout = k;
  double omhk = 0.5 * (1.0 + e);                       /* 1 - hk */
  double rc = d / c, lw = log(c / omhk);
  double xc = (k == 0.0) ? -lw : -expm1(k * lw) / k;
  double a = xc - rc * exp(k * lw);
  double disc = c * d * (9.0 * c * d - 4.0 * (1.0 - 2.0 * k) * (1.0 - 2.0 * h));
  if (!(disc >= 0.0)) return NAN;
  double sq = sqrt(disc), rph;
  if (h >= 0.0) rph = 2.0 * d * (1.0 - 2.0 * h) / (3.0 * c * d + sq) + h;
  else rph = 2.0 * e * omhk / (3.0 * c * d + 2.0 * c * (1.0 - 2.0 * k) * h + sq);
  if (!(rph > 0.0)) return NAN;
  double lb = -log(rph), b = (k == 0.0) ? -lb : -expm1(k * lb) / k;
  return a - b;
}
/* the one k with equal readings at h < 1/2; NaN otherwise. Bisection: the discrepancy is positive
 * towards the frontier (at k = -0.49 for h >= -1) and negative towards k = 1. */
static double bt_locus_k(double h) {
  if (!(h < 0.5) || !isfinite(h)) return NAN;
  double lo, hi, k, D, sgn;
  if (h < -1.0) { lo = -700.0; hi = log1p(-2.0 * h); sgn = 1.0; }     /* positive at lo */
  else          { lo = -36.0;  hi = log(1.49);       sgn = -1.0; }    /* negative at lo */
  D = bt_locus_disc(lo, h, &k);
  if (!(D * sgn > 0.0)) return (h < -1.0) ? 1.0 / (2.0 * h) : NAN;
  for (int it = 0; it < 100; it++) {
    double m = 0.5 * (lo + hi);
    if (m == lo || m == hi) break;
    D = bt_locus_disc(m, h, &k);
    if (D * sgn > 0.0) lo = m; else hi = m;
    if (hi - lo < 1e-16 * (1.0 + fabs(lo))) break;
  }
  bt_locus_disc(0.5 * (lo + hi), h, &k);
  return k;
}
void arcbt_locus_k(const double *h, const int *n, double *k) {
  for (int i = 0; i < *n; i++) k[i] = bt_locus_k(h[i]);
}

/* ========================== the shared variable-projection machinery ========================== */
/* q = (location, log scale, shape...). The kinds:
 *   BT_FULL   kappa, q[2] = t, h = 1/2 - exp(t), k on the locus
 *   BT_TAB    kappa, q[2] = h, k by linear interpolation of a table
 *   BT_FREE   companion, q[2] = log(1+alpha), q[3] = log(1+beta)
 *   BT_CON    companion, q[2] = alpha, beta by linear interpolation of a table */
enum { BT_FULL = 0, BT_TAB = 1, BT_FREE = 2, BT_CON = 3 };
typedef struct {
  const double *x, *y;
  int n, d, kind;
  const double *tx, *ty; int nt;       /* shape table */
  double lo[BT_DMAX], hi[BT_DMAX];     /* box of q */
  double klo, khi;
  double *F;                           /* n: curve of the simplex objective */
  double *W;                           /* (4 + 2d + 5) n: work of the LM step */
} bt_ctx;

static void bt_shape(const bt_ctx *c, const double *q, double *s1, double *s2) {
  switch (c->kind) {
  case BT_FULL: *s2 = 0.5 - exp(q[2]); *s1 = bt_locus_k(*s2); break;
  case BT_TAB:  *s2 = q[2]; *s1 = arcef_interp(c->tx, c->ty, c->nt, q[2], 1); break;
  case BT_FREE: *s1 = expm1(q[2]); *s2 = expm1(q[3]); break;
  default:      *s1 = q[2]; *s2 = arcef_interp(c->tx, c->ty, c->nt, q[2], 1);
  }
}
static void bt_bc_curve(const double *x, int n, double mu, double sr, double a, double b, double *F) {
  ib_par P; arceq2_ib_init(&P, a + 1.0, b + 1.0);
  double s = NAN;
  for (int i = 0; i < n; i++) {
    double z = (x[i] - mu) / sr;
    if (z <= 0.0) { F[i] = 0.0; continue; }            /* exactly 0 and 1 off the support */
    if (z >= 1.0) { F[i] = 1.0; continue; }
    if (z < 1e-300) z = 1e-300;
    s = arceq2_ib_at_frac(&P, z, s);
    F[i] = exp(s > 0.0 ? -log1p(exp(-s)) : s - log1p(exp(s)));
    if (!isfinite(F[i])) { F[i] = NAN; s = NAN; }
  }
}
/* the standard curve at q; 0 if q is inadmissible */
static int bt_curve(const bt_ctx *c, const double *q, double *F) {
  for (int j = 0; j < c->d; j++) if (!(q[j] >= c->lo[j] && q[j] <= c->hi[j])) return 0;
  double s1, s2, sc = exp(q[1]);
  if (!isfinite(sc) || !(sc > 0.0)) return 0;
  bt_shape(c, q, &s1, &s2);
  if (!isfinite(s1) || !isfinite(s2)) return 0;
  if (c->kind <= BT_TAB) {
    if (s1 < c->klo || s1 > c->khi) return 0;
    for (int i = 0; i < c->n; i++) F[i] = bt_k4F((c->x[i] - q[0]) / sc, s1, s2, NULL);
  } else bt_bc_curve(c->x, c->n, q[0], sc, s1, s2, F);
  for (int i = 0; i < c->n; i++) if (!isfinite(F[i])) return 0;
  return 1;
}

typedef struct { double mx, mF, my, sxx, gx, dF, g1, m; } bt_lin;
/* the linear part by centring and one Gram-Schmidt step; 0 when F is not a third direction */
static int bt_linear(const bt_ctx *c, const double *F, bt_lin *L) {
  int n = c->n;
  double a = 0, b = 0, d = 0, lo = HUGE_VAL, hi = -HUGE_VAL;
  for (int i = 0; i < n; i++) {
    a += c->x[i]; b += F[i]; d += c->y[i];
    if (F[i] < lo) lo = F[i];
    if (F[i] > hi) hi = F[i];
  }
  if (!(hi - lo >= 1e-8)) return 0;
  a /= n; b /= n; d /= n;
  double xx = 0, xF = 0, FF = 0, xy = 0;
  for (int i = 0; i < n; i++) {
    double u = c->x[i] - a, v = F[i] - b;
    xx += u * u; xF += u * v; FF += v * v; xy += u * (c->y[i] - d);
  }
  if (!(xx > 0.0)) return 0;
  double g = xF / xx, D = 0, e = 0, sy = xy / xx;
  for (int i = 0; i < n; i++) {                        /* the part of F off (1, x), and y on it */
    double u = c->x[i] - a, w = F[i] - b - g * u;
    D += w * w; e += w * (c->y[i] - d - sy * u);
  }
  if (!(D > 1e-14 * FF)) return 0;
  L->mx = a; L->mF = b; L->my = d; L->sxx = xx; L->gx = g; L->dF = D;
  L->g1 = e / D; L->m = sy - L->g1 * g;
  return isfinite(L->g1);
}
static double bt_res(const bt_ctx *c, const double *F, const bt_lin *L, double *r) {
  double v = 0;
  for (int i = 0; i < c->n; i++) {
    double e = c->y[i] - L->my - L->m * (c->x[i] - L->mx) - L->g1 * (F[i] - L->mF);
    if (r) r[i] = e;
    v += e * e;
  }
  return v;
}
static double bt_rss(const bt_ctx *c, const double *q, double *F) {
  bt_lin L;
  if (!bt_curve(c, q, F) || !bt_linear(c, F, &L)) return BT_BAD;
  double v = bt_res(c, F, &L, NULL);
  return isfinite(v) ? v : BT_BAD;
}
static double bt_obj(const double *q, const void *vc) {
  const bt_ctx *c = (const bt_ctx *)vc;
  return bt_rss(c, q, c->F);
}
/* dF/dq in J (d columns, the first dj of them); F in F; P, M scratch. The shape column is a central
 * difference inside the box */
static int bt_jac(const bt_ctx *c, const double *q, int dj, double *F, double *J, double *P, double *D5) {
  int n = c->n;
  double s1, s2, sc = exp(q[1]);
  bt_shape(c, q, &s1, &s2);
  if (c->kind <= BT_TAB) {
    for (int i = 0; i < n; i++) {
      double t = (c->x[i] - q[0]) / sc, f;
      F[i] = bt_k4F(t, s1, s2, &f);
      J[i] = -f / sc; J[n + i] = -t * f;
    }
    if (dj > 2) {
      double *C = J + 2 * (size_t)n, a[3] = {q[0], q[1], q[2]}, b[3] = {q[0], q[1], q[2]};
      double st = (c->kind == BT_FULL) ? 1e-5 : 1e-6 * (c->hi[2] - c->lo[2]);
      a[2] = fmin(a[2] + st, c->hi[2]); b[2] = fmax(b[2] - st, c->lo[2]);
      double den = a[2] - b[2];
      if (!(den > 0.0)) return 0;
      if (!bt_curve(c, a, C) || !bt_curve(c, b, P)) return 0;
      for (int i = 0; i < n; i++) C[i] = (C[i] - P[i]) / den;
    }
    return 1;
  }
  arcef_jac(c->x, n, q[0], sc, s1, s2, 0.0, D5);
  double sl = 0.0;
  if (c->kind == BT_CON) {
    double st = 1e-7 * (c->hi[2] - c->lo[2]), a = fmin(q[2] + st, c->hi[2]), b = fmax(q[2] - st, c->lo[2]);
    if (a > b) sl = (arcef_interp(c->tx, c->ty, c->nt, a, 1) - arcef_interp(c->tx, c->ty, c->nt, b, 1)) / (a - b);
  }
  for (int i = 0; i < n; i++) {
    F[i] = D5[i];
    J[i] = D5[n + i]; J[n + i] = D5[2 * (size_t)n + i] * sc;
    if (dj > 2) J[2 * (size_t)n + i] = (c->kind == BT_CON)
      ? D5[3 * (size_t)n + i] + sl * D5[4 * (size_t)n + i] : D5[3 * (size_t)n + i] * (1.0 + s1);
    if (dj > 3) J[3 * (size_t)n + i] = D5[4 * (size_t)n + i] * (1.0 + s2);
  }
  return 1;
}
/* Levenberg-Marquardt from q (in place) over the first dj coordinates; returns the sum of squares.
 * A coordinate at its bound whose step points outwards is frozen for that iteration. */
static double bt_lm(const bt_ctx *c, double *q, int dj, int maxit) {
  int n = c->n;
  double *F = c->W, *r = c->W + n, *P = c->W + 2 * (size_t)n, *J = c->W + 3 * (size_t)n;
  double *D5 = J + (size_t)c->d * n;
  for (int j = 0; j < c->d; j++) { if (q[j] < c->lo[j]) q[j] = c->lo[j]; if (q[j] > c->hi[j]) q[j] = c->hi[j]; }
  double f = bt_rss(c, q, F), lam = 1e-2;
  int small = 0;
  if (f >= BT_BAD) return f;
  for (int it = 0; it < maxit; it++) {
    bt_lin L;
    if (!bt_jac(c, q, dj, F, J, P, D5) || !bt_linear(c, F, &L)) break;
    bt_res(c, F, &L, r);
    double A[BT_DMAX][BT_DMAX], g[BT_DMAX];
    int frz[BT_DMAX];
    for (int k = 0; k < dj; k++) {
      double *C = J + (size_t)k * n, cm = 0, cx = 0, cF = 0;
      for (int i = 0; i < n; i++) cm += C[i];
      cm /= n;
      for (int i = 0; i < n; i++) {
        double z = C[i] - cm, u = c->x[i] - L.mx;
        cx += z * u; cF += z * (F[i] - L.mF - L.gx * u);
      }
      double p1 = cF / L.dF, p0 = cx / L.sxx - p1 * L.gx;
      g[k] = 0;
      for (int i = 0; i < n; i++) {
        C[i] = L.g1 * (C[i] - cm - p0 * (c->x[i] - L.mx) - p1 * (F[i] - L.mF));
        g[k] += C[i] * r[i];
      }
      for (int l = 0; l <= k; l++) {
        const double *Cl = J + (size_t)l * n;
        double s = 0;
        for (int i = 0; i < n; i++) s += C[i] * Cl[i];
        A[k][l] = s;
      }
      frz[k] = (q[k] <= c->lo[k] && g[k] < 0.0) || (q[k] >= c->hi[k] && g[k] > 0.0);
    }
    for (int k = 0; k < dj; k++) if (frz[k]) {
      g[k] = 0;
      for (int l = 0; l < k; l++) A[k][l] = 0;
      for (int l = k + 1; l < dj; l++) A[l][k] = 0;
    }
    int ok = 0;
    for (int tr = 0; tr < 14; tr++) {
      double L3[BT_DMAX][BT_DMAX], z[BT_DMAX], st[BT_DMAX];
      int bad = 0;
      for (int k = 0; k < dj && !bad; k++) for (int l = 0; l <= k; l++) {
        double t = A[k][l];
        if (k == l) t += lam * (A[k][k] > 0.0 ? A[k][k] : 1e-12);
        for (int h = 0; h < l; h++) t -= L3[k][h] * L3[l][h];
        if (k == l) { if (!(t > 0.0)) { bad = 1; break; } L3[k][k] = sqrt(t); }
        else L3[k][l] = t / L3[l][l];
      }
      if (bad) { lam *= 10.0; continue; }
      for (int k = 0; k < dj; k++) {
        double t = g[k];
        for (int h = 0; h < k; h++) t -= L3[k][h] * z[h];
        z[k] = t / L3[k][k];
      }
      for (int k = dj - 1; k >= 0; k--) {
        double t = z[k];
        for (int h = k + 1; h < dj; h++) t -= L3[h][k] * st[h];
        st[k] = t / L3[k][k];
      }
      double qn[BT_DMAX] = {0};
      for (int k = 0; k < c->d; k++) qn[k] = q[k];
      for (int k = 0; k < dj; k++) {
        qn[k] += st[k];
        if (qn[k] < c->lo[k]) qn[k] = c->lo[k];
        if (qn[k] > c->hi[k]) qn[k] = c->hi[k];
      }
      double fn = bt_rss(c, qn, P);
      if (fn < f) {
        small = (f - fn < 1e-12 * f) ? small + 1 : 0;
        for (int k = 0; k < c->d; k++) q[k] = qn[k];
        f = fn; ok = 1;
        lam = lam * 0.3 > 1e-12 ? lam * 0.3 : 1e-12;
        break;
      }
      lam *= 4.0;
    }
    if (!ok || small >= 3) break;
  }
  return f;
}
static double bt_polish(const bt_ctx *c, double *q, double f, int passes) {
  double p[BT_DMAX], g;
  for (int s = 0; s < passes; s++) {
    arck4_nmd(bt_obj, c, c->d, q, 600, p, &g);
    if (!(g < f)) break;
    int tiny = (f - g <= 1e-13 * f);
    f = g; memcpy(q, p, (size_t)c->d * sizeof(double));
    if (tiny) break;
  }
  return f;
}
/* The profile over the shape coordinate: location and scale refitted at every node by continuation
 * from q, both ways, and the fit resumed from the three lowest local minima of that profile. This
 * reaches basins at very large scales, which a grid of starts misses. With a linearly interpolated
 * table the nodes are those of the table, where the objective has kinks that neither LM nor a
 * simplex leaves when they are local minima: there the two cells beside a node are searched
 * separately. For the computed locus the nodes are BT_NG equally spaced values of the coordinate. */
#define BT_NG 48
static double bt_nodes(const bt_ctx *c, double *q, double f, int it) {
  int j = c->d - 1, tab = (c->kind != BT_FULL), i0 = 0, m = BT_NG;
  double gr[BT_NG];
  const double *tx = c->tx;
  if (tab) {
    int i1 = c->nt - 1;
    while (i0 < c->nt && tx[i0] < c->lo[j]) i0++;
    while (i1 >= 0 && tx[i1] > c->hi[j]) i1--;
    m = i1 - i0 + 1;
  } else {
    for (int i = 0; i < m; i++) gr[i] = c->lo[j] + (c->hi[j] - c->lo[j]) * i / (m - 1.0);
    gr[m - 1] = c->hi[j];
    tx = gr;
  }
  if (m < 2) return f;
  double *pf = (double *)malloc((size_t)m * sizeof(double));
  double *pq = (double *)malloc((size_t)m * 2 * sizeof(double));
  if (!pf || !pq) { free(pf); free(pq); return f; }
  int is = 0;
  while (is < m - 1 && tx[i0 + is + 1] <= q[j]) is++;
  for (int dir = -1; dir <= 1; dir += 2) {
    double p[BT_DMAX] = {q[0], q[1], q[2], q[3]};
    for (int i = (dir < 0) ? is : is + 1; i >= 0 && i < m; i += dir) {
      p[j] = tx[i0 + i];
      pf[i] = bt_lm(c, p, j, 25);
      pq[2 * i] = p[0]; pq[2 * i + 1] = p[1];
      if (pf[i] >= BT_BAD) { p[0] = q[0]; p[1] = q[1]; }
    }
  }
  double last = -HUGE_VAL;
  for (int rank = 0; rank < 3; rank++) {
    int ib = -1;
    for (int i = 0; i < m; i++) {
      if (pf[i] >= BT_BAD || (i > 0 && pf[i-1] < pf[i]) || (i < m - 1 && pf[i+1] < pf[i])) continue;
      if (pf[i] <= last) continue;
      if (ib < 0 || pf[i] < pf[ib]) ib = i;
    }
    if (ib < 0) break;
    last = pf[ib];
    for (int cell = ib - 1; cell <= ib; cell++) {
      if (tab && (cell < 0 || cell >= m - 1)) continue;
      if (!tab && cell < ib) continue;
      bt_ctx cc = *c;
      double p[BT_DMAX] = {pq[2 * ib], pq[2 * ib + 1], tx[i0 + ib], 0};
      if (tab) {
        cc.lo[j] = tx[i0 + cell]; cc.hi[j] = tx[i0 + cell + 1];
        p[j] = 0.5 * (cc.lo[j] + cc.hi[j]);
      }
      double g = bt_lm(&cc, p, c->d, it);
      if (!tab || g < f * 1.001) g = bt_polish(&cc, p, g, tab ? 1 : 2);
      if (g < f) { f = g; memcpy(q, p, (size_t)c->d * sizeof(double)); }
    }
    if (pf[ib] < f) { f = pf[ib]; q[0] = pq[2 * ib]; q[1] = pq[2 * ib + 1]; q[j] = tx[i0 + ib]; }
  }
  free(pf); free(pq);
  return f;
}
static double bt_converge(const bt_ctx *c, double *q, int it, int passes) {
  double f = bt_lm(c, q, c->d, it);
  return (f >= BT_BAD) ? f : bt_polish(c, q, f, passes);
}

/* empirical quantile, R type 7, of sorted v */
static double bt_q7(const double *v, int n, double p) {
  double t = (n - 1) * p; int i = (int)floor(t);
  if (i >= n - 1) return v[n - 1];
  return v[i] + (t - i) * (v[i + 1] - v[i]);
}
static int bt_cmp(const void *a, const void *b) {
  double x = *(const double *)a, y = *(const double *)b;
  return (x > y) - (x < y);
}

/* ---- a batch: T tasks (run x kind), each with its own starts ---- */
typedef struct { int run, kind, piece; int s0, ns; } bt_task;

static size_t bt_worklen(int n) { return (size_t)n * (1 + 3 + BT_DMAX + 5); }
static double *bt_grow(double *W, int *cap, int n) {
  if (n <= *cap) return W;
  free(W);
  W = (double *)malloc(bt_worklen(n) * sizeof(double));
  *cap = W ? n : 0;
  return W;
}

/* the K best distinct first-stage points of a task, best first; returns how many */
static int bt_select(const bt_task *t, int d, const double *v, const double *Q, int K, int *idx) {
  int m = 0, used = 0;
  char *dead = (char *)calloc((size_t)t->ns, 1);
  if (!dead) return 0;
  while (m < K && used < t->ns) {
    int b = -1;
    for (int s = 0; s < t->ns; s++)
      if (!dead[s] && (b < 0 || v[t->s0 + s] < v[t->s0 + b])) b = s;
    if (b < 0) break;
    dead[b] = 1; used++;
    if (v[t->s0 + b] >= BT_BAD) break;
    const double *qb = Q + BT_DMAX * (size_t)(t->s0 + b);
    int dup = 0;
    for (int e = 0; e < m && !dup; e++) {
      const double *qe = Q + BT_DMAX * (size_t)idx[e];
      dup = fabs(v[t->s0 + b] - v[idx[e]]) <= 1e-4 * v[idx[e]] &&
            fabs(qb[0] - qe[0]) <= 1e-2 * exp(qe[1]);
      for (int k = 1; k < d && dup; k++) dup = fabs(qb[k] - qe[k]) <= 1e-2;
    }
    if (!dup) idx[m++] = t->s0 + b;
  }
  free(dead);
  return m;
}

/* multi-start over all tasks. mk fills the context of a task (its F and W are set here). bv, bq get
 * the best value and point of every task; from (or -1) names the earlier task whose best point,
 * mapped by remap, is a further start of this one. */
typedef void (*bt_mk)(const void *data, const bt_task *t, bt_ctx *c);
static int bt_run(const bt_task *task, int T, const double *st, int M, bt_mk mk, const void *data,
                  const int *off, int it1, int it2, int K, double *bv, double *bq) {
  double *v = (double *)malloc((size_t)M * sizeof(double));
  double *Q = (double *)malloc((size_t)M * BT_DMAX * sizeof(double));
  double *cv = (double *)malloc((size_t)T * K * sizeof(double));
  double *cq = (double *)malloc((size_t)T * K * BT_DMAX * sizeof(double));
  int *idx = (int *)malloc((size_t)T * K * sizeof(int)), *ow = (int *)malloc((size_t)M * sizeof(int));
  int ok = v && Q && cv && cq && idx && ow;
  if (ok) {
    for (int t = 0; t < T; t++) for (int s = 0; s < task[t].ns; s++) ow[task[t].s0 + s] = t;
    for (int i = 0; i < M; i++) v[i] = BT_BAD;
    for (int i = 0; i < T * K; i++) { cv[i] = BT_BAD; idx[i] = -1; }
#ifdef _OPENMP
#pragma omp parallel num_threads(arc_threads())
#endif
    {
      int cap = 0; double *W = NULL;
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 1)
#endif
      for (int i = 0; i < M; i++) {
        const bt_task *t = task + ow[i];
        int n = off[t->run + 1] - off[t->run];
        W = bt_grow(W, &cap, n);
        double *q = Q + BT_DMAX * (size_t)i;
        memcpy(q, st + BT_DMAX * (size_t)i, BT_DMAX * sizeof(double));
        if (!W) continue;
        bt_ctx c; mk(data, t, &c); c.F = W; c.W = W + n;
        v[i] = bt_lm(&c, q, c.d, it1);
      }
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 1)
#endif
      for (int t = 0; t < T; t++) {
        bt_ctx c; mk(data, task + t, &c);
        int m = bt_select(task + t, c.d, v, Q, K, idx + (size_t)t * K);
        for (int k = m; k < K; k++) idx[(size_t)t * K + k] = -1;
      }
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 1)
#endif
      for (int i = 0; i < T * K; i++) {
        if (idx[i] < 0) continue;
        const bt_task *t = task + i / K;
        int n = off[t->run + 1] - off[t->run];
        W = bt_grow(W, &cap, n);
        if (!W) continue;
        bt_ctx c; mk(data, t, &c); c.F = W; c.W = W + n;
        double *q = cq + BT_DMAX * (size_t)i;
        memcpy(q, Q + BT_DMAX * (size_t)idx[i], BT_DMAX * sizeof(double));
        cv[i] = bt_converge(&c, q, it2, 3);
      }
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 1)
#endif
      for (int t = 0; t < T; t++) {
        int b = -1;
        for (int k = 0; k < K; k++) { int i = t * K + k; if (cv[i] < BT_BAD && (b < 0 || cv[i] < cv[b])) b = i; }
        bv[t] = HUGE_VAL;
        for (int k = 0; k < BT_DMAX; k++) bq[BT_DMAX * (size_t)t + k] = NAN;
        if (b < 0) continue;
        double q[BT_DMAX];
        memcpy(q, cq + BT_DMAX * (size_t)b, sizeof q);
        double f = cv[b];
        if (task[t].kind != BT_FREE) {
          int n = off[task[t].run + 1] - off[task[t].run];
          W = bt_grow(W, &cap, n);
          if (W) {
            bt_ctx c; mk(data, task + t, &c); c.F = W; c.W = W + n;
            f = bt_nodes(&c, q, f, it2);
          }
        }
        bv[t] = f;
        memcpy(bq + BT_DMAX * (size_t)t, q, sizeof q);
      }
      free(W);
    }
  }
  free(v); free(Q); free(cv); free(cq); free(idx); free(ow);
  return ok;
}

/* (g0, m, g1, location, scale, shape1, shape2, rss) of a point */
static void bt_theta(const bt_ctx *c, const double *q, double *o) {
  bt_lin L;
  for (int i = 0; i < 8; i++) o[i] = NAN;
  if (!bt_curve(c, q, c->F) || !bt_linear(c, c->F, &L)) return;
  o[7] = bt_res(c, c->F, &L, NULL);
  o[0] = L.my - L.m * L.mx - L.g1 * L.mF; o[1] = L.m; o[2] = L.g1;
  o[3] = q[0]; o[4] = exp(q[1]);
  bt_shape(c, q, &o[5], &o[6]);
}

/* ==================================== kappa on the locus ====================================== */
typedef struct {
  const double *X, *Y; const int *off; int np; const double *TH, *TK; const int *toff;
  const double *bnd; double tlo, thi, tmid;   /* t of h_hi, h_lo, and of h = 0 held between them */
} bt_kdat;

static void bt_kdat_init(bt_kdat *D) {
  if (D->np > 0) return;
  /* the part of the locus inside the box: k increases in h up to the corner, so the upper end is
   * h_hi or the h at which k = k_hi, whichever comes first */
  double hhi = fmin(D->bnd[3], 0.5 - 1e-9), hlo = D->bnd[2];
  if (bt_locus_k(hhi) > D->bnd[1]) {
    double a = fmin(hlo, 0.0), b = hhi;
    for (int i = 0; i < 80; i++) { double m = 0.5 * (a + b); if (bt_locus_k(m) > D->bnd[1]) b = m; else a = m; }
    hhi = a;
  }
  D->tlo = log(0.5 - hhi); D->thi = log(0.5 - hlo);
  D->tmid = fmin(fmax(log(0.5), D->tlo), D->thi);
}
static void bt_mk_k(const void *data, const bt_task *t, bt_ctx *c) {
  const bt_kdat *D = (const bt_kdat *)data;
  memset(c, 0, sizeof *c);
  c->x = D->X + D->off[t->run]; c->y = D->Y + D->off[t->run];
  c->n = D->off[t->run + 1] - D->off[t->run];
  c->d = 3; c->kind = t->kind;
  c->lo[0] = c->lo[1] = -HUGE_VAL; c->hi[0] = c->hi[1] = HUGE_VAL;
  c->klo = D->bnd[0]; c->khi = D->bnd[1];
  if (t->kind == BT_FULL) {              /* piece 0: h <= 0, 1: h >= 0, 2: all of it */
    c->lo[2] = (t->piece == 0) ? D->tmid : D->tlo; c->hi[2] = (t->piece == 1) ? D->tmid : D->thi;
  }
  else {
    c->tx = D->TH + D->toff[t->piece]; c->ty = D->TK + D->toff[t->piece];
    c->nt = D->toff[t->piece + 1] - D->toff[t->piece];
    c->lo[2] = fmax(c->tx[0], D->bnd[2]); c->hi[2] = fmin(c->tx[c->nt - 1], D->bnd[3]);
  }
}
static void bt_k_readings(const double *o, double *ab) {
  int g = 0; double r[3], th[6] = {0.0, 1.0, o[3], o[4], o[5], o[6]};   /* of the shape, not of g1 */
  ab[0] = ab[1] = NAN;
  if (!isfinite(o[7])) return;
  arck4_readings(th, &g, r);
  ab[0] = r[0]; ab[1] = r[1];
}
/* npiece = 0: the whole locus, computed here, inside bounds = (k_lo, k_hi, h_lo, h_hi); otherwise
 * npiece tables (TH increasing within each, split by toff), each interpolated linearly and never
 * extrapolated. out: 11 per run (g0, m, g1, mu, sg, k, h, rss, a, b, piece); outp: the
 * same 11 for the best fit on every piece, pieces within runs; the whole locus is searched as two, h <= 0 and h >= 0. */
void arcbt_k4eq(const double *X, const double *Y, const int *off, const int *nr,
                const int *npiece, const double *TH, const double *TK, const int *toff,
                const double *bounds, const int *it1, const int *it2, const int *K,
                double *out, double *outp) {
  int R = *nr, np = *npiece > 0 ? *npiece : 2, T = R * np;
  static const double pm[5] = {0.15, 0.3, 0.5, 0.7, 0.9}, ps[3] = {0.04, 0.12, 0.35};
  int nh = 8, ns = 5 * 3 * nh;
  bt_kdat D = {X, Y, off, *npiece, TH, TK, toff, bounds, 0, 0, 0};
  bt_kdat_init(&D);
  bt_task *task = (bt_task *)malloc((size_t)T * sizeof(bt_task));
  double *st = (double *)malloc((size_t)T * ns * BT_DMAX * sizeof(double));
  double *bq = (double *)malloc((size_t)T * BT_DMAX * sizeof(double));
  for (int i = 0; i < 11 * R; i++) out[i] = NAN;
  for (int i = 0; i < 11 * T; i++) outp[i] = NAN;
  double *pv = (double *)malloc((size_t)T * sizeof(double));
  if (task && st && bq && pv) {
    for (int t = 0; t < T; t++) {
      bt_task tk = {t / np, *npiece > 0 ? BT_TAB : BT_FULL, t % np, t * ns, ns};
      task[t] = tk;
      bt_ctx c; bt_mk_k(&D, task + t, &c);
      double *v = (double *)malloc((size_t)c.n * sizeof(double));
      double qm[5] = {0}, rg = 0;
      if (v) {
        memcpy(v, c.x, (size_t)c.n * sizeof(double));
        qsort(v, (size_t)c.n, sizeof(double), bt_cmp);
        for (int i = 0; i < 5; i++) qm[i] = bt_q7(v, c.n, pm[i]);
        rg = v[c.n - 1] - v[0];
        free(v);
      }
      double *s = st + (size_t)t * ns * BT_DMAX;
      for (int i = 0; i < 5; i++) for (int j = 0; j < 3; j++) for (int l = 0; l < nh; l++) {
        double w = c.hi[2] - c.lo[2];
        s[0] = qm[i]; s[1] = log(rg * ps[j]);
        s[2] = c.lo[2] + w * (0.03 + 0.94 * l / (nh - 1.0)); s[3] = 0.0;
        s += BT_DMAX;
      }
    }
    if (bt_run(task, T, st, T * ns, bt_mk_k, &D, off, *it1, *it2, *K, pv, bq)) {
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1) num_threads(arc_threads())
#endif
      for (int r = 0; r < R; r++) {
        int b = -1, cap = 0;
        for (int p = 0; p < np; p++) if (isfinite(pv[r * np + p]) && (b < 0 || pv[r * np + p] < pv[r * np + b])) b = p;
        if (b < 0) continue;
        bt_task tk = task[r * np + b];
        if (tk.kind == BT_FULL) tk.piece = 2;             /* the halves are only for the search */
        bt_ctx c; bt_mk_k(&D, &tk, &c);
        double *W = bt_grow(NULL, &cap, c.n);
        if (!W) continue;
        c.F = W; c.W = W + c.n;
        double *q = bq + BT_DMAX * (size_t)(r * np + b), *o = out + 11 * (size_t)r;
        if (tk.kind == BT_FULL) {
          double p[BT_DMAX] = {q[0], q[1], q[2], q[3]};
          double f = bt_converge(&c, p, *it2, 1);
          if (f < pv[r * np + b]) { pv[r * np + b] = f; memcpy(q, p, sizeof p); }
        }
        bt_theta(&c, q, o);
        bt_k_readings(o, o + 8);
        o[10] = b + 1;
        for (int p = 0; p < np; p++) {
          if (!isfinite(pv[r * np + p])) continue;
          double *op = outp + 11 * (size_t)(r * np + p);
          bt_ctx cp; bt_mk_k(&D, task + r * np + p, &cp);
          if (tk.kind == BT_FULL) { bt_mk_k(&D, &tk, &cp); }
          cp.F = W; cp.W = W + cp.n;
          bt_theta(&cp, bq + BT_DMAX * (size_t)(r * np + p), op);
          bt_k_readings(op, op + 8);
          op[10] = p + 1;
        }
        free(W);
      }
    }
  }
  free(pv);
  free(task); free(st); free(bq);
}

/* ==================================== the companion host ====================================== */
typedef struct { const double *X, *Y; const int *off; const double *cal, *cbe; int nc; } bt_hdat;
static void bt_mk_h(const void *data, const bt_task *t, bt_ctx *c) {
  const bt_hdat *D = (const bt_hdat *)data;
  memset(c, 0, sizeof *c);
  c->x = D->X + D->off[t->run]; c->y = D->Y + D->off[t->run];
  c->n = D->off[t->run + 1] - D->off[t->run];
  c->kind = t->kind;
  c->lo[0] = c->lo[1] = -HUGE_VAL; c->hi[0] = c->hi[1] = HUGE_VAL;
  if (t->kind == BT_CON) {
    c->d = 3; c->tx = D->cal; c->ty = D->cbe; c->nt = D->nc;
    c->lo[2] = D->cal[0]; c->hi[2] = D->cal[D->nc - 1];
  } else {                                             /* alpha, beta in [-1 + 1e-6, 20] */
    c->d = 4;
    c->lo[2] = c->lo[3] = log(1e-6); c->hi[2] = c->hi[3] = log(21.0);
  }
}
/* outf, outc: 10 per run (g0, m, g1, s, w, alpha, beta, rss, a, b), the free and the constrained fit.
 * The constrained optimum is a point of the free family, so the free search also starts there. */
void arcbt_host(const double *X, const double *Y, const int *off, const int *nr,
                const double *cal, const double *cbe, const int *nc,
                const int *it1, const int *it2, const int *K, double *outf, double *outc) {
  int R = *nr, T = 2 * R;
  static const double ps[3] = {0.02, 0.12, 0.3}, pw[3] = {0.5, 0.8, 1.2};
  static const double fa[4] = {-0.85, -0.6, -0.35, -0.1}, fb[3] = {-0.85, -0.5, -0.15};
  static const double ca[4] = {0.05, 0.3, 0.6, 0.9};
  int nf = 108, nk = 36, M = R * (nf + nk);
  bt_hdat D = {X, Y, off, cal, cbe, *nc};
  bt_task *task = (bt_task *)malloc((size_t)T * sizeof(bt_task));
  double *st = (double *)calloc((size_t)M * BT_DMAX, sizeof(double));
  double *bq = (double *)malloc((size_t)T * BT_DMAX * sizeof(double));
  double *bv = (double *)malloc((size_t)T * sizeof(double));
  for (int i = 0; i < 10 * R; i++) outf[i] = outc[i] = NAN;
  if (task && st && bq && bv) {
    int s = 0;
    for (int r = 0; r < R; r++) {
      bt_task tf = {r, BT_FREE, 0, s, nf}, tc = {r, BT_CON, 0, s + nf, nk};
      task[2 * r] = tf; task[2 * r + 1] = tc;
      int n = off[r + 1] - off[r];
      double qs[3] = {0}, rg = 1;
      double *v = (double *)malloc((size_t)n * sizeof(double));
      if (v) {
        memcpy(v, X + off[r], (size_t)n * sizeof(double));
        qsort(v, (size_t)n, sizeof(double), bt_cmp);
        for (int i = 0; i < 3; i++) qs[i] = bt_q7(v, n, ps[i]);
        rg = v[n - 1] - v[0];
        free(v);
      }
      double *p = st + BT_DMAX * (size_t)s;
      for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) for (int a = 0; a < 4; a++)
        for (int b = 0; b < 3; b++) {
          p[0] = qs[i]; p[1] = log(rg * pw[j]); p[2] = log1p(fa[a]); p[3] = log1p(fb[b]); p += BT_DMAX;
        }
      for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) for (int a = 0; a < 4; a++) {
        p[0] = qs[i]; p[1] = log(rg * pw[j]); p[2] = cal[0] + ca[a] * (cal[*nc - 1] - cal[0]); p += BT_DMAX;
      }
      s += nf + nk;
    }
    if (bt_run(task, T, st, M, bt_mk_h, &D, off, *it1, *it2, *K, bv, bq)) {
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1) num_threads(arc_threads())
#endif
      for (int r = 0; r < R; r++) {
        int n = off[r + 1] - off[r], cap = 0;
        double *W = bt_grow(NULL, &cap, n);
        if (!W) continue;
        bt_ctx cf, cc;
        bt_mk_h(&D, task + 2 * r, &cf); bt_mk_h(&D, task + 2 * r + 1, &cc);
        cf.F = cc.F = W; cf.W = cc.W = W + n;
        double *qf = bq + BT_DMAX * (size_t)(2 * r), *qc = qf + BT_DMAX;
        if (isfinite(bv[2 * r + 1])) {
          double a, b;
          bt_shape(&cc, qc, &a, &b);
          double q[BT_DMAX] = {qc[0], qc[1], log1p(a), log1p(b)};
          double f = bt_converge(&cf, q, *it2, 3);
          if (f < bv[2 * r] || !isfinite(bv[2 * r])) { bv[2 * r] = f; memcpy(qf, q, sizeof q); }
          double *o = outc + 10 * (size_t)r;
          bt_theta(&cc, qc, o);
          double th[6] = {o[0], o[2], o[3], o[4], o[5], o[6]};
          if (isfinite(o[7])) arcef_bc_ip(th, o + 8);
        }
        if (isfinite(bv[2 * r]) && bv[2 * r] < BT_BAD) {
          double *o = outf + 10 * (size_t)r;
          bt_theta(&cf, qf, o);
          double th[6] = {o[0], o[2], o[3], o[4], o[5], o[6]};
          if (isfinite(o[7])) arcef_bc_ip(th, o + 8);
        }
        free(W);
      }
    }
  }
  free(task); free(st); free(bq); free(bv);
}

/* ======================================== bootstraps ========================================== */
static int bt_bootlen(const double *r, int n, int block) {
  int b = block > 0 ? block : arcef_blocklen(r, n);
  return b > n ? n : b;
}
/* Moving-block residual bootstrap of a fitted kappa curve with drift, many runs in one call.
 * fit: 7 per run (g0, m, g1, mu, sg, k, h); box: 4 per run, the box of that run's own point fit;
 * gev[r]: the point fit is the h = 0 submodel, so h is held at 0 in the replicates; block: 0 reads
 * the block length off the residual autocorrelation, 1 resamples singly. Every replicate is refitted
 * by arck4_fit_varpro from the fitted shape. out: 5 x B per run (a, b, k, h, sel): the readings, the
 * fitted shape, and whether n log(rss_0/rss) - 2 < 0 would choose h = 0 (NaN where h is held). bl:
 * the block length used. seed: one per run; replicate b of a run depends on (seed, b) only. */
void arcbt_k4_boot(const double *X, const double *Y, const int *off, const int *nr,
                   const double *fit, const double *box, const int *gev, const int *B,
                   const int *seed, const int *block, const int *maxit, int *bl, double *out) {
  int R = *nr, NB = *B, N = off[R];
  for (long i = 0; i < 5L * NB * R; i++) out[i] = NAN;
  double *m0 = (double *)malloc((size_t)N * sizeof(double));
  double *r0 = (double *)malloc((size_t)N * sizeof(double));
  if (m0 && r0) {
    for (int r = 0; r < R; r++) {
      const double *f = fit + 7 * (size_t)r;
      for (int i = off[r]; i < off[r + 1]; i++) {
        m0[i] = f[0] + f[1] * X[i] + f[2] * k4_F(X[i], f[3], f[4], f[5], f[6]);
        r0[i] = Y[i] - m0[i];
      }
      bl[r] = bt_bootlen(r0 + off[r], off[r + 1] - off[r], *block);
    }
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1) num_threads(arc_threads())
#endif
    for (int t = 0; t < R * NB; t++) {
      int r = t / NB, b = t % NB, n = off[r + 1] - off[r], one = 1, g = 0;
      const double *f = fit + 7 * (size_t)r, *x = X + off[r];
      int ok = isfinite(f[4]) && f[4] > 0.0 && isfinite(f[5]) && isfinite(f[6]);
      double *yb = ok ? (double *)malloc((size_t)n * sizeof(double)) : NULL;
      if (!yb) continue;
      unsigned long long s = arcef_seed(seed[r], b);
      arcef_blockdraw(m0 + off[r], r0 + off[r], n, bl[r], &s, yb);
      double p0[4] = {f[3], log(f[4]), f[5], gev[r] ? 0.0 : f[6]};
      double hold = gev[r] ? 0.0 : NAN, zero = 0.0, o[8], og[8], rd[3];
      double *res = out + 5 * ((size_t)r * NB + b);
      arck4_fit_varpro(x, yb, &n, p0, &one, maxit, box + 4 * (size_t)r, &hold, o);
      if (isfinite(o[7]) && o[7] < 1e11) {
        double th[6] = {o[0], o[2], o[3], o[4], o[5], o[6]};
        arck4_readings(th, &g, rd);
        res[0] = rd[0]; res[1] = rd[1]; res[2] = o[5]; res[3] = o[6];
        if (!gev[r]) {
          p0[3] = 0.0;
          arck4_fit_varpro(x, yb, &n, p0, &one, maxit, box + 4 * (size_t)r, &zero, og);
          if (isfinite(og[7])) res[4] = (n * log(og[7] / o[7]) - 2.0 < 0.0) ? 1.0 : 0.0;
        }
      }
      free(yb);
    }
  }
  free(m0); free(r0);
}
/* The same for the constrained fit: fit as returned by arcbt_k4eq (its first seven numbers), piece
 * the table piece it lies on (1-based; 1 for the full locus). Replicates are refitted on that piece,
 * from the fitted point. out: 4 x B per run (a, b, k, h). */
void arcbt_k4eq_boot(const double *X, const double *Y, const int *off, const int *nr,
                     const double *fit, const int *piece, const int *npiece, const double *TH,
                     const double *TK, const int *toff, const double *bounds, const int *B,
                     const int *seed, const int *block, const int *it2, int *bl, double *out) {
  int R = *nr, NB = *B, N = off[R];
  for (long i = 0; i < 4L * NB * R; i++) out[i] = NAN;
  double *m0 = (double *)malloc((size_t)N * sizeof(double));
  double *r0 = (double *)malloc((size_t)N * sizeof(double));
  if (m0 && r0) {
    bt_kdat D = {X, Y, off, *npiece, TH, TK, toff, bounds, 0, 0, 0};
    bt_kdat_init(&D);
    for (int r = 0; r < R; r++) {
      const double *f = fit + 7 * (size_t)r;
      for (int i = off[r]; i < off[r + 1]; i++) {
        m0[i] = f[0] + f[1] * X[i] + f[2] * bt_k4F((X[i] - f[3]) / f[4], f[5], f[6], NULL);
        r0[i] = Y[i] - m0[i];
      }
      bl[r] = bt_bootlen(r0 + off[r], off[r + 1] - off[r], *block);
    }
#ifdef _OPENMP
#pragma omp parallel num_threads(arc_threads())
#endif
    {
      int cap = 0; double *W = NULL;
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 1)
#endif
      for (int t = 0; t < R * NB; t++) {
        int r = t / NB, b = t % NB, n = off[r + 1] - off[r];
        const double *f = fit + 7 * (size_t)r;
        double *res = out + 4 * ((size_t)r * NB + b);
        if (!(isfinite(f[4]) && f[4] > 0.0 && isfinite(f[6]))) continue;
        W = bt_grow(W, &cap, n);
        double *yb = (double *)malloc((size_t)n * sizeof(double));
        if (W && yb) {
          unsigned long long s = arcef_seed(seed[r], b);
          arcef_blockdraw(m0 + off[r], r0 + off[r], n, bl[r], &s, yb);
          bt_task tk = {r, *npiece > 0 ? BT_TAB : BT_FULL, *npiece > 0 ? piece[r] - 1 : 2, 0, 0};
          bt_ctx c; bt_mk_k(&D, &tk, &c);
          c.y = yb; c.F = W; c.W = W + n;
          double q[BT_DMAX] = {f[3], log(f[4]), *npiece > 0 ? f[6] : log(0.5 - f[6]), 0};
          double v = bt_converge(&c, q, *it2, 1);
          if (v < BT_BAD && *npiece > 0) v = bt_nodes(&c, q, v, *it2);
          if (v < 1e11) {
            double o[8];
            bt_theta(&c, q, o);
            bt_k_readings(o, res);
            res[2] = o[5]; res[3] = o[6];
          }
        }
        free(yb);
      }
      free(W);
    }
  }
  free(m0); free(r0);
}
