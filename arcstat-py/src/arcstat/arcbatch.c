/* arcbatch.c -- fits and bootstraps that take MANY runs in one call and share one parallel pool.
 *
 *   arcbt_k4eq       kappa curve with drift, shape held on the equivalence locus
 *   arcbt_host       beta-companion curve with drift, free and on its equivalence curve
 *   arcbt_k4_boot    residual bootstrap of a free (or h = 0) kappa fit
 *   arcbt_k4eq_boot  residual bootstrap of the constrained kappa fit
 *   arcbt_wald       inverse Gaussian curve with drift, with its readings
 *   arcbt_lmomw      L-moment matching over the observed window
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
double arcef_k4F(double t, double k, double h, double *dFdt) {
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
 *   BT_K4     kappa, q[2] = k, q[3] = h;  BT_K4H the same with h held (h0)
 *   BT_K4E    kappa with k < 0 in the coordinates of its lower end, q = (xi, log tau, k, h) with
 *             xi = mu + sg/k, tau = -sg/k; for k below -1 location and scale are nearly collinear
 *             and these are not;  BT_K4EH the same with h held
 *   BT_FREE   companion, q[2] = log(1+alpha), q[3] = log(1+beta)
 *   BT_CON    companion, q[2] = alpha, beta by linear interpolation of a table, or, without a table,
 *             the exact equivalence curve
 *   BT_WALD   inverse Gaussian, q = (log mu, log lambda) */
enum { BT_FULL = 0, BT_TAB = 1, BT_K4 = 2, BT_K4H = 3, BT_K4E = 4, BT_K4EH = 5,
       BT_FREE = 6, BT_CON = 7, BT_WALD = 8 };
typedef struct {
  const double *x, *y;
  int n, d, kind;
  const double *tx, *ty; int nt;       /* shape table */
  double lo[BT_DMAX], hi[BT_DMAX];     /* box of q */
  double klo, khi, h0;
  int nodrift;                         /* hold m at 0 */
  double e0, e1;                       /* kinds FULL and TAB: the ends of the searched locus */
  int phys;                            /* admit only fits that are physical, see bt_physical */
  double *F;                           /* n: curve of the simplex objective */
  double *W;                           /* (4 + 2d + 5) n: work of the LM step */
} bt_ctx;

static void bt_shape(const bt_ctx *c, const double *q, double *s1, double *s2) {
  switch (c->kind) {
  case BT_FULL: *s2 = 0.5 - exp(q[2]); *s1 = bt_locus_k(*s2); break;
  case BT_TAB:  *s2 = q[2]; *s1 = arcef_interp(c->tx, c->ty, c->nt, q[2], 1); break;
  case BT_K4E:
  case BT_K4:   *s1 = q[2]; *s2 = q[3]; break;
  case BT_K4EH:
  case BT_K4H:  *s1 = q[2]; *s2 = c->h0; break;
  case BT_FREE: *s1 = expm1(q[2]); *s2 = expm1(q[3]); break;
  case BT_WALD: *s1 = *s2 = 0.0; break;
  default:      *s1 = q[2];
                if (c->nt > 0) *s2 = arcef_interp(c->tx, c->ty, c->nt, q[2], 1);
                else arceq2_bstar(&q[2], s2);
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
/* the kappa distribution function at u = (x - xi)/tau, k < 0: y = u^(1/k) */
double arcef_k4Fe(double u, double k, double h, double *dFdu) {
  if (dFdu) *dFdu = 0.0;
  if (!(u > 0.0)) return 0.0;
  double y = exp(log(u) / k);
  if (!(y < HUGE_VAL)) return 0.0;
  double hy = h * y, F;
  if (hy >= 1.0) return 0.0;
  F = (h == 0.0) ? exp(-y) : exp(log1p(-hy) / h);
  if (dFdu && F > 0.0) {
    double d = -F * y / ((1.0 - hy) * k * u);
    if (isfinite(d)) *dFdu = d;
  }
  return F;
}
/* exp(t^2) erfc(t) for t >= 0 (continued fraction beyond 4) */
static double bt_erfcx(double t) {
  if (t < 4.0) return exp(t * t) * erfc(t);
  double f = 0.0;
  for (int k = 25; k >= 1; k--) f = 0.5 * k / (t + f);
  return 0.56418958354775628695 / (t + f);
}
/* Wald (inverse Gaussian) distribution function at x = mu t, phi = lambda/mu: with z = sqrt(phi/t),
 * Phi(z(t-1)) + exp(2 phi) Phi(-z(t+1)), the product written as exp(-z^2 (t-1)^2/2) erfcx/2, which
 * neither overflows nor cancels */
static double bt_waldF(double t, double phi) {
  if (!(t > 0.0)) return 0.0;
  double z = sqrt(phi / t), a = z * (t - 1.0), b = z * (t + 1.0);
  return 0.5 * erfc(-a * M_SQRT1_2) + 0.5 * exp(-0.5 * a * a) * bt_erfcx(b * M_SQRT1_2);
}
/* (tangent, second, mode) in units of mu. The mode is the root of phi t^2 + 3 t - phi = 0; the second
 * reading the inflection of the density below it, where P^2 + 1.5 t^2 - phi t = 0 with
 * P = (phi (1 - t^2) - 3 t)/2; both exist for every phi > 0, so there is no grid to miss them */
static void bt_wald_read(double phi, double *o) {
  o[0] = o[1] = o[2] = NAN;
  if (!(phi > 0.0) || !isfinite(phi)) return;
  double t = 2.0 * phi / (3.0 + sqrt(9.0 + 4.0 * phi * phi));
  double lf = 0.5 * log(phi / (6.283185307179586476925286766559 * t * t * t)) - phi * (t - 1.0) * (t - 1.0) / (2.0 * t);
  o[2] = t;
  o[0] = t - bt_waldF(t, phi) / exp(lf);
  double lo = 0.0, hi = 1.0;                              /* s = t/mode in (0, 1) */
  for (int i = 0; i < 100; i++) {
    double s = 0.5 * (lo + hi), v = s * t;
    /* P = phi (1 - v^2)/2 - 1.5 v, zero at the mode, so it is written from there */
    double P = 0.5 * (1.0 - s) * t * (phi * t * (1.0 + s) + 3.0);
    if (P * P + v * (1.5 * v - phi) > 0.0) lo = s; else hi = s;
    if (hi - lo < 1e-16) break;
  }
  o[1] = 0.5 * (lo + hi) * t;
}
/* the standard curve at q; 0 if q is inadmissible */
static int bt_curve(const bt_ctx *c, const double *q, double *F) {
  for (int j = 0; j < c->d; j++) if (!(q[j] >= c->lo[j] && q[j] <= c->hi[j])) return 0;
  if (c->kind == BT_WALD) {
    double mu = exp(q[0]), phi = exp(q[1] - q[0]);
    for (int i = 0; i < c->n; i++) if (!isfinite(F[i] = bt_waldF(c->x[i] / mu, phi))) return 0;
    return 1;
  }
  double s1, s2, sc = exp(q[1]);
  if (!isfinite(sc) || !(sc > 0.0)) return 0;
  bt_shape(c, q, &s1, &s2);
  if (!isfinite(s1) || !isfinite(s2)) return 0;
  if (c->kind <= BT_K4EH) {
    if (s1 < c->klo || s1 > c->khi) return 0;
    if (c->kind >= BT_K4E) {
      for (int i = 0; i < c->n; i++) F[i] = arcef_k4Fe((c->x[i] - q[0]) / sc, s1, s2, NULL);
      return 1;
    }
    for (int i = 0; i < c->n; i++) F[i] = arcef_k4F((c->x[i] - q[0]) / sc, s1, s2, NULL);
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
  if (!c->nodrift && !(xx > 0.0)) return 0;
  if (c->nodrift) { xF = xy = 0.0; xx = HUGE_VAL; }
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
  if (c->kind == BT_WALD) {
    if (!bt_curve(c, q, F)) return 0;
    for (int j = 0; j < dj; j++) {
      double a[BT_DMAX] = {q[0], q[1]}, b[BT_DMAX] = {q[0], q[1]}, *C = J + (size_t)j * n;
      a[j] = fmin(a[j] + 1e-5, c->hi[j]); b[j] = fmax(b[j] - 1e-5, c->lo[j]);
      if (!bt_curve(c, a, C) || !bt_curve(c, b, P)) return 0;
      for (int i = 0; i < n; i++) C[i] = (C[i] - P[i]) / (a[j] - b[j]);
    }
    return 1;
  }
  bt_shape(c, q, &s1, &s2);
  if (c->kind <= BT_K4EH) {
    for (int i = 0; i < n; i++) {
      double t = (c->x[i] - q[0]) / sc, f;
      F[i] = (c->kind >= BT_K4E) ? arcef_k4Fe(t, s1, s2, &f) : arcef_k4F(t, s1, s2, &f);
      J[i] = -f / sc; J[n + i] = -t * f;
    }
    for (int j = 2; j < dj; j++) {
      double *C = J + (size_t)j * n, a[BT_DMAX] = {0}, b[BT_DMAX] = {0};
      for (int l = 0; l < c->d; l++) a[l] = b[l] = q[l];
      double st = (c->kind == BT_FULL) ? 1e-5 : (c->kind == BT_TAB) ? 1e-6 * (c->hi[2] - c->lo[2])
                                                                   : 1e-6 * (1.0 + fabs(q[j]));
      a[j] = fmin(a[j] + st, c->hi[j]); b[j] = fmax(b[j] - st, c->lo[j]);
      double den = a[j] - b[j];
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
    double qa[BT_DMAX] = {q[0], q[1], a, 0}, qb[BT_DMAX] = {q[0], q[1], b, 0}, s0, sa, sb;
    if (a > b) { bt_shape(c, qa, &s0, &sa); bt_shape(c, qb, &s0, &sb); sl = (sa - sb) / (a - b); }
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
      if (c->kind == BT_WALD) {                        /* at most two log units */
        double big = fmax(fabs(st[0]), fabs(st[1]));
        if (big > 2.0) { st[0] *= 2.0 / big; st[1] *= 2.0 / big; }
      }
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
  int j = c->d - 1, tab = (c->nt > 0), i0 = 0, m = BT_NG;
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
            fabs(qb[0] - qe[0]) <= 1e-2 * (d == 2 ? 1.0 : exp(qe[1]));
      for (int k = 1; k < d && dup; k++) dup = fabs(qb[k] - qe[k]) <= 1e-2;
    }
    if (!dup) idx[m++] = t->s0 + b;
  }
  free(dead);
  return m;
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
  if (c->kind == BT_K4E || c->kind == BT_K4EH) { o[3] = q[0] + o[4]; o[4] *= -o[5]; }
}
/* a kappa fit is physical if it rises (g1 > 0) and its tangent reading lies in [0, max x] */
static int bt_physical(const double *o, const double *x, int n) {
  double m = x[0];
  for (int i = 1; i < n; i++) if (x[i] > m) m = x[i];
  return isfinite(o[7]) && o[2] > 0.0 && isfinite(o[8]) && o[8] >= 0.0 && o[8] <= m;
}
static int bt_admit(const bt_ctx *c, const double *q) {
  if (!c->phys) return 1;
  double o[10], rd[3];
  int g = 0;
  bt_theta(c, q, o);
  if (!isfinite(o[7])) return 0;
  double th[6] = {o[0], o[2], o[3], o[4], o[5], o[6]};
  arck4_readings(th, &g, rd);
  o[8] = rd[0];
  return bt_physical(o, c->x, c->n);
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
        if (v[i] < BT_BAD && !bt_admit(&c, q)) v[i] = BT_BAD;
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
        if (cv[i] < BT_BAD && !bt_admit(&c, q)) cv[i] = BT_BAD;
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
        if (task[t].kind <= BT_TAB || task[t].kind == BT_CON) {
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

/* ---- the physical set: the points at which bt_physical holds, an open condition ----
 * Two different things can be asked of it, and they are kept apart. (i) An interior minimum: a local
 * minimum of the sum of squares in the box that lies in the set. (ii) The infimum over the set, which
 * may be approached only at its boundary (the reading running to 0 or to the last time, or the
 * steepest point being lost), and is then not a fit. */
static double bt_bobj(const double *q, const void *vc) {          /* the objective, BT_BAD outside */
  const bt_ctx *c = (const bt_ctx *)vc;
  double f = bt_rss(c, q, c->F);
  return (f < BT_BAD && bt_admit(c, q)) ? f : BT_BAD;
}
/* a simplex whose size is 0.02 of the scale, resp. of a log unit, so that it stays in its basin */
typedef struct { const bt_ctx *c; const double *q0; double d0; } bt_near;
static double bt_nobj(const double *z, const void *vn) {
  const bt_near *w = (const bt_near *)vn;
  double q[BT_DMAX] = {0};
  for (int j = 0; j < w->c->d; j++) q[j] = w->q0[j] + z[j] * (j == 0 ? w->d0 : 0.2);
  return bt_bobj(q, w->c);
}
static double bt_nearpolish(const bt_ctx *c, double *q, double f) {
  double q0[BT_DMAX] = {0}, z0[BT_DMAX] = {0}, z[BT_DMAX], g;
  memcpy(q0, q, (size_t)c->d * sizeof(double));
  bt_near w = {c, q0, 0.2 * exp(q[1])};
  arck4_nmd(bt_nobj, &w, c->d, z0, 400, z, &g);
  if (!(g < f)) return f;
  for (int j = 0; j < c->d; j++) q[j] = q0[j] + z[j] * (j == 0 ? w.d0 : 0.2);
  return g;
}
/* is q in the interior of the set, by a margin of 1e-4 in each of its conditions? On the locus the
 * ends of the searched part are edges too: they are where the search stops, not where the curve does. */
static int bt_inside(const bt_ctx *c, const double *q) {
  double o[10], rd[3], m = c->x[0];
  int g = 0;
  if (c->kind <= BT_TAB && !(q[2] > c->e0 + 1e-6 && q[2] < c->e1 - 1e-6)) return 0;
  bt_theta(c, q, o);
  if (!isfinite(o[7]) || !(o[2] > 0.0)) return 0;
  double th[6] = {o[0], 1.0, o[3], o[4], o[5], o[6]};
  arck4_readings(th, &g, rd);
  for (int i = 1; i < c->n; i++) if (c->x[i] > m) m = c->x[i];
  return isfinite(rd[0]) && rd[0] > 1e-4 * m && rd[0] < (1.0 - 1e-4) * m && o[6] < 1.0 - 1e-4 &&
         o[5] < 1.0 - 1e-4 && o[5] * o[6] < 1.0 - 1e-4;
}
/* descend without the conditions; a minimum counts if every descent ends inside the set. */
static double bt_interior(const bt_ctx *c, double *q, int it) {
  double f = BT_BAD;
  for (int s = 0; s < 8; s++) {
    double g = bt_lm(c, q, c->d, it);
    if (g >= BT_BAD || !bt_admit(c, q)) return BT_BAD;
    double h = bt_nearpolish(c, q, g);
    int done = (g - h <= 1e-10 * g) && (f >= BT_BAD || f - h <= 1e-10 * h || s > 0);
    f = h;
    if (done) return bt_inside(c, q) ? f : BT_BAD;
  }
  return BT_BAD;                                       /* still descending: not a minimum */
}
/* the infimum over the set from q, by simplices that may not leave it */
static double bt_infimum(const bt_ctx *c, double *q) {
  double f = bt_bobj(q, c), p[BT_DMAX], g;
  if (f >= BT_BAD) return f;
  for (int s = 0; s < 12; s++) {
    double f0 = f;
    arck4_nmd(bt_bobj, c, c->d, q, 500, p, &g);
    if (g < f) { f = g; memcpy(q, p, (size_t)c->d * sizeof(double)); }
    f = bt_nearpolish(c, q, f);
    if (f0 - f <= 1e-8 * f0) break;
  }
  return f;
}
/* the lowest interior physical minimum reached from the profile over the shape coordinate (nodes as
 * in bt_nodes), continued from q both ways; BT_BAD if none. q is replaced when one is found. */
static double bt_physnodes(const bt_ctx *c, double *q, int it) {
  int j = c->d - 1, tab = (c->nt > 0), i0 = 0, m = 4 * BT_NG;
  double gr[4 * BT_NG], best = BT_BAD;
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
  if (m < 2) return best;
  double *pf = (double *)malloc((size_t)m * sizeof(double)), *pq = (double *)malloc((size_t)m * 2 * sizeof(double));
  if (pf && pq) {
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
    for (int i = 0; i < m; i++) {
      if (pf[i] >= BT_BAD || (i > 0 && pf[i - 1] < pf[i]) || (i < m - 1 && pf[i + 1] < pf[i])) continue;
      for (int cell = i - 1; cell <= i; cell++) {
        if (tab && (cell < 0 || cell >= m - 1)) continue;
        if (!tab && cell < i) continue;
        bt_ctx cc = *c;
        double p[BT_DMAX] = {pq[2 * i], pq[2 * i + 1], tx[i0 + i], 0};
        if (tab) { cc.lo[j] = tx[i0 + cell]; cc.hi[j] = tx[i0 + cell + 1]; p[j] = 0.5 * (cc.lo[j] + cc.hi[j]); }
        double g = bt_interior(&cc, p, it);
        if (g < best) { best = g; memcpy(q, p, (size_t)c->d * sizeof(double)); }
      }
    }
  }
  free(pf); free(pq);
  return best;
}
/* T tasks with their starts, as bt_run, all of kinds BT_FULL, BT_TAB, BT_K4 or BT_K4H and with
 * c.phys set by mk. iv, iq: the lowest interior physical minimum (HUGE_VAL if none); fv, fq: the
 * lowest value found over the physical set (HUGE_VAL if it has no point). wantmin = 0 skips the
 * first, which then is HUGE_VAL. */
static int bt_physrun(const bt_task *task, int T, const double *st, int M, bt_mk mk, const void *data,
                      const int *off, int it1, int it2, int K, int wantmin,
                      double *iv, double *iq, double *fv, double *fq) {
  double *v = (double *)malloc((size_t)M * sizeof(double)), *v0 = (double *)malloc((size_t)M * sizeof(double));
  double *Q = (double *)malloc((size_t)M * BT_DMAX * sizeof(double));
  size_t N = (size_t)T * 2 * K;
  double *mv = (double *)malloc(N * sizeof(double)), *mq = (double *)malloc(N * BT_DMAX * sizeof(double));
  double *bv = (double *)malloc(N * sizeof(double)), *bq = (double *)malloc(N * BT_DMAX * sizeof(double));
  int *idx = (int *)malloc(N * sizeof(int)), *ow = (int *)malloc((size_t)M * sizeof(int));
  int ok = v && v0 && Q && mv && mq && bv && bq && idx && ow;
  if (ok) {
    for (int t = 0; t < T; t++) for (int s = 0; s < task[t].ns; s++) ow[task[t].s0 + s] = t;
    for (int i = 0; i < M; i++) v[i] = v0[i] = BT_BAD;
    for (size_t i = 0; i < N; i++) { mv[i] = bv[i] = BT_BAD; idx[i] = -1; }
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
        v0[i] = bt_bobj(q, &c);
        v[i] = bt_lm(&c, q, c.d, it1);
        if (v[i] < BT_BAD && !bt_admit(&c, q)) v[i] = BT_BAD;
      }
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 1)
#endif
      for (int t = 0; t < T; t++) {        /* K after some steps, K of the starts themselves */
        bt_ctx c; mk(data, task + t, &c);
        int *id = idx + (size_t)t * 2 * K;
        int m = bt_select(task + t, c.d, v, Q, K, id);
        for (int k = m; k < K; k++) id[k] = -1;
        m = bt_select(task + t, c.d, v0, st, K, id + K);
        for (int k = m; k < K; k++) id[K + k] = -1;
      }
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 1)
#endif
      for (int i = 0; i < (int)N; i++) {
        if (idx[i] < 0) continue;
        int k = i % (2 * K);
        const bt_task *t = task + i / (2 * K);
        int n = off[t->run + 1] - off[t->run];
        W = bt_grow(W, &cap, n);
        if (!W) continue;
        bt_ctx c; mk(data, t, &c); c.F = W; c.W = W + n;
        const double *from = (k < K ? Q : st) + BT_DMAX * (size_t)idx[i];
        double *q = mq + BT_DMAX * (size_t)i, *b = bq + BT_DMAX * (size_t)i;
        memcpy(q, from, BT_DMAX * sizeof(double)); memcpy(b, from, BT_DMAX * sizeof(double));
        if (wantmin) mv[i] = bt_interior(&c, q, it2);
        bv[i] = bt_infimum(&c, b);
      }
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 1)
#endif
      for (int t = 0; t < T; t++) {
        int b = -1, f = -1;
        for (int i = t * 2 * K; i < (t + 1) * 2 * K; i++) {
          if (mv[i] < BT_BAD && (b < 0 || mv[i] < mv[b])) b = i;
          if (bv[i] < BT_BAD && (f < 0 || bv[i] < bv[f])) f = i;
        }
        iv[t] = fv[t] = HUGE_VAL;
        for (int k = 0; k < BT_DMAX; k++) iq[BT_DMAX * (size_t)t + k] = fq[BT_DMAX * (size_t)t + k] = NAN;
        if (f >= 0) { fv[t] = bv[f]; memcpy(fq + BT_DMAX * (size_t)t, bq + BT_DMAX * (size_t)f, BT_DMAX * sizeof(double)); }
        if (b >= 0) { iv[t] = mv[b]; memcpy(iq + BT_DMAX * (size_t)t, mq + BT_DMAX * (size_t)b, BT_DMAX * sizeof(double)); }
        if (!wantmin || task[t].kind > BT_TAB || (f < 0 && b < 0)) continue;
        /* the profile along the locus, from the best point of either kind */
        int n = off[task[t].run + 1] - off[task[t].run];
        W = bt_grow(W, &cap, n);
        if (!W) continue;
        bt_ctx c; mk(data, task + t, &c); c.F = W; c.W = W + n;
        for (int from = 0; from < 2; from++) {
          if ((from == 0 && b < 0) || (from == 1 && f < 0)) continue;
          double q[BT_DMAX];
          memcpy(q, (from == 0 ? iq : fq) + BT_DMAX * (size_t)t, sizeof q);
          double g = bt_physnodes(&c, q, it2);
          if (g < iv[t]) { iv[t] = g; memcpy(iq + BT_DMAX * (size_t)t, q, sizeof q); }
        }
      }
      free(W);
    }
    for (int t = 0; t < T; t++) if (iv[t] < fv[t]) {
      fv[t] = iv[t]; memcpy(fq + BT_DMAX * (size_t)t, iq + BT_DMAX * (size_t)t, BT_DMAX * sizeof(double));
    }
  }
  free(v); free(v0); free(Q); free(mv); free(mq); free(bv); free(bq); free(idx); free(ow);
  return ok;
}
/* 1 interior minimum, the lowest value found; 2 interior minimum, but lower values towards the
 * boundary; 3 physical points but no interior minimum; 0 no physical point. o: status, infimum rss, k,
 * h, tangent reading there */
static void bt_status(const bt_ctx *c, double iv, double fv, const double *fq, double *o) {
  for (int i = 0; i < 5; i++) o[i] = NAN;
  o[0] = !(fv < BT_BAD) ? 0.0 : !(iv < BT_BAD) ? 3.0 : (fv < iv * (1.0 - 1e-6)) ? 2.0 : 1.0;
  if (!(fv < BT_BAD)) return;
  double t[8], rd[3];
  int g = 0;
  bt_theta(c, fq, t);
  double th[6] = {0.0, 1.0, t[3], t[4], t[5], t[6]};
  arck4_readings(th, &g, rd);
  o[1] = fv; o[2] = t[5]; o[3] = t[6]; o[4] = rd[0];
}

/* ==================================== kappa on the locus ====================================== */
typedef struct {
  int nodrift, phys;
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
  c->d = 3; c->kind = t->kind; c->nodrift = D->nodrift; c->phys = D->phys;
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
  c->e0 = c->lo[2]; c->e1 = c->hi[2];
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
                const double *bounds, const int *drift, const int *it1, const int *it2,
                const int *K, const int *gate, double *out, double *outp, double *stat) {
  int R = *nr, np = *npiece > 0 ? *npiece : 2, T = R * np;
  static const double pm[5] = {0.15, 0.3, 0.5, 0.7, 0.9}, ps[3] = {0.04, 0.12, 0.35};
  int nh = 8, ns = 5 * 3 * nh;
  bt_kdat D = {!*drift, *gate, X, Y, off, *npiece, TH, TK, toff, bounds, 0, 0, 0};
  bt_kdat_init(&D);
  bt_task *task = (bt_task *)malloc((size_t)T * sizeof(bt_task));
  double *st = (double *)malloc((size_t)T * ns * BT_DMAX * sizeof(double));
  double *bq = (double *)malloc((size_t)T * BT_DMAX * sizeof(double));
  for (int i = 0; i < 11 * R; i++) out[i] = NAN;
  for (int i = 0; i < 11 * T; i++) outp[i] = NAN;
  for (int i = 0; i < 5 * R; i++) stat[i] = NAN;
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
    if (*gate) {
      /* both halves of the whole locus are searched in the whole box: h = 0 is not a boundary */
      double *fv = (double *)malloc((size_t)T * sizeof(double)), *fq = (double *)malloc((size_t)T * BT_DMAX * sizeof(double));
      if (*npiece <= 0) for (int t = 0; t < T; t++) task[t].piece = 2;
      if (fv && fq && bt_physrun(task, T, st, T * ns, bt_mk_k, &D, off, *it1, *it2, *K, 1, pv, bq, fv, fq)) {
        for (int r = 0; r < R; r++) {
          int b = -1, f = -1, n = off[r + 1] - off[r];
          for (int p = r * np; p < (r + 1) * np; p++) {
            if (pv[p] < BT_BAD && (b < 0 || pv[p] < pv[b])) b = p;
            if (fv[p] < BT_BAD && (f < 0 || fv[p] < fv[f])) f = p;
          }
          double *F = (double *)malloc((size_t)n * sizeof(double));
          if (!F) continue;
          bt_ctx c; bt_mk_k(&D, task + (f >= 0 ? f : r * np), &c); c.F = F;
          bt_status(&c, b >= 0 ? pv[b] : HUGE_VAL, f >= 0 ? fv[f] : HUGE_VAL,
                    fq + BT_DMAX * (size_t)(f >= 0 ? f : 0), stat + 5 * (size_t)r);
          for (int p = 0; p < np; p++) {
            int t = r * np + p;
            if (!(pv[t] < BT_BAD)) continue;
            bt_mk_k(&D, task + t, &c); c.F = F;
            double *op = outp + 11 * (size_t)t;
            bt_theta(&c, bq + BT_DMAX * (size_t)t, op);
            bt_k_readings(op, op + 8);
            op[10] = p + 1;
            if (t == b) {
              memcpy(out + 11 * (size_t)r, op, 11 * sizeof(double));
              if (*npiece <= 0) out[11 * (size_t)r + 10] = op[6] < 0.0 ? 1.0 : 2.0;
            }
          }
          free(F);
        }
      }
      free(fv); free(fq);
    } else if (bt_run(task, T, st, T * ns, bt_mk_k, &D, off, *it1, *it2, *K, pv, bq)) {
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

/* ================================== the free kappa fit ======================================== */
typedef struct { const double *X, *Y; const int *off; double box[4]; } bt_fdat;
static void bt_mk_f(const void *data, const bt_task *t, bt_ctx *c) {
  const bt_fdat *D = (const bt_fdat *)data;
  memset(c, 0, sizeof *c);
  c->x = D->X + D->off[t->run]; c->y = D->Y + D->off[t->run];
  c->n = D->off[t->run + 1] - D->off[t->run];
  c->kind = t->kind; c->d = (t->kind == BT_K4 || t->kind == BT_K4E) ? 4 : 3;
  c->lo[0] = c->lo[1] = -HUGE_VAL; c->hi[0] = c->hi[1] = HUGE_VAL;
  c->lo[2] = c->klo = D->box[0]; c->hi[2] = c->khi = D->box[1];
  c->lo[3] = D->box[2]; c->hi[3] = D->box[3]; c->h0 = 0.0;
  if (t->kind >= BT_K4E) c->hi[2] = fmin(c->hi[2], -0.03);
}
/* the fits of the runs in act[0..na-1], all in one pool; kind BT_K4, or BT_K4H with h = h0 */
typedef struct { bt_fdat f; double h0; int phys, nodrift; } bt_fh;
static void bt_mk_fh(const void *data, const bt_task *t, bt_ctx *c) {
  const bt_fh *D = (const bt_fh *)data;
  bt_mk_f(&D->f, t, c); c->h0 = D->h0; c->phys = D->phys; c->nodrift = D->nodrift;
}
/* starts off the paper's grid, at most 320: scales from a twenty-fifth of the window to three windows
 * and shapes up to the corner (1/2, 1) and beyond h = 1. Some of the lowest minima (falling curves,
 * ends of the support far outside the run) are only reached from here. */
static int bt_wide(const double *x, int n, const double *box, int held, double *st) {
  static const double pm[4] = {0.3, 0.5, 0.7, 0.9}, ps[4] = {0.04, 0.12, 0.35, 3.0};
  static const double pk[4] = {-0.6, -0.02, 0.4, 0.85}, ph[5] = {-20.0, -2.0, 0.4, 0.75, 1.15};
  double *v = (double *)malloc((size_t)n * sizeof(double));
  if (!v) return 0;
  memcpy(v, x, (size_t)n * sizeof(double));
  qsort(v, (size_t)n, sizeof(double), bt_cmp);
  int m = 0;
  for (int a = 0; a < 4; a++) for (int b = 0; b < 4; b++) for (int c = 0; c < 4; c++)
    for (int e = 0; e < (held ? 1 : 5); e++) {
      if (!(pk[c] >= box[0] && pk[c] <= box[1])) continue;
      if (!held && !(ph[e] > box[2] && ph[e] <= box[3])) continue;
      st[0] = bt_q7(v, n, pm[a]); st[1] = log((v[n - 1] - v[0]) * ps[b]); st[2] = pk[c];
      st[3] = held ? 0.0 : ph[e];
      st += BT_DMAX; m++;
    }
  free(v);
  return m;
}
/* with stat (5 na), also the status and infimum of the physical set, against the fit returned */
static int bt_estarts(const double *x, const double *y, int n, const double *box, int held, int dense,
                      double *st) {
  static const double xs[6] = {0.55, 0.75, 0.92, 0.3, 0.65, 1.1}, ts[6] = {0.5, 3.0, 20.0, 0.15, 1.2, 120.0};
  static const double ks[15] = {-28, -14, -9, -6, -4.5, -3.2, -2.5, -1.8, -1.4, -0.95, -0.5, -20, -2, -0.7, -0.25};
  static const double hs[13] = {-12, -1.5, -0.2, 0.15, 0.5, 0.8, -40, -4, -0.6, 0.02, 0.3, 0.65, 0.95};
  if (!(box[0] < -0.15)) return 0;
  double *ys = (double *)malloc((size_t)n * sizeof(double));
  if (!ys) return 0;
  int w = 21, m = 0;
  if (n < w) w = n % 2 ? n : n - 1;
  arck4_runmed(y, &n, &w, ys);
  double lo = ys[0], hi = ys[0], q[3] = {0.25, 0.5, 0.75}, t[3];
  for (int i = 1; i < n; i++) { if (ys[i] < lo) lo = ys[i]; if (ys[i] > hi) hi = ys[i]; }
  for (int j = 0; j < 3; j++) { int i = 0; while (i < n - 1 && ys[i] < lo + q[j] * (hi - lo)) i++; t[j] = x[i]; }
  free(ys);
  double m0 = t[1], s0 = fmax((t[2] - t[0]) / 1.5, 0.02 * (x[n - 1] - x[0]));
  int nx = dense ? 6 : 3, nt = dense ? 6 : 3, nk = dense ? 15 : 11, nh = held ? 1 : (dense ? 13 : 6);
  for (int c = 0; c < nk; c++) {
    if (!(ks[c] >= box[0] && ks[c] <= box[1])) continue;
    for (int e = 0; e < nh; e++) {
      if (!held && !(hs[e] >= box[2] && hs[e] <= box[3])) continue;
      for (int a = 0; a < nx; a++) for (int b = 0; b < nt; b++) for (int w = 0; w < 2; w++) {
        /* w = 1: the lower end of the support well before the data, as on a long scale */
        st[1] = log(s0 * ts[b]); st[0] = x[0] + (m0 - x[0]) * xs[a] - (w ? 0.6 * exp(st[1]) : 0.0);
        st[2] = ks[c]; st[3] = held ? 0.0 : hs[e];
        st += BT_DMAX; m++;
      }
    }
  }
  return m;
}
/* with stat (5 na), also the status and infimum of the physical set, against the fit returned */
static void bt_freefit(const double *X, const double *Y, const int *off, const int *act, int na,
                       const double *box, double hfix, int phys, int drift, int dense, int it1, int it2, int K,
                       double *res /* 10 na */, double *stat) {
  int held = isfinite(hfix), T = 2 * na, per = 700 + (dense ? 72 * 15 * 13 : 1188);
  bt_fh D; D.f.X = X; D.f.Y = Y; D.f.off = off; memcpy(D.f.box, box, sizeof D.f.box); D.h0 = held ? hfix : 0.0; D.phys = phys; D.nodrift = !drift;
  bt_task *task = (bt_task *)malloc((size_t)T * sizeof(bt_task));
  double *st = (double *)malloc((size_t)na * per * BT_DMAX * sizeof(double));
  double *bv = (double *)malloc((size_t)T * sizeof(double));
  double *bq = (double *)malloc((size_t)T * BT_DMAX * sizeof(double));
  for (int i = 0; i < 10 * na; i++) res[i] = NAN;
  if (task && st && bv && bq) {
    int M = 0;
    for (int i = 0; i < na; i++) {
      int r = act[i], n = off[r + 1] - off[r];
      bt_task t = {r, held ? BT_K4H : BT_K4, 0, M, 0};
      t.ns = arck4_starts(X + off[r], Y + off[r], n, box, held, st + BT_DMAX * (size_t)M);
      t.ns += bt_wide(X + off[r], n, box, held, st + BT_DMAX * (size_t)(M + t.ns));
      task[2 * i] = t; M += t.ns;
      bt_task e = {r, held ? BT_K4EH : BT_K4E, 0, M, 0};
      e.ns = bt_estarts(X + off[r], Y + off[r], n, box, held, dense, st + BT_DMAX * (size_t)M);
      task[2 * i + 1] = e; M += e.ns;
    }
    int done = M > 0 && bt_run(task, T, st, M, bt_mk_fh, &D, off, it1, it2, K, bv, bq);
    if (done && stat) {
      double *w = (double *)malloc((size_t)T * (2 + 2 * BT_DMAX) * sizeof(double));
      D.phys = 1;
      double *fv = w ? w + T : NULL, *fq = w ? w + (2 + BT_DMAX) * (size_t)T : NULL;
      if (w && bt_physrun(task, T, st, M, bt_mk_fh, &D, off, it1, it2, K, 0, w, w + 2 * T, fv, fq)) {
        for (int i = 0; i < na; i++) {
          int b = (bv[2 * i + 1] < bv[2 * i]) ? 2 * i + 1 : 2 * i, f = (fv[2 * i + 1] < fv[2 * i]) ? 2 * i + 1 : 2 * i;
          bt_ctx c, cf; bt_mk_fh(&D, task + b, &c); bt_mk_fh(&D, task + f, &cf);
          c.F = cf.F = (double *)malloc((size_t)c.n * sizeof(double));
          if (!c.F) continue;
          /* the fit counts as an interior minimum if it is one; it also bounds the infimum */
          const double *q = bq + BT_DMAX * (size_t)b, *g = fq + BT_DMAX * (size_t)f;
          int in = bv[b] < BT_BAD && bt_inside(&c, q);
          double v = fv[f];
          if (bv[b] < BT_BAD && bt_admit(&c, q) && bv[b] < v)
            bt_status(&c, in ? bv[b] : HUGE_VAL, bv[b], q, stat + 5 * (size_t)i);
          else bt_status(&cf, in ? bv[b] : HUGE_VAL, v, g, stat + 5 * (size_t)i);
          free(c.F);
        }
      }
      free(w);
      D.phys = phys;
    }
    if (done) {
      for (int i = 0; i < na; i++) {
        int b = (bv[2 * i + 1] < bv[2 * i]) ? 2 * i + 1 : 2 * i;
        if (!isfinite(bv[b]) || bv[b] >= BT_BAD) continue;
        bt_ctx c; bt_mk_fh(&D, task + b, &c);
        c.F = (double *)malloc((size_t)c.n * sizeof(double));
        if (!c.F) continue;
        double *o = res + 10 * (size_t)i, rd[3];
        bt_theta(&c, bq + BT_DMAX * (size_t)b, o);
        double th[6] = {o[0], o[2], o[3], o[4], o[5], o[6]};
        int g = 0;
        if (isfinite(o[7])) { arck4_readings(th, &g, rd); o[8] = rd[0]; o[9] = rd[1]; }
        free(c.F);
      }
    }
  }
  free(task); free(st); free(bv); free(bq);
}
/* Free fits with drift of all runs. floors: nfl lower bounds on h, deepest first. gate = 0: the fit
 * in the box with the first floor. gate = 1: the first floor whose best fit is physical, the last
 * floor's fit if none is. gate = 2: the lowest PHYSICAL minimum in the box with the first floor; NaN
 * where none was found. hfix finite: h is held there and floors and gate are not used. gevsel: the
 * h = 0 submodel is fitted too and chosen when n log(rss0/rss) - 2 < 0 (and, with gate, when it is
 * physical). out: 14 per run (g0, m, g1, mu, sg, k, h, rss, a, b, gev, dAIC, floor, physical). */
void arcbt_k4fit(const double *X, const double *Y, const int *off, const int *nr, const double *bounds,
                 const double *floors, const int *nfl, const int *gate, const int *gevsel,
                 const double *hfix, const int *drift, const int *it1, const int *it2, const int *K,
                 const int *inf, const int *dense, double *out, double *stat) {
  int R = *nr, held = isfinite(*hfix);
  int *act = (int *)malloc((size_t)R * sizeof(int)), na = R;
  double *res = (double *)malloc((size_t)R * 10 * sizeof(double));
  for (int i = 0; i < 14 * R; i++) out[i] = NAN;
  for (int i = 0; i < 5 * R; i++) stat[i] = NAN;
  if (act && res) {
    for (int r = 0; r < R; r++) act[r] = r;
    int nf = (held || *gate != 1) ? 1 : *nfl;
    for (int f = 0; f < nf && na > 0; f++) {
      double box[4] = {bounds[0], bounds[1], held ? bounds[2] : floors[f], bounds[3]};
      bt_freefit(X, Y, off, act, na, box, *hfix, *gate == 2, *drift, *dense, *it1, *it2, *K, res, *gate == 2 && *inf ? stat : NULL);
      int nn = 0;
      for (int i = 0; i < na; i++) {
        int r = act[i];
        double *o = out + 14 * (size_t)r;
        memcpy(o, res + 10 * (size_t)i, 10 * sizeof(double));
        o[10] = 0.0; o[12] = held ? bounds[2] : floors[f];
        o[13] = bt_physical(o, X + off[r], off[r + 1] - off[r]);
        if (o[13] == 0.0) act[nn++] = r;
      }
      na = *gate ? nn : 0;
    }
    if (*gevsel && !held) {
      for (int r = 0; r < R; r++) act[r] = r;
      bt_freefit(X, Y, off, act, R, bounds, 0.0, 0, *drift, *dense, *it1, *it2, *K, res, NULL);
      for (int r = 0; r < R; r++) {
        double *o = out + 14 * (size_t)r, *g = res + 10 * (size_t)r;
        int n = off[r + 1] - off[r], ph = bt_physical(g, X + off[r], n);
        if (!isfinite(g[7])) continue;
        if (isfinite(o[7])) o[11] = n * log(g[7] / o[7]) - 2.0;
        if (!isfinite(o[7]) || (o[11] < 0.0 && (ph || !*gate))) {
          memcpy(o, g, 10 * sizeof(double));
          o[10] = 1.0; o[13] = ph;
        }
      }
    }
  }
  free(act); free(res);
}
/* one run, for arcmv.c: the lowest physical minimum (the lowest of all if none is physical) from the
 * start grid, in the layout of arck4_fit_varpro */
int arcbt_k4one(const double *x, const double *y, int n, const double *box, double hfix, double *o8) {
  int off[2] = {0, n}, act[1] = {0};
  double r[10];
  bt_freefit(x, y, off, act, 1, box, hfix, 1, 1, 0, 12, 200, 8, r, NULL);
  if (!isfinite(r[7])) bt_freefit(x, y, off, act, 1, box, hfix, 0, 1, 0, 12, 200, 8, r, NULL);
  for (int i = 0; i < 8; i++) o8[i] = r[i];
  return isfinite(r[7]) && r[7] < 1e11;
}

/* ==================================== the companion host ====================================== */
typedef struct { const double *X, *Y; const int *off; const double *cal, *cbe; int nc; double eps, alo, ahi; } bt_hdat;
static void bt_mk_h(const void *data, const bt_task *t, bt_ctx *c) {
  const bt_hdat *D = (const bt_hdat *)data;
  memset(c, 0, sizeof *c);
  c->x = D->X + D->off[t->run]; c->y = D->Y + D->off[t->run];
  c->n = D->off[t->run + 1] - D->off[t->run];
  c->kind = t->kind;
  c->lo[0] = c->lo[1] = -HUGE_VAL; c->hi[0] = c->hi[1] = HUGE_VAL;
  if (t->kind == BT_CON) {
    c->d = 3; c->tx = D->cal; c->ty = D->cbe; c->nt = D->nc;
    c->lo[2] = D->alo; c->hi[2] = D->ahi;
  } else {                                             /* alpha, beta in [-1 + eps, 20] */
    c->d = 4;
    c->lo[2] = c->lo[3] = log(D->eps); c->hi[2] = c->hi[3] = log(21.0);
  }
}
/* outf, outc: 10 per run (g0, m, g1, s, w, alpha, beta, rss, a, b), the free and the constrained fit.
 * The constrained optimum is a point of the free family, so the free search also starts there. */
void arcbt_host(const double *X, const double *Y, const int *off, const int *nr,
                const double *cal, const double *cbe, const int *nc, const double *eps,
                const int *dense, const int *it1, const int *it2, const int *K,
                double *outf, double *outc) {
  int R = *nr, T = 2 * R;
  static const double ps[6] = {0.02, 0.12, 0.3, 0.005, 0.06, 0.5}, pw[6] = {0.5, 0.8, 1.2, 3, 30, 1000};
  static const double fa[8] = {-0.85, -0.6, -0.35, -0.1, -0.95, -0.5, -0.2, 0.3};
  static const double fb[6] = {-0.85, -0.5, -0.15, -0.99, -0.9999, 0.5};
  static const double ca[8] = {0.05, 0.3, 0.6, 0.9, 0.0, 0.15, 0.45, 0.75};
  int dn = *dense ? 2 : 1, nf = 108 * dn * dn * dn * dn, nk = 36 * dn * dn * dn, M = R * (nf + nk);
  bt_hdat D = {X, Y, off, cal, cbe, *nc, *eps, 0, 0};
  if (*nc > 0) { D.alo = cal[0]; D.ahi = cal[*nc - 1]; }
  else {
    /* the whole curve: from the alpha whose beta is -1 + eps to just short of alpha = -1/2, where
     * beta -> 0; found by bisection, beta increasing in alpha */
    double a = -0.645, b = -0.6, v;
    for (int i = 0; i < 60; i++) {
      double m = 0.5 * (a + b);
      arceq2_bstar(&m, &v);
      if (isfinite(v) && v >= -1.0 + *eps) b = m; else a = m;
    }
    D.alo = b; D.ahi = -0.5 - 1e-6;
  }
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
      double qs[6] = {0}, rg = 1;
      double *v = (double *)malloc((size_t)n * sizeof(double));
      if (v) {
        memcpy(v, X + off[r], (size_t)n * sizeof(double));
        qsort(v, (size_t)n, sizeof(double), bt_cmp);
        for (int i = 0; i < 6; i++) qs[i] = bt_q7(v, n, ps[i]);
        rg = v[n - 1] - v[0];
        free(v);
      }
      double *p = st + BT_DMAX * (size_t)s;
      for (int i = 0; i < 3 * dn; i++) for (int j = 0; j < 3 * dn; j++) for (int a = 0; a < 4 * dn; a++)
        for (int b = 0; b < 3 * dn; b++) {
          p[0] = qs[i]; p[1] = log(rg * pw[j]); p[2] = log1p(fa[a]); p[3] = log1p(fb[b]); p += BT_DMAX;
        }
      for (int i = 0; i < 3 * dn; i++) for (int j = 0; j < 3 * dn; j++) for (int a = 0; a < 4 * dn; a++) {
        p[0] = qs[i]; p[1] = log(rg * pw[j]); p[2] = D.alo + ca[a] * (D.ahi - D.alo); p += BT_DMAX;
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

/* ========================================= the Wald curve ====================================== */
typedef struct { const double *X, *Y; const int *off; int nodrift; } bt_wdat;
static void bt_mk_w(const void *data, const bt_task *t, bt_ctx *c) {
  const bt_wdat *D = (const bt_wdat *)data;
  memset(c, 0, sizeof *c);
  c->x = D->X + D->off[t->run]; c->y = D->Y + D->off[t->run];
  c->n = D->off[t->run + 1] - D->off[t->run];
  c->kind = BT_WALD; c->d = 2; c->nodrift = D->nodrift;
  c->lo[0] = c->lo[1] = log(1e-8); c->hi[0] = log(1e8); c->hi[1] = log(1e12);
}
/* out: 10 per run (g0, m, g1, mu, lambda, phi, rss, a, b, mode); the readings are those of the fitted
 * shape, whatever the sign of g1 and wherever they lie */
void arcbt_wald(const double *X, const double *Y, const int *off, const int *nr, const int *drift,
                const int *it1, const int *it2, const int *K, double *out) {
  int R = *nr;
  static const double m0[7] = {2, 4, 6, 9, 13, 20, 35}, l0[7] = {2, 10, 50, 200, 800, 3000, 12000};
  static const double fr[5] = {0.1, 0.25, 0.5, 1.0, 2.5};
  bt_wdat D = {X, Y, off, !*drift};
  bt_task *task = (bt_task *)malloc((size_t)R * sizeof(bt_task));
  double *st = (double *)calloc((size_t)R * 154 * BT_DMAX, sizeof(double));
  double *bv = (double *)malloc((size_t)R * sizeof(double));
  double *bq = (double *)malloc((size_t)R * BT_DMAX * sizeof(double));
  for (int i = 0; i < 10 * R; i++) out[i] = NAN;
  if (task && st && bv && bq) {
    int M = 0;
    for (int r = 0; r < R; r++) {
      const double *x = X + off[r], *y = Y + off[r];
      int n = off[r + 1] - off[r], w = n < 21 ? (n % 2 ? n : n - 1) : 21, nm = 0, k = 0;
      double mm[22], *ys = (double *)malloc((size_t)n * sizeof(double)), lo = y[0], hi = y[0], xm = x[0], md = 0.0;
      for (int i = 1; i < n; i++) { if (y[i] < lo) lo = y[i]; if (y[i] > hi) hi = y[i]; if (x[i] > xm) xm = x[i]; }
      if (ys) {
        arck4_runmed(y, &n, &w, ys);
        for (int i = 0; i < n; i++) if (ys[i] >= lo + 0.5 * (hi - lo)) { md = x[i]; break; }
        free(ys);
      }
      for (int i = 0; i < 7; i++) mm[nm++] = m0[i];       /* the grid of the induction papers, in hours */
      mm[nm++] = md; mm[nm++] = 0.5 * md; mm[nm++] = 2.0 * md;
      for (int i = 0; i < 5; i++) mm[nm++] = fr[i] * xm;  /* and one that does not depend on the unit */
      bt_task t = {r, BT_WALD, 0, M, 0};
      for (int i = 0; i < nm; i++) for (int j = 0; j < 7; j++) {
        if (!(mm[i] > 0.0)) continue;
        double *p = st + BT_DMAX * (size_t)(M + k);
        p[0] = log(mm[i]); p[1] = log(i < 10 ? l0[j] : mm[i] * (j < 6 ? 0.1 * pow(10.0, j) : 0.5)); k++;
      }
      t.ns = k; task[r] = t; M += k;
    }
    if (bt_run(task, R, st, M, bt_mk_w, &D, off, *it1, *it2, *K, bv, bq)) {
      for (int r = 0; r < R; r++) {
        if (!isfinite(bv[r]) || bv[r] >= BT_BAD) continue;
        bt_ctx c; bt_mk_w(&D, task + r, &c);
        c.F = (double *)malloc((size_t)c.n * sizeof(double));
        if (!c.F) continue;
        double t[8], rd[3], *o = out + 10 * (size_t)r, mu = exp(bq[BT_DMAX * (size_t)r]);
        bt_theta(&c, bq + BT_DMAX * (size_t)r, t);
        o[0] = t[0]; o[1] = t[1]; o[2] = t[2]; o[3] = mu; o[4] = t[4]; o[5] = t[4] / mu; o[6] = t[7];
        bt_wald_read(o[5], rd);
        for (int i = 0; i < 3; i++) o[7 + i] = mu * rd[i];
        free(c.F);
      }
    }
  }
  free(task); free(st); free(bv); free(bq);
}
void arcbt_wald_cdf(const double *x, const int *n, const double *mu, const double *lam, double *out) {
  for (int i = 0; i < *n; i++) out[i] = isnan(x[i]) ? NAN : bt_waldF(x[i] / *mu, *lam / *mu);
}
/* out: 3 per pair (a, b, mode) */
void arcbt_wald_readings(const double *mu, const double *lam, const int *n, double *out) {
  for (int i = 0; i < *n; i++) {
    bt_wald_read(lam[i] / mu[i], out + 3 * (size_t)i);
    for (int j = 0; j < 3; j++) out[3 * (size_t)i + j] *= mu[i];
  }
}

/* ============================ windowed L-moment matching of a curve ============================= */
/* ratios of the L-moments of the quantile values q at levels u, both restricted to the window the
 * levels span and centred on it; left Riemann sums */
static void bt_tauw(const double *u, const double *q, int n, double *t) {
  t[0] = t[1] = NAN;
  if (n < 50) return;
  double sw = 0, sq = 0, c2 = 0, c3 = 0, c4 = 0;
  for (int i = 0; i < n; i++) {
    double w = (i < n - 1) ? u[i + 1] - u[i] : u[n - 1] - u[n - 2];
    sw += w; sq += w * q[i];
  }
  double m = sq / sw;
  for (int i = 0; i < n; i++) {
    double w = (i < n - 1) ? u[i + 1] - u[i] : u[n - 1] - u[n - 2], d = w * (q[i] - m), v = u[i];
    c2 += d * (2 * v - 1); c3 += d * (6 * v * v - 6 * v + 1); c4 += d * (20 * v * v * v - 30 * v * v + 12 * v - 1);
  }
  if (!isfinite(c2) || fabs(c2) < 1e-12) return;
  t[0] = c3 / c2; t[1] = c4 / c2;
}
typedef struct { double t3, t4; const double *u; int n; const double *bnd; double *q, *v; } bt_lctx;
static void bt_model_tau(const bt_lctx *c, double k, double h, double *t) {
  double mu = 0.0, sg = 1.0;
  int m = 0;
  arck4_q(c->u, &c->n, &mu, &sg, &k, &h, c->q);
  for (int i = 0; i < c->n; i++) if (isfinite(c->q[i])) { c->q[m] = c->q[i]; c->v[m] = c->u[i]; m++; }
  bt_tauw(c->v, c->q, m, t);
}
static double bt_lobj(const double *p, const void *vc) {
  const bt_lctx *c = (const bt_lctx *)vc;
  double t[2];
  if (!(p[0] >= c->bnd[0] && p[0] <= c->bnd[1] && p[1] >= c->bnd[2] && p[1] <= c->bnd[3])) return 1e12;
  bt_model_tau(c, p[0], p[1], t);
  if (!isfinite(t[0]) || !isfinite(t[1])) return 1e12;
  return (t[0] - c->t3) * (t[0] - c->t3) + (t[1] - c->t4) * (t[1] - c->t4);
}
/* damped Gauss-Newton on the two equations, p in place; its end point can lie on the box */
static double bt_lgn(const bt_lctx *c, double *p) {
  double f = bt_lobj(p, c), t0[2], lam = 1e-3;
  if (f >= 1e12) return f;
  for (int it = 0; it < 40 && f > 1e-24; it++) {
    double J[2][2], r[2];
    bt_model_tau(c, p[0], p[1], t0);
    r[0] = c->t3 - t0[0]; r[1] = c->t4 - t0[1];
    for (int j = 0; j < 2; j++) {
      double h = 1e-6, q[2] = {p[0], p[1]}, t1[2];
      if (q[j] + h > c->bnd[2 * j + 1]) h = -h;
      q[j] += h;
      bt_model_tau(c, q[0], q[1], t1);
      J[0][j] = (t1[0] - t0[0]) / h; J[1][j] = (t1[1] - t0[1]) / h;
    }
    double a = J[0][0] * J[0][0] + J[1][0] * J[1][0], b = J[0][0] * J[0][1] + J[1][0] * J[1][1];
    double d = J[0][1] * J[0][1] + J[1][1] * J[1][1], g0 = J[0][0] * r[0] + J[1][0] * r[1], g1 = J[0][1] * r[0] + J[1][1] * r[1];
    int ok = 0;
    for (int tr = 0; tr < 10 && !ok; tr++) {
      double A = a * (1 + lam), D = d * (1 + lam), det = A * D - b * b;
      double q[2] = {p[0] + (D * g0 - b * g1) / det, p[1] + (A * g1 - b * g0) / det};
      if (!isfinite(q[0]) || !isfinite(q[1])) break;
      for (int j = 0; j < 2; j++) { q[j] = fmin(fmax(q[j], c->bnd[2 * j]), c->bnd[2 * j + 1]); }
      double fn = bt_lobj(q, c);
      if (fn < f) { p[0] = q[0]; p[1] = q[1]; ok = (f - fn > 1e-9 * f) ? 1 : 2; f = fn; lam *= 0.2; }
      else lam *= 8.0;
    }
    if (ok != 1) break;
  }
  return f;
}
static double bt_med(double *b, int m) { qsort(b, (size_t)m, sizeof(double), bt_cmp); return b[m / 2]; }
/* running median with the end rule of R's runmed (endrule "median") */
static void bt_runmed(const double *y, int n, int w, double *o, double *b) {
  int h = w / 2;
  for (int i = 0; i < n; i++) {
    if (i < h || i >= n - h) { o[i] = y[i]; continue; }
    memcpy(b, y + i - h, (size_t)w * sizeof(double));
    o[i] = bt_med(b, w);
  }
  if (h < 1 || n < w) return;
  double *t = b + w;                                     /* the ends are smoothed from o itself */
  memcpy(t, o, (size_t)n * sizeof(double));
  for (int j = 1; j < h; j++) {
    memcpy(b, t, (size_t)(2 * j + 1) * sizeof(double)); o[j] = bt_med(b, 2 * j + 1);
    memcpy(b, t + n - 2 * j - 1, (size_t)(2 * j + 1) * sizeof(double)); o[n - 1 - j] = bt_med(b, 2 * j + 1);
  }
  double e[3] = {t[0], o[1], 3 * o[1] - 2 * o[2]}, f[3] = {t[n - 1], o[n - 2], 3 * o[n - 2] - 2 * o[n - 3]};
  o[0] = bt_med(e, 3); o[n - 1] = bt_med(f, 3);
}
/* Each run's curve is put on the probability scale by its linear part lin = (g0, m, g1) only: running
 * median, (y - g0 - m x)/g1 held in [0, 1] and made nondecreasing. Its third and fourth L-moment
 * ratios over the observed window [a, b] of levels (inside [0.005, 0.995]) are matched by those of
 * the kappa quantile function over the same window. ng, nm: levels for the data and the model. out: 10
 * per run (k, h, objective, t3, t4, a, b, number of exact matches, their smallest and largest h); NaN where fewer than 100 points are inside. */
void arcbt_lmomw(const double *X, const double *Y, const int *off, const int *nr, const double *lin,
                 const double *bounds, const int *window, const int *ng, const int *nm, double *out) {
  int R = *nr;
  static const double k0[4] = {-0.8, -0.4, 0, 0.3}, h0[6] = {-1.5, -0.4, -0.05, 0.05, 0.4, 1};
  for (int i = 0; i < 10 * R; i++) out[i] = NAN;
  double *um = (double *)malloc((size_t)R * *nm * sizeof(double));
  double *fv = (double *)malloc((size_t)R * 24 * sizeof(double)), *fp = (double *)malloc((size_t)R * 48 * sizeof(double));
  double *gv = (double *)malloc((size_t)R * 6 * sizeof(double)), *gp = (double *)malloc((size_t)R * 12 * sizeof(double));
  if (um && fv && fp && gv && gp) {
    for (int i = 0; i < 24 * R; i++) fv[i] = BT_BAD;
#ifdef _OPENMP
#pragma omp parallel num_threads(arc_threads())
#endif
    {
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 1)
#endif
      for (int r = 0; r < R; r++) {
        int n = off[r + 1] - off[r], m = 0, w = *window, G = *ng;
        const double *x = X + off[r], *L = lin + 3 * (size_t)r;
        double *o = out + 10 * (size_t)r;
        double *ys = (double *)malloc((size_t)(4 * n + w + 2 * G) * sizeof(double));
        if (!ys) continue;
        double *F = ys + n, *b = ys + 2 * n, *u = ys + 4 * n + w, *q = u + G, top = 0.0;
        bt_runmed(Y + off[r], n, w, ys, b);
        double *xx = b;                                  /* reused: x of the kept levels */
        for (int i = 0; i < n; i++) {
          double f = (ys[i] - L[0] - L[1] * x[i]) / L[2];
          f = f < 0.0 ? 0.0 : (f > 1.0 ? 1.0 : f);
          if (f > top) top = f;
          if (!(top > 1e-6 && top < 1.0 - 1e-6)) continue;
          F[m] = top; xx[m] = x[i]; m++;
        }
        if (m >= 100 && F[m - 1] > F[0]) {               /* interpolation as approx(ties = "ordered") */
          double a = fmax(0.005, F[0]), e = fmin(0.995, F[m - 1]);
          int j = 0, g = 0;
          for (int i = 0; i < G; i++) {
            double v = (i == G - 1) ? e : a + (e - a) * i / (G - 1.0);
            if (v < F[0] || v > F[m - 1]) continue;
            while (j < m - 2 && F[j + 1] <= v) j++;
            u[g] = v;
            q[g] = (v >= F[j + 1]) ? xx[j + 1] : xx[j] + (xx[j + 1] - xx[j]) * (v - F[j]) / (F[j + 1] - F[j]);
            g++;
          }
          bt_tauw(u, q, g, o + 3);
          o[5] = a; o[6] = e;
          for (int i = 0; i < *nm; i++) um[(size_t)r * *nm + i] = a + (e - a) * i / (*nm - 1.0);
        }
        free(ys);
      }
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 1)
#endif
      for (int t = 0; t < 24 * R; t++) {
        int r = t / 24, s = t % 24;
        const double *o = out + 10 * (size_t)r;
        if (!isfinite(o[3]) || !isfinite(o[4])) continue;
        double *q = (double *)malloc((size_t)*nm * 2 * sizeof(double));
        if (!q) continue;
        bt_lctx c = {o[3], o[4], um + (size_t)r * *nm, *nm, bounds, q, q + *nm};
        double st[2] = {k0[s / 6], h0[s % 6]};
        fv[t] = bt_lgn(&c, st);
        memcpy(fp + 2 * (size_t)t, st, sizeof st);
        free(q);
      }
      /* the six lowest of each run are carried to convergence */
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 1)
#endif
      for (int t = 0; t < 6 * R; t++) {
        int r = t / 6, rank = t % 6, b = -1;
        const double *o = out + 10 * (size_t)r, *v = fv + 24 * (size_t)r;
        for (int s = 0; s < 24 && b < 0; s++) {
          int c = 0;
          for (int l = 0; l < 24; l++) if (v[l] < v[s] || (v[l] == v[s] && l < s)) c++;
          if (c == rank) b = s;
        }
        gv[t] = BT_BAD;
        if (b < 0 || !(v[b] < 1e11)) continue;
        double *q = (double *)malloc((size_t)*nm * 2 * sizeof(double));
        if (!q) continue;
        bt_lctx c = {o[3], o[4], um + (size_t)r * *nm, *nm, bounds, q, q + *nm};
        double p1[2], f1, p2[2], f2;
        arck4_nmd(bt_lobj, &c, 2, fp + 2 * (size_t)(24 * r + b), 300, p1, &f1);
        arck4_nmd(bt_lobj, &c, 2, p1, 300, p2, &f2);
        gv[t] = f2 < f1 ? f2 : f1;
        memcpy(gp + 2 * (size_t)t, f2 < f1 ? p2 : p1, sizeof p1);
        free(q);
      }
    }
    for (int r = 0; r < R; r++) {
      /* the 24 Gauss-Newton points and the six polished ones; two ratios do not always determine the
       * shape, so every distinct exact match is counted and the one with the smallest |h| returned */
      double V[30], P[60], sol[60], *o = out + 10 * (size_t)r;
      int b = -1, ns = 0;
      for (int s = 0; s < 30; s++) {
        V[s] = s < 24 ? fv[24 * r + s] : gv[6 * r + s - 24];
        for (int j = 0; j < 2; j++) P[2 * s + j] = s < 24 ? fp[2 * (24 * r + s) + j] : gp[2 * (6 * r + s - 24) + j];
        if (V[s] < 1e11 && (b < 0 || V[s] < V[b])) b = s;
      }
      if (b < 0) continue;
      for (int s = 0; s < 30; s++) {
        if (!(V[s] < 1e-14)) continue;
        int dup = 0;
        for (int e = 0; e < ns && !dup; e++) dup = fabs(P[2 * s] - sol[2 * e]) < 1e-3 && fabs(P[2 * s + 1] - sol[2 * e + 1]) < 1e-3;
        if (!dup) { sol[2 * ns] = P[2 * s]; sol[2 * ns + 1] = P[2 * s + 1]; ns++; }
        if (s != b && V[b] < 1e-14 && fabs(P[2 * s + 1]) < fabs(P[2 * b + 1]) - 1e-3) b = s;
      }
      o[0] = P[2 * b]; o[1] = P[2 * b + 1]; o[2] = V[b]; o[7] = ns; o[8] = o[9] = o[1];
      for (int e = 0; e < ns; e++) { o[8] = fmin(o[8], sol[2 * e + 1]); o[9] = fmax(o[9], sol[2 * e + 1]); }
    }
  }
  free(um); free(fv); free(fp); free(gv); free(gp);
}

/* ======================================== bootstraps ========================================== */
static int bt_bootlen(const double *r, int n, int block) {
  int b = block > 0 ? block : arcef_blocklen(r, n);
  return b > n ? n : b;
}
/* Moving-block residual bootstrap of a fitted kappa curve, many runs in one call.
 * fit: 7 per run (g0, m, g1, mu, sg, k, h); box: 4 per run, the box of that run's own point fit;
 * gev[r]: the point fit is the h = 0 submodel, so h is held at 0 in the replicates; drift[r]: the
 * point fit carries the linear drift m x (0: m is held at 0 in the mean the residuals are taken from
 * and in every refit, as in the point fit); gate: 1 keeps a replicate only if its refit is physical
 * (g1 > 0, tangent reading in [0, max x]), the rule of a point fit made with gate = "best"; block: 0
 * reads the block length off the residual autocorrelation, 1 resamples singly. Every replicate is
 * refitted by arck4_fit_varpro_dr from the fitted shape. out: 6 x B per run (a, b, k, h, sel, m): the
 * readings, the fitted shape, whether n log(rss_0/rss) - 2 < 0 would choose h = 0 (NaN where h
 * is held), and the fitted slope (exactly 0 without drift). bl: the block length used. seed: one per run; replicate b of a run depends on (seed, b)
 * only. Every setting of the point fit reaches the replicates: box, submodel, drift and gate. */
void arcbt_k4_boot(const double *X, const double *Y, const int *off, const int *nr,
                   const double *fit, const double *box, const int *gev, const int *drift,
                   const int *gate, const int *B, const int *seed, const int *block, const int *maxit,
                   int *bl, double *out) {
  int R = *nr, NB = *B, N = off[R];
  for (long i = 0; i < 6L * NB * R; i++) out[i] = NAN;
  double *m0 = (double *)malloc((size_t)N * sizeof(double));
  double *r0 = (double *)malloc((size_t)N * sizeof(double));
  if (m0 && r0) {
    for (int r = 0; r < R; r++) {
      const double *f = fit + 7 * (size_t)r;
      for (int i = off[r]; i < off[r + 1]; i++) {
        m0[i] = f[0] + (drift[r] ? f[1] * X[i] : 0.0) + f[2] * k4_F(X[i], f[3], f[4], f[5], f[6]);
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
      double *res = out + 6 * ((size_t)r * NB + b);
      arck4_fit_varpro_dr(x, yb, &n, p0, &one, maxit, box + 4 * (size_t)r, &hold, drift[r], o);
      if (isfinite(o[7]) && o[7] < 1e11) {
        double th[6] = {o[0], o[2], o[3], o[4], o[5], o[6]};
        arck4_readings(th, &g, rd);
        double ph[10] = {o[0], o[1], o[2], o[3], o[4], o[5], o[6], o[7], rd[0], rd[1]};
        if (*gate && !bt_physical(ph, x, n)) { free(yb); continue; }
        res[0] = rd[0]; res[1] = rd[1]; res[2] = o[5]; res[3] = o[6]; res[5] = o[1];
        if (!gev[r]) {
          p0[3] = 0.0;
          arck4_fit_varpro_dr(x, yb, &n, p0, &one, maxit, box + 4 * (size_t)r, &zero, drift[r], og);
          if (isfinite(og[7])) res[4] = (n * log(og[7] / o[7]) - 2.0 < 0.0) ? 1.0 : 0.0;
        }
      }
      free(yb);
    }
  }
  free(m0); free(r0);
}
/* The same for the constrained fit: fit as returned by arcbt_k4eq (its first seven numbers), piece
 * the table piece it lies on (1-based; 1 for the full locus), drift[r] whether the point fit carries
 * the drift m x (0: m is held at 0 in the mean and in every refit; until 0.4.0 every replicate was
 * refitted WITH drift whatever the point fit was). Replicates are refitted on that piece, from the
 * fitted point, under the point fit's bounds, gate and drift. out: 5 x B per run (a, b, k, h, m), m the
 * fitted slope (exactly 0 without drift). */
void arcbt_k4eq_boot(const double *X, const double *Y, const int *off, const int *nr,
                     const double *fit, const int *piece, const int *npiece, const double *TH,
                     const double *TK, const int *toff, const double *bounds, const int *B,
                     const int *seed, const int *block, const int *it2, const int *gate,
                     const int *drift, int *bl, double *out) {
  int R = *nr, NB = *B, N = off[R];
  for (long i = 0; i < 5L * NB * R; i++) out[i] = NAN;
  double *m0 = (double *)malloc((size_t)N * sizeof(double));
  double *r0 = (double *)malloc((size_t)N * sizeof(double));
  if (m0 && r0) {
    bt_kdat D = {0, *gate, X, Y, off, *npiece, TH, TK, toff, bounds, 0, 0, 0};   /* nodrift: per run, below */
    bt_kdat_init(&D);
    for (int r = 0; r < R; r++) {
      const double *f = fit + 7 * (size_t)r;
      for (int i = off[r]; i < off[r + 1]; i++) {
        m0[i] = f[0] + (drift[r] ? f[1] * X[i] : 0.0) + f[2] * arcef_k4F((X[i] - f[3]) / f[4], f[5], f[6], NULL);
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
        double *res = out + 5 * ((size_t)r * NB + b);
        if (!(isfinite(f[4]) && f[4] > 0.0 && isfinite(f[6]))) continue;
        W = bt_grow(W, &cap, n);
        double *yb = (double *)malloc((size_t)n * sizeof(double));
        if (W && yb) {
          unsigned long long s = arcef_seed(seed[r], b);
          arcef_blockdraw(m0 + off[r], r0 + off[r], n, bl[r], &s, yb);
          bt_task tk = {r, *npiece > 0 ? BT_TAB : BT_FULL, *npiece > 0 ? piece[r] - 1 : 2, 0, 0};
          bt_ctx c; bt_mk_k(&D, &tk, &c);
          c.nodrift = !drift[r];
          c.y = yb; c.F = W; c.W = W + n;
          double q[BT_DMAX] = {f[3], log(f[4]), *npiece > 0 ? f[6] : log(0.5 - f[6]), 0};
          double v;
          if (*gate) v = (*npiece > 0) ? bt_physnodes(&c, q, *it2) : bt_interior(&c, q, *it2);
          else {
            v = bt_converge(&c, q, *it2, 1);
            if (v < BT_BAD && *npiece > 0) v = bt_nodes(&c, q, v, *it2);
          }
          if (v < 1e11) {
            double o[8];
            bt_theta(&c, q, o);
            bt_k_readings(o, res);
            res[2] = o[5]; res[3] = o[6]; res[4] = o[1];
          }
        }
        free(yb);
      }
      free(W);
    }
  }
  free(m0); free(r0);
}

/* Variance components of the pipelines of a multiverse: A is B x M (rows replicates). Column-wise the
 * pipelines with at least frac*B finite values are kept, then the replicates finite in all of them.
 * With S the M x M covariance of the kept columns over the B' kept rows, v = mean diagonal, w = mean
 * diagonal minus mean off-diagonal of S (the method-by-replicate variance, because the replicates are
 * shared), and m the variance of the M column means:
 *   E m = delta^2 + w/B',  delta^2 = the variance of the pipelines' differences GIVEN THE DATA,
 * so delta^2 = max(0, m - w/B'). The differences given the data contain the data's own sampling error
 * once more: E delta^2 = tau^2 + w, so tau^2 = max(0, m - w (1 + 1/B')).
 * out: M kept, B', v, w, delta^2, delta^2/(delta^2 + v), tau^2, tau^2/(tau^2 + v). */
void arcbt_mvshare(const double *A, const int *nb, const int *nm, const double *frac, double *out) {
  int B = *nb, M = *nm, k = 0, b = 0;
  for (int i = 0; i < 8; i++) out[i] = NAN;
  int *keep = (int *)malloc((size_t)M * sizeof(int)), *row = (int *)malloc((size_t)B * sizeof(int));
  if (keep && row) {
    for (int m = 0; m < M; m++) {
      int c = 0;
      for (int i = 0; i < B; i++) c += isfinite(A[i * M + m]) ? 1 : 0;
      if (c >= *frac * B) keep[k++] = m;
    }
    for (int i = 0; i < B; i++) {
      int ok = 1;
      for (int j = 0; j < k; j++) ok = ok && isfinite(A[i * M + keep[j]]);
      if (ok) row[b++] = i;
    }
    if (k >= 3 && b >= 3) {
      /* sum of S over all pairs is the variance of the row sums */
      double v = 0, g = 0, st = 0, gm = 0, m2 = 0;
      double *mu = (double *)calloc((size_t)k, sizeof(double));
      if (mu) {
        for (int j = 0; j < k; j++) { for (int i = 0; i < b; i++) mu[j] += A[row[i] * M + keep[j]]; mu[j] /= b; g += mu[j]; }
        for (int j = 0; j < k; j++) { double d = mu[j] - g / k; gm += d * d; }
        for (int i = 0; i < b; i++) {
          double t = 0;
          for (int j = 0; j < k; j++) { double d = A[row[i] * M + keep[j]] - mu[j]; v += d * d; t += d; }
          st += t * t;
        }
        v /= (b - 1.0) * k; st /= b - 1.0; m2 = gm / (k - 1.0);
        double c = (st - k * v) / (k * (k - 1.0)), w = v - c;
        double de = fmax(0.0, m2 - w / b), ta = fmax(0.0, m2 - w * (1.0 + 1.0 / b));
        out[0] = k; out[1] = b; out[2] = v; out[3] = w;
        out[4] = de; out[5] = de / (de + v); out[6] = ta; out[7] = ta / (ta + v);
        free(mu);
      }
    }
  }
  free(keep); free(row);
}
