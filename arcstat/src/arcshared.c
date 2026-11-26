/* arcshared.c -- kappa curves that share one shape, and the permutation test built on them.
 *
 * Every run j of a set has its own baseline, drift, plateau, location and scale; the shape (k, h) is
 * common:  y = g0_j + m_j x + g1_j F((x - mu_j)/sg_j; k, h).  The sum of squares is minimised over all
 * of them, the linear part exactly.
 *
 * At a fixed shape the runs separate, so the profile P_j(k, h) = min over (mu_j, sg_j) is tabulated
 * once per run on a grid of shapes, each cell from its neighbours and from cold starts. The sum over a
 * set is read off the tables; its lowest local minima are the starts of a joint Levenberg-Marquardt
 * search in all the parameters, whose normal equations are block-diagonal in the runs and are solved
 * through the Schur complement of the shape. The best end points are checked run by run from cold
 * starts and polished by a simplex over the shape with every run re-solved. Many sets of the same runs
 * share the tables. For k < -0.15 a run is carried by the lower end of its support and the distance
 * scale there, (xi, log tau) = (mu + sg/k, log(-sg/k)), and by (mu, log sg) otherwise.
 *
 * Parallel over runs, then over sets; every loop writes its own slots, so the result does not depend
 * on the number of threads.
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
#include "arcef.h"
#include "arcshared.h"
#include "arcthreads.h"

#define SH_BAD 1e12
#define SH_EP (-0.15)

typedef struct { const double *x; int n; double *cx, *ry, sxx, m0, s0, x0, rg; } sh_run;   /* ry: y off (1, x) */

/* curve and its two local derivatives (D1, D2 may be NULL) */
static void sh_curve(const sh_run *r, const double *q, double k, double h, int ep, double *F, double *D1,
                     double *D2) {
  double sc = exp(q[1]);
  for (int i = 0; i < r->n; i++) {
    double t = (r->x[i] - q[0]) / sc, d;
    F[i] = ep ? arcef_k4Fe(t, k, h, D1 ? &d : NULL) : arcef_k4F(t, k, h, D1 ? &d : NULL);
    if (D1) { D1[i] = -d / sc; D2[i] = -t * d; }
  }
}
/* put a column off (1, x) in place (off 1 only without drift: cx is then NULL) */
static void sh_project(const sh_run *r, double *f) {
  double m = 0, c = 0;
  for (int i = 0; i < r->n; i++) m += f[i];
  m /= r->n;
  for (int i = 0; i < r->n; i++) f[i] -= m;
  if (!r->cx) return;
  for (int i = 0; i < r->n; i++) c += r->cx[i] * f[i];
  c /= r->sxx;
  for (int i = 0; i < r->n; i++) f[i] -= c * r->cx[i];
}
/* F projected in place; its coefficient in *g; the sum of squares */
static double sh_rss(const sh_run *r, double *F, double *g) {
  double ff = 0, lo = HUGE_VAL, hi = -HUGE_VAL;
  for (int i = 0; i < r->n; i++) {
    if (!isfinite(F[i])) return SH_BAD;
    ff += F[i] * F[i];
    if (F[i] < lo) lo = F[i];
    if (F[i] > hi) hi = F[i];
  }
  if (!(hi - lo >= 1e-8)) return SH_BAD;
  sh_project(r, F);
  double d = 0, e = 0, v = 0;
  for (int i = 0; i < r->n; i++) { d += F[i] * F[i]; e += F[i] * r->ry[i]; }
  if (!(d > 1e-14 * ff)) return SH_BAD;
  *g = e / d;
  for (int i = 0; i < r->n; i++) { double u = r->ry[i] - *g * F[i]; v += u * u; }
  return isfinite(v) ? v : SH_BAD;
}
static double sh_val(const sh_run *r, const double *q, double k, double h, int ep, double *F) {
  double g;
  if (!isfinite(q[0]) || !(fabs(q[1]) < 60.0)) return SH_BAD;
  sh_curve(r, q, k, h, ep, F, NULL, NULL);
  return sh_rss(r, F, &g);
}
/* local columns times g, projected, off F as well (F projected, dd = sum F^2) */
static void sh_column(const sh_run *r, double *C, const double *F, double g) {
  sh_project(r, C);
  double d = 0, e = 0;
  for (int i = 0; i < r->n; i++) { d += F[i] * F[i]; e += F[i] * C[i]; }
  e /= d;
  for (int i = 0; i < r->n; i++) C[i] = g * (C[i] - e * F[i]);
}
/* Levenberg-Marquardt in the two local parameters at a fixed shape; W holds 4 n */
static double sh_inner(const sh_run *r, double *q, double k, double h, int ep, int maxit, double *W) {
  int n = r->n;
  double *F = W, *A = W + n, *B = W + 2 * (size_t)n, *P = W + 3 * (size_t)n, g, lam = 1e-2;
  double f = sh_val(r, q, k, h, ep, F);
  int small = 0;
  if (f >= SH_BAD) return f;
  for (int it = 0; it < maxit; it++) {
    sh_curve(r, q, k, h, ep, F, A, B);
    if (sh_rss(r, F, &g) >= SH_BAD) break;
    sh_column(r, A, F, g); sh_column(r, B, F, g);
    double aa = 0, ab = 0, bb = 0, ga = 0, gb = 0;
    for (int i = 0; i < n; i++) {
      double u = r->ry[i] - g * F[i];
      aa += A[i] * A[i]; ab += A[i] * B[i]; bb += B[i] * B[i]; ga += A[i] * u; gb += B[i] * u;
    }
    int ok = 0;
    for (int tr = 0; tr < 10; tr++) {
      double a = aa * (1 + lam) + 1e-300, d = bb * (1 + lam) + 1e-300, det = a * d - ab * ab;
      double s0 = (d * ga - ab * gb) / det, s1 = (a * gb - ab * ga) / det;
      if (!isfinite(s0) || !isfinite(s1)) { lam *= 10.0; continue; }
      if (fabs(s1) > 2.0) { s0 *= 2.0 / fabs(s1); s1 *= 2.0 / fabs(s1); }
      double qn[2] = {q[0] + s0, q[1] + s1}, fn = sh_val(r, qn, k, h, ep, P);
      if (fn < f) {
        small = (f - fn < 1e-10 * f) ? small + 1 : 0;
        q[0] = qn[0]; q[1] = qn[1]; f = fn; ok = 1;
        lam = lam * 0.3 > 1e-12 ? lam * 0.3 : 1e-12;
        break;
      }
      lam *= 5.0;
    }
    if (!ok || small >= 2) break;
  }
  return f;
}
static void sh_convert(double *q, double k, int from_ep, int to_ep) {
  if (from_ep == to_ep) return;
  double s = exp(q[1]);
  if (to_ep) { q[0] = q[0] + s / k; q[1] = log(s / -k); }       /* (mu, sg) to (xi, tau) */
  else       { q[0] = q[0] + s;     q[1] = log(-k * s); }
}
/* the best of q (if it is a start: warm) and a set of cold starts, all solved at (k, h) */
static double sh_cold(const sh_run *r, double *q, double k, double h, int warm, int dense, double *W) {
  /* In the end-point coordinates the lower end may lie far before the data when the scale is long:
   the curve then tends to a power of a linear function of x, along a ridge whose direction is the
   ratio of the distance to the lower end and the scale, so the long scales are started at several
   such ratios. */
  int ep = k < SH_EP, ni = dense ? 6 : 3, nj = dense ? 12 : 6, nw = ep ? 2 : 1;
  /* the four best distinct end points after a few steps are carried to convergence */
  double bv[4] = {SH_BAD, SH_BAD, SH_BAD, SH_BAD}, bp[4][2] = {{0, 0}}, b[2] = {q[0], q[1]}, bf = SH_BAD;
  for (int s = warm ? -1 : 0; s < ni * nj * nw; s++) {
    double p[2];
    if (s < 0) { p[0] = q[0]; p[1] = q[1]; }
    else {
      int iw = s % nw, ij = (s / nw) % nj, ii = s / (nw * nj);
      double u = ii / (ni - 1.0), v = ij / (nj - 1.0);
      p[1] = log(r->s0) + (ep ? log(0.3) + v * log(3e5) : log(0.3) + v * log(300.0));
      p[0] = r->x0 + (r->m0 - r->x0) * (ep ? 0.30 + 0.70 * u : 0.65 + 0.80 * u);
      if (iw > 0) {
        /* a long scale: the valley in the ratio w = (x0 - xi)/tau is narrow, so it is scanned */
        if (ii > 0 || v < 0.4) continue;
        double bw = SH_BAD, tau = exp(p[1]), wb = 0.0;
        for (int i = 0; i <= 80; i++) {
          double z[2] = {r->x0 - (i / 40.0) * tau, p[1]}, g = sh_val(r, z, k, h, 1, W);
          if (g < bw) { bw = g; wb = i / 40.0; }
        }
        for (int i = -8; i <= 8; i++) {
          double w = wb + i / 400.0, z[2] = {r->x0 - w * tau, p[1]}, g = sh_val(r, z, k, h, 1, W);
          if (g < bw) { bw = g; p[0] = z[0]; }
        }
        if (!(bw < SH_BAD)) continue;
        if (p[0] == r->x0 + (r->m0 - r->x0) * 0.30) p[0] = r->x0 - wb * tau;
      }
    }
    double f = sh_inner(r, p, k, h, ep, 6, W);
    /* distinct end points only: a basin reached from many starts takes one place */
    int dup = -1;
    for (int z = 0; z < 4 && dup < 0; z++)
      if (bv[z] < SH_BAD && fabs(p[1] - bp[z][1]) < 0.25 && fabs(p[0] - bp[z][0]) < 0.1 * exp(fmin(p[1], bp[z][1]))) dup = z;
    if (dup >= 0) {
      if (!(f < bv[dup])) continue;
      for (int z = dup; z < 3; z++) { bv[z] = bv[z + 1]; bp[z][0] = bp[z + 1][0]; bp[z][1] = bp[z + 1][1]; }
      bv[3] = SH_BAD;
    }
    int at = 4;
    while (at > 0 && f < bv[at - 1]) at--;
    if (at < 4) {
      for (int z = 3; z > at; z--) { bv[z] = bv[z - 1]; bp[z][0] = bp[z - 1][0]; bp[z][1] = bp[z - 1][1]; }
      bv[at] = f; bp[at][0] = p[0]; bp[at][1] = p[1];
    }
  }
  for (int z = 0; z < 4; z++) {
    if (!(bv[z] < SH_BAD)) continue;
    double f = sh_inner(r, bp[z], k, h, ep, 200, W);
    if (f < bf) { bf = f; b[0] = bp[z][0]; b[1] = bp[z][1]; }
  }
  if (bf >= SH_BAD) return bf;
  /* a long scale may sit on a ridge that falls slowly towards the power-of-a-linear-function limit;
   it is followed by doubling the scale with the ratio w held and re-solving */
  if (ep && b[1] > log(30.0 * r->s0)) {
    for (int t = 0; t < 12; t++) {
      double tau = exp(b[1]), w = (r->x0 - b[0]) / tau, z[2] = {r->x0 - w * 2.0 * tau, b[1] + log(2.0)};
      double g = sh_inner(r, z, k, h, ep, 60, W);
      if (!(g < bf * (1 - 1e-9))) break;
      bf = g; b[0] = z[0]; b[1] = z[1];
    }
  }
  q[0] = b[0]; q[1] = b[1];
  return bf;
}

/* ---- the grid of shapes inside a box ---- */
typedef struct { int nk, nh; double k[64], h[64]; } sh_grid;
static int sh_cmp(const void *a, const void *b) {
  double x = *(const double *)a, y = *(const double *)b;
  return (x > y) - (x < y);
}
static int sh_axis(const double *v, int n, double lo, double hi, int dense, double *o) {
  int m = 0;
  for (int i = 0; i < n; i++) if (v[i] > lo && v[i] < hi) o[m++] = v[i];
  o[m++] = lo; o[m++] = hi;
  qsort(o, (size_t)m, sizeof(double), sh_cmp);
  if (dense) {
    double t[64];
    int c = 0;
    for (int i = 0; i < m; i++) { t[c++] = o[i]; if (i + 1 < m) t[c++] = 0.5 * (o[i] + o[i + 1]); }
    memcpy(o, t, (size_t)c * sizeof(double));
    m = c;
  }
  return m;
}
static void sh_mkgrid(const double *bnd, int dense, sh_grid *G) {
  static const double ks[24] = {-28, -18, -12, -8, -5.5, -4, -3, -2.4, -1.9, -1.5, -1.2, -0.97, -0.8, -0.65,
                                -0.5, -0.35, -0.2, -0.05, 0.1, 0.25, 0.4, 0.6, 0.8, 0.94};
  static const double hs[24] = {-50, -25, -12, -6, -4, -2.9, -2.2, -1.6, -1.2, -0.9, -0.6, -0.4, -0.25, -0.1,
                                0.03, 0.15, 0.28, 0.4, 0.5, 0.6, 0.75, 0.9, 1.3, 2.5};
  G->nk = sh_axis(ks, 24, bnd[0], bnd[1], dense, G->k);
  G->nh = sh_axis(hs, 24, bnd[2], bnd[3], dense, G->h);
}
/* tab: 3 per cell (value, two local parameters in the coordinates of that k) */
static void sh_table(const sh_run *r, const sh_grid *G, int dense, double *tab, double *W) {
  for (int a = 0; a < G->nk; a++) for (int b = 0; b < G->nh; b++) {
    double *c = tab + 3 * ((size_t)a * G->nh + b), q[2] = {0, 0}, p[2];
    int ep = G->k[a] < SH_EP, warm = 0;
    if (b > 0 && c[-3] < SH_BAD) { q[0] = c[-2]; q[1] = c[-1]; warm = 1; }
    double f = sh_cold(r, q, G->k[a], G->h[b], warm, dense, W);
    if (a > 0) {
      const double *u = c - 3 * (size_t)G->nh;
      if (u[0] < SH_BAD) {
        p[0] = u[1]; p[1] = u[2];
        sh_convert(p, G->k[a - 1], G->k[a - 1] < SH_EP, ep);
        double g = sh_inner(r, p, G->k[a], G->h[b], ep, 60, W);
        if (g < f) { f = g; q[0] = p[0]; q[1] = p[1]; }
      }
    }
    c[0] = f; c[1] = q[0]; c[2] = q[1];
  }
  /* back over the grid: a good basin found late spreads to the earlier cells */
  for (int a = G->nk - 1; a >= 0; a--) for (int b = G->nh - 1; b >= 0; b--) {
    double *c = tab + 3 * ((size_t)a * G->nh + b), p[2];
    int ep = G->k[a] < SH_EP;
    for (int w = 0; w < 2; w++) {
      const double *u = NULL;
      int ea = ep; double ka = G->k[a];
      if (w == 0 && b + 1 < G->nh) u = c + 3;
      if (w == 1 && a + 1 < G->nk) { u = c + 3 * (size_t)G->nh; ka = G->k[a + 1]; ea = ka < SH_EP; }
      if (!u || !(u[0] < SH_BAD)) continue;
      p[0] = u[1]; p[1] = u[2];
      if (!(ea == 0 && ep == 1 && !(G->k[a] < 0))) sh_convert(p, ep ? G->k[a] : ka, ea, ep);
      double g = sh_inner(r, p, G->k[a], G->h[b], ep, 40, W);
      if (g < c[0]) { c[0] = g; c[1] = p[0]; c[2] = p[1]; }
    }
  }
}

/* ---- one set of runs ---- */
typedef struct {
  const sh_run *run; const int *id; int m;
  double klo, khi, plo, phi;              /* box in (k, p), p = asinh(2h) */
  double *W;                              /* 8 n of the longest run */
} sh_set;
static double sh_h(double p) { return 0.5 * sinh(p); }
static double sh_total(const sh_set *S, const double *sh, const double *loc, int ep, double *each) {
  double t = 0;
  if (ep && !(sh[0] < -0.02)) return SH_BAD;
  for (int j = 0; j < S->m; j++) {
    double f = sh_val(S->run + S->id[j], loc + 2 * j, sh[0], sh_h(sh[1]), ep, S->W);
    if (f >= SH_BAD) return SH_BAD;
    if (each) each[j] = f;
    t += f;
  }
  return t;
}
/* joint search; sh = (k, p), loc 2 per run in the coordinates of ep = sh[0] < SH_EP; scr: 24 m */
static double sh_joint(const sh_set *S, double *sh, double *loc, int maxit, double *scr) {
  int m = S->m, ep = sh[0] < SH_EP, small = 0;
  double lam = 1e-2, f = sh_total(S, sh, loc, ep, NULL);
  if (f >= SH_BAD) return f;
  double *L = scr, *E = scr + 3 * m, *gl = scr + 7 * m, *st = scr + 9 * m, *ln = scr + 11 * m;
  for (int it = 0; it < maxit; it++) {
    double As[3] = {0, 0, 0}, gs[2] = {0, 0};
    int bad = 0;
    for (int j = 0; j < m && !bad; j++) {
      const sh_run *r = S->run + S->id[j];
      int n = r->n;
      double *F = S->W, *C[4], *P = S->W + 5 * (size_t)n, g;
      for (int c = 0; c < 4; c++) C[c] = S->W + (1 + c) * (size_t)n;
      const double *q = loc + 2 * j;
      for (int c = 0; c < 2; c++) {                   /* the shape columns, inside the box */
        double a[2] = {sh[0], sh[1]}, b[2] = {sh[0], sh[1]}, d = 1e-6;
        a[c] = fmin(a[c] + d, c ? S->phi : S->khi); b[c] = fmax(b[c] - d, c ? S->plo : S->klo);
        if (!(a[c] > b[c])) { for (int i = 0; i < n; i++) C[2 + c][i] = 0.0; continue; }
        sh_curve(r, q, a[0], sh_h(a[1]), ep, C[2 + c], NULL, NULL);
        sh_curve(r, q, b[0], sh_h(b[1]), ep, P, NULL, NULL);
        for (int i = 0; i < n; i++) C[2 + c][i] = (C[2 + c][i] - P[i]) / (a[c] - b[c]);
      }
      sh_curve(r, q, sh[0], sh_h(sh[1]), ep, F, C[0], C[1]);
      if (sh_rss(r, F, &g) >= SH_BAD) { bad = 1; break; }
      double M[4][4], v[4];
      for (int c = 0; c < 4; c++) sh_column(r, C[c], F, g);
      for (int c = 0; c < 4; c++) {
        v[c] = 0;
        for (int i = 0; i < n; i++) v[c] += C[c][i] * (r->ry[i] - g * F[i]);
        for (int d = 0; d <= c; d++) { M[c][d] = 0; for (int i = 0; i < n; i++) M[c][d] += C[c][i] * C[d][i]; }
      }
      L[3 * j] = M[0][0]; L[3 * j + 1] = M[1][0]; L[3 * j + 2] = M[1][1];
      E[4 * j] = M[2][0]; E[4 * j + 1] = M[2][1]; E[4 * j + 2] = M[3][0]; E[4 * j + 3] = M[3][1];
      gl[2 * j] = v[0]; gl[2 * j + 1] = v[1];
      As[0] += M[2][2]; As[1] += M[3][2]; As[2] += M[3][3]; gs[0] += v[2]; gs[1] += v[3];
    }
    if (bad) break;
    int fz0 = (sh[0] <= S->klo && gs[0] < 0) || (sh[0] >= S->khi && gs[0] > 0);
    int fz1 = (sh[1] <= S->plo && gs[1] < 0) || (sh[1] >= S->phi && gs[1] > 0);
    int ok = 0;
    for (int tr = 0; tr < 12; tr++) {
      double a = As[0] * (1 + lam) + 1e-300, b = As[1], d = As[2] * (1 + lam) + 1e-300, u = gs[0], w = gs[1];
      for (int j = 0; j < m; j++) {                    /* Schur complement */
        double p = L[3 * j] * (1 + lam) + 1e-300, q = L[3 * j + 1], r = L[3 * j + 2] * (1 + lam) + 1e-300, det = p * r - q * q;
        const double *e = E + 4 * j, *g = gl + 2 * j;
        double x0 = (r * e[0] - q * e[1]) / det, x1 = (p * e[1] - q * e[0]) / det;
        double y0 = (r * e[2] - q * e[3]) / det, y1 = (p * e[3] - q * e[2]) / det;
        a -= x0 * e[0] + x1 * e[1]; b -= x0 * e[2] + x1 * e[3]; d -= y0 * e[2] + y1 * e[3];
        u -= x0 * g[0] + x1 * g[1]; w -= y0 * g[0] + y1 * g[1];
      }
      double s0, s1;
      if (fz0 && fz1) s0 = s1 = 0;
      else if (fz0) { s0 = 0; s1 = w / d; }
      else if (fz1) { s1 = 0; s0 = u / a; }
      else { double det = a * d - b * b; s0 = (d * u - b * w) / det; s1 = (a * w - b * u) / det; }
      if (!isfinite(s0) || !isfinite(s1)) { lam *= 10; continue; }
      double lim = fmax(1.0, fmax(fabs(s0) / (0.5 + 0.5 * fabs(sh[0])), fabs(s1)));
      double sn[2] = {sh[0] + s0 / lim, sh[1] + s1 / lim};
      sn[0] = fmin(fmax(sn[0], S->klo), S->khi); sn[1] = fmin(fmax(sn[1], S->plo), S->phi);
      double t0 = sn[0] - sh[0], t1 = sn[1] - sh[1];
      for (int j = 0; j < m; j++) {
        double p = L[3 * j] * (1 + lam) + 1e-300, q = L[3 * j + 1], r = L[3 * j + 2] * (1 + lam) + 1e-300, det = p * r - q * q;
        const double *e = E + 4 * j;
        double g0 = gl[2 * j] - e[0] * t0 - e[2] * t1, g1 = gl[2 * j + 1] - e[1] * t0 - e[3] * t1;
        st[2 * j] = (r * g0 - q * g1) / det; st[2 * j + 1] = (p * g1 - q * g0) / det;
        if (fabs(st[2 * j + 1]) > 2.0) { double z = 2.0 / fabs(st[2 * j + 1]); st[2 * j] *= z; st[2 * j + 1] *= z; }
        ln[2 * j] = loc[2 * j] + st[2 * j]; ln[2 * j + 1] = loc[2 * j + 1] + st[2 * j + 1];
      }
      double fn = sh_total(S, sn, ln, ep, NULL);
      if (fn < f) {
        small = (f - fn < 1e-10 * f) ? small + 1 : 0;
        f = fn; sh[0] = sn[0]; sh[1] = sn[1];
        memcpy(loc, ln, 2 * (size_t)m * sizeof(double));
        lam = lam * 0.3 > 1e-12 ? lam * 0.3 : 1e-12;
        ok = 1;
        int en = sh[0] < SH_EP;
        if (en != ep) { for (int j = 0; j < m; j++) sh_convert(loc + 2 * j, sh[0], ep, en); ep = en; }
        break;
      }
      lam *= 5.0;
    }
    if (!ok || small >= 2) break;
  }
  return f;
}
/* the profile over the shape, every run re-solved from where it last was at its best value */
typedef struct { const sh_set *S; double *loc, *tmp; double best; int ep; } sh_prof;
static double sh_pobj(const double *p, const void *vp) {
  sh_prof *c = (sh_prof *)vp;
  const sh_set *S = c->S;
  double k = p[0], h = sh_h(p[1]), t = 0;
  if (!(k >= S->klo && k <= S->khi && p[1] >= S->plo && p[1] <= S->phi)) return SH_BAD;
  int ep = k < SH_EP;
  for (int j = 0; j < S->m; j++) {
    double q[2] = {c->loc[2 * j], c->loc[2 * j + 1]};
    sh_convert(q, k < 0 ? k : -0.1, c->ep, ep);
    double f = sh_inner(S->run + S->id[j], q, k, h, ep, 25, S->W);
    if (f >= SH_BAD) return SH_BAD;
    t += f;
    c->tmp[2 * j] = q[0]; c->tmp[2 * j + 1] = q[1];
  }
  if (t < c->best) { c->best = t; c->ep = ep; memcpy(c->loc, c->tmp, 2 * (size_t)S->m * sizeof(double)); }
  return t;
}
/* out: rss, k, h; loc: 2 per run, in the coordinates of k. tab: the table of every run */
static void sh_fit1(const sh_set *S, const sh_grid *G, double *const *tab, int dense, double *out, double *loc) {
  int m = S->m, nc = G->nk * G->nh, KC = dense ? 16 : 6;
  double *tot = (double *)malloc((size_t)nc * sizeof(double));
  double *cl = (double *)malloc((size_t)KC * 2 * m * sizeof(double)), *scr = (double *)malloc((size_t)(24 * m + 8) * sizeof(double));
  double cf[16], cs[16][2];
  int ci[16], n = 0;
  out[0] = HUGE_VAL; out[1] = out[2] = NAN;
  if (tot && cl && scr) {
    for (int i = 0; i < nc; i++) {
      tot[i] = 0;
      for (int j = 0; j < m; j++) tot[i] += tab[S->id[j]][3 * (size_t)i];
    }
    for (int i = 0; i < nc; i++) {                      /* the local minima of the table, lowest first */
      int a = i / G->nh, b = i % G->nh, mn = tot[i] < SH_BAD;
      for (int da = -1; da <= 1 && mn; da++) for (int db = -1; db <= 1; db++) {
        int a2 = a + da, b2 = b + db;
        if (a2 < 0 || a2 >= G->nk || b2 < 0 || b2 >= G->nh || (!da && !db)) continue;
        if (tot[a2 * G->nh + b2] < tot[i]) { mn = 0; break; }
      }
      if (!mn) continue;
      int p = n < KC ? n++ : KC;
      if (p == KC) { if (!(tot[i] < tot[ci[KC - 1]])) continue; p = KC - 1; }
      while (p > 0 && tot[ci[p - 1]] > tot[i]) { ci[p] = ci[p - 1]; p--; }
      ci[p] = i;
    }
    for (int c = 0; c < n; c++) {
      double *l = cl + (size_t)c * 2 * m;
      cs[c][0] = G->k[ci[c] / G->nh]; cs[c][1] = asinh(2.0 * G->h[ci[c] % G->nh]);
      for (int j = 0; j < m; j++) { const double *t = tab[S->id[j]] + 3 * (size_t)ci[c]; l[2 * j] = t[1]; l[2 * j + 1] = t[2]; }
      cf[c] = sh_joint(S, cs[c], l, 60, scr);
    }
    /* the two best: every run from cold starts at that shape, then the joint search again */
    for (int rank = 0; rank < (dense ? 6 : 2) && rank < n; rank++) {
      int b = rank;
      for (int c = rank; c < n; c++) if (cf[c] < cf[b]) b = c;
      if (b != rank) {
        double t = cf[b]; cf[b] = cf[rank]; cf[rank] = t;
        for (int j = 0; j < 2; j++) { t = cs[b][j]; cs[b][j] = cs[rank][j]; cs[rank][j] = t; }
        for (int j = 0; j < 2 * m; j++) { t = cl[(size_t)b * 2 * m + j]; cl[(size_t)b * 2 * m + j] = cl[(size_t)rank * 2 * m + j]; cl[(size_t)rank * 2 * m + j] = t; }
      }
      double *l = cl + (size_t)rank * 2 * m;
      for (int pass = 0; pass < 3; pass++) {
        double f = 0;
        for (int j = 0; j < m; j++) f += sh_cold(S->run + S->id[j], l + 2 * j, cs[rank][0], sh_h(cs[rank][1]), 1, dense, S->W);
        double g = sh_joint(S, cs[rank], l, 60, scr);
        int done = !(g < cf[rank] * (1 - 1e-7));
        cf[rank] = fmin(cf[rank], g);
        (void)f;
        if (done) break;
      }
      double *q = scr, *t = scr + 2 * m;
      memcpy(q, l, 2 * (size_t)m * sizeof(double));
      sh_prof P = {S, q, t, HUGE_VAL, cs[rank][0] < SH_EP};
      double z[2] = {cs[rank][0], cs[rank][1]}, zn[2], g;
      sh_pobj(z, &P);
      arck4_nmd(sh_pobj, &P, 2, z, dense ? 120 : 40, zn, &g);
      if (P.best < cf[rank]) {
        /* the simplex may end beside the shape of its best value only by rounding */
        double s2[2] = {zn[0], zn[1]}, *l2 = scr + 4 * m;
        memcpy(l2, P.loc, 2 * (size_t)m * sizeof(double));
        double f2 = sh_total(S, s2, l2, s2[0] < SH_EP, NULL);
        if (f2 < cf[rank]) {
          double *scr2 = (double *)malloc((size_t)(24 * m) * sizeof(double));
          if (scr2) { f2 = sh_joint(S, s2, l2, 40, scr2); free(scr2); }
          cf[rank] = f2; cs[rank][0] = s2[0]; cs[rank][1] = s2[1];
          memcpy(l, l2, 2 * (size_t)m * sizeof(double));
        }
      }
      if (cf[rank] < out[0]) {
        out[0] = cf[rank]; out[1] = cs[rank][0]; out[2] = sh_h(cs[rank][1]);
        memcpy(loc, l, 2 * (size_t)m * sizeof(double));
      }
    }
  }
  free(tot); free(cl); free(scr);
}

/* all sets; loc receives 2 per member, in order */
static int sh_core(const double *X, const double *Y, const int *off, int R, int nsub, const int *sub,
                   const int *soff, const double *bnd, int drift, int dense, double *out, double *outrun) {
  int nmax = 0, N = off[R], ok = 0;
  sh_grid G;
  sh_mkgrid(bnd, dense, &G);
  sh_run *run = (sh_run *)calloc((size_t)R, sizeof(sh_run));
  double *buf = (double *)malloc(2 * (size_t)N * sizeof(double));
  double **tab = (double **)calloc((size_t)R, sizeof(double *));
  char *used = (char *)calloc((size_t)R, 1);
  for (int i = 0; i < 3 * nsub; i++) out[i] = NAN;
  if (run && buf && tab && used) {
    ok = 1;
    for (int i = 0; i < soff[nsub]; i++) used[sub[i]] = 1;
    for (int r = 0; r < R; r++) {
      sh_run *q = run + r;
      const double *x = X + off[r], *y = Y + off[r];
      int n = q->n = off[r + 1] - off[r];
      if (n > nmax) nmax = n;
      q->x = x; q->cx = buf + off[r]; q->ry = buf + N + off[r];
      double mx = 0, lo = y[0], hi = y[0], t[3];
      for (int i = 0; i < n; i++) { mx += x[i]; if (y[i] < lo) lo = y[i]; if (y[i] > hi) hi = y[i]; }
      mx /= n;
      for (int i = 0; i < n; i++) { q->cx[i] = x[i] - mx; q->sxx += q->cx[i] * q->cx[i]; q->ry[i] = y[i]; }
      if (!drift) q->cx = NULL;
      sh_project(q, q->ry);
      for (int j = 0; j < 3; j++) { int i = 0; while (i < n - 1 && y[i] < lo + 0.25 * (j + 1) * (hi - lo)) i++; t[j] = x[i]; }
      q->x0 = x[0]; q->rg = x[n - 1] - x[0]; q->m0 = t[1]; q->s0 = fmax((t[2] - t[0]) / 1.5, 0.02 * q->rg);
      if (used[r]) { tab[r] = (double *)malloc(3 * (size_t)G.nk * G.nh * sizeof(double)); if (!tab[r]) ok = 0; }
    }
  }
  if (ok) {
#ifdef _OPENMP
#pragma omp parallel num_threads(arc_threads())
#endif
    {
      double *W = (double *)malloc(8 * (size_t)nmax * sizeof(double));
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 1)
#endif
      for (int r = 0; r < R; r++) if (W && tab[r]) sh_table(run + r, &G, dense, tab[r], W);
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 1)
#endif
      for (int s = 0; s < nsub; s++) {
        int m = soff[s + 1] - soff[s];
        double *loc = (double *)malloc(2 * (size_t)m * sizeof(double));
        if (W && loc) {
          sh_set S = {run, sub + soff[s], m, bnd[0], bnd[1], asinh(2 * bnd[2]), asinh(2 * bnd[3]), W};
          double *o = out + 3 * (size_t)s;
          sh_fit1(&S, &G, tab, dense, o, loc);
          if (outrun && o[0] < SH_BAD) for (int j = 0; j < m; j++) {
            double *q = outrun + 6 * (size_t)(soff[s] + j), g = NAN;   /* g0, m, g1, mu, sg, rss */
            const sh_run *r = run + S.id[j];
            const double *y = Y + off[S.id[j]];
            int ep = o[1] < SH_EP;
            sh_curve(r, loc + 2 * j, o[1], o[2], ep, W, NULL, NULL);
            memcpy(W + r->n, W, (size_t)r->n * sizeof(double));
            q[5] = sh_rss(r, W, &g);
            double mF = 0, my = 0, mx = 0, xF = 0, xy = 0;
            for (int i = 0; i < r->n; i++) { mF += W[r->n + i]; my += y[i]; mx += r->x[i]; }
            mF /= r->n; my /= r->n; mx /= r->n;
            if (r->cx) for (int i = 0; i < r->n; i++) { xF += r->cx[i] * W[r->n + i]; xy += r->cx[i] * y[i]; }
            q[2] = g; q[1] = r->cx ? (xy - g * xF) / r->sxx : 0.0; q[0] = my - q[1] * mx - g * mF;
            double z[2] = {loc[2 * j], loc[2 * j + 1]};
            sh_convert(z, o[1], ep, 0);
            q[3] = z[0]; q[4] = exp(z[1]);
          }
        }
        free(loc);
      }
      free(W);
    }
  }
  if (tab) for (int r = 0; r < R; r++) free(tab[r]);
  free(tab); free(run); free(buf); free(used);
  return ok;
}
void arcsh_fit(const double *X, const double *Y, const int *off, const int *nr, const int *nsub,
               const int *sub, const int *soff, const double *bounds, const int *drift,
               const int *dense, double *out, double *outrun) {
  for (int i = 0; i < 6 * soff[*nsub]; i++) outrun[i] = NAN;
  sh_core(X, Y, off, *nr, *nsub, sub, soff, bounds, *drift, *dense, out, outrun);
}

/* ---- the permutation test ---- */
/* all label assignments with the given group sizes, as the sum of the costs of their groups */
typedef struct { int n, G; const int *size; const double *cost; double obs; const double *sl; int ns; double *cnt, N, mn; } sh_alloc;
static void sh_rec(sh_alloc *A, int g, unsigned left, double sum) {
  if (g == A->G) {
    A->N += 1.0;
    if (sum < A->mn) A->mn = sum;
    for (int s = 0; s < A->ns; s++) if (sum <= A->obs * (1.0 + A->sl[s])) A->cnt[s] += 1.0;
    return;
  }
  /* every subset of left with size[g] members, in increasing order of mask */
  int id[32], m = 0, c[32], k = A->size[g];
  for (int i = 0; i < A->n; i++) if (left >> i & 1u) id[m++] = i;
  for (int i = 0; i < k; i++) c[i] = i;
  for (;;) {
    unsigned u = 0;
    for (int i = 0; i < k; i++) u |= 1u << id[c[i]];
    sh_rec(A, g + 1, left & ~u, sum + A->cost[u]);
    int i = k - 1;
    while (i >= 0 && c[i] == m - k + i) i--;
    if (i < 0) break;
    c[i]++;
    for (int j = i + 1; j < k; j++) c[j] = c[j - 1] + 1;
  }
}
/* Problem q consists of the runs poff[q] .. poff[q+1]-1 with group labels grp (0, 1, ...). RSS_TR is
 * the sum of squares of one shape for all of them, RSS_W that of a shape per group; every set of runs
 * of a group's size is fitted once, and all label assignments are formed from the fits. p[s] is the
 * share of assignments, the observed one included, with RSS_W <= observed (1 + slack[s]): slack is an
 * allowance for the numerical error of the fits and counts near-ties against rejection. out: 5 per
 * problem (RSS_TR, RSS_W, Lambda_T = RSS_TR/RSS_W, number of assignments, smallest RSS_W); p: nslack
 * per problem, NaN when there are more than maxalloc assignments or more than 14 runs. */
void arcsh_perm(const double *X, const double *Y, const int *off, const int *nprob, const int *poff,
                const int *grp, const double *bounds, const int *drift, const int *dense,
                const double *slack, const int *nslack, const int *maxalloc, double *out, double *p) {
  int Q = *nprob, R = poff[Q], nsub = 0, tot = 0;
  for (int i = 0; i < 5 * Q; i++) out[i] = NAN;
  for (int i = 0; i < Q * *nslack; i++) p[i] = NAN;
  int *full = (int *)calloc((size_t)Q, sizeof(int));
  int (*size)[32] = (int (*)[32])calloc((size_t)Q, sizeof(*size));
  int *G = (int *)calloc((size_t)Q, sizeof(int)), *first = (int *)calloc((size_t)Q + 1, sizeof(int));
  if (!full || !size || !G || !first) { free(full); free(size); free(G); free(first); return; }
  for (int pass = 0; pass < 2; pass++) {
    int *sub = NULL, *soff = NULL;
    double *o = NULL;
    if (pass == 1) {
      sub = (int *)malloc((size_t)tot * sizeof(int)); soff = (int *)malloc((size_t)(nsub + 1) * sizeof(int));
      o = (double *)malloc(3 * (size_t)nsub * sizeof(double));
      if (!sub || !soff || !o) { free(sub); free(soff); free(o); break; }
      soff[0] = 0;
    }
    int s = 0, t = 0;
    for (int q = 0; q < Q; q++) {
      int n = poff[q + 1] - poff[q];
      if (pass == 0) {
        G[q] = 0;
        for (int i = 0; i < n; i++) { int g = grp[poff[q] + i]; if (g >= 0 && g < 32) { size[q][g]++; if (g + 1 > G[q]) G[q] = g + 1; } }
        double N = 1, left = n;
        for (int g = 0; g < G[q]; g++) for (int i = 1; i <= size[q][g]; i++) N *= left-- / i;
        full[q] = n <= 14 && N <= *maxalloc;
        first[q] = s;
      }
      /* set 0 of a problem is all its runs; then the observed groups, or every mask of a group's size */
      if (pass == 1) { for (int i = 0; i < n; i++) sub[t + i] = poff[q] + i; soff[s + 1] = t + n; }
      t += n; s++;
      if (!full[q]) {
        for (int g = 0; g < G[q]; g++) {
          if (pass == 1) { int c = 0; for (int i = 0; i < n; i++) if (grp[poff[q] + i] == g) sub[t + c++] = poff[q] + i; soff[s + 1] = t + size[q][g]; }
          t += size[q][g]; s++;
        }
      } else for (unsigned u = 1; u < (1u << n); u++) {
        int c = 0, want = 0;
        for (int i = 0; i < n; i++) c += u >> i & 1u;
        for (int g = 0; g < G[q]; g++) want |= size[q][g] == c;
        if (!want || c == n) continue;
        if (pass == 1) { int d = 0; for (int i = 0; i < n; i++) if (u >> i & 1u) sub[t + d++] = poff[q] + i; soff[s + 1] = t + c; }
        t += c; s++;
      }
    }
    if (pass == 0) { nsub = s; tot = t; first[Q] = s; continue; }
    if (sh_core(X, Y, off, R, nsub, sub, soff, bounds, *drift, *dense, o, NULL)) {
      for (int q = 0; q < Q; q++) {
        int n = poff[q + 1] - poff[q], a = first[q];
        double w = 0;
        if (!full[q]) for (int g = 0; g < G[q]; g++) w += o[3 * (a + 1 + g)];
        else {
          double *cost = (double *)malloc(((size_t)1 << n) * sizeof(double)), *cnt = (double *)calloc((size_t)*nslack, sizeof(double));
          if (cost && cnt) {
            int c = a + 1;
            for (unsigned u = 1; u < (1u << n); u++) {
              int b = 0, want = 0;
              for (int i = 0; i < n; i++) b += u >> i & 1u;
              for (int g = 0; g < G[q]; g++) want |= size[q][g] == b;
              cost[u] = NAN;
              if (b == n && want) cost[u] = o[3 * a];
              else if (want) cost[u] = o[3 * c++];
            }
            for (int g = 0; g < G[q]; g++) { unsigned u = 0; for (int i = 0; i < n; i++) if (grp[poff[q] + i] == g) u |= 1u << i; w += u ? cost[u] : 0.0; }
            sh_alloc A = {n, G[q], size[q], cost, w, slack, *nslack, cnt, 0.0, HUGE_VAL};
            sh_rec(&A, 0, (1u << n) - 1u, 0.0);
            out[5 * q + 3] = A.N; out[5 * q + 4] = A.mn;
            for (int i = 0; i < *nslack; i++) p[q * *nslack + i] = cnt[i] / A.N;
          }
          free(cost); free(cnt);
        }
        out[5 * q] = o[3 * a]; out[5 * q + 1] = w; out[5 * q + 2] = o[3 * a] / w;
      }
    }
    free(sub); free(soff); free(o);
  }
  free(full); free(size); free(G); free(first);
}
