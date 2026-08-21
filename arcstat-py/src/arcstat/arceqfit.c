/* arceqfit.c -- every fitting and simulation loop of the equivalence-system paper, in the
 * back end: the beta-companion curve fits (free and equivalence-constrained), their residual
 * bootstraps, the kappa locus-constrained fit and bootstrap, the exact-score manifold test
 * simulation of the paper's Table 4, the three-estimator study of its Table 1, and two thin
 * grid sweeps (the kappa wall scan and the E-in-beta sections). The R and Python fronts are
 * orchestration only.
 *
 * Randomness never touches a host RNG: every replicate draws from its own splitmix64 stream
 * seeded by (seed, replicate index), with normals by Box-Muller from those uniforms, so all
 * outputs are byte-identical across thread counts and between the fronts.
 */

/* Floating-point contraction is pinned OFF so that the R and Python fronts, which compile these
   sources with different flags, cannot differ in whether multiply-add pairs are fused. Without
   this the two fronts agree to about fifteen digits and disagree in the last bit, which is enough
   to fail the parity harnesses -- and it surfaces unpredictably, because whether a fusion happens
   depends on register allocation, so an unrelated edit can expose it. */
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
#include "arceq.h"
#include "arceq2.h"

#ifdef _OPENMP
static int ef_threads(void) {
  int t = omp_get_max_threads();
  return t > 2 ? 2 : (t < 1 ? 1 : t);
}
#endif

/* ---- deterministic RNG: splitmix64 uniforms, Box-Muller normals ---- */
static double ef_unif(unsigned long long *s) {
  unsigned long long z = (*s += 0x9E3779B97F4A7C15ULL);
  z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
  z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
  z ^= (z >> 31);
  return (double)(z >> 11) * (1.0 / 9007199254740992.0);
}
static double ef_norm(unsigned long long *s) {
  double u1 = ef_unif(s), u2 = ef_unif(s);
  if (u1 < 1e-300) u1 = 1e-300;
  return sqrt(-2.0 * log(u1)) * cos(6.283185307179586476925286766559 * u2);
}
static unsigned long long ef_seed(int seed, int b) {
  /* The starting state must be SCRAMBLED, not offset: splitmix64 advances its state by the
   * golden-ratio constant per draw, so states offset by multiples of it are the same sequence
   * shifted -- replicate streams built that way overlap almost entirely, and an earlier
   * version of this helper did exactly that. Two rounds of the output scrambler place the
   * 2^64 starting states effectively at random, and overlap has birthday-bound probability. */
  unsigned long long z = (unsigned long long)seed * 0x9E3779B97F4A7C15ULL
                         ^ ((unsigned long long)(b + 1) * 0xBF58476D1CE4E5B9ULL);
  z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
  z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
  z ^= (z >> 31);
  return z;
}

static double ef_logis(double p) { return 1.0 / (1.0 + exp(-p)); }
static double ef_lbeta(double a, double b) { return lgamma(a) + lgamma(b) - lgamma(a + b); }

/* linear interpolation on a table; rule 1 returns NAN outside, rule 2 clamps */
static double ef_interp(const double *xs, const double *ys, int n, double x, int rule) {
  if (n < 2 || !isfinite(x)) return NAN;
  if (x <= xs[0])   return (rule == 2) ? ys[0]   : (x == xs[0] ? ys[0] : NAN);
  if (x >= xs[n-1]) return (rule == 2) ? ys[n-1] : (x == xs[n-1] ? ys[n-1] : NAN);
  int lo = 0, hi = n - 1;
  while (hi - lo > 1) { int mid = (lo + hi) / 2; if (xs[mid] <= x) lo = mid; else hi = mid; }
  double w = (x - xs[lo]) / (xs[lo+1] - xs[lo]);
  return ys[lo] + w * (ys[lo+1] - ys[lo]);
}

/* moving-block machinery shared by the residual bootstraps: instrument-log residuals are
 * strongly serially dependent (lag-one autocorrelation up to 0.996 on the application's fits),
 * and iid resampling would understate every bootstrap spread; the block length is read from
 * the residual autocorrelation, as everywhere else in this programme's software. */
static int ef_blocklen(const double *r, int n) {
  double rb = 0.0;
  for (int i = 0; i < n; i++) rb += r[i];
  rb /= n;
  double den = 0.0;
  for (int i = 0; i < n; i++) den += (r[i] - rb) * (r[i] - rb);
  int L = n - 2 < 200 ? n - 2 : 200, bl = 0;
  double band = 2.0 / sqrt((double)n);
  for (int l = 1; l <= L; l++) {
    double sacf = 0.0;
    for (int i = 0; i + l < n; i++) sacf += (r[i] - rb) * (r[i + l] - rb);
    if (fabs(sacf / den) < band) { bl = l; break; }
  }
  if (!bl) bl = (int)(pow((double)n, 1.0 / 3.0) + 0.5);
  return bl < 2 ? 2 : bl;
}
static void ef_blockdraw(const double *m0, const double *r0, int n, int bl,
                         unsigned long long *s, double *yb) {
  int st_ok = n - bl + 1, pos = 0;
  while (pos < n) {
    int i0 = (int)(ef_unif(s) * st_ok);
    if (i0 >= st_ok) i0 = st_ok - 1;
    for (int i = 0; i < bl && pos < n; i++, pos++) yb[pos] = m0[pos] + r0[i0 + i];
  }
}

/* ---- the beta-companion mean curve m(x) = g0 + g1 * F(z; alpha, beta) ----
 * F inverts the scaled incomplete beta. Point-wise numerical inversion (arceq2_bc_cdf) costs
 * thirty-odd forward evaluations per point and dominated every fit; the curve is instead
 * inverted through one forward quantile grid per evaluation -- the quantile is strictly
 * increasing, so monotone interpolation of the inverse is exact to grid resolution, far below
 * the fitting noise. */
#define EF_NG 2049
/* The inverse-grid depends ONLY on the two shape exponents, never on the location, the scale or
 * the affine layer. Building it is the dominant cost of every mean-curve evaluation (2049
 * incomplete-beta inversions), so it is split out here and reused: a score evaluation perturbs
 * mu and sigma for four of its six Jacobian columns, and those four share the base grid with
 * the curve itself. Five builds per evaluation instead of nine, with the grid values, and
 * therefore every downstream number, bit for bit unchanged. */
static void ef_grid(double al, double be, double *ug, double *Qg, double *Bab) {
  int ng = EF_NG;
  for (int i = 0; i < ng; i++) ug[i] = 1e-9 + (1.0 - 2e-9) * (double)i / (double)(ng - 1);
  arceq2_bc_q(&ng, ug, &al, &be, Qg);
  *Bab = exp(ef_lbeta(al + 1.0, be + 1.0));
}

/* mean curve at a precomputed shape grid */
static void ef_mbc_g(const double *x, int n, const double *th, double *out,
                     const double *ug, const double *Qg, double Bab) {
  double g0 = th[0], g1 = th[1], mu = th[2], sr = th[3];
  int ng = EF_NG;
  for (int i = 0; i < n; i++) {
    double z = (x[i] - mu) / sr;
    if (z < 1e-9) z = 1e-9;
    if (z > 1.0 - 1e-9) z = 1.0 - 1e-9;
    double t = z * Bab, u;
    if (t <= Qg[0]) u = ug[0];
    else if (t >= Qg[ng - 1]) u = ug[ng - 1];
    else {
      int lo = 0, hi = ng - 1;
      while (hi - lo > 1) { int mid = (lo + hi) / 2; if (Qg[mid] <= t) lo = mid; else hi = mid; }
      double w = (t - Qg[lo]) / (Qg[lo + 1] - Qg[lo]);
      u = ug[lo] + w * (ug[lo + 1] - ug[lo]);
    }
    out[i] = g0 + g1 * u;
  }
}

static void ef_mbc(const double *x, int n, const double *th, double *out, double *zbuf) {
  (void)zbuf;
  double g0 = th[0], g1 = th[1], mu = th[2], sr = th[3], al = th[4], be = th[5];
  double Bab = exp(ef_lbeta(al + 1.0, be + 1.0));
  double ug[EF_NG], Qg[EF_NG];
  int ng = EF_NG;
  for (int i = 0; i < ng; i++) ug[i] = 1e-9 + (1.0 - 2e-9) * (double)i / (double)(ng - 1);
  arceq2_bc_q(&ng, ug, &al, &be, Qg);
  for (int i = 0; i < n; i++) {
    double z = (x[i] - mu) / sr;
    if (z < 1e-9) z = 1e-9;
    if (z > 1.0 - 1e-9) z = 1.0 - 1e-9;
    double t = z * Bab, u;
    if (t <= Qg[0]) u = ug[0];
    else if (t >= Qg[ng - 1]) u = ug[ng - 1];
    else {
      int lo = 0, hi = ng - 1;
      while (hi - lo > 1) { int mid = (lo + hi) / 2; if (Qg[mid] <= t) lo = mid; else hi = mid; }
      double w = (t - Qg[lo]) / (Qg[lo + 1] - Qg[lo]);
      u = ug[lo] + w * (ug[lo + 1] - ug[lo]);
    }
    out[i] = g0 + g1 * u;
  }
}

/* ---- beta-companion SSE objective, free (d=6) or constrained to a curve table (d=5) ---- */
typedef struct {
  const double *x, *y; int n;
  int eqcon;                       /* 1: beta from the curve table */
  const double *cal, *cbe; int nc; /* equivalence curve table */
  double alo, ahi;                 /* constrained alpha window */
  double *work, *zbuf;             /* n-length scratch */
} bc_ctx;

static void bc_theta(const double *p, const bc_ctx *c, double *th) {
  double al = c->eqcon ? c->alo + ef_logis(p[4]) * (c->ahi - c->alo)
                       : -1.0 + ef_logis(p[4]) * 0.98;
  double be = c->eqcon ? ef_interp(c->cal, c->cbe, c->nc, al, 1)
                       : -1.0 + ef_logis(p[5]) * 1.48;
  th[0] = p[0]; th[1] = exp(p[1]); th[2] = p[2]; th[3] = exp(p[3]); th[4] = al; th[5] = be;
}
static double bc_obj(const double *p, const void *vc) {
  const bc_ctx *c = (const bc_ctx *)vc;
  double th[6]; bc_theta(p, c, th);
  if (!isfinite(th[5])) return 1e12;
  ef_mbc(c->x, c->n, th, c->work, c->zbuf);
  double s = 0.0;
  for (int i = 0; i < c->n; i++) { double r = c->y[i] - c->work[i]; s += r * r; }
  return isfinite(s) ? s : 1e12;
}

/* readings of a fitted companion curve, in observation units: out = (a, b) */
static void bc_ip(const double *th, double *out) {
  out[0] = out[1] = NAN;
  if (!isfinite(th[4]) || !isfinite(th[5])) return;
  int nt = 0, ng = 20001; double r[5];
  arceq_readings(&th[4], &th[5], &nt, NULL, &ng, r);
  double Bab = exp(ef_lbeta(th[4] + 1.0, th[5] + 1.0));
  if (isfinite(r[1])) out[0] = th[2] + th[3] * r[1] / Bab;
  if (isfinite(r[2])) out[1] = th[2] + th[3] * r[2] / Bab;
}

/* entry: fit from ns starts (row-major, np columns each; np = 6 free, 5 constrained),
 * two Nelder-Mead passes per start, best kept. out = (theta[6], sse, a, b). */
void arceqfit_bc(const double *x, const double *y, const int *n,
                 const double *starts, const int *ns, const int *np, const int *maxit,
                 const double *cal, const double *cbe, const int *nc,
                 const double *alo, const double *ahi, double *out) {
  int N = *n, S = *ns, d = *np;
  double *work = (double *)malloc((size_t)N * sizeof(double));
  double *zbuf = (double *)malloc((size_t)N * sizeof(double));
  for (int i = 0; i < 9; i++) out[i] = NAN;
  if (!work || !zbuf) { free(work); free(zbuf); return; }
  bc_ctx c = { x, y, N, d == 5, cal, cbe, *nc, *alo, *ahi, work, zbuf };
  double bestp[6], bestf = HUGE_VAL;
  for (int s = 0; s < S; s++) {
    double p1[6], f1, p2[6], f2;
    arck4_nmd(bc_obj, &c, d, starts + (size_t)s * d, *maxit, p1, &f1);
    arck4_nmd(bc_obj, &c, d, p1, *maxit, p2, &f2);
    double f = (f2 < f1) ? f2 : f1;
    const double *p = (f2 < f1) ? p2 : p1;
    if (f < bestf) { bestf = f; memcpy(bestp, p, (size_t)d * sizeof(double)); }
  }
  if (bestf < 1e12) {
    double th[6]; bc_theta(bestp, &c, th);
    memcpy(out, th, 6 * sizeof(double));
    out[6] = bestf;
    bc_ip(th, out + 7);
  }
  free(work); free(zbuf);
}

/* entry: iid residual bootstrap of the CONSTRAINED companion fit, warm-started, OpenMP over
 * replicates; out_a receives the B tangent readings. */
void arceqfit_bc_boot(const double *x, const double *y, const int *n, const double *th,
                      const int *B, const int *seed, const int *maxit,
                      const double *cal, const double *cbe, const int *nc,
                      const double *alo, const double *ahi, double *out_a) {
  int N = *n, NB = *B;
  double *m0 = (double *)malloc((size_t)N * sizeof(double));
  double *z0 = (double *)malloc((size_t)N * sizeof(double));
  double *r0 = (double *)malloc((size_t)N * sizeof(double));
  for (int b = 0; b < NB; b++) out_a[b] = NAN;
  if (!m0 || !z0 || !r0) { free(m0); free(z0); free(r0); return; }
  ef_mbc(x, N, th, m0, z0);
  for (int i = 0; i < N; i++) r0[i] = y[i] - m0[i];
  int bl = ef_blocklen(r0, N);
  double w0 = (th[4] - *alo) / (*ahi - *alo);
  if (w0 < 1e-6) w0 = 1e-6;
  if (w0 > 1.0 - 1e-6) w0 = 1.0 - 1e-6;
  double p0[5] = { th[0], log(th[1]), th[2], log(th[3]), log(w0 / (1.0 - w0)) };
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(ef_threads())
#endif
  for (int b = 0; b < NB; b++) {
    double *yb = (double *)malloc((size_t)N * sizeof(double));
    double *wk = (double *)malloc((size_t)N * sizeof(double));
    double *zb = (double *)malloc((size_t)N * sizeof(double));
    if (yb && wk && zb) {
      unsigned long long s = ef_seed(*seed, b);
      ef_blockdraw(m0, r0, N, bl, &s, yb);
      bc_ctx c = { x, yb, N, 1, cal, cbe, *nc, *alo, *ahi, wk, zb };
      double p1[5], f1, p2[5], f2;
      arck4_nmd(bc_obj, &c, 5, p0, *maxit, p1, &f1);
      arck4_nmd(bc_obj, &c, 5, p1, *maxit, p2, &f2);
      double thb[6], ip[2];
      bc_theta(f2 < f1 ? p2 : p1, &c, thb);
      bc_ip(thb, ip);
      out_a[b] = ip[0];
    }
    free(yb); free(wk); free(zb);
  }
  free(m0); free(z0); free(r0);
}

/* ---- kappa fit constrained to the equivalence locus table h -> k ---- */
typedef struct {
  const double *x, *y; int n;
  const double *lh, *lk; int nl;   /* locus table */
  double hlo, hhi;
} kq_ctx;
static void kq_theta(const double *p, const kq_ctx *c, double *th) {
  double h = c->hlo + ef_logis(p[4]) * (c->hhi - c->hlo);
  th[0] = p[0]; th[1] = exp(p[1]); th[2] = p[2]; th[3] = exp(p[3]);
  th[4] = ef_interp(c->lh, c->lk, c->nl, h, 2); th[5] = h;
}
static double kq_obj(const double *p, const void *vc) {
  const kq_ctx *c = (const kq_ctx *)vc;
  double th[6]; kq_theta(p, c, th);
  if (!isfinite(th[4])) return 1e12;
  double s = 0.0;
  for (int i = 0; i < c->n; i++) {
    double m = th[0] + th[1] * k4_F(c->x[i], th[2], th[3], th[4], th[5]);
    double r = c->y[i] - m; s += r * r;
  }
  return isfinite(s) ? s : 1e12;
}
void arceqfit_k4eq(const double *x, const double *y, const int *n,
                   const double *starts, const int *ns, const int *maxit,
                   const double *lh, const double *lk, const int *nl,
                   const double *hlo, const double *hhi, double *out) {
  kq_ctx c = { x, y, *n, lh, lk, *nl, *hlo, *hhi };
  double bestp[5], bestf = HUGE_VAL;
  for (int i = 0; i < 8; i++) out[i] = NAN;
  for (int s = 0; s < *ns; s++) {
    double p1[5], f1, p2[5], f2;
    arck4_nmd(kq_obj, &c, 5, starts + (size_t)s * 5, *maxit, p1, &f1);
    arck4_nmd(kq_obj, &c, 5, p1, *maxit, p2, &f2);
    double f = (f2 < f1) ? f2 : f1;
    const double *p = (f2 < f1) ? p2 : p1;
    if (f < bestf) { bestf = f; memcpy(bestp, p, 5 * sizeof(double)); }
  }
  if (bestf < 1e12) {
    double th[6]; kq_theta(bestp, &c, th);
    memcpy(out, th, 6 * sizeof(double));
    out[6] = bestf;
    int g = 4000; double r[3];
    arck4_readings(th, &g, r);
    out[7] = r[0];
  }
}
void arceqfit_k4eq_boot(const double *x, const double *y, const int *n, const double *th,
                        const double *p0in, const int *B, const int *seed, const int *maxit,
                        const double *lh, const double *lk, const int *nl,
                        const double *hlo, const double *hhi, double *out_a) {
  int N = *n, NB = *B;
  double *m0 = (double *)malloc((size_t)N * sizeof(double));
  double *r0 = (double *)malloc((size_t)N * sizeof(double));
  for (int b = 0; b < NB; b++) out_a[b] = NAN;
  if (!m0 || !r0) { free(m0); free(r0); return; }
  for (int i = 0; i < N; i++) {
    m0[i] = th[0] + th[1] * k4_F(x[i], th[2], th[3], th[4], th[5]);
    r0[i] = y[i] - m0[i];
  }
  int bl = ef_blocklen(r0, N);
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(ef_threads())
#endif
  for (int b = 0; b < NB; b++) {
    double *yb = (double *)malloc((size_t)N * sizeof(double));
    if (yb) {
      unsigned long long s = ef_seed(*seed, b);
      ef_blockdraw(m0, r0, N, bl, &s, yb);
      kq_ctx c = { x, yb, N, lh, lk, *nl, *hlo, *hhi };
      double p1[5], f1, p2[5], f2;
      arck4_nmd(kq_obj, &c, 5, p0in, *maxit, p1, &f1);
      arck4_nmd(kq_obj, &c, 5, p1, *maxit, p2, &f2);
      double thb[6]; kq_theta(f2 < f1 ? p2 : p1, &c, thb);
      int g = 4000; double r[3];
      arck4_readings(thb, &g, r);
      out_a[b] = r[0];
      free(yb);
    }
  }
  free(m0); free(r0);
}

/* ---- the exact score statistic of the manifold test: F of the residual regression on the
 * six numerically differentiated Jacobian columns of the companion mean curve ---- */
static double ef_scoreF6(const double *x, const double *y, int n, const double *th,
                         double *wk /* 9n scratch */) {
  double *e0 = wk, *zb = wk + n, *mp = wk + 2 * n, *mm = wk + 3 * n, *J = wk + 4 * n; /* J: 6n? */
  /* layout: e0, zb, mp, mm each n; J needs 6n -> total 10n supplied by caller */
  ef_mbc(x, n, th, e0, zb);
  for (int i = 0; i < n; i++) e0[i] = y[i] - e0[i];
  for (int j = 0; j < 6; j++) {
    double tp[6], tm[6];
    memcpy(tp, th, sizeof(tp)); memcpy(tm, th, sizeof(tm));
    double h = 1e-5 * (fabs(th[j]) > 1.0 ? fabs(th[j]) : 1.0);
    tp[j] += h; tm[j] -= h;
    ef_mbc(x, n, tp, mp, zb);
    ef_mbc(x, n, tm, mm, zb);
    for (int i = 0; i < n; i++) J[(size_t)j * n + i] = (mp[i] - mm[i]) / (2.0 * h);
  }
  /* normal equations with a relative ridge; e1 = e0 - J beta */
  double G[36], v[6];
  for (int a = 0; a < 6; a++) {
    v[a] = 0.0;
    for (int i = 0; i < n; i++) v[a] += J[(size_t)a * n + i] * e0[i];
    for (int b2 = a; b2 < 6; b2++) {
      double s = 0.0;
      for (int i = 0; i < n; i++) s += J[(size_t)a * n + i] * J[(size_t)b2 * n + i];
      G[a * 6 + b2] = G[b2 * 6 + a] = s;
    }
  }
  /* Degeneracy guard: at a parameter point whose mean curve is (near-)flat over the design
   * the Jacobian columns vanish, the ridge keeps the solve finite, and the statistic returns
   * ~0 -- an escape hatch a profiling minimiser will find and live in. The statistic is only
   * meaningful where the model is locally identifiable, so refuse rank-deficient Jacobians:
   * the intercept column always has ||.||^2 = n, hence the relative floor on the smallest
   * Gram diagonal detects any vanished column. */
  double dmx = 0.0, dmn = 0.0;
  for (int a = 0; a < 6; a++) {
    double d = G[a * 6 + a];
    if (a == 0 || d > dmx) dmx = d;
    if (a == 0 || d < dmn) dmn = d;
  }
  if (!(dmx > 0.0) || dmn < 1e-10 * dmx) return NAN;
  double tr = 0.0; for (int a = 0; a < 6; a++) tr += G[a * 6 + a];
  for (int a = 0; a < 6; a++) G[a * 6 + a] += 1e-12 * (tr > 0 ? tr : 1.0);
  /* Cholesky solve */
  double L[36]; memset(L, 0, sizeof(L));
  for (int a = 0; a < 6; a++) {
    for (int b2 = 0; b2 <= a; b2++) {
      double s = G[a * 6 + b2];
      for (int k2 = 0; k2 < b2; k2++) s -= L[a * 6 + k2] * L[b2 * 6 + k2];
      if (a == b2) { if (s <= 0) return NAN; L[a * 6 + a] = sqrt(s); }
      else L[a * 6 + b2] = s / L[b2 * 6 + b2];
    }
  }
  double w1[6], beta[6];
  for (int a = 0; a < 6; a++) {
    double s = v[a];
    for (int k2 = 0; k2 < a; k2++) s -= L[a * 6 + k2] * w1[k2];
    w1[a] = s / L[a * 6 + a];
  }
  for (int a = 5; a >= 0; a--) {
    double s = w1[a];
    for (int k2 = a + 1; k2 < 6; k2++) s -= L[k2 * 6 + a] * beta[k2];
    beta[a] = s / L[a * 6 + a];
  }
  double sse0 = 0.0, sse1 = 0.0;
  for (int i = 0; i < n; i++) {
    double fit = 0.0;
    for (int j = 0; j < 6; j++) fit += J[(size_t)j * n + i] * beta[j];
    double e1 = e0[i] - fit;
    sse0 += e0[i] * e0[i]; sse1 += e1 * e1;
  }
  if (sse1 <= 0) return NAN;
  double num = sse0 - sse1;
  if (num < 0.0) num = 0.0;                /* roundoff can push the difference epsilon-negative */
  return (num / 6.0) / (sse1 / (double)(n - 6));
}

/* profiled minimum of the score statistic over the manifold interpolant.
 * The nuisance location and scale are profiled inside compact windows around the design --
 * mu in [xlo - R/2, xlo + R], sigma in [R/5, 3R] for design range R -- so the search cannot
 * drift to configurations whose support misses the data, where the score statistic is
 * degenerate. The containment argument of the level guarantee needs only that the truth lie
 * inside the searched set, which these windows ensure by a wide margin. */
typedef struct {
  const double *x, *y; int n;
  const double *cal, *cbe; int nc;
  double alo, ahi;                 /* alpha window of the profiling transform */
  double xlo, xrng;                /* design left edge and range, for the nuisance windows */
  double *wk;
} mf_ctx;
/* T at (mu, sg, al, be) with the affine layer profiled EXACTLY: the intercept and curve
 * columns [1, F] lie in the Jacobian's span, so the six-column projection residual (sse1) is
 * invariant to (gamma0, gamma1) and the statistic's minimum over them is attained at the
 * least-squares fit of y on [1, F]. The remaining Jacobian columns are used as dF/d(mu, sg,
 * al, be) -- the true columns times gamma1, a scaling that leaves their span unchanged. */
static double ef_scoreF6p(const double *x, const double *y, int n,
                          double mu, double sg, double al, double be, double *wk) {
  /* Workspace layout: [e0 | reserved | Fp | Fm | J], each n long. The second slot is not used
     by this routine; the pointers below carry explicit offsets, so it stays reserved rather
     than being reclaimed, and the caller's allocation is unchanged. */
  double *e0 = wk, *Fp = wk + 2 * n, *Fm = wk + 3 * n, *J = wk + 4 * n;
  double thF[6] = { 0.0, 1.0, mu, sg, al, be };
  double ug[EF_NG], Qg[EF_NG], Bab;
  ef_grid(al, be, ug, Qg, &Bab);            /* the base shape grid, reused four times below */
  /* Identifiability floor, independent of any caller-side window: the six columns must be
   * informed by more than a handful of observations. Design points outside the fitted
   * support are clamped to the endpoints by ef_mbc and contribute nothing to any column,
   * so a support that has slid off the design leaves the statistic supported on the few
   * interior points -- and a statistic supported on fewer points than it has columns is
   * degenerate, whatever its Gram conditioning looks like. */
  int nin = 0;
  for (int i = 0; i < n; i++) {
    double z = (x[i] - mu) / sg;
    if (z > 1e-9 && z < 1.0 - 1e-9) nin++;
  }
  if (nin < 12 || nin < n / 10) return NAN;
  ef_mbc_g(x, n, thF, J + n, ug, Qg, Bab);            /* column 1: F itself */
  for (int i = 0; i < n; i++) J[i] = 1.0;             /* column 0: intercept */
  /* least squares of y on [1, F] */
  double s1 = 0, sF = 0, sFF = 0, sy = 0, sFy = 0;
  for (int i = 0; i < n; i++) {
    double F = J[n + i];
    s1 += 1.0; sF += F; sFF += F * F; sy += y[i]; sFy += F * y[i];
  }
  double det = s1 * sFF - sF * sF;
  if (!(fabs(det) > 1e-12 * (s1 * sFF + 1.0))) return NAN;   /* F ~ constant: degenerate */
  double g1 = (s1 * sFy - sF * sy) / det, g0 = (sy - g1 * sF) / s1;
  double sse0 = 0.0;
  for (int i = 0; i < n; i++) {
    e0[i] = y[i] - g0 - g1 * J[n + i];
    sse0 += e0[i] * e0[i];
  }
  /* central-difference columns dF/d(mu, sg, al, be) */
  static const int idx[4] = { 2, 3, 4, 5 };
  for (int j = 0; j < 4; j++) {
    double tp[6], tm[6];
    memcpy(tp, thF, sizeof(tp)); memcpy(tm, thF, sizeof(tm));
    double v = thF[idx[j]];
    double h = 1e-5 * (fabs(v) > 1.0 ? fabs(v) : 1.0);
    tp[idx[j]] += h; tm[idx[j]] -= h;
    if (idx[j] < 4) {                       /* mu, sigma: shape unchanged, so reuse the grid */
      ef_mbc_g(x, n, tp, Fp, ug, Qg, Bab);
      ef_mbc_g(x, n, tm, Fm, ug, Qg, Bab);
    } else {                                /* alpha, beta: the grid itself moves */
      double ugp[EF_NG], Qgp[EF_NG], Babp;
      ef_grid(tp[4], tp[5], ugp, Qgp, &Babp);
      ef_mbc_g(x, n, tp, Fp, ugp, Qgp, Babp);
      ef_grid(tm[4], tm[5], ugp, Qgp, &Babp);
      ef_mbc_g(x, n, tm, Fm, ugp, Qgp, Babp);
    }
    for (int i = 0; i < n; i++) J[(size_t)(j + 2) * n + i] = (Fp[i] - Fm[i]) / (2.0 * h);
  }
  /* ridge-stabilised projection of e0 on all six columns, with the degeneracy guard of
   * ef_scoreF6 */
  double G[36], v6[6];
  for (int a = 0; a < 6; a++) {
    v6[a] = 0.0;
    for (int i = 0; i < n; i++) v6[a] += J[(size_t)a * n + i] * e0[i];
    for (int b2 = a; b2 < 6; b2++) {
      double s = 0.0;
      for (int i = 0; i < n; i++) s += J[(size_t)a * n + i] * J[(size_t)b2 * n + i];
      G[a * 6 + b2] = G[b2 * 6 + a] = s;
    }
  }
  double dmx = 0.0, dmn = 0.0;
  for (int a = 0; a < 6; a++) {
    double d = G[a * 6 + a];
    if (a == 0 || d > dmx) dmx = d;
    if (a == 0 || d < dmn) dmn = d;
  }
  if (!(dmx > 0.0) || dmn < 1e-10 * dmx) return NAN;
  double tr = 0.0; for (int a = 0; a < 6; a++) tr += G[a * 6 + a];
  for (int a = 0; a < 6; a++) G[a * 6 + a] += 1e-12 * (tr > 0 ? tr : 1.0);
  double L[36]; memset(L, 0, sizeof(L));
  for (int a = 0; a < 6; a++) {
    for (int b2 = 0; b2 <= a; b2++) {
      double s = G[a * 6 + b2];
      for (int k2 = 0; k2 < b2; k2++) s -= L[a * 6 + k2] * L[b2 * 6 + k2];
      if (a == b2) { if (s <= 0) return NAN; L[a * 6 + a] = sqrt(s); }
      else L[a * 6 + b2] = s / L[b2 * 6 + b2];
    }
  }
  double w1[6], beta[6];
  for (int a = 0; a < 6; a++) {
    double s = v6[a];
    for (int k2 = 0; k2 < a; k2++) s -= L[a * 6 + k2] * w1[k2];
    w1[a] = s / L[a * 6 + a];
  }
  for (int a = 5; a >= 0; a--) {
    double s = w1[a];
    for (int k2 = a + 1; k2 < 6; k2++) s -= L[k2 * 6 + a] * beta[k2];
    beta[a] = s / L[a * 6 + a];
  }
  double sse1 = 0.0;
  for (int i = 0; i < n; i++) {
    double fit = 0.0;
    for (int j = 0; j < 6; j++) fit += J[(size_t)j * n + i] * beta[j];
    double e1 = e0[i] - fit;
    sse1 += e1 * e1;
  }
  if (sse1 <= 0) return NAN;
  double num = sse0 - sse1;
  if (num < 0.0) num = 0.0;
  return (num / 6.0) / (sse1 / (double)(n - 6));
}

/* Support-overlap requirement on the searched set. The score statistic is only meaningful
 * where the model is locally identifiable, and identifiability needs the fitted support to
 * carry the design: a member whose support slides off the design's edge until one point
 * remains interior supports all six Jacobian columns on that one observation, passes any
 * Gram-conditioning guard, and returns a zero statistic -- a degenerate hole INSIDE the
 * nuisance windows, which an earlier version of this routine searched freely. The set is
 * therefore defined by the overlap itself: the fitted support [mu, mu+sg] must cover the
 * design's central half, which every parameter near the truth does with room to spare, so
 * the containment argument behind the level guarantee is untouched. */
static int mf_overlap(const mf_ctx *c, double mu, double sg) {
  double q1 = c->xlo + 0.25 * c->xrng, q3 = c->xlo + 0.75 * c->xrng;
  return (mu <= q1) && (mu + sg >= q3);
}

/* profiling objective in the three remaining coordinates: windowed mu and sigma, and the
 * manifold coordinate alpha */
static double mf_obj3(const double *p, const void *vc) {
  const mf_ctx *c = (const mf_ctx *)vc;
  double al = c->alo + ef_logis(p[2]) * (c->ahi - c->alo);
  double be = ef_interp(c->cal, c->cbe, c->nc, al, 1);
  if (!isfinite(be)) return 1e10;
  /* transforms chosen so the whole reachable box is compatible with the overlap condition
   * and the default start sits comfortably inside it: mu in [xlo-R/2, xlo+R/4] (the upper
   * end is the overlap's own bound), sg in [R/2, 3R] */
  double mu = c->xlo - 0.5 * c->xrng + ef_logis(p[0]) * 0.75 * c->xrng;
  double sg = 0.5 * c->xrng * exp(ef_logis(p[1]) * 1.7917594692280550);  /* log(6): up to 3R */
  if (!mf_overlap(c, mu, sg)) return 1e10;
  double v = ef_scoreF6p(c->x, c->y, c->n, mu, sg, al, be, c->wk);
  return isfinite(v) ? v : 1e10;
}

/* Least-squares objective on the same three coordinates and the same set: the residual sum of
 * squares after the affine layer is profiled out. The score statistic's surface is razor sharp
 * at small noise -- a direct descent on it steps over the minimum and returns a value above the
 * statistic at the truth, which the audit then reports as an optimiser failure (it did, at 13
 * and 97 per cent, before this stage existed). The least-squares surface is smooth and its
 * minimiser is close to the score minimiser, so it is used to locate the neighbourhood first. */
static double mf_sse3(const double *p, const void *vc) {
  const mf_ctx *c = (const mf_ctx *)vc;
  double al = c->alo + ef_logis(p[2]) * (c->ahi - c->alo);
  double be = ef_interp(c->cal, c->cbe, c->nc, al, 1);
  if (!isfinite(be)) return 1e30;
  double mu = c->xlo - 0.5 * c->xrng + ef_logis(p[0]) * 0.75 * c->xrng;
  double sg = 0.5 * c->xrng * exp(ef_logis(p[1]) * 1.7917594692280550);
  if (!mf_overlap(c, mu, sg)) return 1e30;
  double thF[6] = { 0.0, 1.0, mu, sg, al, be };
  double *F0 = c->wk + 4 * c->n;
  double ug[EF_NG], Qg[EF_NG], Bab;
  ef_grid(al, be, ug, Qg, &Bab);
  ef_mbc_g(c->x, c->n, thF, F0, ug, Qg, Bab);
  double s1 = (double)c->n, sF = 0, sFF = 0, sy = 0, sFy = 0;
  for (int i = 0; i < c->n; i++) {
    double F = F0[i];
    sF += F; sFF += F * F; sy += c->y[i]; sFy += F * c->y[i];
  }
  double det = s1 * sFF - sF * sF;
  if (!(fabs(det) > 1e-12 * (s1 * sFF + 1.0))) return 1e30;
  double g1 = (s1 * sFy - sF * sy) / det, g0 = (sy - g1 * sF) / s1;
  double sse = 0.0;
  for (int i = 0; i < c->n; i++) { double e = c->y[i] - g0 - g1 * F0[i]; sse += e * e; }
  return sse;
}

/* the profiled minimum: locate by least squares along a multistart in the manifold coordinate
 * (the hard direction), then descend the score statistic from each located neighbourhood.
 * st[0..2] give the (mu, sigma, alpha-centre) starting transforms. */
static double ef_profmin(const mf_ctx *c, const double *st, int maxit) {
  double pas[5] = { st[2] - 2.2, st[2] - 1.1, st[2], st[2] + 1.1, st[2] + 2.2 };
  double fbest = 1e30;
  int it1 = maxit / 2 > 300 ? maxit / 2 : 300;
  for (int a0 = 0; a0 < 5; a0++) {
    double s0[3] = { st[0], st[1], pas[a0] }, pls[3], fls, p1[3], f1, p2[3], f2;
    /* stage 1: least squares locates the neighbourhood */
    arck4_nmd(mf_sse3, c, 3, s0, it1, pls, &fls);
    if (!isfinite(fls) || fls >= 1e29) continue;
    /* stage 2: the score statistic itself, from the located point and from the raw start */
    arck4_nmd(mf_obj3, c, 3, pls, maxit, p1, &f1);
    if (f1 < fbest) fbest = f1;
    arck4_nmd(mf_obj3, c, 3, p1, maxit, p2, &f2);
    if (f2 < fbest) fbest = f2;
  }
  return fbest;
}

/* entry: the whole Table-4 simulation. Four arms: truth, beta+0.20, beta-0.15, beta-0.30.
 * out: Tlev[R], Taud[R], Tp2[R], Tp15[R], Tp4[R], then the three displaced arms' REFERENCE
 * statistics Tref2[R], Tref15[R], Tref4[R] (8R values).
 *
 * The reference statistics are the power certificate. The audit in the truth arm bounds the
 * profiled minimum by the statistic at the truth; off the manifold there is no truth to use,
 * but the nearest manifold member to the arm's noiseless curve -- the mimicry fit, supplied by
 * the caller in thref as (mu, sigma, alpha, beta) per arm -- is a point of the searched set, so
 * the statistic there is an upper bound the profiled minimum must not exceed. An optimiser that
 * fails off the manifold returns too large a minimum and inflates power, so this is the
 * certificate power claims need and level claims already have. */
void arceqfit_msim(const int *nM, const double *sde, const int *Rm, const int *seed,
                   const double *th0, const double *st5, const int *maxit,
                   const double *cal, const double *cbe, const int *nc,
                   const double *alo, const double *ahi, const double *thref, double *out) {
  int N = *nM, R = *Rm;
  double shifts[4] = { 0.0, 0.20, -0.15, -0.30 };
  double *x = (double *)malloc((size_t)N * sizeof(double));
  if (!x) return;
  for (int i = 0; i < N; i++) x[i] = th0[2] + th0[3] * (double)i / (double)(N - 1);
  for (long i = 0; i < 8L * R; i++) out[i] = NAN;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(ef_threads())
#endif
  for (int t = 0; t < 4 * R; t++) {
    int arm = t / R, r = t % R;
    double *y  = (double *)malloc((size_t)N * sizeof(double));
    double *m  = (double *)malloc((size_t)N * sizeof(double));
    double *zb = (double *)malloc((size_t)N * sizeof(double));
    double *wk = (double *)malloc((size_t)N * 10 * sizeof(double));
    if (y && m && zb && wk) {
      double th2[6]; memcpy(th2, th0, sizeof(th2));
      th2[5] += shifts[arm];
      ef_mbc(x, N, th2, m, zb);
      unsigned long long s = ef_seed(*seed + 1000 * arm, r);
      for (int i = 0; i < N; i++) y[i] = m[i] + (*sde) * ef_norm(&s);
      mf_ctx c = { x, y, N, cal, cbe, *nc, *alo, *ahi, x[0], x[N - 1] - x[0], wk };
      double T = ef_profmin(&c, st5, *maxit);
      if (arm == 0) {
        out[r] = T;
        out[(size_t)R + r] = ef_scoreF6(x, y, N, th0, wk);   /* audit certificate */
      } else {
        out[(size_t)(arm + 1) * R + r] = T;
        /* reference statistic at this arm's supplied manifold member */
        const double *tr = thref + 4 * (arm - 1);
        out[(size_t)(arm + 4) * R + r] = ef_scoreF6p(x, y, N, tr[0], tr[1], tr[2], tr[3], wk);
      }
    }
    free(y); free(m); free(zb); free(wk);
  }
  free(x);
}

/* entry: the profiled score statistic at one parameter point, exposed so the paper can
 * EXERCISE the identifiability floor rather than merely assert it -- a guard that is never
 * triggered by the simulation is a guard that has not been tested. Returns NaN exactly where
 * the statistic is degenerate. */
void arceqfit_score_at(const double *x, const double *y, const int *n, const double *mu,
                       const double *sg, const double *al, const double *be, double *out) {
  double *wk = (double *)malloc((size_t)(*n) * 10 * sizeof(double));
  if (!wk) { *out = NAN; return; }
  *out = ef_scoreF6p(x, y, *n, *mu, *sg, *al, *be, wk);
  free(wk);
}

/* entry: the moving-block length rule applied to a residual vector, exposed so the text can
 * quote the block length the bootstraps actually used. */
void arceqfit_blocklen(const double *r, const int *n, int *out) {
  *out = ef_blocklen(r, *n);
}

/* entry: null draws of the profiled statistic at one manifold point -- the truth arm of the
 * simulation above, at a caller-chosen member, for the no-sharper-calibration demonstration. */
void arceqfit_nullT(const int *nM, const double *sde, const int *R2, const int *seed,
                    const double *th0, const double *st5, const int *maxit,
                    const double *cal, const double *cbe, const int *nc,
                    const double *alo, const double *ahi, double *out) {
  int N = *nM, R = *R2;
  double *x = (double *)malloc((size_t)N * sizeof(double));
  if (!x) return;
  for (int i = 0; i < N; i++) x[i] = th0[2] + th0[3] * (double)i / (double)(N - 1);
  for (int r = 0; r < R; r++) out[r] = NAN;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(ef_threads())
#endif
  for (int r = 0; r < R; r++) {
    double *y  = (double *)malloc((size_t)N * sizeof(double));
    double *m  = (double *)malloc((size_t)N * sizeof(double));
    double *zb = (double *)malloc((size_t)N * sizeof(double));
    double *wk = (double *)malloc((size_t)N * 10 * sizeof(double));
    if (y && m && zb && wk) {
      ef_mbc(x, N, th0, m, zb);
      unsigned long long s = ef_seed(*seed, r);
      for (int i = 0; i < N; i++) y[i] = m[i] + (*sde) * ef_norm(&s);
      mf_ctx c = { x, y, N, cal, cbe, *nc, *alo, *ahi, x[0], x[N - 1] - x[0], wk };
      out[r] = ef_profmin(&c, st5, *maxit);
    }
    free(y); free(m); free(zb); free(wk);
  }
  free(x);
}

/* ---- the three-estimator study of Table 1 ---- */
/* regularised incomplete beta and its inverse through the arceq2 kernel machinery */
static double ef_pbeta(double u, double a, double b) {
  double al = a - 1.0, be = b - 1.0, out;
  int n1 = 1;
  arceq2_bc_q(&n1, &u, &al, &be, &out);            /* unscaled B(u; a, b) */
  return out / exp(ef_lbeta(a, b));
}
static double ef_qbeta(double xreg, double a, double b) {
  double al = a - 1.0, be = b - 1.0, out;
  double xs = xreg * exp(ef_lbeta(a, b));
  int n1 = 1;
  arceq2_bc_cdf(&n1, &xs, &al, &be, &out);         /* inverse of the unscaled B */
  return out;
}
typedef struct { const double *x; int n; int which; const double *targ; const double *g200; } es_ctx;
static double es_obj(const double *p, const void *vc) {
  const es_ctx *c = (const es_ctx *)vc;
  double al = p[0], be = p[1];
  if (al <= -0.99 || al >= -0.01 || be <= -0.99 || be >= 0.5) return 1e9;
  double a = al + 1.0, b = be + 1.0;
  if (c->which == 0) {                     /* negative log-likelihood */
    double s1 = 0.0, s2 = 0.0;
    for (int i = 0; i < c->n; i++) {
      double w = ef_qbeta(c->x[i], a, b);
      if (w < 1e-12) w = 1e-12;
      if (w > 1.0 - 1e-12) w = 1.0 - 1e-12;
      s1 += log(w); s2 += log(1.0 - w);
    }
    double v = -((double)c->n * ef_lbeta(a, b) - al * s1 - be * s2);
    return isfinite(v) ? v : 1e9;
  }
  if (c->which == 1) {                     /* method of moments */
    double m1 = 1.0 - a / (a + b);
    double em2 = 0.0;
    for (int i = 0; i < 200; i++) em2 += ef_qbeta(c->g200[i], a, b) * c->g200[i];
    em2 = 2.0 * em2 / 200.0;
    double m2c = (1.0 - em2) - m1 * m1;
    double d1 = m1 - c->targ[0], d2 = m2c - c->targ[1];
    return d1 * d1 + d2 * d2;
  }
  /* L-moment ratios, closed form */
  double m1 = a / (a + b);
  double m2 = m1 * (a + 1.0) / (a + b + 1.0);
  double m3 = m2 * (a + 2.0) / (a + b + 2.0);
  double m4 = m3 * (a + 3.0) / (a + b + 3.0);
  double l2 = m1 - m2, l3 = 3.0 * m2 - 2.0 * m3 - m1, l4 = m1 - 6.0 * m2 + 10.0 * m3 - 5.0 * m4;
  double d1 = l3 / l2 - c->targ[0], d2 = l4 / l2 - c->targ[1];
  double v = d1 * d1 + d2 * d2;
  return isfinite(v) ? v : 1e9;
}
static int ef_cmp(const void *a, const void *b) {
  double d = *(const double *)a - *(const double *)b;
  return d < 0 ? -1 : (d > 0 ? 1 : 0);
}
/* entry: R replicates x |ns| sample sizes x 3 estimators x 2 parameters.
 * out is R x 3 x 2 x nn, index r + R*(e + 3*(j + 2*ni)). */
void arceqfit_estsim(const double *a0, const double *b0, const int *R, const int *nsizes,
                     const int *nn, const int *seed, const int *maxit, double *out) {
  int Rr = *R, NN = *nn;
  int maxn = 0;
  for (int i = 0; i < NN; i++) if (nsizes[i] > maxn) maxn = nsizes[i];
  double g200[200];
  for (int i = 0; i < 200; i++) g200[i] = ((double)i + 0.5) / 200.0;
  long tot = (long)Rr * 3 * 2 * NN;
  for (long i = 0; i < tot; i++) out[i] = NAN;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(ef_threads())
#endif
  for (int t = 0; t < Rr * NN; t++) {
    int ni = t / Rr, r = t % Rr, n = nsizes[ni];
    double *xs = (double *)malloc((size_t)n * sizeof(double));
    if (xs) {
      unsigned long long s = ef_seed(*seed + 7000 * ni, r);
      for (int i = 0; i < n; i++) xs[i] = ef_pbeta(ef_unif(&s), *a0, *b0);
      double st[2] = { -0.5, -0.3 };
      /* MLE */
      es_ctx c0 = { xs, n, 0, NULL, g200 };
      double p1[2], f1;
      arck4_nmd(es_obj, &c0, 2, st, *maxit, p1, &f1);
      if (f1 < 1e9) {
        out[r + (size_t)Rr * (0 + 3 * (0 + 2 * ni))] = p1[0];
        out[r + (size_t)Rr * (0 + 3 * (1 + 2 * ni))] = p1[1];
      }
      /* method of moments */
      double mo[2]; double mx = 0.0;
      for (int i = 0; i < n; i++) mx += xs[i];
      mx /= n;
      double vv = 0.0;
      for (int i = 0; i < n; i++) vv += (xs[i] - mx) * (xs[i] - mx);
      mo[0] = mx; mo[1] = vv / (n - 1);
      es_ctx c1 = { xs, n, 1, mo, g200 };
      arck4_nmd(es_obj, &c1, 2, st, *maxit, p1, &f1);
      if (f1 < 1e9) {
        out[r + (size_t)Rr * (1 + 3 * (0 + 2 * ni))] = p1[0];
        out[r + (size_t)Rr * (1 + 3 * (1 + 2 * ni))] = p1[1];
      }
      /* sample L-moment ratios */
      qsort(xs, (size_t)n, sizeof(double), ef_cmp);
      double b0s = 0.0, b1s = 0.0, b2s = 0.0, b3s = 0.0;
      for (int i = 0; i < n; i++) {
        double j = (double)(i + 1);
        b0s += xs[i];
        b1s += (j - 1.0) / (n - 1.0) * xs[i];
        b2s += (j - 1.0) * (j - 2.0) / ((n - 1.0) * (n - 2.0)) * xs[i];
        b3s += (j - 1.0) * (j - 2.0) * (j - 3.0) / ((n - 1.0) * (n - 2.0) * (n - 3.0)) * xs[i];
      }
      b0s /= n; b1s /= n; b2s /= n; b3s /= n;
      double l2s = 2.0 * b1s - b0s;
      double tg[2] = { (6.0 * b2s - 6.0 * b1s + b0s) / l2s,
                       (20.0 * b3s - 30.0 * b2s + 12.0 * b1s - b0s) / l2s };
      es_ctx c2 = { xs, n, 2, tg, g200 };
      arck4_nmd(es_obj, &c2, 2, st, *maxit, p1, &f1);
      if (f1 < 1e9) {
        out[r + (size_t)Rr * (2 + 3 * (0 + 2 * ni))] = p1[0];
        out[r + (size_t)Rr * (2 + 3 * (1 + 2 * ni))] = p1[1];
      }
      free(xs);
    }
  }
}

/* closed-form (tau3, tau4) of the kernel member at (alpha, beta) */
static void es_tau(double al, double be, double *out) {
  double a = al + 1.0, b = be + 1.0;
  double m1 = a / (a + b);
  double m2 = m1 * (a + 1.0) / (a + b + 1.0);
  double m3 = m2 * (a + 2.0) / (a + b + 2.0);
  double m4 = m3 * (a + 3.0) / (a + b + 3.0);
  double l2 = m1 - m2;
  out[0] = (3.0 * m2 - 2.0 * m3 - m1) / l2;
  out[1] = (m1 - 6.0 * m2 + 10.0 * m3 - 5.0 * m4) / l2;
}

/* entry: sample L-moment ratio pairs (t3, t4) from R independent companion samples of size n,
 * for the delta-method sandwich check, together with the exact-root inversion of each
 * draw through the rational ratio map; per-replicate splitmix64 streams. out is R x 4:
 * (t3, t4, alpha-hat, beta-hat). */
void arceqfit_taus(const double *a0, const double *b0, const int *R, const int *n,
                   const int *seed, double *out) {
  int Rr = *R, N = *n;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(ef_threads())
#endif
  for (int r = 0; r < Rr; r++) {
    double *xs = (double *)malloc((size_t)N * sizeof(double));
    if (xs) {
      unsigned long long s = ef_seed(*seed, r);
      for (int i = 0; i < N; i++) xs[i] = ef_pbeta(ef_unif(&s), *a0, *b0);
      qsort(xs, (size_t)N, sizeof(double), ef_cmp);
      double b0s = 0, b1s = 0, b2s = 0, b3s = 0;
      for (int i = 0; i < N; i++) {
        double j = (double)(i + 1);
        b0s += xs[i];
        b1s += (j - 1.0) / (N - 1.0) * xs[i];
        b2s += (j - 1.0) * (j - 2.0) / ((N - 1.0) * (N - 2.0)) * xs[i];
        b3s += (j - 1.0) * (j - 2.0) * (j - 3.0) / ((N - 1.0) * (N - 2.0) * (N - 3.0)) * xs[i];
      }
      b0s /= N; b1s /= N; b2s /= N; b3s /= N;
      double l2s = 2.0 * b1s - b0s;
      double t3 = (6.0 * b2s - 6.0 * b1s + b0s) / l2s;
      double t4 = (20.0 * b3s - 30.0 * b2s + 12.0 * b1s - b0s) / l2s;
      out[r] = t3;
      out[(size_t)Rr + r] = t4;
      /* exact-root inversion of the rational ratio map by damped Newton from the truth:
       * the estimator the delta-method proposition actually describes */
      double al = *a0 - 1.0, be = *b0 - 1.0;
      for (int it = 0; it < 12; it++) {
        double f0[2], fp[2], fm[2], Jl[4];
        es_tau(al, be, f0);
        double hh = 1e-7;
        es_tau(al + hh, be, fp); es_tau(al - hh, be, fm);
        Jl[0] = (fp[0] - fm[0]) / (2 * hh); Jl[2] = (fp[1] - fm[1]) / (2 * hh);
        es_tau(al, be + hh, fp); es_tau(al, be - hh, fm);
        Jl[1] = (fp[0] - fm[0]) / (2 * hh); Jl[3] = (fp[1] - fm[1]) / (2 * hh);
        double det = Jl[0] * Jl[3] - Jl[1] * Jl[2];
        if (fabs(det) < 1e-14) break;
        double e1 = f0[0] - t3, e2 = f0[1] - t4;
        double da = (Jl[3] * e1 - Jl[1] * e2) / det;
        double db = (-Jl[2] * e1 + Jl[0] * e2) / det;
        al -= da; be -= db;
        if (fabs(da) + fabs(db) < 1e-12) break;
      }
      out[2 * (size_t)Rr + r] = al;
      out[3 * (size_t)Rr + r] = be;
      free(xs);
    } else { out[r] = NAN; out[(size_t)Rr + r] = NAN;
             out[2 * (size_t)Rr + r] = NAN; out[3 * (size_t)Rr + r] = NAN; }
  }
}

/* ---- thin sweeps ---- */
/* derivative reading b of the standardised kappa member (0,1,0,1,k,h) over a k grid */
void arceqfit_bsweep(const double *ks, const int *nk, const double *h, double *out_b) {
  for (int i = 0; i < *nk; i++) {
    double th[6] = { 0.0, 1.0, 0.0, 1.0, ks[i], *h };
    int g = 4000; double r[3];
    arck4_readings(th, &g, r);
    out_b[i] = r[1];
  }
}
/* both readings of the standardised kappa member over a k grid at fixed h */
void arceqfit_absweep(const double *ks, const int *nk, const double *h,
                      double *out_a, double *out_b) {
  for (int i = 0; i < *nk; i++) {
    double th[6] = { 0.0, 1.0, 0.0, 1.0, ks[i], *h };
    int g = 4000; double r[3];
    arck4_readings(th, &g, r);
    out_a[i] = r[0]; out_b[i] = r[1];
  }
}
/* closed-form discrepancy E over a beta grid at fixed alpha */
void arceqfit_Esweep(const double *alpha, const double *bg, const int *nb, double *out) {
  for (int i = 0; i < *nb; i++) arceq2_E(alpha, &bg[i], &out[i]);
}
