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

#include "arcthreads.h"
#include "arcef.h"

/* ---- deterministic RNG: splitmix64 uniforms, Box-Muller normals ---- */
double arcef_unif(unsigned long long *s) {
  unsigned long long z = (*s += 0x9E3779B97F4A7C15ULL);
  z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
  z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
  z ^= (z >> 31);
  return (double)(z >> 11) * (1.0 / 9007199254740992.0);
}
static double ef_norm(unsigned long long *s) {
  double u1 = arcef_unif(s), u2 = arcef_unif(s);
  if (u1 < 1e-300) u1 = 1e-300;
  return sqrt(-2.0 * log(u1)) * cos(6.283185307179586476925286766559 * u2);
}
unsigned long long arcef_seed(int seed, int b) {
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
double arcef_interp(const double *xs, const double *ys, int n, double x, int rule) {
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
int arcef_blocklen(const double *r, int n) {
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
void arcef_blockdraw(const double *m0, const double *r0, int n, int bl,
                         unsigned long long *s, double *yb) {
  int st_ok = n - bl + 1, pos = 0;
  while (pos < n) {
    int i0 = (int)(arcef_unif(s) * st_ok);
    if (i0 >= st_ok) i0 = st_ok - 1;
    for (int i = 0; i < bl && pos < n; i++, pos++) yb[pos] = m0[pos] + r0[i0 + i];
  }
}

/* ---- the beta-companion mean curve m(x) = g0 + g1 * F(z; alpha, beta) ----
 * F inverts the scaled incomplete beta. Each point is solved to machine precision by a cubic
 * iteration on the logit scale that starts from its neighbour's solution (two to three incomplete
 * betas per point on sorted data). Earlier versions interpolated the inverse linearly on 2049
 * equally spaced quantile levels, which was up to 2e-4 wrong in F on the first cells and, worse,
 * made every numerical derivative of the curve a step function. */
#define EF_ZMIN 1e-9
static void ef_mbc(const double *x, int n, const double *th, double *out, double *zbuf) {
  (void)zbuf;
  double g0 = th[0], g1 = th[1], mu = th[2], sr = th[3];
  ib_par P; arceq2_ib_init(&P, th[4] + 1.0, th[5] + 1.0);
  double s = NAN;
  for (int i = 0; i < n; i++) {
    double z = (x[i] - mu) / sr;
    if (z < EF_ZMIN) z = EF_ZMIN;
    if (z > 1.0 - EF_ZMIN) z = 1.0 - EF_ZMIN;
    s = arceq2_ib_at_frac(&P, z, s);
    out[i] = g0 + g1 / (1.0 + exp(-s));
    if (!isfinite(out[i])) { out[i] = NAN; s = NAN; }
  }
}

/* F and its four derivatives dF/d(mu, sigma, alpha, beta) at the design; columns of length n in D
 * (F first). At a point outside the support, where the curve is clamped, the derivatives are zero;
 * the clamp is at zmin and 1 - zmin, and zmin = 0 means F exactly 0 and 1 there.
 * With I(u; a, b) the regularised incomplete beta, F solves I(F) = z, so dF/dz = f = 1/I'(F) and
 * dF/da = -f dI/da at fixed u; the parameter derivative is a central difference of the smaller
 * tail of I, whose error is about 1e-10. Returns the number of points inside the support. */
int arcef_jac(const double *x, int n, double mu, double sg, double al, double be, double zmin,
              double *D) {
  double a = al + 1.0, b = be + 1.0, h = 1e-5;
  if (h > 0.2 * a) h = 0.2 * a;
  if (h > 0.2 * b) h = 0.2 * b;
  ib_par P, Pap, Pam, Pbp, Pbm;
  arceq2_ib_init(&P, a, b);
  arceq2_ib_init(&Pap, a + h, b); arceq2_ib_init(&Pam, a - h, b);
  arceq2_ib_init(&Pbp, a, b + h); arceq2_ib_init(&Pbm, a, b - h);
  int nin = 0;
  double s = NAN;
  for (int i = 0; i < n; i++) {
    double z = (x[i] - mu) / sg;
    int inside = (z > zmin && z < 1.0 - zmin);
    if (zmin == 0.0 && !inside) {
      D[i] = (z <= 0.0) ? 0.0 : (z >= 1.0 ? 1.0 : NAN);
      D[n + i] = D[2 * (size_t)n + i] = D[3 * (size_t)n + i] = D[4 * (size_t)n + i] = 0.0;
      continue;
    }
    if (z < zmin) z = zmin;
    if (z > 1.0 - zmin) z = 1.0 - zmin;
    if (z < 1e-300) z = 1e-300;
    s = arceq2_ib_at_frac(&P, z, s);
    double lu, lv;
    if (s > 0.0) { lu = -log1p(exp(-s)); lv = lu - s; } else { lv = -log1p(exp(s)); lu = lv + s; }
    D[i] = exp(lu);
    if (!isfinite(D[i])) { D[i] = NAN; s = NAN; }
    D[n + i] = D[2 * (size_t)n + i] = D[3 * (size_t)n + i] = D[4 * (size_t)n + i] = 0.0;
    if (!inside) continue;
    nin++;
    double f = P.T * exp(-al * lu - be * lv), dIa, dIb;
    if (z <= 0.5) {
      dIa = arceq2_ib_lower(&Pap, s) / Pap.T - arceq2_ib_lower(&Pam, s) / Pam.T;
      dIb = arceq2_ib_lower(&Pbp, s) / Pbp.T - arceq2_ib_lower(&Pbm, s) / Pbm.T;
    } else {
      dIa = arceq2_ib_upper(&Pam, s) / Pam.T - arceq2_ib_upper(&Pap, s) / Pap.T;
      dIb = arceq2_ib_upper(&Pbm, s) / Pbm.T - arceq2_ib_upper(&Pbp, s) / Pbp.T;
    }
    D[n + i] = -f / sg;
    D[2 * (size_t)n + i] = -z * f / sg;
    D[3 * (size_t)n + i] = -f * dIa / (2.0 * h);
    D[4 * (size_t)n + i] = -f * dIb / (2.0 * h);
  }
  return nin;
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
  double be = c->eqcon ? arcef_interp(c->cal, c->cbe, c->nc, al, 1)
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
void arcef_bc_ip(const double *th, double *out) {
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
    arcef_bc_ip(th, out + 7);
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
  int bl = arcef_blocklen(r0, N);
  double w0 = (th[4] - *alo) / (*ahi - *alo);
  if (w0 < 1e-6) w0 = 1e-6;
  if (w0 > 1.0 - 1e-6) w0 = 1.0 - 1e-6;
  double p0[5] = { th[0], log(th[1]), th[2], log(th[3]), log(w0 / (1.0 - w0)) };
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(arc_threads())
#endif
  for (int b = 0; b < NB; b++) {
    double *yb = (double *)malloc((size_t)N * sizeof(double));
    double *wk = (double *)malloc((size_t)N * sizeof(double));
    double *zb = (double *)malloc((size_t)N * sizeof(double));
    if (yb && wk && zb) {
      unsigned long long s = arcef_seed(*seed, b);
      arcef_blockdraw(m0, r0, N, bl, &s, yb);
      bc_ctx c = { x, yb, N, 1, cal, cbe, *nc, *alo, *ahi, wk, zb };
      double p1[5], f1, p2[5], f2;
      arck4_nmd(bc_obj, &c, 5, p0, *maxit, p1, &f1);
      arck4_nmd(bc_obj, &c, 5, p1, *maxit, p2, &f2);
      double thb[6], ip[2];
      bc_theta(f2 < f1 ? p2 : p1, &c, thb);
      arcef_bc_ip(thb, ip);
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
  th[4] = arcef_interp(c->lh, c->lk, c->nl, h, 2); th[5] = h;
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
  int bl = arcef_blocklen(r0, N);
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(arc_threads())
#endif
  for (int b = 0; b < NB; b++) {
    double *yb = (double *)malloc((size_t)N * sizeof(double));
    if (yb) {
      unsigned long long s = arcef_seed(*seed, b);
      arcef_blockdraw(m0, r0, N, bl, &s, yb);
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

/* ---- the exact score statistic of the manifold test: F of the residual regression on the six
 * Jacobian columns of the companion mean curve ---- */
/* T from the residual e0 (length n) and the five columns D of arcef_jac; the intercept is added here.
 * The projection is an orthogonalisation (Gram-Schmidt applied twice) of unit-scaled columns, so
 * the value does not depend on the scale of the columns and carries no ridge. Rank-deficient
 * Jacobians are refused: at a point whose mean curve is (near-)flat over the design the columns
 * vanish or coincide, the statistic tends to zero, and a profiling minimiser would live there. */
static double ef_score_core(const double *e0, int n, const double *D, double *Qm /* 6n */) {
  double dmx = (double)n, dmn = (double)n, nrm[6];
  nrm[0] = sqrt((double)n);
  for (int a = 1; a < 6; a++) {
    const double *c = D + (size_t)(a - 1) * n;
    double d = 0.0;
    for (int i = 0; i < n; i++) d += c[i] * c[i];
    if (!isfinite(d)) return NAN;
    if (d > dmx) dmx = d;
    if (d < dmn) dmn = d;
    nrm[a] = sqrt(d);
  }
  if (!(dmx > 0.0) || dmn < 1e-10 * dmx) return NAN;
  double sse0 = 0.0, num = 0.0;
  for (int i = 0; i < n; i++) sse0 += e0[i] * e0[i];
  for (int a = 0; a < 6; a++) {
    double *q = Qm + (size_t)a * n;
    if (a == 0) for (int i = 0; i < n; i++) q[i] = 1.0 / nrm[0];
    else        for (int i = 0; i < n; i++) q[i] = D[(size_t)(a - 1) * n + i] / nrm[a];
    for (int pass = 0; pass < 2; pass++)
      for (int c = 0; c < a; c++) {
        const double *p = Qm + (size_t)c * n;
        double dt = 0.0;
        for (int i = 0; i < n; i++) dt += p[i] * q[i];
        for (int i = 0; i < n; i++) q[i] -= dt * p[i];
      }
    double nn = 0.0;
    for (int i = 0; i < n; i++) nn += q[i] * q[i];
    nn = sqrt(nn);
    if (!(nn > 1e-8)) return NAN;
    double pr = 0.0;
    for (int i = 0; i < n; i++) { q[i] /= nn; pr += q[i] * e0[i]; }
    num += pr * pr;
  }
  double sse1 = sse0 - num;
  if (!(sse1 > 0.0)) return NAN;
  return (num / 6.0) / (sse1 / (double)(n - 6));
}
static double ef_scoreF6(const double *x, const double *y, int n, const double *th,
                         double *wk /* 12n scratch */) {
  double *e0 = wk, *D = wk + n, *Qm = wk + 6 * (size_t)n;
  arcef_jac(x, n, th[2], th[3], th[4], th[5], EF_ZMIN, D);
  for (int i = 0; i < n; i++) e0[i] = y[i] - th[0] - th[1] * D[i];
  return ef_score_core(e0, n, D, Qm);
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
static double ef_affine_resid(const double *y, const double *F, int n, double *e0) {
  double sF = 0.0, sFF = 0.0, sy = 0.0, sFy = 0.0;
  for (int i = 0; i < n; i++) { sF += F[i]; sy += y[i]; }
  sF /= n; sy /= n;
  for (int i = 0; i < n; i++) { sFF += (F[i] - sF) * (F[i] - sF); sFy += (F[i] - sF) * (y[i] - sy); }
  if (!(sFF > 1e-12 * (n * 1.0))) return NAN;          /* F ~ constant: degenerate */
  double g1 = sFy / sFF, sse = 0.0;
  for (int i = 0; i < n; i++) { e0[i] = (y[i] - sy) - g1 * (F[i] - sF); sse += e0[i] * e0[i]; }
  return sse;
}
static double ef_scoreF6p(const double *x, const double *y, int n,
                          double mu, double sg, double al, double be, double *wk /* 12n */) {
  double *e0 = wk, *D = wk + n, *Qm = wk + 6 * (size_t)n;
  /* Identifiability floor, independent of any caller-side window: the six columns must be
   * informed by more than a handful of observations. Design points outside the fitted
   * support are clamped to the endpoints and contribute nothing to any column,
   * so a support that has slid off the design leaves the statistic supported on the few
   * interior points -- and a statistic supported on fewer points than it has columns is
   * degenerate, whatever its conditioning looks like. */
  int nin = 0;
  for (int i = 0; i < n; i++) {
    double z = (x[i] - mu) / sg;
    if (z > EF_ZMIN && z < 1.0 - EF_ZMIN) nin++;
  }
  if (nin < 12 || nin < n / 10) return NAN;
  arcef_jac(x, n, mu, sg, al, be, EF_ZMIN, D);
  if (!isfinite(ef_affine_resid(y, D, n, e0))) return NAN;
  return ef_score_core(e0, n, D, Qm);
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

/* profiling objective at a natural point: mu and sigma inside the windows, alpha in its window
 * (ends included), beta on the tabulated arc. 1e10 outside the searched set or where the
 * statistic is degenerate. */
static double mf_T(const mf_ctx *c, double mu, double sg, double al) {
  double R = c->xrng, tol = 1e-12 * R;
  if (!(al >= c->alo) || !(al <= c->ahi)) return 1e10;
  if (!(mu >= c->xlo - 0.5 * R - tol) || !(mu <= c->xlo + 0.25 * R + tol)) return 1e10;
  if (!(sg >= 0.5 * R - tol) || !(sg <= 3.0 * R + tol)) return 1e10;
  if (!mf_overlap(c, mu, sg)) return 1e10;
  double be = arcef_interp(c->cal, c->cbe, c->nc, al, 1);
  if (!isfinite(be)) return 1e10;
  double v = ef_scoreF6p(c->x, c->y, c->n, mu, sg, al, be, c->wk);
  return isfinite(v) ? v : 1e10;
}
/* residual sum of squares at the same point, after the affine layer is profiled out */
static double mf_sse(const mf_ctx *c, double mu, double sg, double al) {
  double R = c->xrng;
  if (!(al >= c->alo) || !(al <= c->ahi)) return 1e30;
  if (!(mu >= c->xlo - 0.5 * R) || !(mu <= c->xlo + 0.25 * R)) return 1e30;
  if (!(sg >= 0.5 * R) || !(sg <= 3.0 * R) || !mf_overlap(c, mu, sg)) return 1e30;
  double th[6] = { 0.0, 1.0, mu, sg, al, arcef_interp(c->cal, c->cbe, c->nc, al, 1) };
  if (!isfinite(th[5])) return 1e30;
  ef_mbc(c->x, c->n, th, c->wk + c->n, NULL);
  double v = ef_affine_resid(c->y, c->wk + c->n, c->n, c->wk);
  return isfinite(v) ? v : 1e30;
}
/* the search coordinates: p = (mu, e, alpha), e = mu + sigma the right end of the support. The
 * right end is the coordinate the statistic is least regular in, so it is a coordinate. */
typedef struct { const mf_ctx *c; int kind, d, id[3]; double p[3]; } pm_ctx;
static double pm_eval(const pm_ctx *q, const double *v) {
  double p[3] = { q->p[0], q->p[1], q->p[2] };
  for (int j = 0; j < q->d; j++) p[q->id[j]] = v[j];
  return q->kind ? mf_sse(q->c, p[0], p[1] - p[0], p[2]) : mf_T(q->c, p[0], p[1] - p[0], p[2]);
}
/* Nelder-Mead over the free coordinates of q from the simplex (v, v + step), until the simplex is
 * below xtol in units of the steps. Returns the value; v is replaced by the best vertex. */
static double pm_nm(const pm_ctx *q, double *v, const double *step, int maxit, double xtol) {
  int d = q->d;
  double V[4][3], f[4];
  for (int i = 0; i <= d; i++) {
    for (int j = 0; j < d; j++) V[i][j] = v[j];
    if (i) V[i][i - 1] += step[i - 1];
    f[i] = pm_eval(q, V[i]);
    if (i && f[i] >= 1e10) { V[i][i - 1] = v[i - 1] - step[i - 1]; f[i] = pm_eval(q, V[i]); }
  }
  for (int it = 0; it < maxit; it++) {
    int l = 0, h = 0, s2 = -1;
    for (int i = 1; i <= d; i++) { if (f[i] < f[l]) l = i; if (f[i] > f[h]) h = i; }
    for (int i = 0; i <= d; i++) if (i != h && (s2 < 0 || f[i] > f[s2])) s2 = i;
    double ext = 0.0;
    for (int i = 0; i <= d; i++) for (int j = 0; j < d; j++) {
      double e = fabs((V[i][j] - V[l][j]) / step[j]);
      if (e > ext) ext = e;
    }
    if (ext < xtol) break;
    double cen[3] = { 0, 0, 0 }, r[3], t[3];
    for (int i = 0; i <= d; i++) if (i != h) for (int j = 0; j < d; j++) cen[j] += V[i][j] / d;
    for (int j = 0; j < d; j++) r[j] = 2.0 * cen[j] - V[h][j];
    double fr = pm_eval(q, r);
    if (fr < f[l]) {
      for (int j = 0; j < d; j++) t[j] = 3.0 * cen[j] - 2.0 * V[h][j];
      double ft = pm_eval(q, t);
      if (ft < fr) { memcpy(V[h], t, sizeof(t)); f[h] = ft; }
      else         { memcpy(V[h], r, sizeof(r)); f[h] = fr; }
    } else if (fr < f[s2]) { memcpy(V[h], r, sizeof(r)); f[h] = fr; }
    else {
      int out = fr < f[h];
      for (int j = 0; j < d; j++) t[j] = 0.5 * (cen[j] + (out ? r[j] : V[h][j]));
      double ft = pm_eval(q, t);
      if (ft < (out ? fr : f[h])) { memcpy(V[h], t, sizeof(t)); f[h] = ft; }
      else for (int i = 0; i <= d; i++) if (i != l) {
        for (int j = 0; j < d; j++) V[i][j] = 0.5 * (V[i][j] + V[l][j]);
        f[i] = pm_eval(q, V[i]);
      }
    }
  }
  int l = 0;
  for (int i = 1; i <= d; i++) if (f[i] < f[l]) l = i;
  memcpy(v, V[l], (size_t)d * sizeof(double));
  return f[l];
}
typedef struct { double p[3], f; } pm_pt;
static void pm_steps(const mf_ctx *c, const pm_ctx *q, double sc, double *st) {
  const double u[3] = { 0.02 * c->xrng, 0.02 * c->xrng, 0.04 * (c->ahi - c->alo) };
  for (int j = 0; j < q->d; j++) st[j] = sc * u[q->id[j]];
}
/* minimise over the coordinates in `mask` (bit j: coordinate j free) from b, keeping b if no gain */
static void pm_polish(const mf_ctx *c, pm_pt *b, int mask, double sc, int maxit, double xtol) {
  pm_ctx q = { c, 0, 0, { 0, 0, 0 }, { b->p[0], b->p[1], b->p[2] } };
  double v[3], st[3];
  for (int j = 0; j < 3; j++) if (mask & (1 << j)) { v[q.d] = b->p[j]; q.id[q.d++] = j; }
  pm_steps(c, &q, sc, st);
  double f = pm_nm(&q, v, st, maxit, xtol);
  if (f < b->f) { b->f = f; for (int j = 0; j < q.d; j++) b->p[q.id[j]] = v[j]; }
}
/* shrinking restarts with alternating orientation, until two in a row gain nothing */
static void pm_tight(const mf_ctx *c, pm_pt *b, int mask) {
  static const double sc[8] = { 0.5, -0.2, 0.08, -0.03, 0.3, -0.1, 0.01, -0.003 };
  int idle = 0;
  for (int k = 0; k < 8 && idle < 3; k++) {
    double f0 = b->f;
    pm_polish(c, b, mask, sc[k], 250, 1e-8);
    idle = (f0 - b->f > 1e-9) ? 0 : idle + 1;
  }
}
/* minimise with e pinned at the faces of the cell structure: just outside and just inside every
 * design point within `w` of b, and over the cells between them */
static void pm_faces(const mf_ctx *c, pm_pt *best, const pm_pt *b, int w) {
  int n = c->n, j0 = 0;
  for (int i = 1; i < n; i++) if (fabs(c->x[i] - b->p[1]) < fabs(c->x[j0] - b->p[1])) j0 = i;
  for (int j = j0 - w; j <= j0 + w; j++) {
    if (j < 0 || j >= n) continue;
    for (int side = 0; side < 3; side++) {
      pm_pt t = *b;
      if (side == 0) t.p[1] = c->x[j];                                   /* point j just outside */
      else if (side == 1) t.p[1] = c->x[j] + (c->x[j] - t.p[0]) * 1.001 * EF_ZMIN; /* just inside */
      else t.p[1] = (j + 1 < n) ? 0.5 * (c->x[j] + c->x[j + 1]) : c->x[j] + 0.5 * (c->x[j] - c->x[j - 1]);
      t.f = mf_T(c, t.p[0], t.p[1] - t.p[0], t.p[2]);
      if (t.f >= 1e10) continue;
      if (side == 2) pm_polish(c, &t, 7, 0.1, 120, 1e-5);
      else {
        pm_polish(c, &t, 5, 0.3, 60, 1e-4);
        if (side == 1) t.p[1] = c->x[j] + (c->x[j] - t.p[0]) * 1.001 * EF_ZMIN;
        double f = mf_T(c, t.p[0], t.p[1] - t.p[0], t.p[2]);
        if (side == 1) t.f = f;
      }
      if (t.f < best->f + 1e-4) {                       /* worth finishing */
        pm_pt u = t;
        pm_polish(c, &u, side == 2 ? 7 : 5, 0.05, 200, 1e-9);
        pm_polish(c, &u, side == 2 ? 7 : 5, -0.01, 200, 1e-9);
        if (u.f < t.f) t = u;
      }
      if (t.f < best->f) *best = t;
    }
  }
}

/* The profiled minimum.
 *
 * On the searched set the statistic has several local minima: a main one near the least-squares
 * fit, stationary points of the residual surface elsewhere, and, because a design point leaves the
 * fit altogether when the right end of the support passes it, a separate minimum in each cell
 * between neighbouring design points, often ON the face of the cell. A descent from a few starts,
 * as earlier versions had, reaches one of them. The search is therefore global first:
 *   1  a dense grid over (mu, e, alpha), e = mu + sigma, clustered where supports fit the design;
 *   2  least-squares descents from the best residual sums of squares on the grid (the surface is
 *      smooth) and from the caller's start, followed by score descents;
 *   3  score descents from every grid-local minimum (the best PM_NL) and the best grid points;
 *   4  the best distinct results are polished to a tolerance of 1e-9 in the statistic, with the
 *      faces of the cells and the kinks and ends of the tabulated arc searched explicitly.
 * Nothing is random, and the result is the smallest value met. */
#define PM_M 11
#define PM_E 15
#define PM_A 11
#define PM_NL 10
static int pm_cmp(const void *a, const void *b) {
  double d = ((const pm_pt *)a)->f - ((const pm_pt *)b)->f;
  return d < 0 ? -1 : (d > 0 ? 1 : 0);
}
static double ef_profmin(const mf_ctx *c, const double *st, int maxit) {
  static const double eg[PM_E] = { 0.755, 0.80, 0.86, 0.91, 0.95, 0.98, 0.998, 1.012, 1.035, 1.07,
                                   1.13, 1.27, 1.55, 2.1, 2.9 };
  (void)maxit;
  double R = c->xrng, *V = (double *)malloc(PM_M * PM_E * PM_A * 2 * sizeof(double));
  if (!V) return NAN;
  double *W = V + PM_M * PM_E * PM_A;                         /* residual sums of squares */
#define PM_AT(i, j, k) (((i) * PM_E + (j)) * PM_A + (k))
  pm_pt G[PM_NL + 8 + 6 + 6], cand;
  double mus[PM_M], als[PM_A];
  for (int i = 0; i < PM_M; i++) {
    double s = -1.0 + 2.0 * i / (PM_M - 1);
    mus[i] = c->xlo + (s < 0.0 ? -0.5 : 0.25) * R * s * s;
  }
  for (int k = 0; k < PM_A; k++) als[k] = c->alo + (c->ahi - c->alo) * (k + 0.5) / PM_A;
  for (int i = 0; i < PM_M; i++) for (int j = 0; j < PM_E; j++) for (int k = 0; k < PM_A; k++) {
    double e = c->xlo + R * eg[j];
    double v = mf_T(c, mus[i], e - mus[i], als[k]);
    V[PM_AT(i, j, k)] = v;
    double ss = 1e30;
    if (v < 1e10) { ss = 0.0; for (int m = 0; m < c->n; m++) ss += c->wk[m] * c->wk[m]; }
    W[PM_AT(i, j, k)] = ss;                                   /* e0 of the evaluation just made */
  }
  int ng = 0;
  /* grid-local minima of the statistic, best PM_NL kept; then the eight best points outright */
  for (int pass = 0; pass < 2; pass++) {
    int n0 = ng, cap = pass ? 8 : PM_NL;
    for (int i = 0; i < PM_M; i++) for (int j = 0; j < PM_E; j++) for (int k = 0; k < PM_A; k++) {
      double v = V[PM_AT(i, j, k)];
      if (v >= 1e10) continue;
      int ok = 1;
      for (int a = -1; a <= 1 && ok; a++) for (int b2 = -1; b2 <= 1 && ok; b2++)
        for (int d = -1; d <= 1; d++) {
          int ii = i + a, jj = j + b2, kk = k + d;
          if (ii < 0 || jj < 0 || kk < 0 || ii >= PM_M || jj >= PM_E || kk >= PM_A) continue;
          if (V[PM_AT(ii, jj, kk)] < v) { ok = 0; break; }
        }
      if (ok == pass) continue;                               /* pass 0: minima; pass 1: the rest */
      pm_pt t = { { mus[i], c->xlo + R * eg[j], als[k] }, v };
      if (ng - n0 < cap) { G[ng++] = t; }
      else {
        int hi = n0;
        for (int m = n0 + 1; m < ng; m++) if (G[m].f > G[hi].f) hi = m;
        if (v < G[hi].f) G[hi] = t;
      }
    }
  }
  /* least squares: the five best residual sums on the grid, and the caller's start */
  for (int m = 0; m < 6; m++) {
    double p[3];
    if (m < 5) {
      int b = 0;
      for (int q = 1; q < PM_M * PM_E * PM_A; q++) if (W[q] < W[b]) b = q;
      if (W[b] >= 1e30) continue;
      W[b] = 1e30;
      p[0] = mus[b / (PM_E * PM_A)]; p[1] = c->xlo + R * eg[(b / PM_A) % PM_E]; p[2] = als[b % PM_A];
    } else {
      p[0] = c->xlo - 0.5 * R + ef_logis(st[0]) * 0.75 * R;
      p[1] = p[0] + 0.5 * R * exp(ef_logis(st[1]) * 1.7917594692280550);
      p[2] = c->alo + ef_logis(st[2]) * (c->ahi - c->alo);
    }
    pm_ctx q = { c, 1, 3, { 0, 1, 2 }, { 0, 0, 0 } };
    double stp[3];
    pm_steps(c, &q, 1.0, stp);
    pm_nm(&q, p, stp, 200, 1e-5);
    pm_steps(c, &q, -0.1, stp);
    pm_nm(&q, p, stp, 200, 1e-7);
    cand.p[0] = p[0]; cand.p[1] = p[1]; cand.p[2] = p[2];
    cand.f = mf_T(c, p[0], p[1] - p[0], p[2]);
    if (cand.f < 1e10) G[ng++] = cand;
  }
  free(V);
  if (ng == 0) return 1e10;
  for (int m = 0; m < ng; m++) pm_polish(c, &G[m], 7, 0.7, 120, 1e-3);
  qsort(G, (size_t)ng, sizeof(pm_pt), pm_cmp);
  /* best distinct four */
  pm_pt D[4]; int nd = 0;
  for (int m = 0; m < ng && nd < 4; m++) {
    int dup = 0;
    for (int q = 0; q < nd; q++)
      if (fabs(G[m].p[0] - D[q].p[0]) < 0.01 * R && fabs(G[m].p[1] - D[q].p[1]) < 0.002 * R &&
          fabs(G[m].p[2] - D[q].p[2]) < 0.02 * (c->ahi - c->alo)) dup = 1;
    if (!dup) D[nd++] = G[m];
  }
  pm_pt best = D[0];
  for (int m = 0; m < nd; m++) {
    pm_tight(c, &D[m], 7);
    if (D[m].f < best.f) best = D[m];
  }
  qsort(D, (size_t)nd, sizeof(pm_pt), pm_cmp);
  for (int round = 0; round < 2; round++) {
    double f0 = best.f;
    if (round == 0) for (int m = 0; m < nd && m < 2; m++) pm_faces(c, &best, &D[m], m ? 2 : 4);
    else { pm_pt b0 = best; pm_faces(c, &best, &b0, 2); }
    /* kinks and ends of the arc near the best point */
    int k0 = 0;
    for (int k = 1; k < c->nc; k++) if (fabs(c->cal[k] - best.p[2]) < fabs(c->cal[k0] - best.p[2])) k0 = k;
    for (int k = k0 - 1; k <= k0 + 3; k++) {
      double al = k == k0 + 2 ? c->alo : (k == k0 + 3 ? c->ahi : (k >= 0 && k < c->nc ? c->cal[k] : NAN));
      if (!(al >= c->alo) || !(al <= c->ahi)) continue;
      pm_pt t = best;
      t.p[2] = al; t.f = mf_T(c, t.p[0], t.p[1] - t.p[0], al);
      if (t.f >= 1e10) continue;
      pm_polish(c, &t, 3, 0.3, 80, 1e-4);
      if (t.f < best.f + 1e-4) { pm_polish(c, &t, 3, -0.05, 200, 1e-9); if (t.f < best.f) best = t; }
    }
    pm_tight(c, &best, 7);
    pm_polish(c, &best, 5, 0.02, 200, 1e-9);
    if (f0 - best.f < 1e-9) break;
  }
  return best.f;
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
#pragma omp parallel for schedule(dynamic) num_threads(arc_threads())
#endif
  for (int t = 0; t < 4 * R; t++) {
    int arm = t / R, r = t % R;
    double *y  = (double *)malloc((size_t)N * sizeof(double));
    double *m  = (double *)malloc((size_t)N * sizeof(double));
    double *zb = (double *)malloc((size_t)N * sizeof(double));
    double *wk = (double *)malloc((size_t)N * 12 * sizeof(double));
    if (y && m && zb && wk) {
      double th2[6]; memcpy(th2, th0, sizeof(th2));
      th2[5] += shifts[arm];
      ef_mbc(x, N, th2, m, zb);
      unsigned long long s = arcef_seed(*seed + 1000 * arm, r);
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
  double *wk = (double *)malloc((size_t)(*n) * 12 * sizeof(double));
  if (!wk) { *out = NAN; return; }
  *out = ef_scoreF6p(x, y, *n, *mu, *sg, *al, *be, wk);
  free(wk);
}

/* entry: the moving-block length rule applied to a residual vector, exposed so the text can
 * quote the block length the bootstraps actually used. */
void arceqfit_blocklen(const double *r, const int *n, int *out) {
  *out = arcef_blocklen(r, *n);
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
#pragma omp parallel for schedule(dynamic) num_threads(arc_threads())
#endif
  for (int r = 0; r < R; r++) {
    double *y  = (double *)malloc((size_t)N * sizeof(double));
    double *m  = (double *)malloc((size_t)N * sizeof(double));
    double *zb = (double *)malloc((size_t)N * sizeof(double));
    double *wk = (double *)malloc((size_t)N * 12 * sizeof(double));
    if (y && m && zb && wk) {
      ef_mbc(x, N, th0, m, zb);
      unsigned long long s = arcef_seed(*seed, r);
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
#pragma omp parallel for schedule(dynamic) num_threads(arc_threads())
#endif
  for (int t = 0; t < Rr * NN; t++) {
    int ni = t / Rr, r = t % Rr, n = nsizes[ni];
    double *xs = (double *)malloc((size_t)n * sizeof(double));
    if (xs) {
      unsigned long long s = arcef_seed(*seed + 7000 * ni, r);
      for (int i = 0; i < n; i++) xs[i] = ef_pbeta(arcef_unif(&s), *a0, *b0);
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
#pragma omp parallel for schedule(dynamic) num_threads(arc_threads())
#endif
  for (int r = 0; r < Rr; r++) {
    double *xs = (double *)malloc((size_t)N * sizeof(double));
    if (xs) {
      unsigned long long s = arcef_seed(*seed, r);
      for (int i = 0; i < N; i++) xs[i] = ef_pbeta(arcef_unif(&s), *a0, *b0);
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
