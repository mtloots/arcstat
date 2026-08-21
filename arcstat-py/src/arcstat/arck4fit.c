/* arck4fit.c -- variable-projection fitting of the drifted kappa response, parallel over starts.
 *
 * The mean curve
 *      y = g0 + m x + g1 F(x; mu, sigma, k, h)
 * is LINEAR in (g0, m, g1), so for any shape those three solve exactly by least squares and the
 * search runs over the four shape parameters alone. That is variable projection, and it is what
 * makes the fit reliable: a direct seven-parameter search on these curves routinely converges to
 * shapes whose induction-period reading is negative or does not exist.
 *
 * The remaining cost is the multi-start, which is embarrassingly parallel: each start is an
 * independent Nelder-Mead descent. The starts are run under OpenMP, each writing only its own slot,
 * and the best is chosen afterwards in a SERIAL pass, so the result does not depend on the thread
 * count or on the order of completion. Thread count is capped at two, per CRAN policy.
 */
#include <math.h>
#include <stdlib.h>
#include <string.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "arck4.h"

#ifdef _OPENMP
static int ak_threads(void) {                 /* only referenced from the parallel region */
  int t = omp_get_max_threads();
  return t > 2 ? 2 : (t < 1 ? 1 : t);        /* CRAN allows at most two */
}
#endif

/* Solve the 3x3 normal equations for (g0, m, g1) against columns (1, x, F) and return the residual
 * sum of squares. Symmetric positive definite in the non-degenerate case; a ridge of last resort
 * keeps a collinear design from producing a spurious zero. */
static double varpro_rss(const double *x, const double *F, const double *y, int n, double *beta) {
  double A[6], b[3];                        /* upper triangle of the 3x3, then the right side */
  double s1 = 0, sx = 0, sf = 0, sxx = 0, sxf = 0, sff = 0, sy = 0, sxy = 0, sfy = 0;
  for (int i = 0; i < n; i++) {
    double xi = x[i], fi = F[i], yi = y[i];
    s1 += 1.0; sx += xi; sf += fi;
    sxx += xi * xi; sxf += xi * fi; sff += fi * fi;
    sy += yi; sxy += xi * yi; sfy += fi * yi;
  }
  A[0] = s1;  A[1] = sx;  A[2] = sf;
  A[3] = sxx; A[4] = sxf; A[5] = sff;
  b[0] = sy;  b[1] = sxy; b[2] = sfy;
  double M[3][3] = { { A[0], A[1], A[2] }, { A[1], A[3], A[4] }, { A[2], A[4], A[5] } };
  for (int d = 0; d < 3; d++) M[d][d] += 1e-10 * (M[d][d] > 0 ? M[d][d] : 1.0);
  /* Gaussian elimination with partial pivoting */
  double aug[3][4];
  for (int r = 0; r < 3; r++) { for (int c = 0; c < 3; c++) aug[r][c] = M[r][c]; aug[r][3] = b[r]; }
  for (int c = 0; c < 3; c++) {
    int piv = c; double best = fabs(aug[c][c]);
    for (int r = c + 1; r < 3; r++) if (fabs(aug[r][c]) > best) { best = fabs(aug[r][c]); piv = r; }
    if (best < 1e-300) return HUGE_VAL;
    if (piv != c) for (int q = 0; q < 4; q++) { double t = aug[c][q]; aug[c][q] = aug[piv][q]; aug[piv][q] = t; }
    for (int r = c + 1; r < 3; r++) {
      double f = aug[r][c] / aug[c][c];
      for (int q = c; q < 4; q++) aug[r][q] -= f * aug[c][q];
    }
  }
  for (int r = 2; r >= 0; r--) {
    double v = aug[r][3];
    for (int c = r + 1; c < 3; c++) v -= aug[r][c] * beta[c];
    beta[r] = v / aug[r][r];
  }
  double rss = 0.0;
  for (int i = 0; i < n; i++) {
    double e = y[i] - (beta[0] + beta[1] * x[i] + beta[2] * F[i]);
    rss += e * e;
  }
  return isfinite(rss) ? rss : HUGE_VAL;
}

typedef struct { const double *x, *y; int n; double *Fbuf; const double *bnd;
                 int hfix; double hval; } vp_ctx;   /* hfix: hold h at hval and search the rest */
/* bnd = (k_lo, k_hi, h_lo, h_hi). The admissible shape region is a MODELLING choice, not a
 * numerical one: on the edible oils a wider k lowered the residual sum of squares yet moved
 * the induction periods away from the laboratory values, so the caller must state it. */

/* objective over the shape: p = (mu, log sigma, k, log h) */
static double vp_obj(const double *p, void *vc) {
  vp_ctx *c = (vp_ctx *)vc;
  double mu = p[0], sg = exp(p[1]), k = p[2];
  double h = c->hfix ? c->hval : exp(p[3]);
  if (!isfinite(sg) || sg <= 0.0 ||
      k < c->bnd[0] || k > c->bnd[1] || h < c->bnd[2] || h > c->bnd[3]) return 1e12;
  for (int i = 0; i < c->n; i++) {
    double v = k4_F(c->x[i], mu, sg, k, h);
    if (!isfinite(v)) return 1e12;
    c->Fbuf[i] = v;
  }
  double lo = c->Fbuf[0], hi = c->Fbuf[0];
  for (int i = 1; i < c->n; i++) { if (c->Fbuf[i] < lo) lo = c->Fbuf[i]; if (c->Fbuf[i] > hi) hi = c->Fbuf[i]; }
  if (hi - lo < 1e-8) return 1e12;
  double beta[3];
  double r = varpro_rss(c->x, c->Fbuf, c->y, c->n, beta);
  return isfinite(r) ? r : 1e12;
}

/* a compact Nelder-Mead over four parameters, deterministic */
static void nmd(double (*f)(const double *, void *), void *ctx, const double *start,
                int d, int maxit, double *xout, double *fout) {
  const int n1 = d + 1;
  double S[5][4], fv[5];
  for (int i = 0; i < n1; i++) {
    for (int j = 0; j < 4; j++) S[i][j] = start[j];
    if (i > 0) S[i][i-1] += (fabs(start[i-1]) > 1e-8 ? 0.10 * fabs(start[i-1]) : 0.10);
    fv[i] = f(S[i], ctx);
  }
  for (int it = 0; it < maxit; it++) {
    int lo = 0, hi = 0, nh = 0;
    for (int i = 1; i < n1; i++) { if (fv[i] < fv[lo]) lo = i; if (fv[i] > fv[hi]) hi = i; }
    nh = (hi == 0) ? 1 : 0;
    for (int i = 0; i < n1; i++) if (i != hi && fv[i] > fv[nh]) nh = i;
    if (fabs(fv[hi] - fv[lo]) <= 1e-12 * (fabs(fv[lo]) + 1e-12)) break;
    double cen[4] = {0,0,0,0};
    for (int i = 0; i < n1; i++) if (i != hi) for (int j = 0; j < d; j++) cen[j] += S[i][j] / d;
    double xr[4], xe[4], xc[4];
    for (int j = 0; j < d; j++) xr[j] = cen[j] + (cen[j] - S[hi][j]);
    double fr = f(xr, ctx);
    if (fr < fv[lo]) {
      for (int j = 0; j < d; j++) xe[j] = cen[j] + 2.0 * (cen[j] - S[hi][j]);
      double fe = f(xe, ctx);
      if (fe < fr) { memcpy(S[hi], xe, sizeof(xe)); fv[hi] = fe; }
      else { memcpy(S[hi], xr, sizeof(xr)); fv[hi] = fr; }
    } else if (fr < fv[nh]) { memcpy(S[hi], xr, sizeof(xr)); fv[hi] = fr; }
    else {
      for (int j = 0; j < d; j++) xc[j] = cen[j] + 0.5 * (S[hi][j] - cen[j]);
      double fc = f(xc, ctx);
      if (fc < fv[hi]) { memcpy(S[hi], xc, sizeof(xc)); fv[hi] = fc; }
      else {
        for (int i = 0; i < n1; i++) if (i != lo) {
          for (int j = 0; j < d; j++) S[i][j] = S[lo][j] + 0.5 * (S[i][j] - S[lo][j]);
          fv[i] = f(S[i], ctx);
        }
      }
    }
  }
  int lo = 0;
  for (int i = 1; i < n1; i++) if (fv[i] < fv[lo]) lo = i;
  memcpy(xout, S[lo], 4 * sizeof(double));
  *fout = fv[lo];
}

/* starts is nstart x 4, row major: (mu, log sigma, k, log h).
 * bounds is (k_lo, k_hi, h_lo, h_hi). out receives (g0, m, g1, mu, sigma, k, h, rss). */
/* hfix < 0 asks for the free four-parameter fit. hfix >= 0 holds h at that value and searches the
 * other three, which is how the SUBMODELS of the family are fitted in their own right: h = 0 is the
 * generalised extreme value distribution and, with k = 0, the Gumbel. That distinction cannot be
 * left to the free search, because h is searched on the log scale and the boundary therefore lies at
 * minus infinity: the optimiser creeps toward it and reports a spurious small h rather than naming
 * the submodel it has effectively chosen. Which of the two fits to report is a model-selection
 * question -- the free model is the larger of a nested pair and can never fit worse -- and it is
 * deliberately left to the caller rather than settled by a tolerance buried here. */
void arck4_fit_varpro(const double *x, const double *y, const int *n, const double *starts,
                      const int *nstart, const int *maxit, const double *bounds,
                      const double *hfix, double *out) {
  int N = *n, S = *nstart, MI = *maxit;
  double *bestp = (double *)calloc((size_t)S * 4, sizeof(double));
  double *bestf = (double *)malloc((size_t)S * sizeof(double));
  if (!bestp || !bestf) { free(bestp); free(bestf); return; }
  for (int s = 0; s < S; s++) bestf[s] = HUGE_VAL;

#ifdef _OPENMP
#pragma omp parallel num_threads(ak_threads())
#endif
  {
    double *Fb = (double *)malloc((size_t)N * sizeof(double));
    if (Fb) {
      vp_ctx c; c.x = x; c.y = y; c.n = N; c.Fbuf = Fb; c.bnd = bounds;
      c.hfix = (*hfix >= 0.0); c.hval = *hfix;
      int dsr = c.hfix ? 3 : 4;    /* with h held, only the other three are searched */
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
      for (int s = 0; s < S; s++) {
        double p1[4], f1, p2[4], f2;
        nmd(vp_obj, &c, starts + (size_t)s * 4, dsr, MI, p1, &f1);
        nmd(vp_obj, &c, p1, dsr, MI, p2, &f2);
        if (f2 < f1) { memcpy(bestp + (size_t)s * 4, p2, sizeof(p2)); bestf[s] = f2; }
        else         { memcpy(bestp + (size_t)s * 4, p1, sizeof(p1)); bestf[s] = f1; }
      }
      free(Fb);
    }
  }
  /* serial selection, so the answer does not depend on the thread count */
  int ib = -1; double fb = HUGE_VAL;
  for (int s = 0; s < S; s++) if (bestf[s] < fb) { fb = bestf[s]; ib = s; }
  if (ib < 0) { out[7] = HUGE_VAL; free(bestp); free(bestf); return; }
  double mu = bestp[ib*4+0], sg = exp(bestp[ib*4+1]), k = bestp[ib*4+2];
  double h = (*hfix >= 0.0) ? *hfix : exp(bestp[ib*4+3]);
  double *Fb = (double *)malloc((size_t)N * sizeof(double));
  double beta[3] = {0,0,0}; double rss = HUGE_VAL;
  if (Fb) {
    for (int i = 0; i < N; i++) Fb[i] = k4_F(x[i], mu, sg, k, h);
    rss = varpro_rss(x, Fb, y, N, beta);
    free(Fb);
  }
  out[0] = beta[0]; out[1] = beta[1]; out[2] = beta[2];
  out[3] = mu; out[4] = sg; out[5] = k; out[6] = h; out[7] = rss;
  free(bestp); free(bestf);
}
