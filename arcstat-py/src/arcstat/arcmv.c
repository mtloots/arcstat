/* arcmv.c -- the fitting-method multiverse of the kappa4 induction-period paper, in the back end.
 *
 * One moving-block residual bootstrap per run, with EIGHT admissible fitting pipelines refitted
 * on every replicate, all reporting the standard tangent reading from the same trace:
 *
 *   1 LS                 four-parameter kappa by variable-projection least squares, with the
 *                        GEV submodel fitted beside it and the nested pair decided by the same
 *                        n log(rss ratio) - 2 rule as the paper's primary fit
 *   2 med9-LS            running median (window 9) then transformed-parameter NLS
 *   3 med21-LS           running median (window 21) then transformed-parameter NLS
 *   4 arc                the banded arc-length estimator
 *   5 GEV-LS             the generalised extreme value submodel in its own right (h held at 0)
 *   6 transient-excised  pipeline 1 on the run past its measured transient cutoff
 *   7 grid-odd           pipeline 1 on the odd-indexed half of the analysis grid
 *   8 grid-even          pipeline 1 on the even-indexed half
 *
 * Replicate refits start from the base fit rather than from the full start grid: a resampled
 * curve differs from the original only by noise, so one local search from there suffices --
 * the same policy as the paper's interval bootstrap. Point estimates use the full start grid.
 *
 * Parallelism is OpenMP over replicates. Each replicate draws its resampling indices from its
 * own splitmix64 stream seeded by (seed, replicate index), so the output is byte-identical
 * whatever the thread count and identical between the R and Python fronts, which never touch
 * their host RNGs here at all.
 */

#include <math.h>
#include <stdlib.h>
#include <string.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "arck4.h"

#ifdef _OPENMP
static int mv_threads(void) {                 /* the back end's CRAN-safe cap, as in arck4fit.c */
  int t = omp_get_max_threads();
  return t > 2 ? 2 : (t < 1 ? 1 : t);
}
#endif

#define MV_M 8
#define MV_GRID 4000

/* ---- deterministic per-replicate RNG: splitmix64 ---- */
static double mv_unif(unsigned long long *s) {
  unsigned long long z = (*s += 0x9E3779B97F4A7C15ULL);
  z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
  z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
  z ^= (z >> 31);
  return (double)(z >> 11) * (1.0 / 9007199254740992.0);
}

/* ---- the start grid of the paper's k4_starts, ported ---- */
static int mv_starts(const double *x, const double *y, int n, double *st /* <= 120*4 */) {
  double *ys = (double *)malloc((size_t)n * sizeof(double));
  if (!ys) return 0;
  int w = 21;
  arck4_runmed(y, &n, &w, ys);
  double lo = ys[0], hi = ys[0];
  for (int i = 1; i < n; i++) { if (ys[i] < lo) lo = ys[i]; if (ys[i] > hi) hi = ys[i]; }
  int kc = (int)(0.05 * n + 0.5); if (kc < 3) kc = 3;
  /* steepest point of the smoothed polyline, ignoring the opening kc increments */
  int im = kc; double dm = -HUGE_VAL;
  for (int i = kc - 1; i < n - 1; i++) {
    double d = (ys[i + 1] - ys[i]) / (x[i + 1] - x[i]);
    if (d > dm) { dm = d; im = i; }
  }
  double mus = x[im];
  double mum = x[n - 1], x25 = x[n - 1], x75 = x[n - 1];
  for (int i = 0; i < n; i++) if (ys[i] >= lo + 0.50 * (hi - lo)) { mum = x[i]; break; }
  for (int i = 0; i < n; i++) if (ys[i] >= lo + 0.25 * (hi - lo)) { x25 = x[i]; break; }
  for (int i = 0; i < n; i++) if (ys[i] >= lo + 0.75 * (hi - lo)) { x75 = x[i]; break; }
  double sg0 = (x75 - x25) / 1.5, sgf = 0.02 * (x[n - 1] - x[0]);
  if (sg0 < sgf) sg0 = sgf;
  free(ys);
  const double fz[3] = {0.5, 1, 2};
  const double k0[5] = {-0.9, -0.6, -0.3, 0, 0.3};
  const double h0[4] = {0.02, 0.1, 0.3, 0.8};
  int nm = (fabs(mus - mum) > 0) ? 2 : 1;
  double m0[2] = {mus, mum};
  int s = 0;
  for (int a = 0; a < nm; a++) for (int b = 0; b < 3; b++)
    for (int c = 0; c < 5; c++) for (int d = 0; d < 4; d++) {
      st[s * 4 + 0] = m0[a]; st[s * 4 + 1] = log(sg0 * fz[b]);
      st[s * 4 + 2] = k0[c]; st[s * 4 + 3] = log(h0[d]); s++;
    }
  return s;
}

/* varpro out layout: (g0, m, g1, mu, sigma, k, h, rss); readings want (g0, g1, mu, sg, k, h) */
static double mv_read_a(const double *vout) {
  if (!isfinite(vout[7]) || vout[7] >= HUGE_VAL) return NAN;
  double th[6] = {vout[0], vout[2], vout[3], vout[4], vout[5], vout[6]};
  for (int i = 0; i < 6; i++) if (!isfinite(th[i])) return NAN;
  int grid = MV_GRID; double r[3];
  arck4_readings(th, &grid, r);
  return isfinite(r[0]) ? r[0] : NAN;
}

/* pipeline 1/6/7/8: free fit vs GEV submodel, decided as the paper's primary fit decides */
static double mv_fit_sel(const double *x, const double *y, int n, const double *starts,
                         int ns, int maxit, const double *bounds) {
  double hfree = -1.0, hgev = 0.0, of[8], og[8];
  arck4_fit_varpro(x, y, &n, starts, &ns, &maxit, bounds, &hfree, of);
  arck4_fit_varpro(x, y, &n, starts, &ns, &maxit, bounds, &hgev, og);
  int fok = isfinite(of[7]) && of[7] < HUGE_VAL, gok = isfinite(og[7]) && og[7] < HUGE_VAL;
  if (!fok) return gok ? mv_read_a(og) : NAN;
  if (gok) {
    double dA = n * log(og[7] / of[7]) - 2.0;
    if (isfinite(dA) && dA < 0) return mv_read_a(og);
  }
  return mv_read_a(of);
}

/* transformed six-parameter start for nls/nalr from a varpro base fit */
static void mv_tr6(const double *vout, double *st6) {
  double g1 = vout[2] > 1e-8 ? vout[2] : 1e-8;
  double sg = vout[4] > 1e-8 ? vout[4] : 1e-8;
  double h  = vout[6] > 1e-3 ? vout[6] : 1e-3;
  st6[0] = vout[0]; st6[1] = log(g1); st6[2] = vout[3];
  st6[3] = log(sg); st6[4] = vout[5]; st6[5] = log(h);
}

static double mv_read_th(const double *fit7) {
  if (!isfinite(fit7[6])) return NAN;
  for (int i = 0; i < 6; i++) if (!isfinite(fit7[i])) return NAN;
  int grid = MV_GRID; double r[3];
  arck4_readings(fit7, &grid, r);
  return isfinite(r[0]) ? r[0] : NAN;
}

/* all eight pipelines on one curve */
static void mv_fit8(const double *x, const double *y, int n, double cutx,
                    const double *starts, int ns, int maxit, const double *bounds,
                    const double *base_st6, double *out8) {
  for (int m = 0; m < MV_M; m++) out8[m] = NAN;
  double *yb = (double *)malloc((size_t)n * sizeof(double));
  double *xb = (double *)malloc((size_t)n * sizeof(double));
  if (!yb || !xb) { free(yb); free(xb); return; }

  out8[0] = mv_fit_sel(x, y, n, starts, ns, maxit, bounds);

  int w9 = 9, w21 = 21;
  double fit7[7];
  arck4_runmed(y, &n, &w9, yb);
  arck4_fit_nls(x, yb, &n, base_st6, fit7);  out8[1] = mv_read_th(fit7);
  arck4_runmed(y, &n, &w21, yb);
  arck4_fit_nls(x, yb, &n, base_st6, fit7);  out8[2] = mv_read_th(fit7);

  int J = 12, w = 9; double lambda = 1.0, p = 2.0;
  arck4_fit_nalr(x, y, &n, &J, &lambda, &w, &p, base_st6, fit7);
  out8[3] = mv_read_th(fit7);

  double hgev = 0.0, og[8];
  arck4_fit_varpro(x, y, &n, starts, &ns, &maxit, bounds, &hgev, og);
  out8[4] = mv_read_a(og);

  int ne = 0;
  for (int i = 0; i < n; i++) if (x[i] >= cutx) { xb[ne] = x[i]; yb[ne] = y[i]; ne++; }
  if (ne > 30) out8[5] = mv_fit_sel(xb, yb, ne, starts, ns, maxit, bounds);

  int no = 0;
  for (int i = 0; i < n; i += 2) { xb[no] = x[i]; yb[no] = y[i]; no++; }
  if (no > 30) out8[6] = mv_fit_sel(xb, yb, no, starts, ns, maxit, bounds);
  no = 0;
  for (int i = 1; i < n; i += 2) { xb[no] = x[i]; yb[no] = y[i]; no++; }
  if (no > 30) out8[7] = mv_fit_sel(xb, yb, no, starts, ns, maxit, bounds);

  free(yb); free(xb);
}

/* entry point. a_pt gets the eight point estimates (full start grid); A gets B x 8 replicate
 * readings row-major (single-start refits from the base fit, block-resampled residuals). */
void arck4_mv_boot(const double *x, const double *y, const int *n, const double *cutx,
                   const int *B, const int *seed, const double *bounds, const int *maxit,
                   double *a_pt, double *A) {
  int N = *n, NB = *B, MI = *maxit;
  for (int m = 0; m < MV_M; m++) a_pt[m] = NAN;
  for (long i = 0; i < (long)(NB > 0 ? NB : 0) * MV_M; i++) A[i] = NAN;

  double *starts = (double *)malloc(120 * 4 * sizeof(double));
  if (!starts) return;
  int ns = mv_starts(x, y, N, starts);
  if (!ns) { free(starts); return; }

  /* base fit: the paper's primary least-squares selection on the full grid of starts */
  double hfree = -1.0, hgev = 0.0, of[8], og[8];
  arck4_fit_varpro(x, y, &N, starts, &ns, &MI, bounds, &hfree, of);
  arck4_fit_varpro(x, y, &N, starts, &ns, &MI, bounds, &hgev, og);
  int fok = isfinite(of[7]) && of[7] < HUGE_VAL;
  if (!fok) { free(starts); return; }
  const double *base = of;
  if (isfinite(og[7]) && og[7] < HUGE_VAL) {
    double dA = N * log(og[7] / of[7]) - 2.0;
    if (isfinite(dA) && dA < 0) base = og;
  }
  double base_st6[6]; mv_tr6(base, base_st6);
  double base_vp[4] = {base[3], log(base[4] > 1e-8 ? base[4] : 1e-8), base[5],
                       log(base[6] > 1e-3 ? base[6] : 1e-3)};

  /* point estimates: full multi-start */
  mv_fit8(x, y, N, *cutx, starts, ns, MI, bounds, base_st6, a_pt);

  if (NB <= 0) { free(starts); return; }

  /* residuals of the base fit and the acf-matched block length */
  double *m0 = (double *)malloc((size_t)N * sizeof(double));
  double *r  = (double *)malloc((size_t)N * sizeof(double));
  if (!m0 || !r) { free(starts); free(m0); free(r); return; }
  for (int i = 0; i < N; i++) {
    double F = k4_F(x[i], base[3], base[4], base[5], base[6]);
    m0[i] = base[0] + base[1] * x[i] + base[2] * F;
    r[i] = y[i] - m0[i];
  }
  double rb = 0; for (int i = 0; i < N; i++) rb += r[i];
  rb /= N;
  double den = 0; for (int i = 0; i < N; i++) den += (r[i] - rb) * (r[i] - rb);
  int L = N - 2 < 200 ? N - 2 : 200, bl = 0;
  double band = 2.0 / sqrt((double)N);
  for (int l = 1; l <= L; l++) {
    double s = 0;
    for (int i = 0; i + l < N; i++) s += (r[i] - rb) * (r[i + l] - rb);
    if (fabs(s / den) < band) { bl = l; break; }
  }
  if (!bl) { bl = (int)(pow((double)N, 1.0 / 3.0) + 0.5); }
  if (bl < 2) bl = 2;
  int nb = (N + bl - 1) / bl, st_ok = N - bl + 1;

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(mv_threads())
#endif
  for (int b = 0; b < NB; b++) {
    double *y2 = (double *)malloc((size_t)N * sizeof(double));
    if (!y2) continue;
    unsigned long long s = (unsigned long long)(*seed) * 0x9E3779B97F4A7C15ULL
                           ^ ((unsigned long long)(b + 1) * 0xBF58476D1CE4E5B9ULL);
    s = (s ^ (s >> 30)) * 0xBF58476D1CE4E5B9ULL;
    s = (s ^ (s >> 27)) * 0x94D049BB133111EBULL;
    s ^= (s >> 31);
    /* scrambled start: offset-by-golden states are the same splitmix sequence shifted */
    int pos = 0;
    for (int k = 0; k < nb && pos < N; k++) {
      int i0 = (int)(mv_unif(&s) * st_ok); if (i0 >= st_ok) i0 = st_ok - 1;
      for (int i = 0; i < bl && pos < N; i++, pos++) y2[pos] = m0[pos] + r[i0 + i];
    }
    /* single start from the base fit: the replicate differs from the original only by noise */
    int one = 1;
    mv_fit8(x, y2, N, *cutx, base_vp, one, MI, bounds, base_st6, A + (size_t)b * MV_M);
    free(y2);
  }
  free(starts); free(m0); free(r);
}
