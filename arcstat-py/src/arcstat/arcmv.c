/* arcmv.c -- the fitting-method multiverse of the kappa4 induction-period paper, in the back end.
 *
 * One moving-block residual bootstrap per run, with EIGHT admissible fitting pipelines refitted
 * on every replicate, all reporting the standard tangent reading from the same trace:
 *
 *   1 LS                 four-parameter kappa with drift by variable-projection least squares, with
 *                        the GEV submodel fitted beside it and the nested pair decided by the same
 *                        n log(rss ratio) - 2 rule as the paper's primary fit
 *   2 med9-LS            running median (window 9), then the kappa curve without drift by least squares
 *   3 med21-LS           the same with window 21
 *   4 arc                the banded arc-length estimator
 *   5 GEV-LS             the generalised extreme value submodel in its own right (h held at 0)
 *   6 transient-excised  pipeline 1 on the run past its measured transient cutoff
 *   7 grid-odd           pipeline 1 on the odd-indexed half of the analysis grid
 *   8 grid-even          pipeline 1 on the even-indexed half
 *
 * Every pipeline searches the box the caller passes, h directly and on both sides of zero. Point
 * estimates use the full start grid; replicate refits start from the pipeline's own point fit
 * rather than from the grid: a resampled curve differs from the original only by noise, so one
 * local search from there suffices -- the same policy as the paper's interval bootstrap.
 *
 * Parallelism is OpenMP over replicates. Each replicate draws its resampling indices from its
 * own splitmix64 stream seeded by (seed, replicate index), so the output is byte-identical
 * whatever the thread count and identical between the R and Python fronts, which never touch
 * their host RNGs here at all.
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
#include "arcthreads.h"
#include "arcef.h"
#include "arcbatch.h"


#define MV_M 8
#define MV_GRID 4000                   /* ignored by arck4_readings, which is exact */

/* ---- the start grid: the paper's k4_startsw, which reaches both signs of h (until 0.3.1.9001 it had
 * positive h only, so a curve whose fit has h < 0 was searched from the wrong side). Starts outside
 * the box of this run are dropped, those at its lower bound on h included. With hz the grid is that
 * of the h = 0 submodel, where the h of a start is not used. ---- */
int arck4_starts(const double *x, const double *y, int n, const double *bnd, int hz,
                 double *st /* <= 300*4 */) {
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
  static const double hw[10] = {-30, -12, -4, -1.5, -1, -0.4, 0.02, 0.1, 0.3, 0.8};
  const double h0z = 0.0, *h0 = hz ? &h0z : hw;
  int nh = hz ? 1 : 10;
  int nm = (fabs(mus - mum) > 0) ? 2 : 1;
  double m0[2] = {mus, mum};
  int s = 0;
  for (int a = 0; a < nm; a++) for (int b = 0; b < 3; b++)
    for (int c = 0; c < 5; c++) for (int d = 0; d < nh; d++) {
      if (!hz && !(h0[d] > bnd[2] && h0[d] <= bnd[3])) continue;
      st[s * 4 + 0] = m0[a]; st[s * 4 + 1] = log(sg0 * fz[b]);
      st[s * 4 + 2] = k0[c]; st[s * 4 + 3] = h0[d]; s++;   /* h directly, not log h */
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

/* one start (a replicate): a local search from it; more (a point fit): the start grid, by the batched
 * Levenberg-Marquardt search. o as arck4_fit_varpro */
static int mv_vp(const double *x, const double *y, int n, const double *st, int ns, int maxit,
                 const double *bounds, double hfix, double *o) {
  for (int i = 0; i < 8; i++) o[i] = NAN;
  if (ns < 1) return 0;
  if (ns > 1) return arcbt_k4one(x, y, n, bounds, hfix, o);   /* a point fit: the whole start grid */
  arck4_fit_varpro(x, y, &n, st, &ns, &maxit, bounds, &hfix, o);
  return isfinite(o[7]) && o[7] < 1e11;
}
/* pipeline 1/6/7/8: free fit vs GEV submodel, decided as the paper's primary fit decides. of and og
 * receive the two fits, from which the replicates of this pipeline start. */
static double mv_fit_sel(const double *x, const double *y, int n, const double *sf, int nf,
                         const double *sg, int ng, int maxit, const double *bounds,
                         double *of, double *og) {
  int fok = mv_vp(x, y, n, sf, nf, maxit, bounds, NAN, of);
  int gok = mv_vp(x, y, n, sg, ng, maxit, bounds, 0.0, og);
  if (!fok) return gok ? mv_read_a(og) : NAN;
  if (gok) {
    double dA = n * log(og[7] / of[7]) - 2.0;
    if (isfinite(dA) && dA < 0) return mv_read_a(og);
  }
  return mv_read_a(of);
}
/* the varpro start that repeats a fit */
static void mv_st(const double *o, double *s) {
  s[0] = o[3]; s[1] = log(o[4] > 1e-8 ? o[4] : 1e-8); s[2] = o[5]; s[3] = o[6];
}

/* ---- pipelines 2 and 3: least squares on a presmoothed trace, the curve without drift ----
 * y = g0 + g1 F with g1 > 0 and the two linear coefficients profiled; (mu, log sg, k, h) searched in
 * the run's own box. The objective is that of arck4_fit_nls, which holds h > 0.02 and |k| < 0.9 and
 * therefore cannot be used for a curve with h < 0. */
typedef struct { const double *x, *y; int n; const double *bnd; double *F; } mv_nls;
static double mv_nls_obj(const double *p, const void *vc) {
  const mv_nls *c = (const mv_nls *)vc;
  double sg = exp(p[1]);
  if (!isfinite(sg) || !(sg > 0.0) || !(p[2] >= c->bnd[0] && p[2] <= c->bnd[1] &&
      p[3] >= c->bnd[2] && p[3] <= c->bnd[3])) return 1e12;
  int n = c->n;
  double mF = 0, my = 0;
  for (int i = 0; i < n; i++) {
    c->F[i] = k4_F(c->x[i], p[0], sg, p[2], p[3]);
    mF += c->F[i]; my += c->y[i];
  }
  if (!isfinite(mF)) return 1e12;
  mF /= n; my /= n;
  double FF = 0, Fy = 0;
  for (int i = 0; i < n; i++) { double v = c->F[i] - mF; FF += v * v; Fy += v * (c->y[i] - my); }
  if (!(FF > 1e-16 * n) || !(Fy > 0.0)) return 1e12;
  double g = Fy / FF, s = 0;
  for (int i = 0; i < n; i++) { double r = c->y[i] - my - g * (c->F[i] - mF); s += r * r; }
  return isfinite(s) ? s : 1e12;
}
/* th: (g0, g1, mu, sg, k, h), rss */
static int mv_fit_nls(const double *x, const double *y, int n, const double *st, int ns, int maxit,
                      const double *bnd, double *th) {
  for (int i = 0; i < 7; i++) th[i] = NAN;
  if (ns < 1) return 0;
  double *P = (double *)malloc((size_t)ns * 4 * sizeof(double));
  double *f = (double *)malloc((size_t)ns * sizeof(double));
  double *F = (double *)malloc((size_t)n * sizeof(double));
  int ok = P && f && F;
  if (ok) {
    for (int s = 0; s < ns; s++) f[s] = 1e12;
#ifdef _OPENMP
#pragma omp parallel num_threads(arc_threads())
#endif
    {
      double *Fb = (double *)malloc((size_t)n * sizeof(double));
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 1)
#endif
      for (int s = 0; s < ns; s++) {
        if (!Fb) continue;
        mv_nls c = {x, y, n, bnd, Fb};
        double p1[4], f1, p2[4], f2;
        arck4_nmd(mv_nls_obj, &c, 4, st + (size_t)s * 4, maxit, p1, &f1);
        arck4_nmd(mv_nls_obj, &c, 4, p1, maxit, p2, &f2);
        memcpy(P + (size_t)s * 4, f2 < f1 ? p2 : p1, sizeof p1);
        f[s] = f2 < f1 ? f2 : f1;
      }
      free(Fb);
    }
    int b = 0;
    for (int s = 1; s < ns; s++) if (f[s] < f[b]) b = s;
    ok = f[b] < 1e11;
    if (ok) {
      mv_nls c = {x, y, n, bnd, F};
      const double *p = P + (size_t)b * 4;
      th[6] = mv_nls_obj(p, &c);
      double mF = 0, my = 0, FF = 0, Fy = 0;
      for (int i = 0; i < n; i++) { mF += F[i]; my += y[i]; }
      mF /= n; my /= n;
      for (int i = 0; i < n; i++) { FF += (F[i] - mF) * (F[i] - mF); Fy += (F[i] - mF) * (y[i] - my); }
      th[1] = Fy / FF; th[0] = my - th[1] * mF;
      th[2] = p[0]; th[3] = exp(p[1]); th[4] = p[2]; th[5] = p[3];
    }
  }
  free(P); free(f); free(F);
  return ok;
}
static void mv_st_th(const double *th, double *s) {
  s[0] = th[2]; s[1] = log(th[3] > 1e-8 ? th[3] : 1e-8); s[2] = th[4]; s[3] = th[5];
}

static double mv_read_th(const double *fit7) {
  if (!isfinite(fit7[6])) return NAN;
  for (int i = 0; i < 6; i++) if (!isfinite(fit7[i])) return NAN;
  int grid = MV_GRID; double r[3];
  arck4_readings(fit7, &grid, r);
  return isfinite(r[0]) ? r[0] : NAN;
}

/* the starts of every pipeline. As a point fit: the grid, shared. As a replicate: one start, the
 * pipeline's own point fit, because a resampled curve differs from the original only by noise. */
typedef struct { double f[MV_M][4], g[MV_M][4], arc[7]; int okf[MV_M], okg[MV_M], oka; } mv_warm;

/* all eight pipelines on one curve; w is read when warm is set and written otherwise */
static void mv_fit8(const double *x, const double *y, int n, double cutx, const double *sf, int nf,
                    const double *sg, int ng, int maxit, const double *bounds, int drift, int warm, mv_warm *w,
                    double *out8) {
  for (int m = 0; m < MV_M; m++) out8[m] = NAN;
  double *yb = (double *)malloc((size_t)n * sizeof(double));
  double *xb = (double *)malloc((size_t)n * sizeof(double));
  if (!yb || !xb) { free(yb); free(xb); return; }
  double of[8], og[8], th[8];

  static const int sel[4] = {0, 5, 6, 7};
  for (int q = 0; q < 4; q++) {
    int m = sel[q], nn = 0;
    if (m == 0) { memcpy(xb, x, (size_t)n * sizeof(double)); memcpy(yb, y, (size_t)n * sizeof(double)); nn = n; }
    else for (int i = 0; i < n; i++)
      if (m == 5 ? (x[i] >= cutx) : (i % 2 == m - 6)) { xb[nn] = x[i]; yb[nn] = y[i]; nn++; }
    if (m > 0 && nn <= 30) continue;
    out8[m] = mv_fit_sel(xb, yb, nn, warm ? w->f[m] : sf, warm ? w->okf[m] : nf,
                         warm ? w->g[m] : sg, warm ? w->okg[m] : ng, maxit, bounds, of, og);
    if (!warm) {
      w->okf[m] = isfinite(of[7]) && of[7] < 1e11; w->okg[m] = isfinite(og[7]) && og[7] < 1e11;
      if (w->okf[m]) mv_st(of, w->f[m]);
      if (w->okg[m]) mv_st(og, w->g[m]);
    }
  }
  if (mv_vp(x, y, n, warm ? w->g[4] : sg, warm ? w->okg[4] : ng, maxit, bounds, 0.0, og)) out8[4] = mv_read_a(og);
  if (!warm) { w->okg[4] = isfinite(og[7]) && og[7] < 1e11; if (w->okg[4]) mv_st(og, w->g[4]); }

  for (int m = 1; m <= 2; m++) {
    int win = (m == 1) ? 9 : 21;
    arck4_runmed(y, &n, &win, yb);
    int ok;
    double slope = 0.0;
    if (drift) {
      ok = mv_vp(x, yb, n, warm ? w->f[m] : sf, warm ? w->okf[m] : nf, maxit, bounds, NAN, of) && of[2] > 0.0;
      if (ok) {
        th[0] = of[0]; slope = of[1]; th[1] = of[2]; th[2] = of[3]; th[3] = of[4]; th[4] = of[5]; th[5] = of[6];
        th[6] = of[7];
      }
    } else ok = mv_fit_nls(x, yb, n, warm ? w->f[m] : sf, warm ? w->okf[m] : nf, maxit, bounds, th);
    if (ok) out8[m] = mv_read_th(th);
    if (!warm) {
      w->okf[m] = ok;
      if (ok) mv_st_th(th, w->f[m]);
      if (ok && m == 1) {                /* the banded fit starts from the window-9 fit */
        w->arc[0] = th[0]; w->arc[1] = log(th[1]); w->arc[2] = th[2]; w->arc[3] = log(th[3]);
        w->arc[4] = th[4]; w->arc[5] = th[5]; w->arc[6] = slope;
      }
      if (m == 1) w->oka = ok;
    }
  }
  if (w->oka) {
    int J = 12, w9 = 9; double lambda = 1.0, p = 2.0;
    arck4_fit_nalr_box(x, y, &n, &J, &lambda, &w9, &p, bounds, &drift, w->arc, th);
    if (th[6] < 1e11) {
      out8[3] = mv_read_th(th);
      if (!warm) {
        w->arc[0] = th[0]; w->arc[1] = log(th[1]); w->arc[2] = th[2]; w->arc[3] = log(th[3]);
        w->arc[4] = th[4]; w->arc[5] = th[5]; w->arc[6] = th[7];
      }
    }
  }
  free(yb); free(xb);
}

/* entry point. a_pt gets the eight point estimates (full start grid, in the box the caller passes,
 * which should be that of the run's own point fit); drift: pipelines 2 to 4 carry the linear drift the
 * others carry; A gets B x 8 replicate readings row-major
 * (block-resampled residuals, each pipeline refitted from its own point fit). */
void arck4_mv_boot(const double *x, const double *y, const int *n, const double *cutx,
                   const int *B, const int *seed, const double *bounds, const int *maxit,
                   const int *drift, double *a_pt, double *A) {
  int N = *n, NB = *B, MI = *maxit, DR = *drift;
  for (int m = 0; m < MV_M; m++) a_pt[m] = NAN;
  for (long i = 0; i < (long)(NB > 0 ? NB : 0) * MV_M; i++) A[i] = NAN;

  double *sf = (double *)malloc(300 * 4 * sizeof(double));
  double *sg = (double *)malloc(300 * 4 * sizeof(double));
  if (!sf || !sg) { free(sf); free(sg); return; }
  int nf = arck4_starts(x, y, N, bounds, 0, sf), ng = arck4_starts(x, y, N, bounds, 1, sg);
  mv_warm w;
  memset(&w, 0, sizeof w);
  mv_fit8(x, y, N, *cutx, sf, nf, sg, ng, MI, bounds, DR, 0, &w, a_pt);
  if (NB <= 0 || !w.okf[0]) { free(sf); free(sg); return; }

  /* residuals of the base fit, the paper's primary selection on the whole curve */
  double of[8], og[8];
  const double *base = of;
  mv_vp(x, y, N, w.f[0], 1, MI, bounds, NAN, of);
  if (w.okg[0] && mv_vp(x, y, N, w.g[0], 1, MI, bounds, 0.0, og)) {
    double dA = N * log(og[7] / of[7]) - 2.0;
    if (isfinite(dA) && dA < 0) base = og;
  }
  free(sf); free(sg);

  double *m0 = (double *)malloc((size_t)N * sizeof(double));
  double *r  = (double *)malloc((size_t)N * sizeof(double));
  if (!m0 || !r) { free(m0); free(r); return; }
  for (int i = 0; i < N; i++) {
    double F = k4_F(x[i], base[3], base[4], base[5], base[6]);
    m0[i] = base[0] + base[1] * x[i] + base[2] * F;
    r[i] = y[i] - m0[i];
  }
  int bl = arcef_blocklen(r, N);

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(arc_threads())
#endif
  for (int b = 0; b < NB; b++) {
    double *y2 = (double *)malloc((size_t)N * sizeof(double));
    if (!y2) continue;
    unsigned long long s = arcef_seed(*seed, b);       /* scrambled: no two streams overlap */
    arcef_blockdraw(m0, r, N, bl, &s, y2);
    mv_warm ww = w;
    mv_fit8(x, y2, N, *cutx, NULL, 0, NULL, 0, MI, bounds, DR, 1, &ww, A + (size_t)b * MV_M);
    free(y2);
  }
  free(m0); free(r);
}
