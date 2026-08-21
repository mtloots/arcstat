/* Four-parameter kappa machinery for the induction-curve paper. Shared verbatim by the R and
   Python front ends; every routine is a pure function of its arguments (fixed quadrature, fixed
   grids, deterministic Nelder-Mead), so parity between the fronts is byte-identical. */
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
#include "arck4.h"

/* ---------- small utilities ---------------------------------------------------------------- */
static int cmp_d(const void *a, const void *b) {
  double x = *(const double *)a, y = *(const double *)b;
  return (x > y) - (x < y);
}
static double med_of(double *buf, int m) {          /* buf is scratch, m <= 64 */
  qsort(buf, (size_t)m, sizeof(double), cmp_d);
  return (m % 2) ? buf[m / 2] : 0.5 * (buf[m / 2 - 1] + buf[m / 2]);
}
static void gl01(int n, double *x01, double *w01) { /* Gauss-Legendre on [0,1] */
  for (int i = 0; i < n; i++) {
    double z = cos(M_PI * (i + 0.75) / (n + 0.5)), dz;
    for (int it = 0; it < 100; it++) {
      double p0 = 1.0, p1 = z;
      for (int k = 2; k <= n; k++) { double p2 = ((2*k-1)*z*p1 - (k-1)*p0)/k; p0 = p1; p1 = p2; }
      double dp = n * (z * p1 - p0) / (z * z - 1.0);
      dz = p1 / dp; z -= dz; if (fabs(dz) < 1e-15) break;
    }
    double p0 = 1.0, p1 = z;
    for (int k = 2; k <= n; k++) { double p2 = ((2*k-1)*z*p1 - (k-1)*p0)/k; p0 = p1; p1 = p2; }
    double dp = n * (z * p1 - p0) / (z * z - 1.0);
    x01[i] = 0.5 * (z + 1.0);
    w01[i] = 1.0 / ((1.0 - z * z) * dp * dp);
  }
}

/* ---------- the family ---------------------------------------------------------------------- */
/* The four-parameter kappa is written
 *      F(x) = [ 1 - h y ]^(1/h),      y = [ 1 - k (x-mu)/sg ]^(1/k),
 * and BOTH exponents are singular at zero. The two limits are not degenerate cases to be avoided:
 * h -> 0 is the generalised extreme value distribution and k -> 0 within it is the Gumbel, and they
 * are the members this data actually wants. Evaluated literally the expressions fail silently
 * rather than loudly -- at h = 0 one gets 1 - h y = 1 and pow(1, inf) = 1, so the distribution
 * function comes back identically one for every argument, which is a valid double and a nonsense
 * distribution. The limits are therefore taken explicitly:
 *      k -> 0:   y = exp(-(x-mu)/sg)
 *      h -> 0:   F = exp(-y)
 * Both are approached continuously, so a hard switch inside the tolerance costs relative accuracy of
 * the order of the tolerance itself and nothing else. */
#define K4_EPS 1e-8

/* y = (1 - k (x-mu)/sg)^(1/k), and its base z, which the density needs separately */
static void k4_y(double x, double mu, double sg, double k, double *y, double *z) {
  double t = (x - mu) / sg;
  if (fabs(k) < K4_EPS) { *z = 1.0; *y = exp(-t); return; }
  double zz = 1.0 - k * t;
  *z = zz;
  *y = (zz <= 0.0) ? ((k > 0.0) ? 0.0 : HUGE_VAL) : exp(log(zz) / k);
}
static double k4_q(double u, double mu, double sg, double k, double h) {
  double w;                                   /* w = (1 - u^h)/h, the Gumbel reduced variate at h=0 */
  if (fabs(h) < K4_EPS) { if (u <= 0.0) return -HUGE_VAL; w = -log(u); }
  else                  w = (1.0 - pow(u, h)) / h;
  if (w <= 0.0) return (k > 0.0) ? mu + sg / k : HUGE_VAL;
  if (fabs(k) < K4_EPS) return mu - sg * log(w);
  return mu + (sg / k) * (1.0 - exp(k * log(w)));
}
double k4_F(double x, double mu, double sg, double k, double h) {
  double y, z; k4_y(x, mu, sg, k, &y, &z);
  if (!(y >= 0.0)) return 0.0;                /* catches the NaN as well as the negative */
  if (y == HUGE_VAL) return 0.0;
  if (fabs(h) < K4_EPS) return exp(-y);
  double u = 1.0 - h * y;
  if (u <= 0.0) return (h > 0.0) ? 0.0 : 1.0;
  return exp(log(u) / h);
}
static double k4_f(double x, double mu, double sg, double k, double h) {
  double y, z; k4_y(x, mu, sg, k, &y, &z);
  if (!(y > 0.0) || y == HUGE_VAL || z <= 0.0) return 0.0;
  double F = k4_F(x, mu, sg, k, h);
  if (F <= 0.0) return 0.0;
  /* z^(1/k - 1) = y/z, which stays finite as k -> 0 because z -> 1 there */
  return (1.0 / sg) * (y / z) * pow(F, 1.0 - h);
}

void arck4_q(const double *u, const int *nu, const double *mu, const double *sg,
             const double *k, const double *h, double *out) {
  for (int i = 0; i < *nu; i++) out[i] = k4_q(u[i], *mu, *sg, *k, *h);
}
void arck4_cdf(const double *x, const int *nx, const double *mu, const double *sg,
               const double *k, const double *h, double *out) {
  for (int i = 0; i < *nx; i++) out[i] = k4_F(x[i], *mu, *sg, *k, *h);
}
void arck4_pdf(const double *x, const int *nx, const double *mu, const double *sg,
               const double *k, const double *h, double *out) {
  for (int i = 0; i < *nx; i++) out[i] = k4_f(x[i], *mu, *sg, *k, *h);
}

/* ---------- theoretical L-moment ratios (shifted-Legendre quadrature) ---------------------- */
void arck4_tau34(const double *k, const double *h, const int *nodes, double *out) {
  int n = *nodes;
  double *x = (double *)malloc((size_t)n * sizeof(double));
  double *w = (double *)malloc((size_t)n * sizeof(double));
  gl01(n, x, w);
  double l1 = 0, l2 = 0, l3 = 0, l4 = 0;
  for (int i = 0; i < n; i++) {
    double u = x[i], q = k4_q(u, 0.0, 1.0, *k, *h);
    double P1 = 2.0*u - 1.0;
    double P2 = 6.0*u*u - 6.0*u + 1.0;
    double P3 = 20.0*u*u*u - 30.0*u*u + 12.0*u - 1.0;
    l1 += w[i]*q; l2 += w[i]*q*P1; l3 += w[i]*q*P2; l4 += w[i]*q*P3;
  }
  free(x); free(w);
  out[0] = l3 / l2; out[1] = l4 / l2; out[2] = l1; out[3] = l2;
}

/* ---------- running median (shrinking symmetric window at the edges) ------------------------ */
void arck4_runmed(const double *y, const int *n, const int *w, double *out) {
  int N = *n, half = (*w) / 2;
  double buf[129];
  for (int i = 0; i < N; i++) {
    int lo = i - half, hi = i + half;
    if (lo < 0) lo = 0;
    if (hi > N - 1) hi = N - 1;
    int m = hi - lo + 1; if (m > 129) m = 129;
    memcpy(buf, y + lo, (size_t)m * sizeof(double));
    out[i] = med_of(buf, m);
  }
}

/* ---------- band arc lengths ----------------------------------------------------------------
   theta = (g0, g1, mu, sg, k, h); model curve m(x) = g0 + g1 F(x). */
/* The arc-length element is the norm of the vector (dx, dy). Which norm is a choice:
   p = 1 gives the ordinary sum, p = 2 the Euclidean arc length, and p = infinity the max-plus
   (tropical) form max(|dx|, |dy|). R's .C interface refuses a non-finite argument, so the
   tropical case is signalled by ANY NON-POSITIVE p; the R and Python fronts accept p = Inf and
   translate it here. */
static double pnorm2(double a, double b, double p) {
  a = fabs(a); b = fabs(b);
  if (!isfinite(p) || p <= 0.0) return a > b ? a : b;
  if (p == 2.0)     return sqrt(a * a + b * b);      /* the exact old path, bit for bit */
  if (p == 1.0)     return a + b;
  return pow(pow(a, p) + pow(b, p), 1.0 / p);
}

static double band_model_one_d(const double *th, double m, double a, double b, int nodes,
                               const double *x01, const double *w01, double p) {
  double s = 0.0;
  for (int i = 0; i < nodes; i++) {
    double u = a + (b - a) * x01[i];
    double d = m + th[1] * k4_f(u, th[2], th[3], th[4], th[5]);
    s += w01[i] * (b - a) * pnorm2(1.0, d, p);
  }
  return s;
}
static double band_model_one(const double *th, double a, double b, int nodes,
                             const double *x01, const double *w01, double p) {
  return band_model_one_d(th, 0.0, a, b, nodes, x01, w01, p);
}
void arck4_band_model(const double *theta, const double *breaks, const int *J,
                      const int *nodes, const double *p, double *out) {
  int n = *nodes;
  double *x = (double *)malloc((size_t)n * sizeof(double));
  double *w = (double *)malloc((size_t)n * sizeof(double));
  gl01(n, x, w);
  for (int j = 0; j < *J; j++) out[j] = band_model_one(theta, breaks[j], breaks[j+1], n, x, w, *p);
  free(x); free(w);
}
void arck4_band_sample(const double *x, const double *y, const int *n,
                       const double *breaks, const int *J, const double *p, double *out) {
  for (int j = 0; j < *J; j++) out[j] = 0.0;
  for (int i = 0; i + 1 < *n; i++) {
    double xm = 0.5 * (x[i] + x[i+1]);
    double d  = pnorm2(x[i+1]-x[i], y[i+1]-y[i], *p);
    for (int j = 0; j < *J; j++)
      if (xm >= breaks[j] && xm < breaks[j+1]) { out[j] += d; break; }
  }
}

/* ---------- the two standard readings (dense grid, differenced third derivative) ------------ */
void arck4_readings(const double *theta, const int *grid, double *out) {
  double g0 = theta[0], g1 = theta[1], mu = theta[2], sg = theta[3], k = theta[4], h = theta[5];
  (void)g0;
  if (fabs(k) < 1e-4) k = (k >= 0 ? 1.0 : -1.0) * 1e-4;
  double x0 = k4_q(0.001, mu, sg, k, h), x1 = k4_q(0.999, mu, sg, k, h);
  out[0] = out[1] = out[2] = NAN;
  if (!isfinite(x0) || !isfinite(x1) || x1 <= x0) return;
  int G = *grid;
  double *m1 = (double *)malloc((size_t)G * sizeof(double));
  double dx = (x1 - x0) / (G - 1);
  int ic = 0; double best = -1.0;
  for (int i = 0; i < G; i++) {
    double xx = x0 + dx * i;
    m1[i] = g1 * k4_f(xx, mu, sg, k, h);
    if (m1[i] > best) { best = m1[i]; ic = i; }
  }
  /* The steepest point must be an INTERIOR maximum of the density. For this family both terms
     of (log q)' are negative whenever k >= 1 and h < 1, so q is strictly decreasing, f = 1/q is
     strictly increasing, and the grid argmax lands on the right-hand endpoint. The "mode" is
     then a boundary artefact and both readings computed from it are meaningless -- the tangent
     reading in particular is finite and plausible-looking. Neither reading exists there, so
     report NaN. The generic quantile-domain reader carries the same guard. */
  if (ic <= 0 || ic >= G - 1) { free(m1); return; }
  double c0 = x0 + dx * ic;
  out[2] = c0;
  out[0] = c0 - (g1 * k4_F(c0, mu, sg, k, h)) / m1[ic];   /* tangent meets the baseline */
  free(m1);

  /* The third-derivative reading lies BELOW the mode, so spreading a uniform grid over the whole
     support wastes nearly all of it on the right tail. For the skewed shapes these curves prefer
     that tail can run to thousands of hours while the reading sits at a few, and the root is then
     missed entirely. Refine on [x0, c0] alone, at the same cost. */
  if (c0 > x0) {
    double *m2 = (double *)malloc((size_t)G * sizeof(double));
    if (m2 == NULL) return;
    double dz = (c0 - x0) / (G - 1);
    double m2max = 0.0;
    for (int i = 0; i < G; i++) {
      m2[i] = g1 * k4_f(x0 + dz * i, mu, sg, k, h);
      if (m2[i] > m2max) m2max = m2[i];
    }
    /* A sign change between two below-noise second differences is rounding jitter, not an
       inflection: at (k, h) = (1, 1/2) the density is exactly linear, the third derivative
       vanishes identically, and the unthresholded detector manufactured a reading out of the
       last bits. A genuine inflection at this grid density carries |d2| of order f'' dz^2,
       many orders above the floor; the threshold removes the artefact without touching any
       real reading. */
    double thr = 1e-13 * m2max;
    for (int i = 1; i + 2 < G; i++) {
      double d2a = (m2[i+1] - 2.0*m2[i] + m2[i-1]);
      double d2b = (m2[i+2] - 2.0*m2[i+1] + m2[i]);
      if (fabs(d2a) > thr && fabs(d2b) > thr &&
          ((d2a > 0 && d2b < 0) || (d2a < 0 && d2b > 0))) { out[1] = x0 + dz * i; break; }
    }
    free(m2);
  }
}

/* ---------- deterministic Nelder-Mead -------------------------------------------------------- */
static void nelder_mead(objfn f, const void *ctx, int d, const double *start, int maxit,
                        double *xout, double *fout) {
  double S[8][7], fv[8];                                    /* d <= 6 */
  int n1 = d + 1;
  for (int i = 0; i < n1; i++) {
    for (int j = 0; j < d; j++) S[i][j] = start[j];
    if (i > 0) S[i][i-1] += (start[i-1] != 0.0 ? 0.10 * fabs(start[i-1]) : 0.10);
    fv[i] = f(S[i], ctx);
  }
  for (int it = 0; it < maxit; it++) {
    int lo = 0, hi = 0, nh = 0;
    for (int i = 1; i < n1; i++) { if (fv[i] < fv[lo]) lo = i; if (fv[i] > fv[hi]) hi = i; }
    nh = (hi == 0) ? 1 : 0;
    for (int i = 0; i < n1; i++) if (i != hi && fv[i] > fv[nh]) nh = i;
    if (fabs(fv[hi] - fv[lo]) < 1e-12 * (fabs(fv[lo]) + 1e-12)) break;
    double cen[7] = {0};
    for (int i = 0; i < n1; i++) if (i != hi) for (int j = 0; j < d; j++) cen[j] += S[i][j] / d;
    double xr[7], xe[7], xc[7];
    for (int j = 0; j < d; j++) xr[j] = cen[j] + 1.0 * (cen[j] - S[hi][j]);
    double fr = f(xr, ctx);
    if (fr < fv[lo]) {
      for (int j = 0; j < d; j++) xe[j] = cen[j] + 2.0 * (cen[j] - S[hi][j]);
      double fe = f(xe, ctx);
      if (fe < fr) { memcpy(S[hi], xe, (size_t)d * sizeof(double)); fv[hi] = fe; }
      else         { memcpy(S[hi], xr, (size_t)d * sizeof(double)); fv[hi] = fr; }
    } else if (fr < fv[nh]) {
      memcpy(S[hi], xr, (size_t)d * sizeof(double)); fv[hi] = fr;
    } else {
      for (int j = 0; j < d; j++) xc[j] = cen[j] + 0.5 * (S[hi][j] - cen[j]);
      double fc = f(xc, ctx);
      if (fc < fv[hi]) { memcpy(S[hi], xc, (size_t)d * sizeof(double)); fv[hi] = fc; }
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
  memcpy(xout, S[lo], (size_t)d * sizeof(double));
  *fout = fv[lo];
}

/* ---------- fit 1: (t3,t4) -> (k,h) ---------------------------------------------------------- */
typedef struct { double t3, t4; int nodes; } lm_ctx;
static double lm_obj(const double *p, const void *vc) {
  const lm_ctx *c = (const lm_ctx *)vc;
  double o[4]; double k = p[0], h = p[1]; int nd = c->nodes;
  if (h <= 0.01 || h > 3.0 || k < -0.9 || k > 0.9) return 1e6;
  arck4_tau34(&k, &h, &nd, o);
  if (!isfinite(o[0]) || !isfinite(o[1])) return 1e6;
  return (o[0]-c->t3)*(o[0]-c->t3) + (o[1]-c->t4)*(o[1]-c->t4);
}
void arck4_fit_lmom(const double *t3, const double *t4, const int *nodes, double *out) {
  lm_ctx c = { *t3, *t4, *nodes };
  const double starts[3][2] = { {0.0, 0.3}, {0.2, 0.5}, {-0.2, 0.2} };
  double bx[2], bf = 1e300;
  for (int s = 0; s < 3; s++) {
    double x[2], fv;
    nelder_mead(lm_obj, &c, 2, starts[s], 300, x, &fv);
    if (fv < bf) { bf = fv; bx[0] = x[0]; bx[1] = x[1]; }
  }
  out[0] = bx[0]; out[1] = bx[1]; out[2] = bf;
}

/* ---------- fit 2: quantile-domain arc-length shape fit ------------------------------------- */
typedef struct { const double *A; const double *bands; int J; int nodes;
                 const double *x01; const double *w01; } aq_ctx;
static double aq_obj(const double *p, const void *vc) {
  const aq_ctx *c = (const aq_ctx *)vc;
  double al = exp(p[0]), k = p[1], h = exp(p[2]);
  if (k < -0.9 || k > 0.9 || h > 3.0) return 1e6;
  double s = 0.0;
  for (int j = 0; j < c->J; j++) {
    double a = c->bands[2*j], b = c->bands[2*j+1], S = 0.0;
    for (int i = 0; i < c->nodes; i++) {
      double u = a + (b - a) * c->x01[i];
      double qd = al * pow((1.0 - pow(u, h)) / h, k - 1.0) * pow(u, h - 1.0);
      S += c->w01[i] * (b - a) * sqrt(1.0 + qd * qd);
    }
    s += (S - c->A[j]) * (S - c->A[j]);
  }
  return isfinite(s) ? s : 1e6;
}
void arck4_fit_aleq(const double *ysorted, const int *n, const double *bands,
                    const int *J, const double *start, double *out) {
  int N = *n, nJ = *J, nodes = 40;
  double *A = (double *)malloc((size_t)nJ * sizeof(double));
  for (int j = 0; j < nJ; j++) A[j] = 0.0;
  for (int i = 0; i + 1 < N; i++) {
    double u  = (double)(i + 1) / (N + 1), u2 = (double)(i + 2) / (N + 1);
    double um = 0.5 * (u + u2);
    double d  = sqrt((ysorted[i+1]-ysorted[i])*(ysorted[i+1]-ysorted[i]) + (u2-u)*(u2-u));
    for (int j = 0; j < nJ; j++)
      if (um >= bands[2*j] && um < bands[2*j+1]) { A[j] += d; break; }
  }
  double x01[40], w01[40];
  gl01(nodes, x01, w01);
  aq_ctx c = { A, bands, nJ, nodes, x01, w01 };
  double bx[3], bf = 1e300;
  const double st2[3][3] = { {0.0, 0.0, -1.2039728043259360},
                             {0.0, 0.3, -0.6931471805599453},
                             {0.3, -0.2, -1.6094379124341003} };
  for (int s = 0; s < 3; s++) {
    double x[3], fv;
    nelder_mead(aq_obj, &c, 3, s == 0 ? start : st2[s], 400, x, &fv);
    if (fv < bf) { bf = fv; memcpy(bx, x, 3 * sizeof(double)); }
  }
  out[0] = bx[1]; out[1] = exp(bx[2]); out[2] = exp(bx[0]); out[3] = bf;
  free(A);
}

/* ---------- fits 3 and 4: curve-domain NLS and NALR ------------------------------------------
   parameter vector in the optimiser: (g0, log g1, mu, log sg, k, log h). */
static void untrans(const double *p, double *th) {
  th[0] = p[0]; th[1] = exp(p[1]); th[2] = p[2]; th[3] = exp(p[3]); th[4] = p[4]; th[5] = exp(p[5]);
}
typedef struct { const double *x, *y; int n; } nls_ctx;
static double nls_obj(const double *p, const void *vc) {
  const nls_ctx *c = (const nls_ctx *)vc;
  double th[6]; untrans(p, th);
  if (th[4] < -0.9 || th[4] > 0.9 || th[5] < 0.02 || th[5] > 3.0) return 1e12;
  double s = 0.0;
  for (int i = 0; i < c->n; i++) {
    double r = c->y[i] - (th[0] + th[1] * k4_F(c->x[i], th[2], th[3], th[4], th[5]));
    s += r * r;
  }
  return isfinite(s) ? s : 1e12;
}
/* The Rancimat conductivity signal carries a slow linear rise: the instrument tracks the
   conductivity of the water trap, and impurities stripped from the specimen by the air passing
   through it raise it steadily (Focke, personal communication, 2026). Modelling it as
       y = g0 + m x + g1 F(x)
   leaves BOTH induction-period readings unchanged -- the third derivative annihilates the linear
   term, and in the tangent construction to the sloping baseline the drift cancels -- so the
   equivalence theory is untouched and only the fit improves. Parameter vector with drift:
   (g0, m, log g1, mu, log sg, k, log h). */
typedef struct { const double *x, *y; int n; } nlsd_ctx;
static double nlsd_obj(const double *p, const void *vc) {
  const nlsd_ctx *c = (const nlsd_ctx *)vc;
  double th[6]; th[0]=p[0]; th[1]=exp(p[2]); th[2]=p[3]; th[3]=exp(p[4]); th[4]=p[5]; th[5]=exp(p[6]);
  double m = p[1];
  if (th[4] < -0.9 || th[4] > 0.9 || th[5] < 0.02 || th[5] > 3.0) return 1e12;
  double s = 0.0;
  for (int i = 0; i < c->n; i++) {
    double r = c->y[i] - (th[0] + m * c->x[i] + th[1] * k4_F(c->x[i], th[2], th[3], th[4], th[5]));
    s += r * r;
  }
  return isfinite(s) ? s : 1e12;
}
/* start and out are length 7: (g0, m, log g1, mu, log sg, k, log h) in, untransformed out,
   with out[7] the objective. */
void arck4_fit_nls_drift(const double *x, const double *y, const int *n,
                         const double *start, double *out) {
  nlsd_ctx c = { x, y, *n };
  double p1[7], f1, p2[7], f2;
  nelder_mead(nlsd_obj, &c, 7, start, 3000, p1, &f1);
  nelder_mead(nlsd_obj, &c, 7, p1, 3000, p2, &f2);
  const double *b = (f2 < f1 ? p2 : p1);
  out[0]=b[0]; out[1]=b[1]; out[2]=exp(b[2]); out[3]=b[3]; out[4]=exp(b[4]); out[5]=b[5];
  out[6]=exp(b[6]); out[7]=(f2 < f1 ? f2 : f1);
}
void arck4_fit_nls(const double *x, const double *y, const int *n,
                   const double *start, double *out) {
  nls_ctx c = { x, y, *n };
  double p1[6], f1, p2[6], f2;
  nelder_mead(nls_obj, &c, 6, start, 2000, p1, &f1);
  nelder_mead(nls_obj, &c, 6, p1, 2000, p2, &f2);
  untrans(f2 < f1 ? p2 : p1, out);
  out[6] = (f2 < f1 ? f2 : f1);
}
typedef struct { const double *A; const double *breaks; int J; int nodes;
                 const double *x01, *w01; double x0, y0, lambda, p; } na_ctx;
static double na_obj(const double *p, const void *vc) {
  const na_ctx *c = (const na_ctx *)vc;
  double th[6]; untrans(p, th);
  if (th[4] < -0.9 || th[4] > 0.9 || th[5] < 0.02 || th[5] > 3.0) return 1e12;
  double s = 0.0;
  for (int j = 0; j < c->J; j++) {
    double S = band_model_one(th, c->breaks[j], c->breaks[j+1], c->nodes, c->x01, c->w01, c->p);
    s += (S - c->A[j]) * (S - c->A[j]);
  }
  double anc = c->y0 - (th[0] + th[1] * k4_F(c->x0, th[2], th[3], th[4], th[5]));
  s += c->lambda * anc * anc;
  return isfinite(s) ? s : 1e12;
}
void arck4_fit_nalr(const double *x, const double *y, const int *n, const int *J,
                    const double *lambda, const int *w, const double *p,
                    const double *start, double *out) {
  int N = *n, nJ = *J, nodes = 60;
  double *ys = (double *)malloc((size_t)N * sizeof(double));
  arck4_runmed(y, n, w, ys);
  double *br = (double *)malloc((size_t)(nJ + 1) * sizeof(double));
  for (int j = 0; j <= nJ; j++) br[j] = x[0] + (x[N-1] - x[0]) * (double)j / nJ;
  br[nJ] += 1e-9 * (x[N-1] - x[0]);
  double *A = (double *)malloc((size_t)nJ * sizeof(double));
  arck4_band_sample(x, ys, n, br, J, p, A);
  /* baseline anchor: median of the smoothed early 3 per cent */
  int n0 = (int)(0.03 * N); if (n0 < 3) n0 = 3;
  double buf[129]; int m = n0 > 129 ? 129 : n0;
  memcpy(buf, ys, (size_t)m * sizeof(double));
  double y0 = med_of(buf, m), x0 = x[m - 1];
  double x01[60], w01[60];
  gl01(nodes, x01, w01);
  na_ctx c = { A, br, nJ, nodes, x01, w01, x0, y0, *lambda, *p };
  double p1[6], f1, p2[6], f2;
  nelder_mead(na_obj, &c, 6, start, 2000, p1, &f1);
  nelder_mead(na_obj, &c, 6, p1, 2000, p2, &f2);
  untrans(f2 < f1 ? p2 : p1, out);
  out[6] = (f2 < f1 ? f2 : f1);
  free(ys); free(br); free(A);
}

/* ---------- exact mean of the tropical band arc length under Gaussian error ------------------- */
/* The tropical (max-plus) arc-length element is max(dx, |dy|) rather than sqrt(dx^2 + dy^2), and
 * that substitution is what makes its SAMPLING BEHAVIOUR exactly computable. Under y = m + e with
 * e iid N(0, sigma^2) the increment dy = dm + de has de ~ N(0, 2 sigma^2), and
 *      E max(d, |Z|),  Z ~ N(mu, s),
 * is elementary in Phi and phi. The Euclidean element's expectation E sqrt(d^2 + Z^2) is a
 * confluent hypergeometric function instead, so this is a property of the tropical form and not a
 * convenience of notation: the observed band arc length can be compared with its own mean rather
 * than with a clean-curve quantity it does not estimate.
 *
 * ONLY THE MEAN IS SUPPLIED, DELIBERATELY. The variance of a band SUM is not the sum of the
 * elements' variances: consecutive increments share an observation, so they are 1-dependent, and
 * although the noise increments correlate at -1/2 the absolute value discards the sign and leaves
 * the elements POSITIVELY correlated (about +0.15 in the checked configuration). Summing as if
 * independent understates the variance of the band by roughly a fifth. The correction needs
 * bivariate normal quantities and is not elementary; it is not attempted here, and no variance is
 * returned rather than a wrong one. */
static double k4_norm_cdf(double z) { return 0.5 * erfc(-z * M_SQRT1_2); }
static double k4_norm_pdf(double z) { return exp(-0.5 * z * z) / sqrt(2.0 * M_PI); }

/* E max(d, |Z|) for Z ~ N(mu, s^2), d >= 0, s > 0 */
static double k4_emax(double d, double mu, double s) {
  if (!(s > 0.0)) return fmax(d, fabs(mu));
  double a = (d - mu) / s, b = (-d - mu) / s;
  double Pa = k4_norm_cdf(a), Pb = k4_norm_cdf(b);
  double hi = mu * (1.0 - Pa) + s * k4_norm_pdf(a);   /* E[Z 1{Z >  d}] */
  double lo = mu * Pb - s * k4_norm_pdf(b);           /* E[Z 1{Z < -d}] */
  return d * (Pa - Pb) + hi - lo;
}

/* Expected tropical band arc lengths of the fitted mean curve observed with noise.
 * theta is the six-vector (g0, g1, mu, sg, k, h); breaks holds J+1 band edges; out gets J values. */
void arck4_band_trop_mean(const double *x, const int *n, const double *theta, const double *sigma,
                          const double *breaks, const int *J, double *out) {
  int N = *n, JJ = *J;
  double s = sqrt(2.0) * (*sigma);
  for (int j = 0; j < JJ; j++) out[j] = 0.0;
  for (int i = 0; i + 1 < N; i++) {
    double d  = x[i+1] - x[i];
    double m0 = theta[0] + theta[1] * k4_F(x[i],   theta[2], theta[3], theta[4], theta[5]);
    double m1 = theta[0] + theta[1] * k4_F(x[i+1], theta[2], theta[3], theta[4], theta[5]);
    double e  = k4_emax(d, m1 - m0, s);
    /* the segment is attributed to the band holding its right endpoint, as the empirical
     * band arc length does */
    double xr = x[i+1];
    int j = JJ - 1;
    for (int q = 0; q < JJ; q++) if (xr <= breaks[q+1]) { j = q; break; }
    if (j < 0) j = 0;
    out[j] += e;
  }
}

/* exported handle on the deterministic Nelder-Mead, so sibling modules share one optimiser */
void arck4_nmd(objfn f, const void *ctx, int d, const double *start, int maxit,
               double *xout, double *fout) {
  nelder_mead(f, ctx, d, start, maxit, xout, fout);
}
