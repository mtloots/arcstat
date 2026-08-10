/* arcdistc.c : shared C back-end for the arc-length distributions package (arcdist).
 * Pure libm, no external numerical libraries. Bound identically from R and Python. */
#include "arcdistc.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

/* shape 1 + sum_r coef[r] P*_{r+1}(u), P* shifted Legendre; recurrence
 * r P_r = (2r-1)(2u-1) P_{r-1} - (r-1) P_{r-2}, P_0=1, P_1=2u-1. */
static double arcq_shape(double u, const double *coef, int k) {
  if (k <= 0) return 1.0;
  double t = 2.0 * u - 1.0, p0 = 1.0, p1 = t, s = 1.0 + coef[0] * p1;
  for (int r = 2; r <= k; r++) {
    double p2 = ((2 * r - 1) * t * p1 - (r - 1) * p0) / r;
    s += coef[r - 1] * p2; p0 = p1; p1 = p2;
  }
  return s;
}

void arcq_qd(const double *u, const int *nu, const double *coef, const int *k,
             const double *sigma, double *out) {
  for (int i = 0; i < *nu; i++) {
    double sh = arcq_shape(u[i], coef, *k);
    out[i] = (*sigma) * (sh > 0 ? sh : 0.0);
  }
}

/* Gauss-Legendre nodes/weights on [0,1] (Newton on the Legendre polynomial) */
static void gauss_legendre(int n, double *x01, double *w01) {
  for (int i = 0; i < n; i++) {
    double z = cos(M_PI * (i + 0.75) / (n + 0.5)), dz;
    for (int it = 0; it < 100; it++) {
      double p0 = 1.0, p1 = z;
      for (int k = 2; k <= n; k++) { double p2 = ((2 * k - 1) * z * p1 - (k - 1) * p0) / k; p0 = p1; p1 = p2; }
      double dp = n * (z * p1 - p0) / (z * z - 1.0);
      dz = p1 / dp; z -= dz; if (fabs(dz) < 1e-15) break;
    }
    double p0 = 1.0, p1 = z;
    for (int k = 2; k <= n; k++) { double p2 = ((2 * k - 1) * z * p1 - (k - 1) * p0) / k; p0 = p1; p1 = p2; }
    double dp = n * (z * p1 - p0) / (z * z - 1.0);
    x01[i] = 0.5 * (z + 1.0);
    w01[i] = 1.0 / ((1.0 - z * z) * dp * dp);   /* = (2/((1-z^2)dp^2))/2 */
  }
}

void arcq_arclength(const double *coef, const int *k, const double *sigma,
                    const int *nodes, double *out) {
  int n = *nodes;
  double *x = (double *)malloc((size_t)n * sizeof(double));
  double *w = (double *)malloc((size_t)n * sizeof(double));
  gauss_legendre(n, x, w);
  double s = 0.0;
  for (int i = 0; i < n; i++) {
    double sh = arcq_shape(x[i], coef, *k); if (sh < 0) sh = 0;
    double qd = (*sigma) * sh;
    s += w[i] * sqrt(1.0 + qd * qd);
  }
  free(x); free(w);
  *out = s;
}

static int cmp_double(const void *a, const void *b) {
  double x = *(const double *)a, y = *(const double *)b;
  return (x > y) - (x < y);
}

/* Closed-form order-two L-moment estimator for the arcq family.  The family's L-moments are an
 * exact linear function of its shape coefficients -- lambda_1 = 1/2 - c1/6, lambda_2 = 1/6 - c2/30,
 * lambda_3 = c1/30, lambda_4 = c2/70 -- so matching inverts without any numerical search:
 *   c2 = 35 t4 / (3 + 7 t4),  c1 = t3 (5 - c2),  sigma = l2/lambda_2,  mu = l1 - sigma lambda_1.
 * The denominator 3 + 7 t4 is bounded below by 5/2 over the family's attainable range t4 >= -1/14,
 * so no guard is needed there; out is set to NA only when the sample L-scale vanishes.
 *
 * The inversion is derived on the ADMISSIBLE set, where the bracket 1 + c1 P1 + c2 P2 stays
 * non-negative so that the positive part is inactive and the L-moments really are the Legendre
 * coefficients.  Sample ratios can land outside it (about 1 per cent of samples at n = 100, and
 * essentially never by n = 400), and there the fitted law's own L-moments are NOT the values matched.
 * out[4] flags this: 1 if the estimate is admissible, 0 if not.
 * out[0..4] = c1, c2, mu, sigma, admissible. */
static int arcq_admissible(double c1, double c2) {
  double A = 6.0 * c2, B = 2.0 * c1 - 6.0 * c2, C = 1.0 - c1 + c2;
  if (C < -1e-12 || A + B + C < -1e-12) return 0;
  if (A > 0.0) {
    double us = -B / (2.0 * A);
    if (us > 0.0 && us < 1.0) return (C - B * B / (4.0 * A) >= -1e-12);
  }
  return 1;
}

void arcq_fit_cf_c(const double *x, const int *n, double *out) {
  int four = 4;
  double L[4];
  sample_lmoments_c(x, n, &four, L);
  if (!(L[1] > 0.0)) { out[0] = out[1] = out[2] = out[3] = out[4] = NAN; return; }
  double t3 = L[2] / L[1], t4 = L[3] / L[1];
  double den = 3.0 + 7.0 * t4;
  if (den == 0.0) { out[0] = out[1] = out[2] = out[3] = out[4] = NAN; return; }
  double c2 = 35.0 * t4 / den;
  double c1 = t3 * (5.0 - c2);
  double lam1 = 0.5 - c1 / 6.0, lam2 = 1.0 / 6.0 - c2 / 30.0;
  double sigma = L[1] / lam2, mu = L[0] - sigma * lam1;
  out[0] = c1; out[1] = c2; out[2] = mu; out[3] = sigma;
  out[4] = (double) arcq_admissible(c1, c2);
}

void sample_lmoments_c(const double *x, const int *n, const int *nmom, double *out) {
  int N = *n, m = *nmom;
  double *v = (double *)malloc((size_t)N * sizeof(double));
  for (int i = 0; i < N; i++) v[i] = x[i];
  qsort(v, (size_t)N, sizeof(double), cmp_double);
  double *b = (double *)malloc((size_t)m * sizeof(double));
  for (int r = 0; r < m; r++) {
    if (r == 0) { double s = 0; for (int i = 0; i < N; i++) s += v[i]; b[0] = s / N; }
    else {
      double s = 0;
      for (int j = 0; j < N; j++) {
        /* weight = C(j, r) / C(N-1, r) as a product (0-based j) */
        double wgt = 1.0;
        for (int t = 0; t < r; t++) wgt *= (double)(j - t) / (double)((N - 1) - t);
        s += wgt * v[j];
      }
      b[r] = s / N;
    }
  }
  for (int r = 0; r < m; r++) {           /* L_{r+1} = sum_k (-1)^{r-k} C(r,k) C(r+k,k) b_k */
    double L = 0.0;
    for (int kk = 0; kk <= r; kk++) {
      double crk = 1.0, crkk = 1.0;
      for (int t = 0; t < kk; t++) { crk *= (double)(r - t) / (double)(t + 1); crkk *= (double)(r + kk - t) / (double)(t + 1); }
      double sign = ((r - kk) % 2 == 0) ? 1.0 : -1.0;
      L += sign * crk * crkk * b[kk];
    }
    out[r] = L;
  }
  free(v); free(b);
}

/* ---- band arc length and the standardised matching estimator ---- */

static double std_normal_pdf(double u) { return exp(-0.5 * u * u) / sqrt(2.0 * M_PI); }

/* inverse standard normal CDF (Acklam's rational approximation, refined by one Halley step) */
static double std_normal_qf(double p) {
  static const double a[6] = {-3.969683028665376e+01, 2.209460984245205e+02, -2.759285104469687e+02,
                               1.383577518672690e+02, -3.066479806614716e+01, 2.506628277459239e+00};
  static const double b[5] = {-5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02,
                               6.680131188771972e+01, -1.328068155288572e+01};
  static const double c_[6] = {-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00,
                               -2.549732539343734e+00, 4.374664141464968e+00, 2.938163982698783e+00};
  static const double d_[4] = {7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00,
                               3.754408661907416e+00};
  double q, r, x;
  if (p <= 0.0 || p >= 1.0) return NAN;
  if (p < 0.02425) { q = sqrt(-2.0 * log(p));
    x = (((((c_[0]*q+c_[1])*q+c_[2])*q+c_[3])*q+c_[4])*q+c_[5]) / ((((d_[0]*q+d_[1])*q+d_[2])*q+d_[3])*q+1.0);
  } else if (p > 1.0 - 0.02425) { q = sqrt(-2.0 * log(1.0 - p));
    x = -(((((c_[0]*q+c_[1])*q+c_[2])*q+c_[3])*q+c_[4])*q+c_[5]) / ((((d_[0]*q+d_[1])*q+d_[2])*q+d_[3])*q+1.0);
  } else { q = p - 0.5; r = q * q;
    x = (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*q / (((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1.0);
  }
  double e = 0.5 * erfc(-x / sqrt(2.0)) - p, u = e * sqrt(2.0 * M_PI) * exp(0.5 * x * x);
  return x - u / (1.0 + 0.5 * x * u);
}

void al_band_model(const double *sigma, const double *a, const double *b, const int *nodes,
                   double *out) {
  double za = std_normal_qf(*a), zb = std_normal_qf(*b);
  int N = (*nodes > 8) ? *nodes : 8;
  /* composite Simpson on [za, zb]: the integrand is analytic and this is ample */
  if (N % 2) N++;
  double h = (zb - za) / N, s = 0.0;
  for (int i = 0; i <= N; i++) {
    double u = za + i * h, f = std_normal_pdf(u);
    double g = sqrt((*sigma) * (*sigma) + f * f);
    double wt = (i == 0 || i == N) ? 1.0 : (i % 2 ? 4.0 : 2.0);
    s += wt * g;
  }
  *out = s * h / 3.0;
}

static int cmp_d(const void *p, const void *q) {
  double x = *(const double *)p, y = *(const double *)q; return (x < y) ? -1 : (x > y);
}

/* type-7 sample quantile, matching R's default */
static double quantile7(const double *v, int n, double p) {
  double hh = (n - 1) * p; int lo = (int)floor(hh);
  if (lo < 0) lo = 0; if (lo > n - 1) lo = n - 1;
  int hi2 = (lo + 1 < n) ? lo + 1 : n - 1;
  return v[lo] + (hh - lo) * (v[hi2] - v[lo]);
}

void al_band_sample(const double *x, const int *n, const double *a, const double *b, double *out) {
  int N = *n;
  double *v = (double *)malloc((size_t)N * sizeof(double));
  for (int i = 0; i < N; i++) v[i] = x[i];
  qsort(v, (size_t)N, sizeof(double), cmp_d);
  double qa = quantile7(v, N, *a), qb = quantile7(v, N, *b);
  double c0 = 1.0 / N, s = 0.0, prev = 0.0; int started = 0, cnt = 0;
  for (int i = 0; i < N; i++) {
    if (v[i] < qa || v[i] > qb) continue;
    if (started) { double d = v[i] - prev; s += sqrt(d * d + c0 * c0); }
    prev = v[i]; started = 1; cnt++;
  }
  free(v);
  *out = (cnt < 2) ? NAN : s;
}

void al_scale(const double *x, const int *n, const double *a, const double *b, double *out) {
  int N = *n;
  double *v = (double *)malloc((size_t)N * sizeof(double));
  for (int i = 0; i < N; i++) v[i] = x[i];
  qsort(v, (size_t)N, sizeof(double), cmp_d);
  double med = quantile7(v, N, 0.5);
  for (int i = 0; i < N; i++) v[i] = fabs(v[i] - med);
  qsort(v, (size_t)N, sizeof(double), cmp_d);
  double s0 = 1.4826 * quantile7(v, N, 0.5);
  free(v);
  if (!(s0 > 0.0) || !isfinite(s0)) { *out = NAN; return; }

  double *z = (double *)malloc((size_t)N * sizeof(double));
  for (int i = 0; i < N; i++) z[i] = x[i] / s0;
  double target; al_band_sample(z, n, a, b, &target);
  free(z);
  if (isnan(target)) { *out = NAN; return; }

  int nodes = 400;
  double lo = 1e-4, hi = 30.0, flo, fhi, fm, mid = 0.0;
  al_band_model(&lo, a, b, &nodes, &flo); flo -= target;
  al_band_model(&hi, a, b, &nodes, &fhi); fhi -= target;
  if (flo * fhi > 0.0) { *out = NAN; return; }              /* no root in the bracket */
  for (int it = 0; it < 200; it++) {
    mid = 0.5 * (lo + hi);
    al_band_model(&mid, a, b, &nodes, &fm); fm -= target;
    if (fm == 0.0 || (hi - lo) < 1e-12) break;
    if (flo * fm < 0.0) { hi = mid; fhi = fm; } else { lo = mid; flo = fm; }
  }
  *out = s0 * mid;
}
