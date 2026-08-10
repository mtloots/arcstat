/* Closed-form machinery of the incomplete-beta equivalence family: regularised
   incomplete beta (Lentz continued fraction) and its inverse, the beta-companion
   quantile/distribution/density functions, the quadratic shoulder, the closed-form
   equivalence discrepancy, and the equivalence-curve solver. Deterministic throughout;
   shared verbatim by the R and Python fronts. */
#if defined(__clang__)
#pragma clang fp contract(off)
#elif defined(__GNUC__)
#pragma GCC optimize ("fp-contract=off")
#else
#pragma STDC FP_CONTRACT OFF
#endif
#include <math.h>
#include "arceq2.h"

/* log Beta(a,b) */
static double lbeta_(double a, double b) {
  return lgamma(a) + lgamma(b) - lgamma(a + b);
}

/* regularised incomplete beta I_x(a,b), Lentz continued fraction (NR-style) */
static double betacf_(double a, double b, double x) {
  const double eps = 3e-16, fpmin = 1e-300;
  double qab = a + b, qap = a + 1.0, qam = a - 1.0;
  double c = 1.0, d = 1.0 - qab * x / qap;
  if (fabs(d) < fpmin) d = fpmin;
  d = 1.0 / d;
  double h = d;
  for (int m = 1; m <= 300; m++) {
    int m2 = 2 * m;
    double aa = m * (b - m) * x / ((qam + m2) * (a + m2));
    d = 1.0 + aa * d; if (fabs(d) < fpmin) d = fpmin;
    c = 1.0 + aa / c; if (fabs(c) < fpmin) c = fpmin;
    d = 1.0 / d; h *= d * c;
    aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
    d = 1.0 + aa * d; if (fabs(d) < fpmin) d = fpmin;
    c = 1.0 + aa / c; if (fabs(c) < fpmin) c = fpmin;
    d = 1.0 / d;
    double del = d * c; h *= del;
    if (fabs(del - 1.0) < eps) break;
  }
  return h;
}

static double ibeta_(double a, double b, double x) {
  if (x <= 0.0) return 0.0;
  if (x >= 1.0) return 1.0;
  double bt = exp(a * log(x) + b * log1p(-x) - lbeta_(a, b));
  if (x < (a + 1.0) / (a + b + 2.0)) return bt * betacf_(a, b, x) / a;
  return 1.0 - bt * betacf_(b, a, 1.0 - x) / b;
}

/* inverse of the regularised incomplete beta: bisection then Newton polish */
static double ibeta_inv_(double a, double b, double p) {
  if (p <= 0.0) return 0.0;
  if (p >= 1.0) return 1.0;
  double lo = 0.0, hi = 1.0, x = 0.5;
  for (int i = 0; i < 200; i++) {
    x = 0.5 * (lo + hi);
    if (ibeta_(a, b, x) < p) lo = x; else hi = x;
    if (hi - lo < 1e-15) break;
  }
  x = 0.5 * (lo + hi);
  double lb = lbeta_(a, b);
  for (int i = 0; i < 4; i++) {
    double f = ibeta_(a, b, x) - p;
    double dens = exp((a - 1.0) * log(x) + (b - 1.0) * log1p(-x) - lb);
    if (dens <= 0.0) break;
    double step = f / dens;
    double xn = x - step;
    if (xn <= 0.0 || xn >= 1.0) break;
    x = xn;
  }
  return x;
}

/* beta-companion family, standardised support [0, B(a,b)] with a=alpha+1, b=beta+1:
   Q(u) = B(a,b) I_u(a,b);  F(x) = I^{-1}_{x/B(a,b)}(a,b);  f(x) = 1/q(F(x)). */
void arceq2_bc_q(const int *n, const double *u, const double *alpha,
                 const double *beta, double *out) {
  double a = *alpha + 1.0, b = *beta + 1.0, B = exp(lbeta_(a, b));
  for (int i = 0; i < *n; i++) out[i] = B * ibeta_(a, b, u[i]);
}

void arceq2_bc_cdf(const int *n, const double *x, const double *alpha,
                   const double *beta, double *out) {
  double a = *alpha + 1.0, b = *beta + 1.0, B = exp(lbeta_(a, b));
  for (int i = 0; i < *n; i++) {
    double z = x[i] / B;
    out[i] = ibeta_inv_(a, b, z);
  }
}

void arceq2_bc_pdf(const int *n, const double *x, const double *alpha,
                   const double *beta, double *out) {
  double al = *alpha, be = *beta;
  double a = al + 1.0, b = be + 1.0, B = exp(lbeta_(a, b));
  for (int i = 0; i < *n; i++) {
    double w = ibeta_inv_(a, b, x[i] / B);
    out[i] = pow(w, -al) * pow(1.0 - w, -be);
  }
}

/* mode u_c = alpha/(alpha+beta); quadratic shoulder u_b:
   A u^2 + B u + C = 0, A = 2(al+be)^2 + al + be, B = -2 al (2 al + 2 be + 1),
   C = al (2 al + 1); root selected in (0, u_c). out[2] = u_b, u_c (NaN if absent). */
void arceq2_ub_quad(const double *alpha, const double *beta, double *out) {
  double al = *alpha, be = *beta;
  double uc = al / (al + be);
  double A = 2.0 * (al + be) * (al + be) + al + be;
  double B = -2.0 * al * (2.0 * al + 2.0 * be + 1.0);
  double C = al * (2.0 * al + 1.0);
  double disc = B * B - 4.0 * A * C;
  out[0] = NAN; out[1] = uc;
  if (disc < 0.0 || !isfinite(disc)) return;
  double sq = sqrt(disc);
  double r1 = (-B - sq) / (2.0 * A), r2 = (-B + sq) / (2.0 * A);
  double lo = fmin(r1, r2), hi = fmax(r1, r2);
  if (lo > 0.0 && lo < uc) out[0] = lo;
  else if (hi > 0.0 && hi < uc) out[0] = hi;
}

/* closed-form discrepancy E = B(a,b){I_uc - I_ub} - uc^(1+al) (1-uc)^be
   (adaptive quadrature unnecessary: the integral IS the incomplete beta). */
void arceq2_E(const double *alpha, const double *beta, double *out) {
  double al = *alpha, be = *beta;
  double ubuc[2];
  arceq2_ub_quad(alpha, beta, ubuc);
  if (!isfinite(ubuc[0])) { *out = NAN; return; }
  double a = al + 1.0, b = be + 1.0, B = exp(lbeta_(a, b));
  double I = B * (ibeta_(a, b, ubuc[1]) - ibeta_(a, b, ubuc[0]));
  *out = I - pow(ubuc[1], 1.0 + al) * pow(1.0 - ubuc[1], be);
}

/* equivalence curve: beta*(alpha) by scan-then-bisection on E, deterministic */
void arceq2_bstar(const double *alpha, double *out) {
  double lo = NAN, hi = NAN, Elo = NAN;
  double prevb = NAN, prevE = NAN;
  for (int i = 1; i <= 600; i++) {
    double b = -0.999 + (0.999 - 1e-6) * (double)(i - 1) / 599.0;
    double E; arceq2_E(alpha, &b, &E);
    if (isfinite(E) && isfinite(prevE) && ((E < 0.0) != (prevE < 0.0))) {
      lo = prevb; hi = b; Elo = prevE; break;
    }
    if (isfinite(E)) { prevb = b; prevE = E; }
  }
  if (!isfinite(lo)) { *out = NAN; return; }
  for (int i = 0; i < 100; i++) {
    double mid = 0.5 * (lo + hi), Em;
    arceq2_E(alpha, &mid, &Em);
    if (!isfinite(Em)) { *out = NAN; return; }
    if ((Em < 0.0) == (Elo < 0.0)) { lo = mid; Elo = Em; } else { hi = mid; }
    if (hi - lo < 1e-14) break;
  }
  *out = 0.5 * (lo + hi);
}
