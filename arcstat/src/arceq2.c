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

/* (lgamma(x+b) - lgamma(x)) / b to full RELATIVE accuracy however small b is (x > 0, x + b > 0):
   the argument is shifted above 20 and the Stirling terms are differenced analytically. */
static double dlgamma_(double x, double b) {
  if (fabs(b) < 1e-290) b = 1e-290;
  double s = 0.0;
  while (x < 20.0) { s -= log1p(b / x); x += 1.0; }
  double y = x + b;
  s += (x - 0.5) * log1p(b / x) + b * log(y) - b;
  s += -b / (12.0 * x * y);
  s += b * (x * x + x * y + y * y) / (360.0 * x * x * x * y * y * y);
  s += (pow(y, -5.0) - pow(x, -5.0)) / 1260.0;
  return s / b;
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

/* ---- the UNSCALED incomplete beta ------------------------------------------------------------
   lower(s) = int_0^u t^(a-1) (1-t)^(b-1) dt, upper(s) = the rest, u = 1/(1+exp(-s)). Working
   unscaled and on the logit scale keeps full relative accuracy at both ends and, above all, as
   b -> 0 (beta -> -1): there B(a,b) ~ 1/b and the textbook 1 - I_{1-u}(b,a) cancels to nothing,
   losing eps/b. */
void arceq2_ib_init(ib_par *P, double a, double b) {
  P->a = a; P->b = b;
  P->T = exp(lbeta_(a, b));
  P->thr = (a + 1.0) / (a + b + 2.0);
  /* log(b B(a,b)) / b and log(a B(a,b)) / a, accurate as b (a) -> 0 */
  P->Gb = (b < 0.1) ? dlgamma_(1.0, b) - dlgamma_(a, b) : 0.0;
  P->Ga = (a < 0.1) ? dlgamma_(1.0, a) - dlgamma_(b, a) : 0.0;
}
static void ib_logs_(double s, double *lu, double *lv) {
  if (s > 0.0) { *lu = -log1p(exp(-s)); *lv = *lu - s; }
  else         { *lv = -log1p(exp(s));  *lu = *lv + s; }
}
/* tail below u when u is beyond the continued fraction's good half; G as in ib_init */
static double ib_far_(double a, double b, double lu, double lv, double T, double G) {
  double v = exp(lv);
  if (b >= 0.1) return T - exp(a * lu + b * lv) * betacf_(b, a, v) / b;
  double z = G - lv, e = (b == 0.0) ? z : expm1(b * z) / b;
  double c = 1.0, sum = 0.0;
  for (int n = 1; n <= 400; n++) {
    c *= (n - a) / n * v;
    double t = c / (n + b);
    sum += t;
    if (fabs(t) < 1e-18 * (1.0 + fabs(sum))) break;
  }
  return exp(b * lv) * (e - sum);
}
double arceq2_ib_lower(const ib_par *P, double s) {
  double lu, lv; ib_logs_(s, &lu, &lv);
  double u = exp(lu);
  if (u < P->thr) return exp(P->a * lu + P->b * lv) * betacf_(P->a, P->b, u) / P->a;
  return ib_far_(P->a, P->b, lu, lv, P->T, P->Gb);
}
double arceq2_ib_upper(const ib_par *P, double s) {
  double lu, lv; ib_logs_(s, &lu, &lv);
  if (exp(lu) >= P->thr) return exp(P->a * lu + P->b * lv) * betacf_(P->b, P->a, exp(lv)) / P->b;
  return ib_far_(P->b, P->a, lv, lu, P->T, P->Ga);
}
/* s with lower(s) = x (up == 0) or upper(s) = x (up == 1), by a bracketed Halley iteration on the
   logarithm of the tail, which is close to linear in s at both ends. s0: starting point, or NAN. */
double arceq2_ib_solve(const ib_par *P, double x, int up, double s0) {
  double a = P->a, b = P->b, lo = -HUGE_VAL, hi = HUGE_VAL, lt = log(x), s = s0;
  if (!isfinite(s)) {
    s = up ? -(lt + log(b)) / b : (lt + log(a)) / a;
    if (up ? s < 0.0 : s > 0.0) s = 0.0;
    if (!isfinite(s)) s = 0.0;
  }
  for (int it = 0; it < 200; it++) {
    double lu, lv; ib_logs_(s, &lu, &lv);
    double L = up ? arceq2_ib_upper(P, s) : arceq2_ib_lower(P, s), d;
    if (!(L > 0.0) || !isfinite(L)) {               /* out of range: pull towards the middle */
      if (s > 0.0) hi = s; else lo = s;
      s *= 0.5; continue;
    }
    double f = log(L) - lt, f1 = exp(a * lu + b * lv) / L;
    if (up) f1 = -f1;
    if (f == 0.0) return s;
    if ((f > 0.0) != (up != 0)) hi = s; else lo = s;
    double nw = f / f1, den = 1.0 - 0.5 * nw * (a * exp(lv) - b * exp(lu) - f1);
    d = -nw / ((den > 0.2 && den < 5.0) ? den : 1.0);
    double sn = s + d;
    if (fabs(d) < 1e-6 + 1e-15 * fabs(s)) return sn;          /* the step is cubic: done */
    if (!(sn > lo) || !(sn < hi)) {
      if (isfinite(lo) && isfinite(hi)) sn = 0.5 * (lo + hi);
      else sn = (f > 0.0) != (up != 0) ? s - (1.0 + fabs(s)) : s + (1.0 + fabs(s));
    }
    s = sn;
    if (isfinite(lo) && isfinite(hi) && hi - lo <= 4e-16 * fabs(s)) break;
  }
  return s;
}
/* point of the standardised companion distribution at fraction z of its support (0 < z < 1) */
static double ib_at_frac_(const ib_par *P, double z, double s0) {
  return z <= 0.5 ? arceq2_ib_solve(P, z * P->T, 0, s0) : arceq2_ib_solve(P, (1.0 - z) * P->T, 1, s0);
}
double arceq2_ib_at_frac(const ib_par *P, double z, double s0) { return ib_at_frac_(P, z, s0); }

/* beta-companion family, standardised support [0, B(a,b)] with a=alpha+1, b=beta+1:
   Q(u) = B(a,b) I_u(a,b);  F(x) = I^{-1}_{x/B(a,b)}(a,b);  f(x) = 1/q(F(x)). */
void arceq2_bc_q(const int *n, const double *u, const double *alpha,
                 const double *beta, double *out) {
  ib_par P; arceq2_ib_init(&P, *alpha + 1.0, *beta + 1.0);
  for (int i = 0; i < *n; i++) {
    if (!(u[i] > 0.0)) out[i] = (u[i] <= 0.0) ? 0.0 : NAN;
    else if (u[i] >= 1.0) out[i] = P.T;
    else out[i] = arceq2_ib_lower(&P, log(u[i]) - log1p(-u[i]));
  }
}

/* s = logit F(x); *end is -1 / +1 when x is at or beyond an end of the support, 2 for NaN */
static double bc_s_(const ib_par *P, double x, int *end) {
  *end = 0;
  if (!(x > 0.0)) { *end = (x <= 0.0) ? -1 : 2; return 0.0; }
  if (x >= P->T) { *end = 1; return 0.0; }
  return (x <= 0.5 * P->T) ? arceq2_ib_solve(P, x, 0, NAN) : arceq2_ib_solve(P, P->T - x, 1, NAN);
}
void arceq2_bc_cdf(const int *n, const double *x, const double *alpha,
                   const double *beta, double *out) {
  ib_par P; arceq2_ib_init(&P, *alpha + 1.0, *beta + 1.0);
  for (int i = 0; i < *n; i++) {
    int e; double s = bc_s_(&P, x[i], &e);
    out[i] = e == 0 ? (s > 0 ? 1.0 / (1.0 + exp(-s)) : exp(s) / (1.0 + exp(s)))
                    : (e == -1 ? 0.0 : (e == 1 ? 1.0 : NAN));
  }
}

void arceq2_bc_pdf(const int *n, const double *x, const double *alpha,
                   const double *beta, double *out) {
  double al = *alpha, be = *beta;
  ib_par P; arceq2_ib_init(&P, al + 1.0, be + 1.0);
  for (int i = 0; i < *n; i++) {
    int e; double s = bc_s_(&P, x[i], &e), lu, lv;
    if (e == 0) { ib_logs_(s, &lu, &lv); out[i] = exp(-al * lu - be * lv); }
    else if (e == 2) out[i] = NAN;
    else out[i] = (e == -1) ? pow(0.0, -al) : pow(0.0, -be);
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
  /* both roots without cancellation: A vanishes as alpha + beta -> -1/2, C as alpha -> -1/2 */
  double sq = sqrt(disc), qq = -0.5 * (B + (B < 0.0 ? -sq : sq));
  double r1 = (A != 0.0) ? qq / A : HUGE_VAL, r2 = (qq != 0.0) ? C / qq : HUGE_VAL;
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
  ib_par P; arceq2_ib_init(&P, al + 1.0, be + 1.0);
  double I = arceq2_ib_lower(&P, log(ubuc[1]) - log1p(-ubuc[1]))
           - arceq2_ib_lower(&P, log(ubuc[0]) - log1p(-ubuc[0]));
  *out = I - pow(ubuc[1], 1.0 + al) * pow(1.0 - ubuc[1], be);
}

/* equivalence curve: beta*(alpha) by scan-then-bisection on E, deterministic. The scan is uniform
   on (-0.999, -1e-6) and geometric beyond both ends, so the roots next to beta = -1 (alpha at the
   end of the curve) and next to beta = 0 (alpha near -1/2) are not lost. */
void arceq2_bstar(const double *alpha, double *out) {
  double lo = NAN, hi = NAN, Elo = NAN;
  double prevb = NAN, prevE = NAN;
  for (int i = -10; i <= 606; i++) {
    double b = -0.999 + (0.999 - 1e-6) * (double)(i - 1) / 599.0;
    if (i < 1) b = -1.0 + pow(10.0, -3.0 + (double)(i - 1));        /* -1 + 1e-14 ... -1 + 1e-4 */
    if (i > 600) b = -pow(10.0, -6.0 - (double)(i - 600));          /* -1e-7 ... -1e-12 */
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
