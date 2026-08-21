/* cfarc.c : characteristic-function arc length for arcstat.
 * The empirical characteristic function of a sample traces a curve in the complex plane; this routine
 * returns the arc length of that curve over a window, the transform-domain companion of the arc-length
 * goodness-of-fit statistic. Pure libm (complex.h); bound identically from R and Python. */
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

#include "cfarc.h"
#include <complex.h>
#include <math.h>

/* Windowed arc length of the empirical characteristic function of a standardised sample.
 * x is the length-n data vector; the sample is standardised internally (mean 0, sd 1) so the result
 * is scale free. The curve t -> phi_n(t) = mean(exp(i t x)) is followed over [0, *T] on *ng grid
 * points and its length doubled for the symmetric negative half. phi_n'(t) = mean(i x exp(i t x)),
 * whose modulus is the speed; the integral of the speed is the arc length, taken by the trapezoidal
 * rule. Writes the arc length to *out. */
void cf_arclength_emp(const double *x, const int *n, const double *T, const int *ng, double *out) {
  int N = *n, G = *ng;
  double mu = 0.0, s2 = 0.0;
  for (int i = 0; i < N; i++) mu += x[i];
  mu /= N;
  for (int i = 0; i < N; i++) { double d = x[i] - mu; s2 += d * d; }
  double sd = sqrt(s2 / (N - 1));
  double dt = (*T) / (G - 1);
  double prev = 0.0, sum = 0.0;
  for (int g = 0; g < G; g++) {
    double t = g * dt;
    double complex dphi = 0.0 + 0.0 * I;
    for (int i = 0; i < N; i++) {
      double z = (x[i] - mu) / sd;               /* standardised observation */
      dphi += I * z * cexp(I * t * z);
    }
    double speed = cabs(dphi / N);
    if (g > 0) sum += 0.5 * (prev + speed) * dt; /* trapezoidal step */
    prev = speed;
  }
  *out = 2.0 * sum;                              /* symmetric negative half */
}
