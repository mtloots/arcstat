/* arcdistc: shared C back-end for the arc-length distributions package (arcdist).
 * Pure libm; no external numerical libraries. Bound identically from R and Python. */
#ifndef ARCDISTC_H
#define ARCDISTC_H

/* quantile density Q'(u) = sigma * max(0, 1 + sum_r coef[r] P*_{r+1}(u)) of the arcq family,
 * P* the shifted Legendre polynomials, evaluated at the nu points u; out has length nu. */
void arcq_qd(const double *u, const int *nu, const double *coef, const int *k,
             const double *sigma, double *out);

/* total arc length of the arcq quantile curve, int_0^1 sqrt(1 + Q'(u)^2) du, by
 * Gauss-Legendre quadrature with `nodes` nodes (geometric convergence). out has length 1. */
void arcq_arclength(const double *coef, const int *k, const double *sigma,
                    const int *nodes, double *out);

/* Band arc length of the model distribution curve over the probability band [a,b] under a normal
 * reference of scale sigma: int_{z_a}^{z_b} sqrt(sigma^2 + phi(u)^2) du, by Gauss-Legendre. */
void al_band_model(const double *sigma, const double *a, const double *b, const int *nodes,
                   double *out);

/* Sample band arc length: piecewise-linear empirical-CDF arc length between the a and b sample
 * quantiles of x (length n), the vertical step being 1/n. */
void al_band_sample(const double *x, const int *n, const double *a, const double *b, double *out);

/* Scale estimate by band matching, in the standardised form that is exactly scale equivariant:
 * the data are divided by their MAD, the matching equation is solved by bisection on the
 * standardised scale, and the root is multiplied back. out is NA when the MAD vanishes or the
 * matching equation has no root in the bracket. */
void al_scale(const double *x, const int *n, const double *a, const double *b, double *out);

/* first nmom sample L-moments (Hosking) of x (length n); out has length nmom. */
void sample_lmoments_c(const double *x, const int *n, const int *nmom, double *out);

/* closed-form order-two L-moment fit of the arcq family;
 * out[0..4] = c1, c2, mu, sigma, admissible (1/0). */
void arcq_fit_cf_c(const double *x, const int *n, double *out);

#endif
