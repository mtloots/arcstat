/* cfarc.h : characteristic-function arc length for arcstat.
 * Pure libm; no external numerical libraries. Bound identically from R and Python. */
#ifndef ARCSTAT_CFARC_H
#define ARCSTAT_CFARC_H
#ifdef __cplusplus
extern "C" {
#endif

/* Windowed arc length of the empirical characteristic function of a standardised length-n sample x,
 * over [0, *T] on *ng grid points, doubled for the symmetric half. Writes the arc length to *out. */
void cf_arclength_emp(const double *x, const int *n, const double *T, const int *ng, double *out);

#ifdef __cplusplus
}
#endif
#endif
