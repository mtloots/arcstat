/* arclen.h : shared C back-end for the arc-length goodness-of-fit test.
 *
 * One C ABI bound by both R (.C) and Python (ctypes). Every entry point takes only
 * pointer arguments so the SAME compiled library serves both bindings.
 */
#ifndef ARCLEN_H
#define ARCLEN_H
#ifdef __cplusplus
extern "C" {
#endif

/* Arc-length statistic S of the probability plot from n PIT values u in [0,1] (unsorted ok).
 * Writes S to *out. */
void al_statistic(const int *n, const double *u, double *out);

/* Saddlepoint right-tail probability P(S >= s) at sample size n (conditional double saddlepoint).
 * Writes the p-value to *out; writes NaN if the saddlepoint solve fails. */
void al_pvalue(const double *s, const int *n, double *out);

/* Exact mean, variance and support [lo, hi) of S at sample size n.
 * Writes {mean, var, lo, hi} to out[0..3]. */
void al_moments(const int *n, double *out);

#ifdef __cplusplus
}
#endif
#endif
