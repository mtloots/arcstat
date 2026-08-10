/* bayesarc: shared C back-end for the Bayesian arc-length goodness-of-fit test.
 * Pure libm; no external numerical libraries. Bound identically from R and Python. */
#ifndef BAYESARC_H
#define BAYESARC_H

/* M Bayesian-bootstrap posterior discrepancies of the P-P ogive for the SORTED
 * probability-integral transforms u (length n). out has length 2*M:
 * out[2*m] = arc-length discrepancy, out[2*m+1] = Kolmogorov-Smirnov discrepancy. */
void bb_post(const int *n, const double *u, const int *M,
             const int *seed, double *out);

/* H0 reference: D uniform data sets of size n, m posterior draws each, pooled.
 * out has length 2*D*m (arc, ks interleaved), the null distribution of the discrepancy. */
void bb_ref(const int *n, const int *D, const int *m,
            const int *seed, double *out);

#endif
