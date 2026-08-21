#ifndef ARCVC_H
#define ARCVC_H

void arcvc_icc_oneway(const double *y, const int *g, const int *n, const int *ng, double *out);
void arcvc_locus_dist(const double *h, const double *k, const int *n,
                      const double *lh, const double *lk, const int *nl, double *out);

#endif
