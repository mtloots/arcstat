#ifndef ARCSHARED_H
#define ARCSHARED_H
void arcsh_fit(const double *X, const double *Y, const int *off, const int *nr, const int *nsub,
               const int *sub, const int *soff, const double *bounds, const int *drift,
               const int *dense, double *out, double *outrun);
void arcsh_perm(const double *X, const double *Y, const int *off, const int *nprob, const int *poff,
                const int *grp, const double *bounds, const int *drift, const int *dense,
                const double *slack, const int *nslack, const int *maxalloc, double *out, double *p);
#endif
