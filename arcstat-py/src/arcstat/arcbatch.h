#ifndef ARCBATCH_H
#define ARCBATCH_H
void arcbt_locus_k(const double *h, const int *n, double *k);
void arcbt_k4eq(const double *X, const double *Y, const int *off, const int *nr,
                const int *npiece, const double *TH, const double *TK, const int *toff,
                const double *bounds, const int *drift, const int *it1, const int *it2,
                const int *K, const int *gate, double *out, double *outp, double *stat);
void arcbt_mvshare(const double *A, const int *nb, const int *nm, const double *frac, double *out);
void arcbt_host(const double *X, const double *Y, const int *off, const int *nr,
                const double *cal, const double *cbe, const int *nc, const double *eps,
                const int *dense, const int *it1, const int *it2, const int *K,
                double *outf, double *outc);
void arcbt_k4fit(const double *X, const double *Y, const int *off, const int *nr, const double *bounds,
                 const double *floors, const int *nfl, const int *gate, const int *gevsel,
                 const double *hfix, const int *drift, const int *it1, const int *it2, const int *K,
                 const int *inf, double *out, double *stat);
void arcbt_wald(const double *X, const double *Y, const int *off, const int *nr, const int *drift,
                const int *it1, const int *it2, const int *K, double *out);
void arcbt_wald_cdf(const double *x, const int *n, const double *mu, const double *lam, double *out);
void arcbt_wald_readings(const double *mu, const double *lam, const int *n, double *out);
void arcbt_lmomw(const double *X, const double *Y, const int *off, const int *nr, const double *lin,
                 const double *bounds, const int *window, const int *ng, const int *nm, double *out);
int arcbt_k4one(const double *x, const double *y, int n, const double *box, double hfix, double *o8);
void arcbt_k4_boot(const double *X, const double *Y, const int *off, const int *nr,
                   const double *fit, const double *box, const int *gev, const int *B,
                   const int *seed, const int *block, const int *maxit, int *bl, double *out);
void arcbt_k4eq_boot(const double *X, const double *Y, const int *off, const int *nr,
                     const double *fit, const int *piece, const int *npiece, const double *TH,
                     const double *TK, const int *toff, const double *bounds, const int *B,
                     const int *seed, const int *block, const int *it2, const int *gate, int *bl,
                     double *out);
#endif
