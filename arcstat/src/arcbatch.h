#ifndef ARCBATCH_H
#define ARCBATCH_H
void arcbt_locus_k(const double *h, const int *n, double *k);
void arcbt_k4eq(const double *X, const double *Y, const int *off, const int *nr,
                const int *npiece, const double *TH, const double *TK, const int *toff,
                const double *bounds, const int *it1, const int *it2, const int *K,
                double *out, double *outp);
void arcbt_host(const double *X, const double *Y, const int *off, const int *nr,
                const double *cal, const double *cbe, const int *nc,
                const int *it1, const int *it2, const int *K, double *outf, double *outc);
void arcbt_k4_boot(const double *X, const double *Y, const int *off, const int *nr,
                   const double *fit, const double *box, const int *gev, const int *B,
                   const int *seed, const int *block, const int *maxit, int *bl, double *out);
void arcbt_k4eq_boot(const double *X, const double *Y, const int *off, const int *nr,
                     const double *fit, const int *piece, const int *npiece, const double *TH,
                     const double *TK, const int *toff, const double *bounds, const int *B,
                     const int *seed, const int *block, const int *it2, int *bl, double *out);
#endif
