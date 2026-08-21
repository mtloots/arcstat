#ifndef ARCEQFIT_H
#define ARCEQFIT_H

void arceqfit_bc(const double *x, const double *y, const int *n,
                 const double *starts, const int *ns, const int *np, const int *maxit,
                 const double *cal, const double *cbe, const int *nc,
                 const double *alo, const double *ahi, double *out);
void arceqfit_bc_boot(const double *x, const double *y, const int *n, const double *th,
                      const int *B, const int *seed, const int *maxit,
                      const double *cal, const double *cbe, const int *nc,
                      const double *alo, const double *ahi, double *out_a);
void arceqfit_k4eq(const double *x, const double *y, const int *n,
                   const double *starts, const int *ns, const int *maxit,
                   const double *lh, const double *lk, const int *nl,
                   const double *hlo, const double *hhi, double *out);
void arceqfit_k4eq_boot(const double *x, const double *y, const int *n, const double *th,
                        const double *p0in, const int *B, const int *seed, const int *maxit,
                        const double *lh, const double *lk, const int *nl,
                        const double *hlo, const double *hhi, double *out_a);
void arceqfit_msim(const int *nM, const double *sde, const int *Rm, const int *seed,
                   const double *th0, const double *st5, const int *maxit,
                   const double *cal, const double *cbe, const int *nc,
                   const double *alo, const double *ahi, const double *thref, double *out);
void arceqfit_score_at(const double *x, const double *y, const int *n, const double *mu,
                       const double *sg, const double *al, const double *be, double *out);
void arceqfit_blocklen(const double *r, const int *n, int *out);
void arceqfit_nullT(const int *nM, const double *sde, const int *R2, const int *seed,
                    const double *th0, const double *st5, const int *maxit,
                    const double *cal, const double *cbe, const int *nc,
                    const double *alo, const double *ahi, double *out);
void arceqfit_estsim(const double *a0, const double *b0, const int *R, const int *nsizes,
                     const int *nn, const int *seed, const int *maxit, double *out);
void arceqfit_taus(const double *a0, const double *b0, const int *R, const int *n,
                   const int *seed, double *out);
void arceqfit_bsweep(const double *ks, const int *nk, const double *h, double *out_b);
void arceqfit_absweep(const double *ks, const int *nk, const double *h,
                      double *out_a, double *out_b);
void arceqfit_Esweep(const double *alpha, const double *bg, const int *nb, double *out);

#endif
