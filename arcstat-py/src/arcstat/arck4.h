#ifndef ARCK4_H
#define ARCK4_H

/* Four-parameter kappa (Hosking 1994) machinery for the induction-curve paper:
   quantile / cdf / pdf, theoretical L-moment ratios by Gauss-Legendre, running-median
   presmoother, band arc lengths (model and sample), the two standard induction-period
   readings, and deterministic Nelder-Mead fits (L-moment shape inversion, quantile-domain
   arc-length shape fit, curve-domain NLS and NALR). All routines are pure functions of
   their arguments so the R and Python front ends agree byte for byte. */

void arck4_q(const double *u, const int *nu, const double *mu, const double *sg,
             const double *k, const double *h, double *out);
void arck4_cdf(const double *x, const int *nx, const double *mu, const double *sg,
               const double *k, const double *h, double *out);
void arck4_pdf(const double *x, const int *nx, const double *mu, const double *sg,
               const double *k, const double *h, double *out);
void arck4_tau34(const double *k, const double *h, const int *nodes, double *out);
double k4_F(double x, double mu, double sg, double k, double h);
void arck4_fit_varpro(const double *, const double *, const int *, const double *,
                      const int *, const int *, const double *, const double *, double *);
void arck4_band_trop_mean(const double *, const int *, const double *, const double *,
                          const double *, const int *, double *);
void arck4_runmed(const double *y, const int *n, const int *w, double *out);
void arck4_band_model(const double *theta, const double *breaks, const int *J,
                      const int *nodes, const double *p, double *out);
void arck4_band_sample(const double *x, const double *y, const int *n,
                       const double *breaks, const int *J, const double *p, double *out);
void arck4_readings(const double *theta, const int *grid, double *out);
void arck4_fit_lmom(const double *t3, const double *t4, const int *nodes, double *out);
void arck4_fit_aleq(const double *ysorted, const int *n, const double *bands,
                    const int *J, const double *start, double *out);
void arck4_fit_nls_drift(const double *x, const double *y, const int *n,
                         const double *start, double *out);
void arck4_fit_nls(const double *x, const double *y, const int *n,
                   const double *start, double *out);
void arck4_fit_nalr(const double *x, const double *y, const int *n, const int *J,
                    const double *lambda, const int *w, const double *p,
                    const double *start, double *out);
typedef double (*objfn)(const double *, const void *);
void arck4_nmd(objfn f, const void *ctx, int d, const double *start, int maxit,
               double *xout, double *fout);
void arck4_mv_boot(const double *x, const double *y, const int *n, const double *cutx,
                   const int *B, const int *seed, const double *bounds, const int *maxit,
                   double *a_pt, double *A);

#endif
