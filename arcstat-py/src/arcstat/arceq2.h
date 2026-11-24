#ifndef ARCEQ2_H
#define ARCEQ2_H

void arceq2_bc_q(const int *n, const double *u, const double *alpha,
                 const double *beta, double *out);
void arceq2_bc_cdf(const int *n, const double *x, const double *alpha,
                   const double *beta, double *out);
void arceq2_bc_pdf(const int *n, const double *x, const double *alpha,
                   const double *beta, double *out);
void arceq2_ub_quad(const double *alpha, const double *beta, double *out);
void arceq2_E(const double *alpha, const double *beta, double *out);
void arceq2_bstar(const double *alpha, double *out);

/* unscaled incomplete beta on the logit scale, shared with arceq.c and arceqfit.c */
typedef struct { double a, b, T, thr, Gb, Ga; } ib_par;
void arceq2_ib_init(ib_par *P, double a, double b);
double arceq2_ib_lower(const ib_par *P, double s);
double arceq2_ib_upper(const ib_par *P, double s);
double arceq2_ib_solve(const ib_par *P, double x, int up, double s0);
double arceq2_ib_at_frac(const ib_par *P, double z, double s0);
#endif
