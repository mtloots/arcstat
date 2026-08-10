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

#endif
