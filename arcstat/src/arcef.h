#ifndef ARCEF_H
#define ARCEF_H
/* helpers of arceqfit.c shared with arcbatch.c; not part of any public interface */
double arcef_unif(unsigned long long *s);
unsigned long long arcef_seed(int seed, int b);
double arcef_interp(const double *xs, const double *ys, int n, double x, int rule);
int arcef_blocklen(const double *r, int n);
void arcef_blockdraw(const double *m0, const double *r0, int n, int bl, unsigned long long *s,
                     double *yb);
int arcef_jac(const double *x, int n, double mu, double sg, double al, double be, double zmin,
              double *D);
void arcef_bc_ip(const double *th, double *out);
#endif
