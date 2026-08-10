#ifndef ARCEQ_H
#define ARCEQ_H
void arceq_readings(const double *alpha, const double *beta, const int *ntheta,
                    const double *theta, const int *ngrid, double *out);
void arceq_readings_vsl(const double *lambda, const double *delta,
                        const int *ngrid, double *out);
#endif
