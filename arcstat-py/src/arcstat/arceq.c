/* Equivalence-system machinery: quantile-domain induction readings and the equivalence
   discrepancy for tilted beta-kernel quantile densities
       q(u) = u^alpha (1-u)^beta exp( sum_j theta_j P_j(u) ),
   with P_j the shifted Legendre polynomials (orders 1..3 supported). Deterministic grid
   evaluation; shared verbatim by the R and Python fronts. */
#if defined(__clang__)
#pragma clang fp contract(off)
#elif defined(__GNUC__)
#pragma GCC optimize ("fp-contract=off")
#else
#pragma STDC FP_CONTRACT OFF
#endif
#include <math.h>
#include <stdlib.h>
#include "arceq.h"

/* van Staden-Loots host: q(u) = (1-d) u^(l-1) + d (1-u)^(l-1) (scale-free) */

/* out[5] = D, a, b, u_c, u_b (a and b as quantile values with Q(grid start) = 0) */
void arceq_readings(const double *alpha, const double *beta, const int *ntheta,
                    const double *theta, const int *ngrid, double *out) {
  int N = *ngrid, nt = *ntheta;
  double al = *alpha, be = *beta;
  double *u  = (double *)malloc((size_t)N * sizeof(double));
  double *q  = (double *)malloc((size_t)N * sizeof(double));
  double *gp = (double *)malloc((size_t)N * sizeof(double));
  double *gpp= (double *)malloc((size_t)N * sizeof(double));
  for (int i = 0; i < N; i++) {
    double ui = 1e-5 + (1.0 - 2e-5) * (double)i / (N - 1);
    u[i] = ui;
    double g   = al * log(ui) + be * log(1.0 - ui);
    double g1  = al / ui - be / (1.0 - ui);
    double g2  = -al / (ui * ui) - be / ((1.0 - ui) * (1.0 - ui));
    if (nt >= 1) { g += theta[0] * (2.0*ui - 1.0); g1 += theta[0] * 2.0; }
    if (nt >= 2) { g += theta[1] * (6.0*ui*ui - 6.0*ui + 1.0);
                   g1 += theta[1] * (12.0*ui - 6.0); g2 += theta[1] * 12.0; }
    if (nt >= 3) { g += theta[2] * (20.0*ui*ui*ui - 30.0*ui*ui + 12.0*ui - 1.0);
                   g1 += theta[2] * (60.0*ui*ui - 60.0*ui + 12.0); g2 += theta[2] * (120.0*ui - 60.0); }
    q[i] = exp(g); gp[i] = g1; gpp[i] = g2;
  }
  out[0] = out[1] = out[2] = out[3] = out[4] = NAN;
  int ic = 0; double qmin = q[0];
  for (int i = 1; i < N; i++) if (q[i] < qmin) { qmin = q[i]; ic = i; }
  if (ic <= 10 || ic >= N - 10) { free(u); free(q); free(gp); free(gpp); return; }
  out[3] = u[ic];
  int ib = -1;
  for (int i = 6; i + 1 < ic; i++) {              /* smallest root of 2 g'^2 = g'' below u_c */
    double r0 = 2.0*gp[i]*gp[i] - gpp[i], r1 = 2.0*gp[i+1]*gp[i+1] - gpp[i+1];
    if ((r0 > 0 && r1 < 0) || (r0 < 0 && r1 > 0)) { ib = i; break; }
  }
  if (ib < 0) { free(u); free(q); free(gp); free(gpp); return; }
  out[4] = u[ib];
  /* Q by trapezoid from the grid start; a and b in those units */
  double Qb = 0.0, Qc = 0.0, acc = 0.0;
  for (int i = 0; i + 1 < N && i + 1 <= ic; i++) {
    acc += 0.5 * (q[i] + q[i+1]) * (u[i+1] - u[i]);
    if (i + 1 == ib) Qb = acc;
    if (i + 1 == ic) Qc = acc;
  }
  out[1] = Qc - u[ic] * q[ic];                     /* tangent reading a */
  out[2] = Qb;                                     /* derivative reading b */
  out[0] = out[2] - out[1];                        /* D = b - a */
  free(u); free(q); free(gp); free(gpp);
}

void arceq_readings_vsl(const double *lambda, const double *delta, const int *ngrid, double *out) {
  int N = *ngrid; double l = *lambda, d = *delta;
  double *u = (double *)malloc((size_t)N * sizeof(double));
  double *q = (double *)malloc((size_t)N * sizeof(double));
  double *gp = (double *)malloc((size_t)N * sizeof(double));
  double *gpp= (double *)malloc((size_t)N * sizeof(double));
  for (int i = 0; i < N; i++) {
    double ui = 1e-5 + (1.0 - 2e-5) * (double)i / (N - 1);
    u[i] = ui;
    double t1 = (1.0-d)*pow(ui, l-1.0), t2 = d*pow(1.0-ui, l-1.0);
    double qv = t1 + t2;
    double q1 = (1.0-d)*(l-1.0)*pow(ui, l-2.0) - d*(l-1.0)*pow(1.0-ui, l-2.0);
    double q2 = (1.0-d)*(l-1.0)*(l-2.0)*pow(ui, l-3.0) + d*(l-1.0)*(l-2.0)*pow(1.0-ui, l-3.0);
    q[i] = qv; gp[i] = q1/qv; gpp[i] = q2/qv - (q1/qv)*(q1/qv);
  }
  out[0] = out[1] = out[2] = out[3] = out[4] = NAN;
  int ic = 0; double qmin = q[0];
  for (int i = 1; i < N; i++) if (q[i] < qmin) { qmin = q[i]; ic = i; }
  if (ic <= 10 || ic >= N - 10) { free(u); free(q); free(gp); free(gpp); return; }
  out[3] = u[ic];
  int ib = -1;
  for (int i = 0; i + 1 < ic; i++) {
    if (u[i] < 1e-3) continue;
    double r0 = 2.0*gp[i]*gp[i] - gpp[i], r1 = 2.0*gp[i+1]*gp[i+1] - gpp[i+1];
    if ((r0 > 0 && r1 < 0) || (r0 < 0 && r1 > 0)) { ib = i; break; }
  }
  if (ib < 0) { free(u); free(q); free(gp); free(gpp); return; }
  out[4] = u[ib];
  double Qb = 0.0, Qc = 0.0, acc = 0.0;
  for (int i = 0; i + 1 < N && i + 1 <= ic; i++) {
    acc += 0.5 * (q[i] + q[i+1]) * (u[i+1] - u[i]);
    if (i + 1 == ib) Qb = acc;
    if (i + 1 == ic) Qc = acc;
  }
  out[1] = Qc - u[ic]*q[ic]; out[2] = Qb; out[0] = out[2] - out[1];
  free(u); free(q); free(gp); free(gpp);
}
