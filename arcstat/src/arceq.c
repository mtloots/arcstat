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

/* Exact integral of q(u) = u^e h(u) over one cell [u0, u1] when h is taken linear between its
   node values h0 and h1 (e > -1). The singular factor u^e is integrated exactly, so the rule is
   exact for a pure power law and its error is only that of a linear interpolant of the smooth
   factor, however steep u^e is on the first cells of the grid. */
static double eq_cell_power(double e, double u0, double u1, double h0, double h1) {
  double e1 = e + 1.0, e2 = e + 2.0;
  double p0 = pow(u0, e1), p1 = pow(u1, e1);
  double base = (p1 - p0) / e1;                                     /* int u^e */
  double mom  = (u1 * p1 - u0 * p0) / e2 - u0 * base;               /* int (u - u0) u^e */
  double s = (h1 - h0) / (u1 - u0);
  return h0 * base + s * mom;
}

/* log q and its first two derivatives for the tilted beta kernel at one point */
static void eq_g(double al, double be, int nt, const double *theta, double ui,
                 double *g, double *g1, double *g2) {
  double gg = al * log(ui) + be * log(1.0 - ui);
  double a1 = al / ui - be / (1.0 - ui);
  double a2 = -al / (ui * ui) - be / ((1.0 - ui) * (1.0 - ui));
  if (nt >= 1) { gg += theta[0] * (2.0*ui - 1.0); a1 += theta[0] * 2.0; }
  if (nt >= 2) { gg += theta[1] * (6.0*ui*ui - 6.0*ui + 1.0);
                 a1 += theta[1] * (12.0*ui - 6.0); a2 += theta[1] * 12.0; }
  if (nt >= 3) { gg += theta[2] * (20.0*ui*ui*ui - 30.0*ui*ui + 12.0*ui - 1.0);
                 a1 += theta[2] * (60.0*ui*ui - 60.0*ui + 12.0); a2 += theta[2] * (120.0*ui - 60.0); }
  *g = gg; *g1 = a1; *g2 = a2;
}

/* Gauss-Legendre, five points, of exp(g) over [lo, hi] (signed); used for the short stretch from
   a grid node to a refined root */
static double eq_quad_q(double al, double be, int nt, const double *theta, double lo, double hi) {
  static const double x[5] = { 0.0, 0.5384693101056831, -0.5384693101056831,
                               0.9061798459386640, -0.9061798459386640 };
  static const double w[5] = { 0.5688888888888889, 0.4786286704993665, 0.4786286704993665,
                               0.2369268850561891, 0.2369268850561891 };
  double c = 0.5 * (lo + hi), h = 0.5 * (hi - lo), s = 0.0, g, g1, g2;
  for (int k = 0; k < 5; k++) { eq_g(al, be, nt, theta, c + h * x[k], &g, &g1, &g2); s += w[k] * exp(g); }
  return s * h;
}

/* bisection for the sign change of f (f = g' when kind == 0, 2 g'^2 - g'' when kind == 1) on
   [lo, hi]; returns NAN when the end values do not bracket a root */
static double eq_root(double al, double be, int nt, const double *theta, int kind,
                      double lo, double hi) {
  double g, g1, g2, flo, fhi, fm;
  eq_g(al, be, nt, theta, lo, &g, &g1, &g2); flo = kind ? 2.0*g1*g1 - g2 : g1;
  eq_g(al, be, nt, theta, hi, &g, &g1, &g2); fhi = kind ? 2.0*g1*g1 - g2 : g1;
  if (!((flo > 0 && fhi < 0) || (flo < 0 && fhi > 0))) return NAN;
  for (int it = 0; it < 80; it++) {
    double m = 0.5 * (lo + hi);
    eq_g(al, be, nt, theta, m, &g, &g1, &g2); fm = kind ? 2.0*g1*g1 - g2 : g1;
    if (fm == 0.0) return m;
    if ((fm > 0) == (flo > 0)) { lo = m; flo = fm; } else { hi = m; }
  }
  return 0.5 * (lo + hi);
}

/* out[5] = D, a, b, u_c, u_b. Q(u) is the integral of q from ZERO, so a and b are quantile-domain
   values on the same footing as the closed forms (incomplete beta) for a beta kernel. The grid
   starts at 1e-5, not at zero; for alpha > -1 the mass on (0, 1e-5) is added and the singular
   power factor is integrated exactly on every cell, which a plain trapezoid cannot do for
   alpha < -1/2. u_c and u_b are the roots of g' = 0 and of 2 g'^2 = g'' located on the grid and
   then refined by bisection; Q is carried from the nearest node to the refined point. */
void arceq_readings(const double *alpha, const double *beta, const int *ntheta,
                    const double *theta, const int *ngrid, double *out) {
  int N = *ngrid, nt = *ntheta;
  double al = *alpha, be = *beta;
  double *u  = (double *)malloc((size_t)N * sizeof(double));
  double *q  = (double *)malloc((size_t)N * sizeof(double));
  double *gp = (double *)malloc((size_t)N * sizeof(double));
  double *gpp= (double *)malloc((size_t)N * sizeof(double));
  for (int i = 0; i < N; i++) {
    double ui = 1e-5 + (1.0 - 2e-5) * (double)i / (N - 1), g, g1, g2;
    u[i] = ui;
    eq_g(al, be, nt, theta, ui, &g, &g1, &g2);
    q[i] = exp(g); gp[i] = g1; gpp[i] = g2;
  }
  out[0] = out[1] = out[2] = out[3] = out[4] = NAN;
  int ic = 0; double qmin = q[0];
  for (int i = 1; i < N; i++) if (q[i] < qmin) { qmin = q[i]; ic = i; }
  if (ic <= 10 || ic >= N - 10) { free(u); free(q); free(gp); free(gpp); return; }
  int ib = -1;
  for (int i = 6; i + 1 < ic; i++) {              /* smallest root of 2 g'^2 = g'' below u_c */
    double r0 = 2.0*gp[i]*gp[i] - gpp[i], r1 = 2.0*gp[i+1]*gp[i+1] - gpp[i+1];
    if ((r0 > 0 && r1 < 0) || (r0 < 0 && r1 > 0)) { ib = i; break; }
  }
  if (ib < 0) { free(u); free(q); free(gp); free(gpp); return; }
  /* Q at the nodes: exact power factor on each cell, plus the mass below the first node */
  int power = (al > -1.0) && (al < 1e6);
  double acc = 0.0, Qb = 0.0, Qc = 0.0;
  if (power && q[0] > 0.0) {
    double h0 = q[0] / pow(u[0], al), h1 = q[1] / pow(u[1], al), s = (h1 - h0) / (u[1] - u[0]);
    double e1 = al + 1.0, e2 = al + 2.0, p0 = pow(u[0], e1);
    acc = h0 * p0 / e1 + s * (u[0] * p0 / e2 - u[0] * p0 / e1);   /* int_0^u0 u^al (h0 + s (u-u0)) */
  }
  for (int i = 0; i + 1 < N && i + 1 <= ic; i++) {
    double seg;
    if (power && q[i] > 0.0 && q[i+1] > 0.0)
      seg = eq_cell_power(al, u[i], u[i+1], q[i] / pow(u[i], al), q[i+1] / pow(u[i+1], al));
    else seg = 0.5 * (q[i] + q[i+1]) * (u[i+1] - u[i]);
    acc += seg;
    if (i + 1 == ib) Qb = acc;
    if (i + 1 == ic) Qc = acc;
  }
  /* refine the two roots between neighbouring nodes and carry Q to them */
  double uc = u[ic], ub = u[ib], Qcr = Qc, Qbr = Qb;
  double rc = eq_root(al, be, nt, theta, 0, u[ic-1], u[ic+1]);
  if (isfinite(rc)) { uc = rc; Qcr = Qc + eq_quad_q(al, be, nt, theta, u[ic], rc); }
  double rb = eq_root(al, be, nt, theta, 1, u[ib], u[ib+1]);
  if (isfinite(rb)) { ub = rb; Qbr = Qb + eq_quad_q(al, be, nt, theta, u[ib], rb); }
  double gc, g1c, g2c;
  eq_g(al, be, nt, theta, uc, &gc, &g1c, &g2c);
  out[3] = uc; out[4] = ub;
  out[1] = Qcr - uc * exp(gc);                     /* tangent reading a */
  out[2] = Qbr;                                    /* derivative reading b */
  out[0] = out[2] - out[1];                        /* D = b - a */
  free(u); free(q); free(gp); free(gpp);
}

/* van Staden-Loots host: q(u) = (1-d) u^(l-1) + d (1-u)^(l-1) (scale-free). Its integral is
   closed, Q(u) = [(1-d) u^l + d (1 - (1-u)^l)] / l, which is used directly: exact at every node and
   at the refined roots, and zero at u = 0 for l > 0. For l <= 0 the integral from zero diverges,
   and Q is measured from the first grid node (u = 1e-5) as before; only the gap D = b - a is then
   meaningful. */
static double vsl_Q(double l, double d, double u) {
  if (fabs(l) < 1e-9) return (1.0 - d) * log(u) - d * log1p(-u);
  return ((1.0 - d) * pow(u, l) - d * expm1(l * log1p(-u))) / l;
}
static void vsl_g(double l, double d, double u, double *gp, double *gpp) {
  double t1 = (1.0-d)*pow(u, l-1.0), t2 = d*pow(1.0-u, l-1.0), qv = t1 + t2;
  double q1 = (1.0-d)*(l-1.0)*pow(u, l-2.0) - d*(l-1.0)*pow(1.0-u, l-2.0);
  double q2 = (1.0-d)*(l-1.0)*(l-2.0)*pow(u, l-3.0) + d*(l-1.0)*(l-2.0)*pow(1.0-u, l-3.0);
  *gp = q1/qv; *gpp = q2/qv - (q1/qv)*(q1/qv);
}
static double vsl_root(double l, double d, int kind, double lo, double hi) {
  double g1, g2, flo, fhi, fm;
  vsl_g(l, d, lo, &g1, &g2); flo = kind ? 2.0*g1*g1 - g2 : g1;
  vsl_g(l, d, hi, &g1, &g2); fhi = kind ? 2.0*g1*g1 - g2 : g1;
  if (!((flo > 0 && fhi < 0) || (flo < 0 && fhi > 0))) return NAN;
  for (int it = 0; it < 80; it++) {
    double m = 0.5 * (lo + hi);
    vsl_g(l, d, m, &g1, &g2); fm = kind ? 2.0*g1*g1 - g2 : g1;
    if (fm == 0.0) return m;
    if ((fm > 0) == (flo > 0)) { lo = m; flo = fm; } else { hi = m; }
  }
  return 0.5 * (lo + hi);
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
    q[i] = (1.0-d)*pow(ui, l-1.0) + d*pow(1.0-ui, l-1.0);
    vsl_g(l, d, ui, &gp[i], &gpp[i]);
  }
  out[0] = out[1] = out[2] = out[3] = out[4] = NAN;
  int ic = 0; double qmin = q[0];
  for (int i = 1; i < N; i++) if (q[i] < qmin) { qmin = q[i]; ic = i; }
  if (ic <= 10 || ic >= N - 10) { free(u); free(q); free(gp); free(gpp); return; }
  int ib = -1;
  for (int i = 0; i + 1 < ic; i++) {
    if (u[i] < 1e-3) continue;
    double r0 = 2.0*gp[i]*gp[i] - gpp[i], r1 = 2.0*gp[i+1]*gp[i+1] - gpp[i+1];
    if ((r0 > 0 && r1 < 0) || (r0 < 0 && r1 > 0)) { ib = i; break; }
  }
  if (ib < 0) { free(u); free(q); free(gp); free(gpp); return; }
  double uc = u[ic], ub = u[ib];
  double rc = vsl_root(l, d, 0, u[ic-1], u[ic+1]);
  if (isfinite(rc)) uc = rc;
  double rb = vsl_root(l, d, 1, u[ib], u[ib+1]);
  if (isfinite(rb)) ub = rb;
  double Q0 = (l > 0.0) ? 0.0 : vsl_Q(l, d, u[0]);
  double Qc = vsl_Q(l, d, uc) - Q0, Qb = vsl_Q(l, d, ub) - Q0;
  double qc = (1.0-d)*pow(uc, l-1.0) + d*pow(1.0-uc, l-1.0);
  out[3] = uc; out[4] = ub;
  out[1] = Qc - uc*qc; out[2] = Qb; out[0] = out[2] - out[1];
  free(u); free(q); free(gp); free(gpp);
}
