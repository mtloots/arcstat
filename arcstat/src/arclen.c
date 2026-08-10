/* arclen.c : shared C back-end for the arc-length goodness-of-fit test.
 * Pure libm; no external numerical libraries. Bound identically from R and Python. */
#include "arclen.h"
#include <math.h>
#include <stdlib.h>

static int cmp_double(const void *a, const void *b) {
  double x = *(const double *)a, y = *(const double *)b;
  return (x > y) - (x < y);
}

void al_statistic(const int *n, const double *u, double *out) {
  int m = *n + 1;
  double *v = (double *)malloc((size_t)*n * sizeof(double));
  for (int i = 0; i < *n; i++) v[i] = u[i];
  qsort(v, (size_t)*n, sizeof(double), cmp_double);
  double c = 1.0 / m, prev = 0.0, s = 0.0;
  for (int i = 0; i < *n; i++) { double g = v[i] - prev; s += sqrt(g * g + c * c); prev = v[i]; }
  s += sqrt((1.0 - prev) * (1.0 - prev) + c * c);           /* final spacing to 1 */
  free(v);
  *out = s;
}

/* --- conditional double-saddlepoint machinery --- */
/* The scalar CGF al_kappa that stood here is gone with the finite-difference helpers it
 * served: al_kappa_all computes the value and all five derivatives analytically in one
 * Simpson pass, which is both exact and better conditioned. */


/* Analytic derivatives of the CGF.  With phi(w) = sqrt(w^2+1) and
 *   g(w; a, b) = exp{a (phi(w) - 1) + (b - 1) w},   I(a,b) = int_0^inf g dw,
 * the CGF is kappa = a + log I, and every derivative it needs is the SAME integral with a
 * polynomial weight: I_a = int (phi-1) g, I_b = int w g, I_aa = int (phi-1)^2 g,
 * I_bb = int w^2 g, I_ab = int (phi-1) w g.  Computing the six accumulators in one Simpson
 * pass removes the finite differencing entirely; the second differences it replaces carried
 * a relative noise of order (quadrature error)/h^2, which made the Newton solve, and hence
 * the tail probability, sensitive to the last bits of its input. */
static void al_kappa_all(double a, double b, double *o) {   /* o: kv, ka, kb, kaa, kbb, kab */
  double d = 1.0 - b - a;
  if (d <= 1e-6) { for (int i = 0; i < 6; i++) o[i] = NAN; return; }
  double Wmax = 60.0 / d; if (Wmax > 400.0) Wmax = 400.0;
  int N = 20000;
  double h = Wmax / N;
  double I = 0.0, Ia = 0.0, Ib = 0.0, Iaa = 0.0, Ibb = 0.0, Iab = 0.0;
  for (int i = 0; i <= N; i++) {
    double w = i * h;
    double p1 = sqrt(w * w + 1.0) - 1.0;
    double g = exp(a * p1 + (b - 1.0) * w);
    double wt = (i == 0 || i == N) ? 1.0 : (i % 2 ? 4.0 : 2.0);
    double gw = wt * g;
    I   += gw;
    Ia  += gw * p1;
    Ib  += gw * w;
    Iaa += gw * p1 * p1;
    Ibb += gw * w * w;
    Iab += gw * p1 * w;
  }
  double c = h / 3.0;
  I *= c; Ia *= c; Ib *= c; Iaa *= c; Ibb *= c; Iab *= c;
  if (!(I > 0.0) || isinf(I)) { for (int i = 0; i < 6; i++) o[i] = NAN; return; }
  double ra = Ia / I, rb = Ib / I;
  o[0] = a + log(I);
  o[1] = 1.0 + ra;
  o[2] = rb;
  o[3] = Iaa / I - ra * ra;
  o[4] = Ibb / I - rb * rb;
  o[5] = Iab / I - ra * rb;
}

/* The finite-difference kappa helpers that stood here were superseded by al_kappa_all,
 * which computes kv, ka, kb, kaa, kbb and kab ANALYTICALLY in one Simpson pass; the
 * differenced versions were the source of the ill-conditioning that fix removed. */

void al_pvalue(const double *s, const int *n, double *out) {
  int m = *n + 1; double t = m * (*s);
  /* solve  ka(a,b) = t/m,  kb(a,b) = 1  by 2-D Newton (Jacobian = Hessian of kappa) */
  double a = 0.01, b = 0.0, target_a = t / (double)m;
  double K[6];
  for (int it = 0; it < 60; it++) {
    al_kappa_all(a, b, K);
    if (isnan(K[0])) { *out = NAN; return; }
    double fa = K[1] - target_a, fb = K[2] - 1.0;
    if (fabs(fa) < 1e-12 && fabs(fb) < 1e-12) break;
    double Jaa = K[3], Jbb = K[4], Jab = K[5];                  /* exact symmetric Jacobian */
    double det = Jaa * Jbb - Jab * Jab;
    if (fabs(det) < 1e-14) { *out = NAN; return; }
    double da = (Jbb * fa - Jab * fb) / det;                   /* Newton step */
    double db = (-Jab * fa + Jaa * fb) / det;
    double step = 1.0;                                          /* damp to stay in the CGF domain a<1-b */
    while (a - step * da >= 1.0 - (b - step * db) - 1e-4 && step > 1e-6) step *= 0.5;
    a -= step * da; b -= step * db;
  }
  al_kappa_all(a, b, K);
  double kva = K[0], Kaa = K[3], Kbb = K[4], Kab = K[5];
  if (isnan(kva)) { *out = NAN; return; }
  double w = (a >= 0 ? 1.0 : -1.0) * sqrt(2.0 * (a * t + b * m - m * kva));
  double u = a * sqrt(m * (Kaa * Kbb - Kab * Kab));
  double Phi = 0.5 * erfc(-w / sqrt(2.0));                    /* pnorm(w) */
  double phi = exp(-0.5 * w * w) / sqrt(2.0 * M_PI);
  *out = 1.0 - Phi - phi * (1.0 / w - 1.0 / u);
}

/* --- exact moments via the Dirichlet spacing integrals (composite Simpson) --- */
static double simpson(double (*f)(double, void *), void *p, double lo, double hi, int N) {
  double h = (hi - lo) / N, sum = f(lo, p) + f(hi, p);
  for (int i = 1; i < N; i++) sum += (i % 2 ? 4.0 : 2.0) * f(lo + i * h, p);
  return sum * h / 3.0;
}
typedef struct { int m; double c0; double a; } mpar;
static double marg1(double u, void *pv) { mpar *p = (mpar *)pv; return sqrt(u*u + p->c0*p->c0) * (p->m - 1) * pow(1 - u, p->m - 2); }
static double marg2(double u, void *pv) { mpar *p = (mpar *)pv; return (u*u + p->c0*p->c0) * (p->m - 1) * pow(1 - u, p->m - 2); }
static double inner(double g2, void *pv) { mpar *p = (mpar *)pv; return sqrt(g2*g2 + p->c0*p->c0) * pow(1 - p->a - g2, p->m - 3); }
static double crossf(double g1, void *pv) { mpar *p = (mpar *)pv; p->a = g1;
  return sqrt(g1*g1 + p->c0*p->c0) * (p->m - 1) * (p->m - 2) * simpson(inner, p, 0.0, 1.0 - g1, 2000); }

void al_moments(const int *n, double *out) {
  int m = *n + 1; double c0 = 1.0 / m; mpar p = { m, c0, 0.0 };
  double E1 = simpson(marg1, &p, 0.0, 1.0, 8000);
  double E2 = simpson(marg2, &p, 0.0, 1.0, 8000);
  double Ecr = simpson(crossf, &p, 0.0, 1.0, 2000);
  out[0] = m * E1;                                            /* mean */
  out[1] = m * (E2 - E1 * E1) + (double)m * (m - 1) * (Ecr - E1 * E1);  /* var */
  out[2] = sqrt(2.0);                                         /* support lo */
  out[3] = 2.0 - 1.0 / m;                                     /* support hi */
}
