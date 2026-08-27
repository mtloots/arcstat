/* arcstatAPI.h : the C entry points arcstat makes callable from other packages.
 *
 * arcstat registers its routines with R_useDynamicSymbols(dll, FALSE), which is what CRAN asks
 * for and which also means another package's C cannot reach these symbols by dynamic lookup. The
 * supported route is R_RegisterCCallable on this side and R_GetCCallable on the other; this
 * header wraps that so a dependent package writes an ordinary call.
 *
 * Usage in a dependent package:  LinkingTo: arcstat, then #include <arcstatAPI.h>.
 * Only the quantile-family routines are exported: the four-parameter kappa distribution of
 * Hosking (1994) and the beta-companion (complementary beta) quantile family, each with the
 * pieces needed to fit by L-moments and evaluate. Nothing else in arcstat is part of this API.
 */
#ifndef ARCSTAT_API_H
#define ARCSTAT_API_H
#include <R_ext/Rdynload.h>
#include <Rinternals.h>

/* four-parameter kappa: quantile function, L-moment ratios, and the L-moment fit */
static R_INLINE void arcstat_k4_q(const double *u, const int *nu, const double *mu,
                                  const double *sg, const double *k, const double *h,
                                  double *out) {
  static void (*fn)(const double *, const int *, const double *, const double *,
                    const double *, const double *, double *) = NULL;
  if (fn == NULL) fn = (void (*)(const double *, const int *, const double *, const double *,
                                 const double *, const double *, double *))
                       R_GetCCallable("arcstat", "arck4_q");
  fn(u, nu, mu, sg, k, h, out);
}

static R_INLINE void arcstat_k4_tau34(const double *k, const double *h, const int *nodes,
                                      double *out) {
  static void (*fn)(const double *, const double *, const int *, double *) = NULL;
  if (fn == NULL) fn = (void (*)(const double *, const double *, const int *, double *))
                       R_GetCCallable("arcstat", "arck4_tau34");
  fn(k, h, nodes, out);
}

static R_INLINE void arcstat_k4_fit_lmom(const double *t3, const double *t4, const int *nodes,
                                         double *out) {
  static void (*fn)(const double *, const double *, const int *, double *) = NULL;
  if (fn == NULL) fn = (void (*)(const double *, const double *, const int *, double *))
                       R_GetCCallable("arcstat", "arck4_fit_lmom");
  fn(t3, t4, nodes, out);
}

/* beta-companion (complementary beta) quantile family */
static R_INLINE void arcstat_bc_q(const int *n, const double *u, const double *alpha,
                                  const double *beta, double *out) {
  static void (*fn)(const int *, const double *, const double *, const double *, double *) = NULL;
  if (fn == NULL) fn = (void (*)(const int *, const double *, const double *, const double *,
                                 double *)) R_GetCCallable("arcstat", "arceq2_bc_q");
  fn(n, u, alpha, beta, out);
}

static R_INLINE void arcstat_bc_pdf(const int *n, const double *x, const double *alpha,
                                    const double *beta, double *out) {
  static void (*fn)(const int *, const double *, const double *, const double *, double *) = NULL;
  if (fn == NULL) fn = (void (*)(const int *, const double *, const double *, const double *,
                                 double *)) R_GetCCallable("arcstat", "arceq2_bc_pdf");
  fn(n, x, alpha, beta, out);
}

#endif
