/* ellstat --- the elliptic moment system on the circle.
 *
 * Basis of modulus k in [0,1):
 *     c_p(t) = cn(2 K(k) p t / pi, k),   s_p(t) = sn(2 K(k) p t / pi, k),
 * which reduces to (cos pt, sin pt) at k = 0.  The transfer to trigonometric
 * moments is triangular with non-zero diagonal, so the system characterises
 * every circular law at every modulus.
 *
 * Pure libm.  Pointer arguments throughout for the .C ABI, so one source
 * compiles under R and under Python via ctypes.  Every loop lives here.
 */
/* Floating-point contraction is pinned OFF so that the R and Python fronts, which compile this
   source with different flags, cannot differ in whether multiply-add pairs are fused. Without it
   the fronts agree to about fifteen digits and disagree in the last bit, which fails the parity
   harness -- and it surfaces unpredictably, since whether a fusion happens depends on register
   allocation, so an unrelated edit can expose it. arcstat pins it in every source for the same
   reason. */
#if defined(__clang__)
#pragma clang fp contract(off)
#elif defined(__GNUC__)
#pragma GCC optimize ("fp-contract=off")
#else
#pragma STDC FP_CONTRACT OFF
#endif

#include <math.h>
#include <stdlib.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ---------------------------------------------------------------- streams */
static unsigned long long sm64(unsigned long long *s){
  unsigned long long z = (*s += 0x9E3779B97F4A7C15ULL);
  z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
  z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
  return z ^ (z >> 31);
}
static double es_unif(unsigned long long *s){ return (sm64(s) >> 11) * 0x1.0p-53; }
static unsigned long long es_seed(unsigned long long base, long rep){
  unsigned long long s = base ^ (0x9E3779B97F4A7C15ULL * (unsigned long long)(rep + 1));
  (void)sm64(&s); (void)sm64(&s);
  return s;
}

/* ------------------------------------------- complete elliptic integrals  */
/* AGM.  E = K (1 - sum_{n>=0} 2^{n-1} c_n^2).                              */
static void es_ke(double k, double *K, double *E){
  double a = 1.0, b = sqrt(1.0 - k * k), c = k, p = 0.5, s = p * c * c;
  int i;
  for (i = 0; i < 60; i++){
    double an = 0.5 * (a + b), bn = sqrt(a * b), cn = 0.5 * (a - b);
    a = an; b = bn; c = cn; p *= 2.0; s += p * c * c;
    if (fabs(c) < 1e-18) break;
  }
  *K = M_PI / (2.0 * a);
  *E = (*K) * (1.0 - s);
}
void C_ell_ke(double *k, int *n, double *K, double *E){
  int i; for (i = 0; i < *n; i++) es_ke(k[i], &K[i], &E[i]);
}
/* nome q = exp(-pi K'/K) */
static double es_nome(double k){
  double K, E, Kp, Ep;
  es_ke(k, &K, &E); es_ke(sqrt(1.0 - k * k), &Kp, &Ep);
  return exp(-M_PI * Kp / K);
}
void C_ell_nome(double *k, int *n, double *out){
  int i; for (i = 0; i < *n; i++) out[i] = es_nome(k[i]);
}

/* --------------------------------------------- Jacobi sn, cn, dn (Landen) */
static void es_jac(double u, double k, double *sn, double *cn, double *dn){
  double a[64], c[64], b, phi;
  int i, N = 0;
  if (k < 1e-14){ *sn = sin(u); *cn = cos(u); *dn = 1.0; return; }
  b = sqrt(1.0 - k * k); a[0] = 1.0; c[0] = k;
  for (i = 1; i < 64; i++){
    a[i] = 0.5 * (a[i-1] + b);
    c[i] = 0.5 * (a[i-1] - b);
    b = sqrt(a[i-1] * b);
    N = i;
    if (fabs(c[i]) < 1e-17) break;
  }
  phi = ldexp(a[N] * u, N);
  for (i = N; i > 0; i--){
    double x = (c[i] / a[i]) * sin(phi);
    if (x >  1.0) x =  1.0;
    if (x < -1.0) x = -1.0;
    phi = 0.5 * (phi + asin(x));
  }
  *sn = sin(phi); *cn = cos(phi);
  *dn = sqrt(1.0 - k * k * (*sn) * (*sn));
}
void C_ell_jac(double *u, int *n, double *k, double *sn, double *cn, double *dn){
  int i; for (i = 0; i < *n; i++) es_jac(u[i], *k, &sn[i], &cn[i], &dn[i]);
}

/* -------------------------------------------------------- the basis c_p, s_p */
void C_ell_basis(double *t, int *nt, int *P, double *k, double *C, double *S){
  double K, E; int i, p;
  es_ke(*k, &K, &E);
  for (p = 1; p <= *P; p++)
    for (i = 0; i < *nt; i++){
      double sn, cn, dn;
      es_jac(2.0 * K * p * t[i] / M_PI, *k, &sn, &cn, &dn);
      C[(p - 1) * (*nt) + i] = cn;
      S[(p - 1) * (*nt) + i] = sn;
    }
}

/* ------------------------------------------------- sample elliptic moments */
/* out[2(p-1)] = mean c_p, out[2(p-1)+1] = mean s_p                          */
void C_ell_smom(double *t, int *n, int *P, double *k, double *out){
  double K, E; int i, p;
  es_ke(*k, &K, &E);
  for (p = 0; p < 2 * (*P); p++) out[p] = 0.0;
  for (i = 0; i < *n; i++)
    for (p = 1; p <= *P; p++){
      double sn, cn, dn;
      es_jac(2.0 * K * p * t[i] / M_PI, *k, &sn, &cn, &dn);
      out[2 * (p - 1)]     += cn / (double)(*n);
      out[2 * (p - 1) + 1] += sn / (double)(*n);
    }
}

/* --------------------------------------- transfer trig -> elliptic moments */
/* c_p is supported on harmonics (2m+1)p with coefficient                    */
/*     A_{p,m} = (2 pi / (k K)) q^{m+1/2} / (1 + q^{2m+1}).                   */
/* out is P x M, row p holds A_{p,0..M-1}.  Diagonal entry A_{p,0} != 0.      */
void C_ell_transfer(double *k, int *P, int *M, double *out){
  double K, E, q, pref; int p, m;
  es_ke(*k, &K, &E);
  q = es_nome(*k);
  pref = 2.0 * M_PI / ((*k) * K);
  for (p = 1; p <= *P; p++)
    for (m = 0; m < *M; m++)
      out[(p - 1) * (*M) + m] =
        pref * pow(q, m + 0.5) / (1.0 + pow(q, 2 * m + 1));
}

/* ------------------------------------------- the diagonal family: dn^2/(4E) */
/* On its own clock the ellipse arc-length law is dn^2(u)/(4E), with         */
/*     dn^2(u) = E/K + (2 pi^2 / K^2) sum_{n>=1} n q^n/(1 - q^{2n}) cos(n pi u/K). */
/* out[0] = E/K, out[n] = (2 pi^2/K^2) n q^n/(1 - q^{2n}) for n = 1..P.       */
void C_ell_dn2mom(double *k, int *P, double *out){
  double K, E, q; int n;
  es_ke(*k, &K, &E);
  q = es_nome(*k);
  out[0] = E / K;
  for (n = 1; n <= *P; n++)
    out[n] = (2.0 * M_PI * M_PI / (K * K)) * n * pow(q, n) / (1.0 - pow(q, 2 * n));
}

/* ------------------------------------------------- elliptic von Mises law  */
/* f(t) propto exp(kappa * cn(2K(t-mu)/pi, k)); normalised numerically on a   */
/* grid of ng points, which is exact to machine precision because the         */
/* integrand is analytic and periodic (trapezoid converges geometrically).    */
static double es_evm_logk(double kap, double k, int ng){
  double K, E, h = 2.0 * M_PI / ng, Z = 0.0; int i;
  es_ke(k, &K, &E);
  for (i = 0; i < ng; i++){
    double t = h * i, sn, cn, dn;
    es_jac(2.0 * K * t / M_PI, k, &sn, &cn, &dn);
    Z += exp(kap * (cn - 1.0)) * h;
  }
  return log(Z) + kap;   /* log of the true normalising constant */
}
void C_ell_evm_d(double *t, int *n, double *mu, double *kap, double *k,
               int *ng, double *out){
  double K, E, lZ; int i;
  es_ke(*k, &K, &E);
  lZ = es_evm_logk(*kap, *k, *ng);
  for (i = 0; i < *n; i++){
    double a = t[i] - *mu, sn, cn, dn;
    a = fmod(a, 2.0 * M_PI); if (a < 0.0) a += 2.0 * M_PI;
    es_jac(2.0 * K * a / M_PI, *k, &sn, &cn, &dn);
    out[i] = exp((*kap) * cn - lZ);
  }
}
void C_ell_evm_r(int *n, double *mu, double *kap, double *k,
               double *seed, double *out){
  double K, E; unsigned long long st = (unsigned long long)(*seed); int i;
  es_ke(*k, &K, &E);
  for (i = 0; i < *n; i++){
    for (;;){
      double t = 2.0 * M_PI * es_unif(&st), sn, cn, dn;
      es_jac(2.0 * K * t / M_PI, *k, &sn, &cn, &dn);
      if (log(es_unif(&st)) < (*kap) * (cn - 1.0)){
        double v = t + *mu;
        v = fmod(v, 2.0 * M_PI); if (v < 0.0) v += 2.0 * M_PI;
        out[i] = v; break;
      }
    }
  }
}

/* =========================================================================
 * Estimation and exact inference.
 *
 * The identity cn^2 + sn^2 = 1 makes (c_p(t), s_p(t)) a point on the unit
 * circle for every t, so the sample elliptic moment is a mean of unit
 * vectors.  Hence for every n and every k:
 *     E[hat E_p] = E_p            exactly (no bias, no asymptotics), and
 *     tr Var(hat E_p) = (1 - |E_p|^2)/n   exactly.
 * ========================================================================= */

/* exact first- and second-order sample theory, for the verification harness */
void C_ell_moment_var(double *t, int *n, int *P, double *k, double *out){
  /* out[3(p-1)+0] = |hat E_p|, [1] = sample tr-variance, [2] = (1-|E|^2)/n */
  double K, E; int i, p;
  es_ke(*k, &K, &E);
  for (p = 1; p <= *P; p++){
    double a = 0.0, b = 0.0, q2 = 0.0, m;
    for (i = 0; i < *n; i++){
      double sn, cn, dn;
      es_jac(2.0 * K * p * t[i] / M_PI, *k, &sn, &cn, &dn);
      a += cn; b += sn; q2 += cn * cn + sn * sn;   /* == n identically */
    }
    a /= (double)(*n); b /= (double)(*n);
    m = a * a + b * b;
    out[3 * (p - 1)]     = sqrt(m);
    out[3 * (p - 1) + 1] = (q2 / (double)(*n) - m) / (double)(*n);
    out[3 * (p - 1) + 2] = (1.0 - m) / (double)(*n);
  }
}

/* ------------------------------- population elliptic moments of a density */
/* Generic: caller supplies a density on a uniform grid over [0, 2pi).      */
void C_ell_popmom_grid(double *f, int *ng, int *P, double *k, double *out){
  double K, E, h = 2.0 * M_PI / (*ng), Z = 0.0; int i, p;
  es_ke(*k, &K, &E);
  for (i = 0; i < *ng; i++) Z += f[i] * h;
  for (p = 0; p < 2 * (*P); p++) out[p] = 0.0;
  for (i = 0; i < *ng; i++){
    double t = h * (i + 0.5), w = f[i] * h / Z, sn, cn, dn;
    for (p = 1; p <= *P; p++){
      es_jac(2.0 * K * p * t / M_PI, *k, &sn, &cn, &dn);
      out[2 * (p - 1)]     += w * cn;
      out[2 * (p - 1) + 1] += w * sn;
    }
  }
}
/* population elliptic moments of the elliptic von Mises, on its own grid   */
static void es_evm_pop(double mu, double kap, double k, int P, int ng, double *out){
  double K, E, h = 2.0 * M_PI / ng, Z = 0.0; int i, p;
  es_ke(k, &K, &E);
  for (p = 0; p < 2 * P; p++) out[p] = 0.0;
  for (i = 0; i < ng; i++){
    double t = h * (i + 0.5), sn, cn, dn;
    es_jac(2.0 * K * t / M_PI, k, &sn, &cn, &dn);
    Z += exp(kap * (cn - 1.0)) * h;
  }
  for (i = 0; i < ng; i++){
    double t = h * (i + 0.5), sn, cn, dn, w;
    es_jac(2.0 * K * t / M_PI, k, &sn, &cn, &dn);
    w = exp(kap * (cn - 1.0)) * h / Z;
    { double a = t + mu; int pp;
      for (pp = 1; pp <= P; pp++){
        double sn2, cn2, dn2;
        es_jac(2.0 * K * pp * a / M_PI, k, &sn2, &cn2, &dn2);
        out[2 * (pp - 1)]     += w * cn2;
        out[2 * (pp - 1) + 1] += w * sn2;
      } }
  }
}
void C_ell_evm_pop(double *mu, double *kap, double *k, int *P, int *ng, double *out){
  es_evm_pop(*mu, *kap, *k, *P, *ng, out);
}

/* ------------------------------------- minimum-distance fit of (mu, kappa) */
static void es_evm_fit(const double *t, int n, double k, int P, int ng,
                       double *mu_hat, double *kap_hat){
  double emp[16], best = 1e300, bm = 0.0, bk = 1.0;
  int gm, gk, it, p;
  double K, E;
  es_ke(k, &K, &E);
  for (p = 0; p < 2 * P; p++) emp[p] = 0.0;
  { int i;
    for (i = 0; i < n; i++)
      for (p = 1; p <= P; p++){
        double sn, cn, dn;
        es_jac(2.0 * K * p * t[i] / M_PI, k, &sn, &cn, &dn);
        emp[2 * (p - 1)]     += cn / (double)n;
        emp[2 * (p - 1) + 1] += sn / (double)n;
      } }
  for (gm = 0; gm < 36; gm++)
    for (gk = 0; gk < 32; gk++){
      double mu = 2.0 * M_PI * gm / 36.0, kap = exp(-3.0 + 6.0 * gk / 31.0);
      double pm[16], d = 0.0;
      es_evm_pop(mu, kap, k, P, ng, pm);
      for (p = 0; p < 2 * P; p++){ double e = pm[p] - emp[p]; d += e * e; }
      if (d < best){ best = d; bm = mu; bk = kap; }
    }
  for (it = 0; it < 30; it++){
    double step = 0.15 * pow(0.82, it); int dm, dk;
    for (dm = -1; dm <= 1; dm++)
      for (dk = -1; dk <= 1; dk++){
        double mu = bm + step * dm, kap = bk * exp(step * dk);
        double pm[16], d = 0.0;
        if (kap < 1e-4 || kap > 1e4) continue;
        es_evm_pop(mu, kap, k, P, ng, pm);
        for (p = 0; p < 2 * P; p++){ double e = pm[p] - emp[p]; d += e * e; }
        if (d < best){ best = d; bm = mu; bk = kap; }
      }
  }
  *mu_hat = fmod(bm, 2.0 * M_PI); if (*mu_hat < 0.0) *mu_hat += 2.0 * M_PI;
  *kap_hat = bk;
}
void C_ell_evm_fit(double *t, int *n, double *k, int *P, int *ng, double *out){
  es_evm_fit(t, *n, *k, *P, *ng, &out[0], &out[1]);
}

/* ============================== EXACT inference =========================== */
/* Exact test of circular uniformity.
 *
 * The null is SIMPLE and fully specified, so it can be simulated directly:
 * draw B samples of size n from the circular uniform, and rank the observed
 * statistic among them. The rank p-value is exact for EVERY n and EVERY B
 * (Dwass 1957; Hope 1968).
 *
 * An earlier version calibrated by random rigid ROTATIONS of the observed
 * sample. That is also exactly level correct, the uniform law being rotation
 * invariant, but it is nearly powerless: |hat E_1| is a modulus, and the
 * modulus of a first harmonic is exactly rotation invariant, so only the small
 * non-fundamental part of cn -- about 1.3 per cent of its energy at k = 0.9 --
 * drives the statistic at all. Direct simulation of the null avoids this
 * entirely and costs nothing.                                                */
void C_ell_exact_unif(double *t, int *n, double *k, int *B,
                    double *seed, double *out){
  double K, E, obs, a = 0.0, b = 0.0; long ge = 0; int r; int i;
  es_ke(*k, &K, &E);
  for (i = 0; i < *n; i++){
    double sn, cn, dn;
    es_jac(2.0 * K * t[i] / M_PI, *k, &sn, &cn, &dn);
    a += cn; b += sn;
  }
  a /= (double)(*n); b /= (double)(*n);
  obs = sqrt(a * a + b * b);
#ifdef _OPENMP
#pragma omp parallel for reduction(+:ge) schedule(static) num_threads(2)
#endif
  for (r = 0; r < *B; r++){
    unsigned long long st = es_seed((unsigned long long)(*seed), r);
    double aa = 0.0, bb = 0.0; int j;
    for (j = 0; j < *n; j++){
      double v = 2.0 * M_PI * es_unif(&st), sn, cn, dn;
      es_jac(2.0 * K * v / M_PI, *k, &sn, &cn, &dn);
      aa += cn; bb += sn;
    }
    aa /= (double)(*n); bb /= (double)(*n);
    if (sqrt(aa * aa + bb * bb) >= obs) ge++;
  }
  out[0] = obs;
  out[1] = (double)(ge + 1) / (double)(*B + 1);   /* exact p-value */
}

/* Exact confidence interval for kappa by Monte Carlo test inversion.
 *
 * TWO CORRECTIONS over the naive construction, both found by demonstrating the
 * coverage rather than trusting the proof.
 *
 * (1) THE NULL IS COMPOSITE. Fixing kappa does not fix the law: the statistic
 *     |hat E_1| also depends on the mean direction, because it is NOT rotation
 *     invariant when k > 0. Simulating at one arbitrary mean direction therefore
 *     calibrates the wrong null and coverage fails (measured 0.587 at mu = 1
 *     against a nominal 0.90). The repair is the device used for the manifold
 *     test of the equivalence-system manuscript: the p-value is a SUPREMUM over
 *     a set of nuisance values chosen by the analyst. Provided that set contains
 *     the true mean direction, the reported p-value is at least the true one and
 *     the interval covers at no less than its nominal level.
 *
 * (2) THE GRID TRUNCATES INWARD. Reporting the convex hull of the retained grid
 *     points loses up to one grid step at each end, which on a typical interval
 *     was a quarter of its width. The endpoints are now obtained by linear
 *     interpolation of the p-value curve across the level, so the answer no
 *     longer depends on the caller's grid resolution to first order.
 *
 * out[0..1] interval, out[2] statistic, out[3..] the (supremum) p-value curve. */
void C_ell_exact_ci(double *t, int *n, double *k, double *kgrid, int *ngr,
                    int *B, double *level, double *seed, int *ng, int *nmu,
                    double *out){
  double K, E, obs, a = 0.0, b = 0.0; int i, g;
  es_ke(*k, &K, &E);
  for (i = 0; i < *n; i++){
    double sn, cn, dn;
    es_jac(2.0 * K * t[i] / M_PI, *k, &sn, &cn, &dn);
    a += cn; b += sn;
  }
  a /= (double)(*n); b /= (double)(*n);
  obs = sqrt(a * a + b * b);
  out[2] = obs;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(2)
#endif
  for (g = 0; g < *ngr; g++){
    double best = 0.0; int im;
    for (im = 0; im < *nmu; im++){
      double mu = 2.0 * M_PI * im / (double)(*nmu);
      long ge = 0, le = 0; int r; double pv;
      for (r = 0; r < *B; r++){
        unsigned long long st =
          es_seed((unsigned long long)(*seed)
                  ^ (0x5DEECE66DULL * (unsigned long long)(g + 1))
                  ^ (0x9E3779B9ULL * (unsigned long long)(im + 1)), r);
        double aa = 0.0, bb = 0.0, s2; int j;
        for (j = 0; j < *n; j++){
          double v, sn, cn, dn, w;
          for (;;){
            v = 2.0 * M_PI * es_unif(&st);
            es_jac(2.0 * K * v / M_PI, *k, &sn, &cn, &dn);
            if (log(es_unif(&st)) < kgrid[g] * (cn - 1.0)) break;
          }
          /* shift to the nuisance mean direction, then evaluate the basis there */
          w = v + mu;
          w = fmod(w, 2.0 * M_PI); if (w < 0.0) w += 2.0 * M_PI;
          es_jac(2.0 * K * w / M_PI, *k, &sn, &cn, &dn);
          aa += cn; bb += sn;
        }
        aa /= (double)(*n); bb /= (double)(*n);
        s2 = sqrt(aa * aa + bb * bb);
        if (s2 >= obs) ge++;
        if (s2 <= obs) le++;
      }
      pv = 2.0 * (double)((ge < le ? ge : le) + 1) / (double)(*B + 1);
      if (pv > 1.0) pv = 1.0;
      if (pv > best) best = pv;          /* supremum over the nuisance set */
    }
    out[3 + g] = best;
  }
  /* endpoints by linear interpolation of the p-value curve across the level */
  { int lo = -1, hi = -1;
    for (g = 0; g < *ngr; g++) if (out[3 + g] > *level){ if (lo < 0) lo = g; hi = g; }
    if (lo < 0){ out[0] = 1e300; out[1] = -1e300; return; }
    if (lo == 0) out[0] = kgrid[0];
    else {
      double p0 = out[3 + lo - 1], p1 = out[3 + lo];
      double f = (p1 > p0) ? (*level - p0) / (p1 - p0) : 0.0;
      if (f < 0.0) f = 0.0; if (f > 1.0) f = 1.0;
      out[0] = kgrid[lo - 1] + f * (kgrid[lo] - kgrid[lo - 1]);
    }
    if (hi == *ngr - 1) out[1] = kgrid[*ngr - 1];
    else {
      double p0 = out[3 + hi], p1 = out[3 + hi + 1];
      double f = (p0 > p1) ? (p0 - *level) / (p0 - p1) : 0.0;
      if (f < 0.0) f = 0.0; if (f > 1.0) f = 1.0;
      out[1] = kgrid[hi] + f * (kgrid[hi + 1] - kgrid[hi]);
    }
  }
}

void C_ell_effstudy(double *k, double *kap0, int *n, int *R, int *ng,
                    double *seed, double *out){
  double K, E, h, a1 = 0.0, tot = 0.0;
  double *gk, *gt, *ge;
  int G = 240, g, i;
  int r;
  es_ke(*k, &K, &E);
  /* eta(k): energy of c_1 in its fundamental */
  h = 2.0 * M_PI / (*ng);
  for (i = 0; i < *ng; i++){
    double t = h * (i + 0.5), sn, cn, dn;
    es_jac(2.0 * K * t / M_PI, *k, &sn, &cn, &dn);
    a1 += cn * cos(t) * h / M_PI;
    tot += cn * cn * h / M_PI;
  }
  out[0] = a1 * a1 / tot;
  /* monotone tables kappa -> E[cn], kappa -> E[cos] */
  gk = (double*)malloc(sizeof(double) * (G + 1));
  gt = (double*)malloc(sizeof(double) * (G + 1));
  ge = (double*)malloc(sizeof(double) * (G + 1));
  for (g = 0; g <= G; g++){
    double kap = exp(-3.0 + 6.0 * g / (double)G), Z = 0.0, sc = 0.0, se = 0.0;
    for (i = 0; i < *ng; i++){
      double t = h * (i + 0.5), sn, cn, dn, w;
      es_jac(2.0 * K * t / M_PI, *k, &sn, &cn, &dn);
      w = exp(kap * (cn - 1.0)) * h;
      Z += w; se += w * cn; sc += w * cos(t);
    }
    gk[g] = kap; ge[g] = se / Z; gt[g] = sc / Z;
  }
  { double s1[2] = {0,0}, s2[2] = {0,0}; long u[2] = {0,0};
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(2)
#endif
    for (r = 0; r < *R; r++){
      unsigned long long st = es_seed((unsigned long long)(*seed), r);
      double mc = 0.0, me = 0.0; int j, w;
      for (j = 0; j < *n; j++){
        double t, sn, cn, dn;
        for (;;){
          t = 2.0 * M_PI * es_unif(&st);
          es_jac(2.0 * K * t / M_PI, *k, &sn, &cn, &dn);
          if (log(es_unif(&st)) < (*kap0) * (cn - 1.0)) break;
        }
        mc += cos(t) / (double)(*n); me += cn / (double)(*n);
      }
      for (w = 0; w < 2; w++){
        const double *tab = (w == 0) ? gt : ge;
        double m = (w == 0) ? mc : me, f, kh;
        int lo = 0, hi = G;
        if (m <= tab[0] || m >= tab[G]) continue;
        while (hi - lo > 1){ int mid = (lo + hi) / 2;
          if (tab[mid] < m) lo = mid; else hi = mid; }
        f = (m - tab[lo]) / (tab[hi] - tab[lo]);
        kh = gk[lo] * pow(gk[hi] / gk[lo], f);
#ifdef _OPENMP
#pragma omp critical
#endif
        { s1[w] += kh; s2[w] += kh * kh; u[w]++; }
      }
    }
    out[1] = s2[0]/u[0] - (s1[0]/u[0])*(s1[0]/u[0]);
    out[2] = s2[1]/u[1] - (s1[1]/u[1])*(s1[1]/u[1]);
    /* fraction of replicates each estimator could solve: if these differ the
       comparison is not paired and the variance ratio is not interpretable */
    out[3] = (double)u[0] / (double)(*R);
    out[4] = (double)u[1] / (double)(*R);
  }
  free(gk); free(gt); free(ge);
}

/* =========================================================================
 * The projected family, and the modulus characterisation.
 *
 * THEOREM. For a one-parameter family with score s = d/dtheta log f, the moment
 * estimator built on a statistic g has asymptotic efficiency exactly
 * corr^2(g, s) in L^2(f), because d/dtheta E[g] = Cov(g, s). The optimal
 * modulus is therefore k* = argmax_k corr^2(c_1(.,k), s), computable by
 * quadrature with no simulation.
 *
 * The demonstration family: Y | U uniform on the ellipse with semi-axes
 * U(1,tau) offset by d along the axis, U ~ GG(1,beta+1,beta), theta = arg(Y).
 * Projecting a uniform law is geometry, so the conditional angular density is
 * elementary and the mixture is one quadrature.
 * ========================================================================= */

/* conditional angular density of the offset ellipse, given the scale U */
static double es_pj_cond(double t, double U, double tau, double d){
  double a = U, b = U * tau, ct = cos(t), st = sin(t);
  double A = ct * ct / (a * a) + st * st / (b * b);
  double B = -2.0 * ct * d / (a * a);
  double C = d * d / (a * a) - 1.0;
  double disc = B * B - 4.0 * A * C, sq, r1, r2;
  if (disc <= 0.0) return 0.0;
  sq = sqrt(disc);
  r1 = (-B - sq) / (2.0 * A); r2 = (-B + sq) / (2.0 * A);
  if (r2 <= 0.0) return 0.0;
  if (r1 < 0.0) r1 = 0.0;
  return (r2 * r2 - r1 * r1) / (2.0 * M_PI * a * b);
}
/* mixture over the generalised-gamma scale, midpoint rule on U = s/(1-s) */
static double es_pj_dens(double t, double tau, double d, double be, int ng){
  double lgn = lgamma(1.0 + 1.0 / be), acc = 0.0; int j;
  for (j = 0; j < ng; j++){
    double s = (j + 0.5) / ng, U = s / (1.0 - s);
    double jac = 1.0 / ((1.0 - s) * (1.0 - s));
    double lg = log(be) + be * log(U) - pow(U, be) - lgn;
    acc += (1.0 / ng) * jac * exp(lg) * es_pj_cond(t, U, tau, d);
  }
  return acc;
}
void C_ell_pj_dens(double *t, int *n, double *tau, double *d, double *be,
                   int *ng, double *out){
  int i; for (i = 0; i < *n; i++) out[i] = es_pj_dens(t[i], *tau, *d, *be, *ng);
}
void C_ell_pj_rand(int *n, double *tau, double *d, double *be,
                   double *seed, double *out){
  unsigned long long st = (unsigned long long)(*seed); int i;
  double sh = 1.0 + 1.0 / (*be);
  for (i = 0; i < *n; i++){
    double g, U, zx, zy, q, yx, yy, v;
    /* generalised gamma scale via a gamma variate (Marsaglia-Tsang) */
    { double shape = (sh >= 1.0) ? sh : sh + 1.0;
      double dd = shape - 1.0 / 3.0, c = 1.0 / sqrt(9.0 * dd);
      for (;;){
        double x, w, u1, u2, z;
        do { u1 = es_unif(&st); u2 = es_unif(&st);
             if (u1 < 1e-300) u1 = 1e-300;
             z = sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
             w = 1.0 + c * z; } while (w <= 0.0);
        w = w * w * w; x = z;
        { double u = es_unif(&st);
          if (u < 1.0 - 0.0331 * x * x * x * x) { g = dd * w; break; }
          if (log(u) < 0.5 * x * x + dd * (1.0 - w + log(w))) { g = dd * w; break; } }
      }
      if (sh < 1.0) g *= pow(es_unif(&st), 1.0 / sh);
      U = pow(g, 1.0 / (*be));
    }
    do { zx = 2.0 * es_unif(&st) - 1.0; zy = 2.0 * es_unif(&st) - 1.0;
         q = zx * zx + zy * zy; } while (q > 1.0);
    yx = U * zx + *d; yy = U * (*tau) * zy;
    v = atan2(yy, yx); if (v < 0.0) v += 2.0 * M_PI;
    out[i] = v;
  }
}

/* efficiency curve: out[g] = corr^2(c_1(.,kgrid[g]), score) at the given
 * parameter value, by quadrature. out[ngr] holds the Fisher information and
 * out[ngr+1] the Godambe residual |dE[g]/dd - Cov(g,s)| at the first grid
 * point, which is the harness's check on the theorem itself.               */
void C_ell_pj_effcurve(double *tau, double *d, double *be, double *kgrid,
                       int *ngr, int *nt, int *ng, double *eps, double *out){
  int NT = *nt, i, g;
  double h = 2.0 * M_PI / NT, I = 0.0;
  double *TG = (double*)malloc(sizeof(double) * NT);
  double *FD = (double*)malloc(sizeof(double) * NT);
  double *SC = (double*)malloc(sizeof(double) * NT);
  for (i = 0; i < NT; i++){
    double fp, fm;
    TG[i] = h * (i + 0.5);
    FD[i] = es_pj_dens(TG[i], *tau, *d, *be, *ng);
    fp = es_pj_dens(TG[i], *tau, *d + *eps, *be, *ng);
    fm = es_pj_dens(TG[i], *tau, *d - *eps, *be, *ng);
    SC[i] = (log(fp > 1e-300 ? fp : 1e-300) - log(fm > 1e-300 ? fm : 1e-300))
            / (2.0 * (*eps));
  }
  for (i = 0; i < NT; i++) I += SC[i] * SC[i] * FD[i] * h;
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(2)
#endif
  for (g = 0; g < *ngr; g++){
    double K, E, mg = 0.0, cov = 0.0, vg = 0.0; int j;
    double *gv = (double*)malloc(sizeof(double) * NT);
    es_ke(kgrid[g], &K, &E);
    for (j = 0; j < NT; j++){
      double sn, cn, dn;
      es_jac(2.0 * K * TG[j] / M_PI, kgrid[g], &sn, &cn, &dn);
      gv[j] = cn; mg += cn * FD[j] * h;
    }
    for (j = 0; j < NT; j++){
      cov += (gv[j] - mg) * SC[j] * FD[j] * h;
      vg  += (gv[j] - mg) * (gv[j] - mg) * FD[j] * h;
    }
    out[g] = cov * cov / (vg * I);
    free(gv);
  }
  out[*ngr] = I;
  /* Godambe check at the first grid point: dE[g]/dd against Cov(g, score) */
  { double K, E, cov = 0.0, Ep = 0.0, Em = 0.0, mg = 0.0; int j;
    es_ke(kgrid[0], &K, &E);
    for (j = 0; j < NT; j++){
      double sn, cn, dn;
      es_jac(2.0 * K * TG[j] / M_PI, kgrid[0], &sn, &cn, &dn);
      mg += cn * FD[j] * h; cov += cn * SC[j] * FD[j] * h;
      Ep += cn * es_pj_dens(TG[j], *tau, *d + *eps, *be, *ng) * h;
      Em += cn * es_pj_dens(TG[j], *tau, *d - *eps, *be, *ng) * h;
    }
    out[*ngr + 1] = fabs((Ep - Em) / (2.0 * (*eps)) - cov);
  }
  free(TG); free(FD); free(SC);
}

/* the optimal modulus, by golden section on the efficiency curve */
void C_ell_pj_kopt(double *tau, double *d, double *be, int *nt, int *ng,
                   double *eps, double *out){
  double a = 0.0, c = 0.9999999, gr = 0.6180339887498949;
  double x1 = c - gr * (c - a), x2 = a + gr * (c - a), f1, f2, tmp[3];
  int one = 1, it;
  C_ell_pj_effcurve(tau, d, be, &x1, &one, nt, ng, eps, tmp); f1 = tmp[0];
  C_ell_pj_effcurve(tau, d, be, &x2, &one, nt, ng, eps, tmp); f2 = tmp[0];
  for (it = 0; it < 60; it++){
    if (f1 < f2){ a = x1; x1 = x2; f1 = f2; x2 = a + gr * (c - a);
      C_ell_pj_effcurve(tau, d, be, &x2, &one, nt, ng, eps, tmp); f2 = tmp[0]; }
    else { c = x2; x2 = x1; f2 = f1; x1 = c - gr * (c - a);
      C_ell_pj_effcurve(tau, d, be, &x1, &one, nt, ng, eps, tmp); f1 = tmp[0]; }
    if (c - a < 1e-9) break;
  }
  out[0] = 0.5 * (a + c);
  C_ell_pj_effcurve(tau, d, be, &out[0], &one, nt, ng, eps, tmp);
  out[1] = tmp[0];
}

/* =========================================================================
 * The cited projected families: projected normal and projected Cauchy.
 *
 * For Y elliptical about mu with scatter Sigma, theta = arg(Y), the angular
 * density is the ray integral int_0^inf r f_Y(r u(theta)) dr. Writing
 *   A = u' Sigma^{-1} u,  b = u' Sigma^{-1} mu,  a = mu' Sigma^{-1} mu,
 * both integrals are elementary.
 *
 * PROJECTED NORMAL (Mardia's offset normal; Presnell et al. 1998; Wang and
 * Gelfand 2013):
 *   f = (2 pi |Sigma|^{1/2})^{-1} [ A^{-1} e^{-a/2}
 *       + b A^{-3/2} sqrt(2 pi) e^{(b^2/A - a)/2} Phi(b/sqrt A) ].
 *
 * PROJECTED CAUCHY (Tsagris and Alenazi 2025), bivariate t on one degree of
 * freedom, f_Y propto (1 + (y-mu)'Sigma^{-1}(y-mu))^{-3/2}:
 *   f = (2 pi |Sigma|^{1/2})^{-1} A^{-3/2}
 *       [ (m^2+s^2)^{-1/2} + (m/s^2)(1 + m (m^2+s^2)^{-1/2}) ],
 *   m = b/A,  s^2 = (1+a)/A - b^2/A^2,  m^2+s^2 = (1+a)/A.
 * No special function is needed for the Cauchy case at all.
 * ========================================================================= */

static double es_Phi(double x){ return 0.5 * erfc(-x * M_SQRT1_2); }

/* scatter diag(1, tau^2), mean (d, 0) */
static void es_proj_abA(double t, double tau, double d,
                        double *A, double *b, double *a){
  double ct = cos(t), st = sin(t), it2 = 1.0 / (tau * tau);
  *A = ct * ct + st * st * it2;
  *b = d * ct;
  *a = d * d;
}
static double es_pnorm_dens(double t, double tau, double d){
  double A, b, a, r1, r2;
  es_proj_abA(t, tau, d, &A, &b, &a);
  r1 = exp(-a / 2.0) / A;
  r2 = b * pow(A, -1.5) * sqrt(2.0 * M_PI)
       * exp((b * b / A - a) / 2.0) * es_Phi(b / sqrt(A));
  return (r1 + r2) / (2.0 * M_PI * tau);
}
static double es_pcauchy_dens(double t, double tau, double d){
  double A, b, a, m, s2, rt;
  es_proj_abA(t, tau, d, &A, &b, &a);
  m = b / A; s2 = (1.0 + a) / A - b * b / (A * A);
  rt = sqrt(A / (1.0 + a));                     /* = 1/sqrt(m^2+s^2) */
  return pow(A, -1.5) * (rt + (m / s2) * (1.0 + m * rt)) / (2.0 * M_PI * tau);
}
void C_ell_projfam_dens(double *t, int *n, double *tau, double *d,
                        int *which, double *out){
  int i;
  for (i = 0; i < *n; i++)
    out[i] = (*which == 0) ? es_pnorm_dens(t[i], *tau, *d)
                           : es_pcauchy_dens(t[i], *tau, *d);
}
/* simulate: normal or bivariate Cauchy (normal scaled by a chi-1 ratio) */
void C_ell_projfam_rand(int *n, double *tau, double *d, int *which,
                        double *seed, double *out){
  unsigned long long st = (unsigned long long)(*seed); int i;
  for (i = 0; i < *n; i++){
    double u1 = es_unif(&st), u2 = es_unif(&st), z1, z2, yx, yy, v, sc = 1.0;
    if (u1 < 1e-300) u1 = 1e-300;
    z1 = sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
    z2 = sqrt(-2.0 * log(u1)) * sin(2.0 * M_PI * u2);
    if (*which == 1){                     /* bivariate t on 1 df */
      double w1 = es_unif(&st), w2 = es_unif(&st), g;
      if (w1 < 1e-300) w1 = 1e-300;
      g = sqrt(-2.0 * log(w1)) * cos(2.0 * M_PI * w2);
      sc = 1.0 / sqrt(g * g);             /* chi^2_1 */
    }
    yx = sc * z1 + *d; yy = sc * (*tau) * z2;
    v = atan2(yy, yx); if (v < 0.0) v += 2.0 * M_PI;
    out[i] = v;
  }
}
/* efficiency curve corr^2(c_1(.,k), score) for the cited families */
void C_ell_projfam_eff(double *tau, double *d, int *which, double *kgrid,
                       int *ngr, int *nt, double *eps, double *out){
  int NT = *nt, i, g;
  double h = 2.0 * M_PI / NT, I = 0.0;
  double *TG = (double*)malloc(sizeof(double) * NT);
  double *FD = (double*)malloc(sizeof(double) * NT);
  double *SC = (double*)malloc(sizeof(double) * NT);
  for (i = 0; i < NT; i++){
    double fp, fm;
    TG[i] = h * (i + 0.5);
    FD[i] = (*which == 0) ? es_pnorm_dens(TG[i], *tau, *d)
                          : es_pcauchy_dens(TG[i], *tau, *d);
    fp = (*which == 0) ? es_pnorm_dens(TG[i], *tau, *d + *eps)
                       : es_pcauchy_dens(TG[i], *tau, *d + *eps);
    fm = (*which == 0) ? es_pnorm_dens(TG[i], *tau, *d - *eps)
                       : es_pcauchy_dens(TG[i], *tau, *d - *eps);
    SC[i] = (log(fp > 1e-300 ? fp : 1e-300) - log(fm > 1e-300 ? fm : 1e-300))
            / (2.0 * (*eps));
  }
  for (i = 0; i < NT; i++) I += SC[i] * SC[i] * FD[i] * h;
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(2)
#endif
  for (g = 0; g < *ngr; g++){
    double K, E, mg = 0.0, cov = 0.0, vg = 0.0; int j;
    double *gv = (double*)malloc(sizeof(double) * NT);
    es_ke(kgrid[g], &K, &E);
    for (j = 0; j < NT; j++){
      double sn, cn, dn;
      es_jac(2.0 * K * TG[j] / M_PI, kgrid[g], &sn, &cn, &dn);
      gv[j] = cn; mg += cn * FD[j] * h;
    }
    for (j = 0; j < NT; j++){
      cov += (gv[j] - mg) * SC[j] * FD[j] * h;
      vg  += (gv[j] - mg) * (gv[j] - mg) * FD[j] * h;
    }
    out[g] = cov * cov / (vg * I);
    free(gv);
  }
  { double Z = 0.0; for (i = 0; i < NT; i++) Z += FD[i] * h;
    out[*ngr] = I; out[*ngr + 1] = Z; }      /* Fisher info, and the mass check */
  free(TG); free(FD); free(SC);
}
void C_ell_projfam_kopt(double *tau, double *d, int *which, int *nt,
                        double *eps, double *out){
  double a = 0.0, c = 0.9999999, gr = 0.6180339887498949;
  double x1 = c - gr * (c - a), x2 = a + gr * (c - a), f1, f2, tmp[3];
  int one = 1, it;
  C_ell_projfam_eff(tau, d, which, &x1, &one, nt, eps, tmp); f1 = tmp[0];
  C_ell_projfam_eff(tau, d, which, &x2, &one, nt, eps, tmp); f2 = tmp[0];
  for (it = 0; it < 60; it++){
    if (f1 < f2){ a = x1; x1 = x2; f1 = f2; x2 = a + gr * (c - a);
      C_ell_projfam_eff(tau, d, which, &x2, &one, nt, eps, tmp); f2 = tmp[0]; }
    else { c = x2; x2 = x1; f2 = f1; x1 = c - gr * (c - a);
      C_ell_projfam_eff(tau, d, which, &x1, &one, nt, eps, tmp); f1 = tmp[0]; }
    if (c - a < 1e-9) break;
  }
  out[0] = 0.5 * (a + c);
  C_ell_projfam_eff(tau, d, which, &out[0], &one, nt, eps, tmp);
  out[1] = tmp[0];
}

/* =========================================================================
 * A catalogue of circular distributions, after Hosking's table of L-moments.
 * For each family the concentration parameter is the estimand; the routine
 * returns the efficiency of the first elliptic moment at each modulus, and
 * the optimal modulus, all by quadrature.
 *
 *   0 von Mises            exp(p1 cos t)
 *   1 wrapped normal       (1/2pi)(1 + 2 sum p1^{j^2} cos jt)
 *   2 wrapped Cauchy       (1/2pi)(1-p1^2)/(1+p1^2-2 p1 cos t)
 *   3 cardioid             (1/2pi)(1 + 2 p1 cos t),  |p1| <= 1/2
 *   4 projected normal     offset p1, scatter diag(1, p2^2)
 *   5 projected Cauchy     offset p1, scatter diag(1, p2^2)
 *   6 elliptic von Mises   exp(p1 cn(2Kt/pi, p2))
 *   7 arc-length ellipse   sqrt(1 - p1^2 cos^2 t)/(4E(p1))
 * ========================================================================= */
static double es_cat_dens(int fam, double t, double p1, double p2){
  switch (fam){
    case 0: return exp(p1 * cos(t));                    /* normalised later */
    case 1: { double s = 1.0; int j;
              for (j = 1; j <= 60; j++) s += 2.0 * pow(p1, (double)(j * j)) * cos(j * t);
              return (s > 1e-12 ? s : 1e-12) / (2.0 * M_PI); }
    case 2: return (1.0 - p1 * p1)
                   / (2.0 * M_PI * (1.0 + p1 * p1 - 2.0 * p1 * cos(t)));
    case 3: return (1.0 + 2.0 * p1 * cos(t)) / (2.0 * M_PI);
    case 4: return es_pnorm_dens(t, p2, p1);
    case 5: return es_pcauchy_dens(t, p2, p1);
    case 6: { double K, E, sn, cn, dn; es_ke(p2, &K, &E);
              es_jac(2.0 * K * t / M_PI, p2, &sn, &cn, &dn);
              return exp(p1 * cn); }                    /* normalised later */
    case 7: { double K, E; es_ke(p1, &K, &E);
              return sqrt(1.0 - p1 * p1 * cos(t) * cos(t)) / (4.0 * E); }
    default: return 1.0 / (2.0 * M_PI);
  }
}
/* efficiency curve for a catalogue family; out[ngr]=Fisher, out[ngr+1]=mass */
void C_ell_cat_eff(int *fam, double *p1, double *p2, int *ord, double *kgrid,
                   int *ngr, int *nt, double *eps, double *out){
  int NT = *nt, i, g;
  double h = 2.0 * M_PI / NT, I = 0.0, Z = 0.0, Zp = 0.0, Zm = 0.0;
  double *TG = (double*)malloc(sizeof(double) * NT);
  double *FD = (double*)malloc(sizeof(double) * NT);
  double *SC = (double*)malloc(sizeof(double) * NT);
  for (i = 0; i < NT; i++){
    TG[i] = h * (i + 0.5);
    FD[i] = es_cat_dens(*fam, TG[i], *p1, *p2);
    Z  += FD[i] * h;
    Zp += es_cat_dens(*fam, TG[i], *p1 + *eps, *p2) * h;
    Zm += es_cat_dens(*fam, TG[i], *p1 - *eps, *p2) * h;
  }
  for (i = 0; i < NT; i++){
    double fp = es_cat_dens(*fam, TG[i], *p1 + *eps, *p2) / Zp;
    double fm = es_cat_dens(*fam, TG[i], *p1 - *eps, *p2) / Zm;
    FD[i] /= Z;                                    /* normalise once, here */
    SC[i] = (log(fp > 1e-300 ? fp : 1e-300) - log(fm > 1e-300 ? fm : 1e-300))
            / (2.0 * (*eps));
  }
  for (i = 0; i < NT; i++) I += SC[i] * SC[i] * FD[i] * h;
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(2)
#endif
  for (g = 0; g < *ngr; g++){
    double K, E, mg = 0.0, cov = 0.0, vg = 0.0; int j;
    double *gv = (double*)malloc(sizeof(double) * NT);
    es_ke(kgrid[g], &K, &E);
    for (j = 0; j < NT; j++){
      double sn, cn, dn;
      es_jac(2.0 * K * (*ord) * TG[j] / M_PI, kgrid[g], &sn, &cn, &dn);
      gv[j] = cn; mg += cn * FD[j] * h;
    }
    for (j = 0; j < NT; j++){
      cov += (gv[j] - mg) * SC[j] * FD[j] * h;
      vg  += (gv[j] - mg) * (gv[j] - mg) * FD[j] * h;
    }
    out[g] = (vg > 0.0 && I > 0.0) ? cov * cov / (vg * I) : 0.0;
    free(gv);
  }
  out[*ngr] = I; out[*ngr + 1] = Z;
  free(TG); free(FD); free(SC);
}
void C_ell_cat_kopt(int *fam, double *p1, double *p2, int *ord, int *nt,
                    double *eps, double *out){
  double a = 0.0, c = 0.9999999, gr = 0.6180339887498949;
  double x1 = c - gr * (c - a), x2 = a + gr * (c - a), f1, f2, tmp[3];
  int one = 1, it;
  C_ell_cat_eff(fam, p1, p2, ord, &x1, &one, nt, eps, tmp); f1 = tmp[0];
  C_ell_cat_eff(fam, p1, p2, ord, &x2, &one, nt, eps, tmp); f2 = tmp[0];
  for (it = 0; it < 60; it++){
    if (f1 < f2){ a = x1; x1 = x2; f1 = f2; x2 = a + gr * (c - a);
      C_ell_cat_eff(fam, p1, p2, ord, &x2, &one, nt, eps, tmp); f2 = tmp[0]; }
    else { c = x2; x2 = x1; f2 = f1; x1 = c - gr * (c - a);
      C_ell_cat_eff(fam, p1, p2, ord, &x1, &one, nt, eps, tmp); f1 = tmp[0]; }
    if (c - a < 1e-9) break;
  }
  out[0] = 0.5 * (a + c);
  C_ell_cat_eff(fam, p1, p2, ord, &out[0], &one, nt, eps, tmp);
  out[1] = tmp[0];
}

/* =========================================================================
 * The joint law of the sample elliptic moments (cf. Hosking's joint
 * distribution of sample L-moments, which is asymptotic; this one is exact).
 *
 * Stack the basis as g = (c_1, s_1, ..., c_P, s_P). Then hat E = mean(g) is a
 * sample mean of a BOUNDED vector, so for every n
 *      E[hat E] = E[g]        and      Cov(hat E) = Sigma / n,
 * with Sigma = Cov(g) exactly. Sigma is returned by quadrature against a
 * supplied density on a uniform grid.
 * ========================================================================= */
void C_ell_joint_cov(double *f, int *ng, int *P, double *k,
                     double *mean, double *Sigma){
  int NG = *ng, D = 2 * (*P), i, p, r, c;
  double K, E, h = 2.0 * M_PI / NG, Z = 0.0;
  double *g = (double*)malloc(sizeof(double) * D);
  es_ke(*k, &K, &E);
  for (i = 0; i < NG; i++) Z += f[i] * h;
  for (r = 0; r < D; r++){ mean[r] = 0.0;
    for (c = 0; c < D; c++) Sigma[r * D + c] = 0.0; }
  for (i = 0; i < NG; i++){
    double t = h * (i + 0.5), w = f[i] * h / Z;
    for (p = 1; p <= *P; p++){
      double sn, cn, dn;
      es_jac(2.0 * K * p * t / M_PI, *k, &sn, &cn, &dn);
      g[2 * (p - 1)] = cn; g[2 * (p - 1) + 1] = sn;
    }
    for (r = 0; r < D; r++){
      mean[r] += w * g[r];
      for (c = r; c < D; c++) Sigma[r * D + c] += w * g[r] * g[c];
    }
  }
  for (r = 0; r < D; r++)
    for (c = r; c < D; c++){
      Sigma[r * D + c] -= mean[r] * mean[c];
      Sigma[c * D + r] = Sigma[r * D + c];   /* symmetric to the last bit */
    }
  free(g);
}
/* realised covariance of the basis on a sample, for the verification harness */
void C_ell_joint_samp(double *t, int *n, int *P, double *k,
                      double *mean, double *Sigma){
  int D = 2 * (*P), i, p, r, c;
  double K, E;
  double *g = (double*)malloc(sizeof(double) * D);
  es_ke(*k, &K, &E);
  for (r = 0; r < D; r++){ mean[r] = 0.0;
    for (c = 0; c < D; c++) Sigma[r * D + c] = 0.0; }
  for (i = 0; i < *n; i++){
    for (p = 1; p <= *P; p++){
      double sn, cn, dn;
      es_jac(2.0 * K * p * t[i] / M_PI, *k, &sn, &cn, &dn);
      g[2 * (p - 1)] = cn; g[2 * (p - 1) + 1] = sn;
    }
    for (r = 0; r < D; r++){
      mean[r] += g[r] / (double)(*n);
      for (c = r; c < D; c++)
        Sigma[r * D + c] += g[r] * g[c] / (double)(*n);
    }
  }
  for (r = 0; r < D; r++)
    for (c = r; c < D; c++){
      Sigma[r * D + c] -= mean[r] * mean[c];
      Sigma[c * D + r] = Sigma[r * D + c];
    }
  free(g);
}

/* =========================================================================
 * Simultaneous exact test of uniformity from the JOINT moment vector.
 *
 * Under uniformity the joint covariance Sigma_0 of (c_1,s_1,...,c_P,s_P) is
 * fully determined, so the quadratic form
 *      T = n * hat E' Sigma_0^{-1} hat E
 * has a null distribution free of unknowns and can be simulated exactly. The
 * null being simple, the rank p-value is exact for every n and every B.
 * This is the joint counterpart of the single-moment test, and it uses the
 * covariance of Theorem (joint) rather than only its diagonal.
 * ========================================================================= */
static void es_solve_sym(double *A, int D, double *b, double *x){
  /* Gauss-Jordan with partial pivoting; D is tiny */
  int i, j, p; double *M = (double*)malloc(sizeof(double) * D * (D + 1));
  for (i = 0; i < D; i++){
    for (j = 0; j < D; j++) M[i * (D + 1) + j] = A[i * D + j];
    M[i * (D + 1) + D] = b[i];
  }
  for (p = 0; p < D; p++){
    int piv = p; double mx = fabs(M[p * (D + 1) + p]);
    for (i = p + 1; i < D; i++)
      if (fabs(M[i * (D + 1) + p]) > mx){ mx = fabs(M[i * (D + 1) + p]); piv = i; }
    if (piv != p)
      for (j = 0; j <= D; j++){ double tmp = M[p * (D + 1) + j];
        M[p * (D + 1) + j] = M[piv * (D + 1) + j]; M[piv * (D + 1) + j] = tmp; }
    { double d = M[p * (D + 1) + p];
      if (fabs(d) < 1e-300) d = 1e-300;
      for (j = p; j <= D; j++) M[p * (D + 1) + j] /= d; }
    for (i = 0; i < D; i++) if (i != p){
      double f = M[i * (D + 1) + p];
      for (j = p; j <= D; j++) M[i * (D + 1) + j] -= f * M[p * (D + 1) + j];
    }
  }
  for (i = 0; i < D; i++) x[i] = M[i * (D + 1) + D];
  free(M);
}
static double es_joint_T(const double *t, int n, int P, double k,
                         double *S0, int D){
  double K, E, *m = (double*)calloc(D, sizeof(double));
  double *x = (double*)malloc(sizeof(double) * D);
  double T = 0.0; int i, p;
  es_ke(k, &K, &E);
  for (i = 0; i < n; i++)
    for (p = 1; p <= P; p++){
      double sn, cn, dn;
      es_jac(2.0 * K * p * t[i] / M_PI, k, &sn, &cn, &dn);
      m[2 * (p - 1)] += cn / (double)n; m[2 * (p - 1) + 1] += sn / (double)n;
    }
  es_solve_sym(S0, D, m, x);
  for (i = 0; i < D; i++) T += m[i] * x[i];
  free(m); free(x);
  return (double)n * T;
}
void C_ell_joint_unif(double *t, int *n, int *P, double *k, int *B, int *ng,
                      double *seed, double *out){
  int D = 2 * (*P), i, NG = *ng;
  double *S0 = (double*)malloc(sizeof(double) * D * D);
  double *m0 = (double*)malloc(sizeof(double) * D);
  double *f  = (double*)malloc(sizeof(double) * NG);
  double *S0c = (double*)malloc(sizeof(double) * D * D);
  double obs; long ge = 0; int r;
  for (i = 0; i < NG; i++) f[i] = 1.0;              /* the uniform density */
  C_ell_joint_cov(f, ng, P, k, m0, S0);
  for (i = 0; i < D * D; i++) S0c[i] = S0[i];
  obs = es_joint_T(t, *n, *P, *k, S0c, D);
#ifdef _OPENMP
#pragma omp parallel for reduction(+:ge) schedule(static) num_threads(2)
#endif
  for (r = 0; r < *B; r++){
    unsigned long long st = es_seed((unsigned long long)(*seed), r);
    double *z = (double*)malloc(sizeof(double) * (*n));
    double *Sc = (double*)malloc(sizeof(double) * D * D);
    int j;
    for (j = 0; j < *n; j++) z[j] = 2.0 * M_PI * es_unif(&st);
    for (j = 0; j < D * D; j++) Sc[j] = S0[j];
    if (es_joint_T(z, *n, *P, *k, Sc, D) >= obs) ge++;
    free(z); free(Sc);
  }
  out[0] = obs;
  out[1] = (double)(ge + 1) / (double)(*B + 1);
  free(S0); free(m0); free(f); free(S0c);
}
