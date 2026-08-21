/* arccirc: circular (wrapped) arc-length family. Pure libm; see arccirc.h for the interface. */
/* Floating-point contraction is pinned OFF so that the R and Python fronts, which compile these
   sources with different flags, cannot differ in whether multiply-add pairs are fused. Without
   this the two fronts agree to about fifteen digits and disagree in the last bit, which is enough
   to fail the parity harnesses -- and it surfaces unpredictably, because whether a fusion happens
   depends on register allocation, so an unrelated edit can expose it. */
#if defined(__clang__)
#pragma clang fp contract(off)
#elif defined(__GNUC__)
#pragma GCC optimize ("fp-contract=off")
#else
#pragma STDC FP_CONTRACT OFF
#endif

#include <math.h>
#include <stdlib.h>
#include "arccirc.h"
#ifdef _OPENMP
#include <omp.h>
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#define ARCC_TWOPI 6.28318530717958647692

/* ---- order-three smooth family ------------------------------------------------------------- */

static double arcc_q3_one(double u, double c3) {
    return 1.0 + 10.0 * c3 * u * (2.0 * u - 1.0) * (u - 1.0);
}

static double arcc_Q3_one(double u, double c3) {
    double v = 1.0 - u;
    return u + 5.0 * c3 * u * u * v * v;
}

void arcc_c3max(double *out) {
    out[0] = 3.0 * sqrt(3.0) / 5.0;
}

void arcc_qd3(const double *u, const int *nu, const double *c3, double *out) {
    int i, n = *nu;
    for (i = 0; i < n; i++) out[i] = arcc_q3_one(u[i], *c3);
}

void arcc_Q3(const double *u, const int *nu, const double *c3, double *out) {
    int i, n = *nu;
    for (i = 0; i < n; i++) out[i] = arcc_Q3_one(u[i], *c3);
}

/* wrap x to [0, 1) */
static double arcc_wrap01(double x) {
    double y = fmod(x, 1.0);
    if (y < 0.0) y += 1.0;
    return y;
}

/* wrap an angle to (-pi, pi] */
static double arcc_wrappi(double x) {
    double y = fmod(x + M_PI, ARCC_TWOPI);
    if (y < 0.0) y += ARCC_TWOPI;
    return y - M_PI;
}

/* Q is strictly increasing from 0 to 1 for admissible c3; invert by bisection. */
static double arcc_Qinv(double t, double c3) {
    double lo = 0.0, hi = 1.0, mid, f;
    int it;
    if (t <= 0.0) return 0.0;
    if (t >= 1.0) return 1.0;
    for (it = 0; it < 80; it++) {
        mid = 0.5 * (lo + hi);
        f = arcc_Q3_one(mid, c3);
        if (f < t) lo = mid; else hi = mid;
    }
    return 0.5 * (lo + hi);
}

void arcc_dens3(const double *theta, const int *nth, const double *c3, const double *mu,
                double *out) {
    int i, n = *nth;
    for (i = 0; i < n; i++) {
        double t = arcc_wrap01((theta[i] - *mu) / ARCC_TWOPI);
        double u = arcc_Qinv(t, *c3);
        double q = arcc_q3_one(u, *c3);
        out[i] = (q > 0.0) ? 1.0 / (ARCC_TWOPI * q) : 0.0;
    }
}

/* Rectangle rule. The integrand exp(2 pi i p Q(u)) is C^1 periodic on [0,1] because Q(0)=0,
 * Q(1)=1 and q(0)=q(1)=1, so the rule converges geometrically rather than at a fixed order. */
void arcc_trigmom3(const int *p, const double *c3, const double *mu, const int *nodes,
                   double *out) {
    int j, N = *nodes;
    double sr = 0.0, si = 0.0, ang, w;
    for (j = 0; j < N; j++) {
        double u = ((double) j + 0.5) / (double) N;
        ang = ARCC_TWOPI * (double) (*p) * arcc_Q3_one(u, *c3);
        sr += cos(ang);
        si += sin(ang);
    }
    sr /= (double) N;
    si /= (double) N;
    /* rotate by mu: phi_p -> exp(i p mu) phi_p */
    w = (double) (*p) * (*mu);
    out[0] = sr * cos(w) - si * sin(w);
    out[1] = sr * sin(w) + si * cos(w);
}

void arcc_rand3(const double *unif, const int *n, const double *c3, const double *mu,
                double *out) {
    int i, m = *n;
    for (i = 0; i < m; i++) {
        double th = ARCC_TWOPI * arcc_Q3_one(unif[i], *c3) + *mu;
        double y = fmod(th, ARCC_TWOPI);
        if (y < 0.0) y += ARCC_TWOPI;
        out[i] = y;
    }
}

void arcc_fit3(const double *theta, const int *n, const int *nodes, double *out) {
    int i, it, m = *n, one = 1;
    double ar = 0.0, ai = 0.0, br = 0.0, bi = 0.0, mod1, psi, cmax, lo, hi, mid, c3, sgn;
    double tm[2], modmax, zero = 0.0;

    for (i = 0; i < m; i++) {
        ar += cos(theta[i]);      ai += sin(theta[i]);
        br += cos(2.0 * theta[i]); bi += sin(2.0 * theta[i]);
    }
    ar /= (double) m; ai /= (double) m;
    br /= (double) m; bi /= (double) m;

    mod1 = sqrt(ar * ar + ai * ai);
    psi  = arcc_wrappi(atan2(bi, br) - 2.0 * atan2(ai, ar));
    sgn  = (psi >= 0.0) ? 1.0 : -1.0;

    arcc_c3max(&cmax);
    arcc_trigmom3(&one, &cmax, &zero, nodes, tm);
    modmax = sqrt(tm[0] * tm[0] + tm[1] * tm[1]);

    out[3] = 1.0;
    if (mod1 <= 0.0) {
        c3 = 0.0;
    } else if (mod1 >= modmax) {
        c3 = cmax;                 /* clipped at the admissible boundary */
        out[3] = 0.0;
    } else {
        lo = 0.0; hi = cmax;       /* |phi_1| is strictly increasing in |c3| */
        for (it = 0; it < 80; it++) {
            mid = 0.5 * (lo + hi);
            arcc_trigmom3(&one, &mid, &zero, nodes, tm);
            if (sqrt(tm[0] * tm[0] + tm[1] * tm[1]) < mod1) lo = mid; else hi = mid;
        }
        c3 = 0.5 * (lo + hi);
    }
    c3 *= sgn;

    arcc_trigmom3(&one, &c3, &zero, nodes, tm);
    out[0] = arcc_wrappi(atan2(ai, ar) - atan2(tm[1], tm[0]));
    out[1] = c3;
    out[2] = mod1;
}

/* ---- general smooth family by Fejer-Riesz ---------------------------------------------------- */

/* |P(exp(2 pi i u))|^2 with P(z) = sum_j (pre[j] + i pim[j]) z^j, evaluated by Horner. */
static double arcc_fr_mod2(double u, const double *pre, const double *pim, int np) {
    int j;
    double zr = cos(ARCC_TWOPI * u), zi = sin(ARCC_TWOPI * u);
    double sr = 0.0, si = 0.0, tr;
    for (j = np - 1; j >= 0; j--) {
        tr = sr * zr - si * zi + pre[j];
        si = sr * zi + si * zr + pim[j];
        sr = tr;
    }
    return sr * sr + si * si;
}

/* the normaliser int_0^1 |P|^2 du equals sum_j |p_j|^2 exactly (Parseval) */
static double arcc_fr_norm(const double *pre, const double *pim, int np) {
    int j; double s = 0.0;
    for (j = 0; j < np; j++) s += pre[j] * pre[j] + pim[j] * pim[j];
    return s;
}

void arcc_qd_fr(const double *u, const int *nu, const double *pre, const double *pim,
                const int *np, double *out) {
    int i, n = *nu;
    double z = arcc_fr_norm(pre, pim, *np);
    for (i = 0; i < n; i++)
        out[i] = (z > 0.0) ? arcc_fr_mod2(u[i], pre, pim, *np) / z : 1.0;
}

/* The Fejer-Riesz quantile function has no closed form, so Q is accumulated by composite Simpson
 * while the rectangle rule sweeps the midpoints. Over each cell the quadratic through the two
 * endpoints and the midpoint gives Q at the midpoint as (h/12)(5 q0 + 8 qm - q1) and over the whole
 * cell the usual (h/6)(q0 + 4 qm + q1), so Q carries O(h^4) error while the outer rule stays
 * spectral on the periodic integrand. */
void arcc_trigmom_fr(const int *p, const double *pre, const double *pim, const int *np,
                     const double *mu, const int *nodes, double *out) {
    int j, N = *nodes;
    double z = arcc_fr_norm(pre, pim, *np);
    double Q = 0.0, Qm, sr = 0.0, si = 0.0, ang, w, h = 1.0 / (double) N;
    double q0, qm, q1;
    if (z <= 0.0) z = 1.0;
    q0 = arcc_fr_mod2(0.0, pre, pim, *np) / z;
    for (j = 0; j < N; j++) {
        double um = ((double) j + 0.5) * h;
        double u1 = ((double) j + 1.0) * h;
        qm = arcc_fr_mod2(um, pre, pim, *np) / z;
        q1 = arcc_fr_mod2(u1, pre, pim, *np) / z;
        Qm = Q + h * (5.0 * q0 + 8.0 * qm - q1) / 12.0;
        ang = ARCC_TWOPI * (double) (*p) * Qm;
        sr += cos(ang);
        si += sin(ang);
        Q += h * (q0 + 4.0 * qm + q1) / 6.0;
        q0 = q1;
    }
    sr /= (double) N;
    si /= (double) N;
    w = (double) (*p) * (*mu);
    out[0] = sr * cos(w) - si * sin(w);
    out[1] = sr * sin(w) + si * cos(w);
}

/* ---- Fejer-Riesz family: quantile function, density, simulation ------------------------------- */

/* Table of Q at nodes+1 equally spaced points, by composite Simpson on each cell. */
static void arcc_fr_Qtable(const double *pre, const double *pim, int np, int nodes, double *tab) {
    int j;
    double z = arcc_fr_norm(pre, pim, np), h = 1.0 / (double) nodes, q0, qm, q1;
    if (z <= 0.0) z = 1.0;
    tab[0] = 0.0;
    q0 = arcc_fr_mod2(0.0, pre, pim, np) / z;
    for (j = 0; j < nodes; j++) {
        qm = arcc_fr_mod2(((double) j + 0.5) * h, pre, pim, np) / z;
        q1 = arcc_fr_mod2(((double) j + 1.0) * h, pre, pim, np) / z;
        tab[j + 1] = tab[j] + h * (q0 + 4.0 * qm + q1) / 6.0;
        q0 = q1;
    }
}

/* linear interpolation of the table at u in [0,1] */
static double arcc_tab_at(const double *tab, int nodes, double u) {
    double x, f;
    int j;
    if (u <= 0.0) return tab[0];
    if (u >= 1.0) return tab[nodes];
    x = u * (double) nodes;
    j = (int) x;
    if (j >= nodes) j = nodes - 1;
    f = x - (double) j;
    return tab[j] * (1.0 - f) + tab[j + 1] * f;
}

/* inverse of the table by bisection on the index, then linear interpolation */
static double arcc_tab_inv(const double *tab, int nodes, double t) {
    int lo = 0, hi = nodes, mid;
    double a, b;
    if (t <= tab[0]) return 0.0;
    if (t >= tab[nodes]) return 1.0;
    while (hi - lo > 1) {
        mid = (lo + hi) / 2;
        if (tab[mid] <= t) lo = mid; else hi = mid;
    }
    a = tab[lo]; b = tab[hi];
    if (b <= a) return (double) lo / (double) nodes;
    return ((double) lo + (t - a) / (b - a)) / (double) nodes;
}

void arcc_Q_fr(const double *u, const int *nu, const double *pre, const double *pim,
               const int *np, const int *nodes, double *out) {
    int i, n = *nu, N = *nodes;
    double *tab = (double *) malloc((size_t) (N + 1) * sizeof(double));
    if (!tab) return;
    arcc_fr_Qtable(pre, pim, *np, N, tab);
    for (i = 0; i < n; i++) out[i] = arcc_tab_at(tab, N, u[i]);
    free(tab);
}

void arcc_dens_fr(const double *theta, const int *nth, const double *pre, const double *pim,
                  const int *np, const double *mu, const int *nodes, double *out) {
    int i, n = *nth, N = *nodes;
    double z = arcc_fr_norm(pre, pim, *np), t, uu, q;
    double *tab = (double *) malloc((size_t) (N + 1) * sizeof(double));
    if (!tab) return;
    if (z <= 0.0) z = 1.0;
    arcc_fr_Qtable(pre, pim, *np, N, tab);
    for (i = 0; i < n; i++) {
        t  = arcc_wrap01((theta[i] - *mu) / ARCC_TWOPI);
        uu = arcc_tab_inv(tab, N, t);
        q  = arcc_fr_mod2(uu, pre, pim, *np) / z;
        out[i] = (q > 0.0) ? 1.0 / (ARCC_TWOPI * q) : 0.0;
    }
    free(tab);
}

void arcc_rand_fr(const double *unif, const int *n, const double *pre, const double *pim,
                  const int *np, const double *mu, const int *nodes, double *out) {
    int i, m = *n, N = *nodes;
    double th, y;
    double *tab = (double *) malloc((size_t) (N + 1) * sizeof(double));
    if (!tab) return;
    arcc_fr_Qtable(pre, pim, *np, N, tab);
    for (i = 0; i < m; i++) {
        th = ARCC_TWOPI * arcc_tab_at(tab, N, unif[i]) + *mu;
        y = fmod(th, ARCC_TWOPI);
        if (y < 0.0) y += ARCC_TWOPI;
        out[i] = y;
    }
    free(tab);
}

/* ---- circular L-moments ------------------------------------------------------------------------ */

static int arcc_dcmp(const void *a, const void *b) {
    double x = *(const double *) a, y = *(const double *) b;
    return (x < y) ? -1 : ((x > y) ? 1 : 0);
}

void arcc_lmom(const double *x, const int *n, const int *nm, double *out) {
    int i, m, N = *n, M = *nm;
    double *xs = (double *) malloc((size_t) N * sizeof(double));
    double c, wr, wi, a0, a1, sr, si;
    if (!xs) return;
    for (i = 0; i < N; i++) xs[i] = x[i];
    qsort(xs, (size_t) N, sizeof(double), arcc_dcmp);
    for (m = 1; m <= M; m++) {
        c = ARCC_TWOPI * (double) m;
        sr = 0.0; si = 0.0;
        for (i = 0; i < N; i++) {
            /* weight = (exp(-i c (i)/N) - exp(-i c (i+1)/N)) / (i c), with the order statistic x_(i+1) */
            a0 = -c * (double) i / (double) N;
            a1 = -c * (double) (i + 1) / (double) N;
            /* numerator (cos a0 - cos a1) + i(sin a0 - sin a1); divide by i c  =>  multiply by -i/c */
            wr = (sin(a0) - sin(a1)) / c;
            wi = -(cos(a0) - cos(a1)) / c;
            sr += xs[i] * wr;
            si += xs[i] * wi;
        }
        out[2 * (m - 1)]     = sr;
        out[2 * (m - 1) + 1] = si;
    }
    free(xs);
}

void arcc_rho_from_lmom(const double *ellre, const double *ellim, const int *nm, double *out) {
    int m, M = *nm;
    double c;
    for (m = 1; m <= M; m++) {
        c = ARCC_TWOPI * (double) m;                 /* rho_m = 1 + i c ell_m */
        out[2 * (m - 1)]     = 1.0 - c * ellim[m - 1];
        out[2 * (m - 1) + 1] = c * ellre[m - 1];
    }
}

/* q(u) = 1 + 2 sum_m Re(rho_m exp(2 pi i m u)) */
static double arcc_q_of_rho(double u, const double *rr, const double *ri, int M) {
    int m; double s = 1.0, a;
    for (m = 1; m <= M; m++) {
        a = ARCC_TWOPI * (double) m * u;
        s += 2.0 * (rr[m - 1] * cos(a) - ri[m - 1] * sin(a));
    }
    return s;
}

void arcc_qmin_rho(const double *rhore, const double *rhoim, const int *nm, const int *nodes,
                   double *out) {
    int j, N = *nodes;
    double mn = 1.0e300, v;
    for (j = 0; j < N; j++) {
        v = arcc_q_of_rho((double) j / (double) N, rhore, rhoim, *nm);
        if (v < mn) mn = v;
    }
    out[0] = mn;
}

void arcc_admiss(const double *rhore, const double *rhoim, const int *nm, const double *margin,
                 const int *nodes, double *out) {
    int m, it, M = *nm;
    double lo = 0.0, hi = 1.0, t, mn;
    double *rr = (double *) malloc((size_t) M * sizeof(double));
    double *ri = (double *) malloc((size_t) M * sizeof(double));
    if (!rr || !ri) { free(rr); free(ri); return; }

    for (m = 0; m < M; m++) { rr[m] = rhore[m]; ri[m] = rhoim[m]; }
    arcc_qmin_rho(rr, ri, nm, nodes, &mn);
    if (mn >= *margin) {                       /* already admissible: report t = 1 exactly */
        for (m = 0; m < M; m++) { out[2 * m] = rhore[m]; out[2 * m + 1] = rhoim[m]; }
        out[2 * M] = 1.0;
        free(rr); free(ri);
        return;
    }
    for (it = 0; it < 80; it++) {              /* bisection on the shrink factor */
        t = 0.5 * (lo + hi);
        for (m = 0; m < M; m++) { rr[m] = t * rhore[m]; ri[m] = t * rhoim[m]; }
        arcc_qmin_rho(rr, ri, nm, nodes, &mn);
        if (mn >= *margin) lo = t; else hi = t;
    }
    t = lo;
    for (m = 0; m < M; m++) { out[2 * m] = t * rhore[m]; out[2 * m + 1] = t * rhoim[m]; }
    out[2 * M] = t;
    free(rr); free(ri);
}

/* ---- Fejer-Riesz spectral factorisation and the closed-form fit -------------------------------- */

/* Minimal complex arithmetic on pairs of doubles, so that nothing beyond libm is required. */
typedef struct { double r, i; } arcc_cx;
static arcc_cx arcc_cxm(arcc_cx a, arcc_cx b) {
    arcc_cx z; z.r = a.r * b.r - a.i * b.i; z.i = a.r * b.i + a.i * b.r; return z;
}
static arcc_cx arcc_cxs(arcc_cx a, arcc_cx b) { arcc_cx z; z.r = a.r - b.r; z.i = a.i - b.i; return z; }
static arcc_cx arcc_cxd(arcc_cx a, arcc_cx b) {
    arcc_cx z; double d = b.r * b.r + b.i * b.i;
    if (d == 0.0) { z.r = 0.0; z.i = 0.0; return z; }
    z.r = (a.r * b.r + a.i * b.i) / d; z.i = (a.i * b.r - a.r * b.i) / d; return z;
}
static double arcc_cxmod(arcc_cx a) { return sqrt(a.r * a.r + a.i * a.i); }

/* Horner evaluation of a polynomial given in increasing powers */
static arcc_cx arcc_polyval(const arcc_cx *c, int d, arcc_cx z) {
    int k; arcc_cx s = c[d];
    for (k = d - 1; k >= 0; k--) s = arcc_cxm(s, z), s.r += c[k].r, s.i += c[k].i;
    return s;
}

void arcc_factorise(const double *rhore, const double *rhoim, const int *nm, double *out) {
    int M = *nm, d = 2 * M, k, i, j, it, nin;
    double maxmod = 0.0, tol, nrm;
    arcc_cx *c, *z, *w, num, den, diff, one;
    if (M < 1) { out[0] = 1.0; out[1] = 0.0; return; }

    c = (arcc_cx *) malloc((size_t) (d + 1) * sizeof(arcc_cx));
    z = (arcc_cx *) malloc((size_t) d * sizeof(arcc_cx));
    w = (arcc_cx *) malloc((size_t) (M + 1) * sizeof(arcc_cx));
    if (!c || !z || !w) { free(c); free(z); free(w); return; }

    /* coefficients c_k = rho_{k-M}, with rho_0 = 1 and rho_{-m} = conj(rho_m) */
    for (k = 0; k <= d; k++) {
        if (k < M)       { c[k].r =  rhore[M - k - 1]; c[k].i = -rhoim[M - k - 1]; }
        else if (k == M) { c[k].r =  1.0;              c[k].i =  0.0;              }
        else             { c[k].r =  rhore[k - M - 1]; c[k].i =  rhoim[k - M - 1]; }
    }
    for (k = 0; k <= d; k++) { double m2 = arcc_cxmod(c[k]); if (m2 > maxmod) maxmod = m2; }
    tol = 1.0e-12 * (maxmod > 0.0 ? maxmod : 1.0);
    while (d > 0 && arcc_cxmod(c[d]) <= tol) d--;      /* trim a degenerate leading coefficient */
    if (d < 2) {                                        /* spectrum is (numerically) flat */
        out[0] = 1.0; out[1] = 0.0;
        for (j = 1; j <= M; j++) { out[2 * j] = 0.0; out[2 * j + 1] = 0.0; }
        free(c); free(z); free(w);
        return;
    }

    /* Durand-Kerner from the standard spiral start */
    for (i = 0; i < d; i++) {
        double ang = 2.0 * M_PI * (double) i / (double) d + 0.4;
        z[i].r = 0.9 * cos(ang); z[i].i = 0.9 * sin(ang);
    }
    for (it = 0; it < 500; it++) {
        double move = 0.0;
        for (i = 0; i < d; i++) {
            num = arcc_polyval(c, d, z[i]);
            den = c[d];
            for (j = 0; j < d; j++) if (j != i) den = arcc_cxm(den, arcc_cxs(z[i], z[j]));
            diff = arcc_cxd(num, den);
            z[i] = arcc_cxs(z[i], diff);
            move += arcc_cxmod(diff);
        }
        if (move < 1.0e-14) break;
    }

    /* keep the d/2 roots of smallest modulus (inside the unit disc) */
    for (i = 0; i < d; i++)                                    /* selection sort by modulus */
        for (j = i + 1; j < d; j++)
            if (arcc_cxmod(z[j]) < arcc_cxmod(z[i])) { arcc_cx t = z[i]; z[i] = z[j]; z[j] = t; }
    nin = d / 2;

    /* P(z) = prod (z - r_i), coefficients in increasing powers */
    one.r = 1.0; one.i = 0.0;
    w[0] = one;
    for (i = 1; i <= M; i++) { w[i].r = 0.0; w[i].i = 0.0; }
    for (i = 0; i < nin; i++) {
        for (j = i + 1; j >= 1; j--) w[j] = arcc_cxs(w[j - 1], arcc_cxm(w[j], z[i]));
        w[0] = arcc_cxm(w[0], z[i]);
        w[0].r = -w[0].r; w[0].i = -w[0].i;
    }
    nrm = 0.0;
    for (j = 0; j <= M; j++) nrm += w[j].r * w[j].r + w[j].i * w[j].i;
    nrm = (nrm > 0.0) ? sqrt(nrm) : 1.0;
    for (j = 0; j <= M; j++) { out[2 * j] = w[j].r / nrm; out[2 * j + 1] = w[j].i / nrm; }

    free(c); free(z); free(w);
}

void arcc_fit_fr(const double *theta, const int *n, const int *nm, const double *margin,
                 const int *nodes, double *out) {
    int i, m, N = *n, M = *nm;
    double ar = 0.0, ai = 0.0, mu, tshrink;
    double *x, *ell, *sh, *elr, *eli, *rhor, *rhoi;

    if (M < 1) return;
    x    = (double *) malloc((size_t) N * sizeof(double));
    ell  = (double *) malloc((size_t) (2 * M) * sizeof(double));
    sh   = (double *) malloc((size_t) (2 * M + 1) * sizeof(double));
    elr  = (double *) malloc((size_t) M * sizeof(double));
    eli  = (double *) malloc((size_t) M * sizeof(double));
    rhor = (double *) malloc((size_t) M * sizeof(double));
    rhoi = (double *) malloc((size_t) M * sizeof(double));
    if (!x || !ell || !sh || !elr || !eli || !rhor || !rhoi) {
        free(x); free(ell); free(sh); free(elr); free(eli); free(rhor); free(rhoi);
        return;
    }

    /* The pair (p, mu) is not identified: shifting the cut by delta is absorbed exactly by the
     * phase ramp p_j -> p_j exp(2 pi i j delta), so the same circular law has a one-parameter
     * family of representations. The gauge is fixed here by cutting at the sample mean direction,
     * which makes the fit a single pass with no iteration anywhere. */
    for (i = 0; i < N; i++) { ar += cos(theta[i]); ai += sin(theta[i]); }
    mu = atan2(ai, ar);
    for (i = 0; i < N; i++) x[i] = arcc_wrap01((theta[i] - mu) / ARCC_TWOPI);

    arcc_lmom(x, n, nm, ell);
    for (m = 0; m < M; m++) { elr[m] = ell[2 * m]; eli[m] = ell[2 * m + 1]; }
    arcc_rho_from_lmom(elr, eli, nm, sh);              /* sh temporarily holds rho */
    for (m = 0; m < M; m++) { rhor[m] = sh[2 * m]; rhoi[m] = sh[2 * m + 1]; }
    arcc_admiss(rhor, rhoi, nm, margin, nodes, sh);
    tshrink = sh[2 * M];
    for (m = 0; m < M; m++) { rhor[m] = sh[2 * m]; rhoi[m] = sh[2 * m + 1]; }
    arcc_factorise(rhor, rhoi, nm, out);               /* writes p into out[0 .. 2M+1] */

    out[2 * (M + 1)]     = mu;
    out[2 * (M + 1) + 1] = tshrink;
    out[2 * (M + 1) + 2] = (tshrink >= 1.0) ? 1.0 : 0.0;

    free(x); free(ell); free(sh); free(elr); free(eli); free(rhor); free(rhoi);
}


/* ---- arc-length tempering of a von Mises base -------------------------------------------------- */

/* Modified Bessel function of the first kind, order zero; Abramowitz and Stegun 9.8.1 and 9.8.2.
 * Only libm is used, as everywhere else in this file. */
static double arcc_bessi0(double x) {
    double ax = fabs(x), y, r;
    if (ax < 3.75) {
        y = x / 3.75; y *= y;
        r = 1.0 + y*(3.5156229 + y*(3.0899424 + y*(1.2067492
            + y*(0.2659732 + y*(0.0360768 + y*0.0045813)))));
    } else {
        y = 3.75 / ax;
        r = (exp(ax)/sqrt(ax)) * (0.39894228 + y*(0.01328592 + y*(0.00225319
            + y*(-0.00157565 + y*(0.00916281 + y*(-0.02057706 + y*(0.02635537
            + y*(-0.01647633 + y*0.00392377))))))));
    }
    return r;
}

/* unnormalised tempered density at theta */
static double arcc_temper_raw(double theta, double kappa, double mu, double s, double i0) {
    double f = exp(kappa * cos(theta - mu)) / (ARCC_TWOPI * i0);
    double sf = s * f;
    return sqrt(1.0 + sf * sf);
}

void arcc_temper_vm(const double *theta, const int *nth, const double *kappa, const double *mu,
                    const double *s, const int *nodes, double *out) {
    int i, j, n = *nth, N = *nodes;
    double i0 = arcc_bessi0(*kappa), z = 0.0, h = ARCC_TWOPI / (double) N;
    for (j = 0; j < N; j++)
        z += arcc_temper_raw(((double) j + 0.5) * h, *kappa, *mu, *s, i0);
    z *= h;
    if (z <= 0.0) z = 1.0;
    for (i = 0; i < n; i++)
        out[i] = arcc_temper_raw(theta[i], *kappa, *mu, *s, i0) / z;
}

void arcc_temper_vm_trigmom(const int *p, const double *kappa, const double *mu, const double *s,
                            const int *nodes, double *out) {
    int j, N = *nodes;
    double i0 = arcc_bessi0(*kappa), z = 0.0, sr = 0.0, si = 0.0;
    double h = ARCC_TWOPI / (double) N, th, g;
    for (j = 0; j < N; j++) {
        th = ((double) j + 0.5) * h;
        g  = arcc_temper_raw(th, *kappa, *mu, *s, i0);
        z  += g;
        sr += g * cos((double) (*p) * th);
        si += g * sin((double) (*p) * th);
    }
    if (z <= 0.0) z = 1.0;
    out[0] = sr / z;
    out[1] = si / z;
}

/* ---- exact inference for the shape parameter ------------------------------------------------
 * Maximum likelihood is unavailable for this family: a free rotation places the near-zero of the
 * quantile density onto an observation, so the likelihood is unbounded. Asymptotics are no help
 * either, because the admissible boundary |c3| = 3 sqrt(3) / 5 is attained. What remains, and is
 * exact, is inversion of a Monte Carlo test on a ROTATION-INVARIANT statistic:
 *   - the rotation is eliminated exactly by invariance, so there is no nuisance to profile;
 *   - the level is exact for every n and every B (Dwass 1957; Hope 1968), not asymptotically;
 *   - the parameter space is a compact interval, so a grid over it is exhaustive rather than a
 *     search, and sampling the family is one transform of a uniform.
 * Two statistics are offered. The modulus of the first trigonometric moment is the paper's own
 * closed-form object and is about twice as efficient here, because it grows linearly in c3 where
 * the arc-length functional grows quadratically and is nearly flat near the circular uniform. The
 * arc-length spacings functional is retained because it is the goodness-of-fit statistic of the
 * programme's own test, and because it responds to structure the first moment cannot see. Both are
 * reflection invariants, so what is bracketed is |c3|; the sign needs a reflection-odd statistic.
 * Recorded resolution is handled by grouping the reference draws exactly as the data are grouped. */
static unsigned long long arcc_seed(int seed, int j) {
    unsigned long long z = (unsigned long long)seed * 0x9E3779B97F4A7C15ULL
                           ^ ((unsigned long long)(j + 1) * 0xBF58476D1CE4E5B9ULL);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}
static double arcc_unif(unsigned long long *s) {
    unsigned long long z = (*s += 0x9E3779B97F4A7C15ULL);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    z ^= (z >> 31);
    return (double)(z >> 11) * (1.0 / 9007199254740992.0);
}
/* group to the nearest of k equal cells on the circle; k <= 0 leaves the value alone */
static double arcc_group(double th01, int k) {
    if (k <= 0) return th01;
    double g = floor(th01 * (double)k + 0.5) / (double)k;
    return g - floor(g);
}
/* statistic on angles given as fractions of a turn in [0,1) */
static double arcc_stat(const double *u01, int n, int which, double *work) {
    if (which == 0) {                                  /* modulus of the first trigonometric moment */
        double cr = 0.0, ci = 0.0;
        for (int i = 0; i < n; i++) { cr += cos(ARCC_TWOPI * u01[i]); ci += sin(ARCC_TWOPI * u01[i]); }
        cr /= n; ci /= n;
        return sqrt(cr * cr + ci * ci);
    }
    for (int i = 0; i < n; i++) work[i] = u01[i];      /* arc length of the circular-gap path */
    qsort(work, (size_t)n, sizeof(double), arcc_dcmp);
    double c = 1.0 / (double)n, s = 0.0;
    for (int i = 0; i < n - 1; i++) { double g = work[i + 1] - work[i]; s += sqrt(g * g + c * c); }
    double gl = 1.0 - (work[n - 1] - work[0]);
    return s + sqrt(gl * gl + c * c);
}
/* entry: exact confidence interval for |c3| by test inversion.
 * out[0], out[1] = interval endpoints (NaN if the set is empty, which is itself informative:
 * no admissible member is consistent with the data); out[2] = the observed statistic;
 * out[3 .. 3+ng-1] = the p-value curve over the grid, so the paper can plot what it inverted. */
void arcc_exact_ci(const double *theta, const int *np, const int *Bp, const int *seedp,
                   const double *cgrid, const int *ngp, const int *statp, const int *groupp,
                   const double *levp, double *out) {
    const int n = *np, B = *Bp, ng = *ngp, which = *statp, grp = *groupp;
    const double lev = *levp;
    double *u = (double *)malloc((size_t)n * sizeof(double));
    double *wk = (double *)malloc((size_t)n * sizeof(double));
    if (!u || !wk) { free(u); free(wk); out[0] = out[1] = 0.0 / 0.0; return; }
    for (int i = 0; i < n; i++) {                      /* radians -> fraction of a turn, grouped */
        double v = theta[i] / ARCC_TWOPI; v -= floor(v);
        u[i] = arcc_group(v, grp);
    }
    const double tobs = arcc_stat(u, n, which, wk);
    free(u); free(wk);
    out[2] = tobs;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
    for (int j = 0; j < ng; j++) {
        double *ub = (double *)malloc((size_t)n * sizeof(double));
        double *w2 = (double *)malloc((size_t)n * sizeof(double));
        if (!ub || !w2) { free(ub); free(w2); out[3 + j] = 0.0 / 0.0; continue; }
        unsigned long long st = arcc_seed(*seedp, j);
        int le = 0, ge = 0;
        for (int b = 0; b < B; b++) {
            double rot = arcc_unif(&st);
            for (int i = 0; i < n; i++) {
                double v = arcc_Q3_one(arcc_unif(&st), cgrid[j]) + rot;
                v -= floor(v);
                ub[i] = arcc_group(v, grp);
            }
            double t = arcc_stat(ub, n, which, w2);
            if (t <= tobs) le++;
            if (t >= tobs) ge++;
        }
        double pl = (1.0 + le) / (B + 1.0), pu = (1.0 + ge) / (B + 1.0);
        double p = 2.0 * (pl < pu ? pl : pu);
        out[3 + j] = p > 1.0 ? 1.0 : p;
        free(ub); free(w2);
    }
    double lo = 0.0 / 0.0, hi = 0.0 / 0.0;
    for (int j = 0; j < ng; j++) if (out[3 + j] > lev) {
        if (!(lo == lo)) lo = cgrid[j];
        hi = cgrid[j];
    }
    out[0] = lo; out[1] = hi;
}

/* entry: goodness of fit of a fitted member, by the arc-length spacings statistic of the
 * programme's own test. The probability-integral transform through the fitted member sends the
 * sample to uniform under the hypothesis, so the reference law is DISTRIBUTION-FREE: the reference
 * draws are uniforms, whatever the fitted shape. out[0] = statistic, out[1] = Monte Carlo p-value,
 * exact in level for any B by the same argument as the interval above. */
void arcc_gof(const double *theta, const int *np, const double *c3p, const double *mup,
              const int *Bp, const int *seedp, double *out) {
    const int n = *np, B = *Bp;
    double *u = (double *)malloc((size_t)n * sizeof(double));
    double *w = (double *)malloc((size_t)n * sizeof(double));
    if (!u || !w) { free(u); free(w); out[0] = out[1] = 0.0 / 0.0; return; }
    for (int i = 0; i < n; i++) {
        double v = (theta[i] - *mup) / ARCC_TWOPI;      /* undo the rotation, then invert Q */
        v -= floor(v);
        u[i] = arcc_Qinv(v, *c3p);
    }
    const double tobs = arcc_stat(u, n, 1, w);
    unsigned long long st = arcc_seed(*seedp, 0);
    int ge = 0;
    for (int b = 0; b < B; b++) {
        for (int i = 0; i < n; i++) u[i] = arcc_unif(&st);
        if (arcc_stat(u, n, 1, w) >= tobs) ge++;
    }
    out[0] = tobs;
    out[1] = (1.0 + ge) / (B + 1.0);
    free(u); free(w);
}
