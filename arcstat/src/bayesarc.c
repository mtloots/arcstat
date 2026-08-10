/* bayesarc.c : shared C back-end for the Bayesian arc-length goodness-of-fit test.
 * Pure libm, no external numerical libraries. Bound identically from R and Python. */
#include "bayesarc.h"
#include <math.h>
#include <stdlib.h>

static const double SQRT2 = 1.41421356237309514880;

/* ---- deterministic PRNG: splitmix64 seed -> xoshiro256** (reproducible across languages) ---- */
static unsigned long long sm_state;
static unsigned long long splitmix64(void) {
  unsigned long long z = (sm_state += 0x9E3779B97F4A7C15ULL);
  z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
  z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
  return z ^ (z >> 31);
}
static unsigned long long xs[4];
static void seed_rng(unsigned long long seed) {
  sm_state = seed;
  for (int i = 0; i < 4; i++) xs[i] = splitmix64();
}
static inline unsigned long long rotl(unsigned long long x, int k) {
  return (x << k) | (x >> (64 - k));
}
static unsigned long long next_u64(void) {
  unsigned long long result = rotl(xs[1] * 5, 7) * 9;
  unsigned long long t = xs[1] << 17;
  xs[2] ^= xs[0]; xs[3] ^= xs[1]; xs[1] ^= xs[2]; xs[0] ^= xs[3];
  xs[2] ^= t; xs[3] = rotl(xs[3], 45);
  return result;
}
/* uniform in (0,1) */
static inline double next_unif(void) {
  return ((next_u64() >> 11) + 0.5) * (1.0 / 9007199254740992.0);
}

/* arc-length and KS discrepancy of the ogive for sorted u (length n) and Dirichlet weights w */
static void ogive_disc(int n, const double *u, const double *w, double *arc, double *ks) {
  double prevx = 0.0, prevy = 0.0, C = 0.0, a = 0.0, k = 0.0;
  for (int i = 0; i < n; i++) {
    double dx = u[i] - prevx;
    double dy = w[i];                    /* vertical step = weight */
    a += sqrt(dx * dx + dy * dy);
    C += w[i];
    double d = fabs(C - u[i]); if (d > k) k = d;
    prevx = u[i]; prevy = C;
  }
  a += (1.0 - prevx);                    /* closing segment to (1,1); vertical part is 0 */
  *arc = a - SQRT2;
  *ks  = k;
}

/* one Bayesian-bootstrap posterior draw: Dirichlet(1,..,1) weights = normalised Exp(1) */
static void one_draw(int n, const double *u, double *w, double *arc, double *ks) {
  double s = 0.0;
  for (int i = 0; i < n; i++) { double e = -log(next_unif()); w[i] = e; s += e; }
  for (int i = 0; i < n; i++) w[i] /= s;
  ogive_disc(n, u, w, arc, ks);
}

void bb_post(const int *n, const double *u, const int *M,
             const int *seed, double *out) {
  seed_rng((unsigned long long)(unsigned int)(*seed));
  double *w = (double *)malloc((size_t)(*n) * sizeof(double));
  for (int m = 0; m < *M; m++)
    one_draw(*n, u, w, &out[2 * m], &out[2 * m + 1]);
  free(w);
}

/* insertion sort (n is modest; keeps the core dependency-free) */
static void isort(double *a, int n) {
  for (int i = 1; i < n; i++) { double v = a[i]; int j = i - 1;
    while (j >= 0 && a[j] > v) { a[j + 1] = a[j]; j--; } a[j + 1] = v; }
}

void bb_ref(const int *n, const int *D, const int *m,
            const int *seed, double *out) {
  seed_rng((unsigned long long)(unsigned int)(*seed));
  int nn = *n;
  double *u = (double *)malloc((size_t)nn * sizeof(double));
  double *w = (double *)malloc((size_t)nn * sizeof(double));
  int idx = 0;
  for (int d = 0; d < *D; d++) {
    for (int i = 0; i < nn; i++) u[i] = next_unif();   /* uniform H0 data */
    isort(u, nn);
    for (int j = 0; j < *m; j++) { one_draw(nn, u, w, &out[2 * idx], &out[2 * idx + 1]); idx++; }
  }
  free(u); free(w);
}
