/* arcvc.c -- variance components and locus geometry for the induction-curve application.
 *
 * Both routines exist because the paper's application section now argues at the level of OILS
 * rather than RUNS: replicate runs on one specimen are not independent, so the clustering has to
 * be measured, and the price of the equivalence constraint has to be related to how far a fitted
 * shape sits from the locus. Both quantities are load-bearing for stated results, so they belong
 * in the shared back end with the rest of the numerics rather than in one front.
 *
 * No R headers: this file is compiled into the standalone shared object that the Python front
 * loads by ctypes as well as into the R package.
 */

/* Floating-point contraction is pinned OFF for the same reason as in arck4.c: R and the
   Python front compile with different flags, and whether a multiply-add pair is fused
   changes the last bits of a distance. Removing the difference at its source keeps the
   parity harness able to demand exact equality rather than a tolerance. */
#if defined(__clang__)
#pragma clang fp contract(off)
#elif defined(__GNUC__)
#pragma GCC optimize ("fp-contract=off")
#else
#pragma STDC FP_CONTRACT OFF
#endif

#include <math.h>
#include <stdlib.h>
#include "arcvc.h"

/* One-way random-effects variance components.
 *
 * y[0..n-1] are the observations, g[0..n-1] the zero-based group labels, ng the number of groups.
 * Groups may be of unequal size, which they are here (three to six runs per oil), so the
 * between-group mean square is divided by the UNBALANCED constant
 *
 *     k0 = (N - sum_i n_i^2 / N) / (ng - 1)
 *
 * and not by the mean group size. Using the mean size is the balanced-design shortcut and it
 * biases the variance ratio whenever the groups differ in size.
 *
 * out[0] intraclass correlation, out[1] within-group sd, out[2] between-group sd.
 * A negative between-group variance estimate is truncated at zero, which is the usual convention
 * and is what the corresponding R and Python fronts document.
 */
void arcvc_icc_oneway(const double *y, const int *g, const int *n, const int *ng, double *out)
{
    int N = *n, G = *ng, i;
    double grand = 0.0, ssb = 0.0, ssw = 0.0, sumsq = 0.0, msb, msw, k0, vb;

    if (N <= 0 || G <= 1 || N <= G) { out[0] = out[1] = out[2] = NAN; return; }

    /* Group sums and counts. Allocated rather than taken from a fixed buffer: the number of groups
     * is a property of the data, and a capped buffer turns a larger study into a silent NaN. */
    double *gs = (double *) calloc((size_t) G, sizeof(double));
    int    *gc = (int *)    calloc((size_t) G, sizeof(int));
    if (gs == NULL || gc == NULL) { free(gs); free(gc); out[0] = out[1] = out[2] = NAN; return; }

    for (i = 0; i < N; i++) {
        int gi = g[i];
        if (gi < 0 || gi >= G) { free(gs); free(gc); out[0] = out[1] = out[2] = NAN; return; }
        gs[gi] += y[i]; gc[gi] += 1; grand += y[i];
    }
    grand /= (double) N;

    for (i = 0; i < G; i++) {
        if (gc[i] <= 0) { free(gs); free(gc); out[0] = out[1] = out[2] = NAN; return; }
        double gm = gs[i] / (double) gc[i];
        ssb   += (double) gc[i] * (gm - grand) * (gm - grand);
        sumsq += (double) gc[i] * (double) gc[i];
    }
    for (i = 0; i < N; i++) {
        double gm = gs[g[i]] / (double) gc[g[i]];
        ssw += (y[i] - gm) * (y[i] - gm);
    }
    free(gs); free(gc);

    msb = ssb / (double) (G - 1);
    msw = ssw / (double) (N - G);
    k0  = ((double) N - sumsq / (double) N) / (double) (G - 1);

    vb = (k0 > 0.0) ? (msb - msw) / k0 : 0.0;
    if (vb < 0.0) vb = 0.0;

    out[0] = (vb + msw > 0.0) ? vb / (vb + msw) : NAN;
    out[1] = sqrt(msw);
    out[2] = sqrt(vb);
}

/* Shortest distance from each fitted shape (h,k) to the equivalence locus.
 *
 * The locus arrives as a polyline (lh, lk) of nl vertices. Distance is measured to the SEGMENTS
 * rather than to the vertices: a vertex-only search overstates the distance by up to half the
 * vertex spacing, which on the grid used here is enough to matter when the quantity is being
 * correlated against the price of the constraint.
 */
void arcvc_locus_dist(const double *h, const double *k, const int *n,
                      const double *lh, const double *lk, const int *nl, double *out)
{
    int N = *n, L = *nl, i, j;

    if (L <= 0) { for (i = 0; i < N; i++) out[i] = NAN; return; }

    for (i = 0; i < N; i++) {
        double px = h[i], py = k[i], best;
        if (!isfinite(px) || !isfinite(py)) { out[i] = NAN; continue; }

        best = INFINITY;
        if (L == 1) {
            double dx = px - lh[0], dy = py - lk[0];
            best = sqrt(dx * dx + dy * dy);
        } else {
            for (j = 0; j < L - 1; j++) {
                double ax = lh[j],     ay = lk[j];
                double bx = lh[j + 1], by = lk[j + 1];
                double vx = bx - ax,   vy = by - ay;
                double wx = px - ax,   wy = py - ay;
                double vv = vx * vx + vy * vy;
                double t  = (vv > 0.0) ? (wx * vx + wy * vy) / vv : 0.0;
                double dx, dy, d;
                if (t < 0.0) t = 0.0;
                if (t > 1.0) t = 1.0;
                dx = px - (ax + t * vx);
                dy = py - (ay + t * vy);
                d  = sqrt(dx * dx + dy * dy);
                if (d < best) best = d;
            }
        }
        out[i] = best;
    }
}
