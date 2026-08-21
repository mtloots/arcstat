/* arccirc: shared C back-end for the circular (wrapped) arc-length family.
 * Pure libm; no external numerical libraries. Bound identically from R and Python.
 *
 * The linear quantile arc-length family has Q(0)=0 and Q(1)=1 identically, so theta = 2 pi Q(U)
 * wraps the circle exactly once: the wrapping series has a single term. Smoothness across the join
 * costs two linear constraints on the shifted-Legendre coefficients,
 *     C^0 : sum over odd k of c_k = 0,      C^1 : sum over even k >= 2 of c_k k(k+1) = 0,
 * which at order two admit only the uniform. Order three is the first smooth non-uniform member and
 * leaves one free shape parameter c3, with the closed forms implemented here.
 *
 * Because Q(0)=0, Q(1)=1 and q(0)=q(1), the integrand of the trigonometric moments is C^1 periodic
 * on [0,1], so the uniform rectangle rule converges geometrically; no quadrature table is required. */
#ifndef ARCCIRC_H
#define ARCCIRC_H

/* Largest admissible |c3|: q >= 0 on [0,1] iff |c3| <= 3 sqrt(3) / 5. out has length 1. */
void arcc_c3max(double *out);

/* Quantile density q(u) = 1 + 10 c3 u (2u-1)(u-1) of the order-three smooth family,
 * at the nu points u; out has length nu. Negative values are returned as they stand so that
 * callers can test admissibility rather than have it silently imposed. */
void arcc_qd3(const double *u, const int *nu, const double *c3, double *out);

/* Quantile function Q(u) = u + 5 c3 u^2 (1-u)^2 at the nu points u; out has length nu. */
void arcc_Q3(const double *u, const int *nu, const double *c3, double *out);

/* Circular density at the nth angles theta (radians), for shape c3 and mean direction mu.
 * g(theta) = 1 / (2 pi q(u)) with u = Q^{-1}(((theta - mu)/(2 pi)) mod 1), the inverse taken by
 * bisection (Q is strictly increasing whenever c3 is admissible). out has length nth. */
void arcc_dens3(const double *theta, const int *nth, const double *c3, const double *mu,
                double *out);

/* Trigonometric moment phi_p = E[exp(i p Theta)] for Theta = 2 pi Q(U) + mu, by the rectangle rule
 * with `nodes` points. out has length 2: out[0] = Re phi_p, out[1] = Im phi_p. */
void arcc_trigmom3(const int *p, const double *c3, const double *mu, const int *nodes,
                   double *out);

/* Simulation. The caller supplies n uniforms, so the routine is deterministic and its output is
 * byte-identical from R and from Python. out has length n and holds angles in [0, 2 pi). */
void arcc_rand3(const double *unif, const int *n, const double *c3, const double *mu, double *out);

/* Method-of-moments fit of (mu, c3) from n angles. |phi_1| is even in c3 and strictly increasing in
 * |c3|, so it identifies the magnitude; the rotation-invariant psi = arg(phi_2) - 2 arg(phi_1),
 * reduced to the principal branch, identifies the sign. out has length 4:
 * out[0] = mu, out[1] = c3, out[2] = |phi_1|, out[3] = 1 if |phi_1| lay inside the attainable
 * range and the inversion converged, 0 if it was clipped to the boundary. */
void arcc_fit3(const double *theta, const int *n, const int *nodes, double *out);

/* General smooth family by Fejer-Riesz: q(u) = |P(exp(2 pi i u))|^2 normalised to integrate to one,
 * with P(z) = sum_{j=0}^{np-1} (pre[j] + i pim[j]) z^j. Every such q is non-negative and C-infinity
 * periodic by construction, so no admissibility or smoothness constraint arises. out has length nu. */
void arcc_qd_fr(const double *u, const int *nu, const double *pre, const double *pim,
                const int *np, double *out);

/* Trigonometric moment of the Fejer-Riesz family, rectangle rule with `nodes` points.
 * out has length 2 (real, imaginary). */
void arcc_trigmom_fr(const int *p, const double *pre, const double *pim, const int *np,
                     const double *mu, const int *nodes, double *out);

/* ---- Fejer-Riesz family: quantile function, density, simulation ------------------------------- */

/* Quantile function Q(u) of the Fejer-Riesz family at the nu points u, by composite Simpson on
 * `nodes` cells. Q(0) = 0 and Q(1) = 1 whatever the coefficients, since q integrates to one. */
void arcc_Q_fr(const double *u, const int *nu, const double *pre, const double *pim,
               const int *np, const int *nodes, double *out);

/* Circular density of the Fejer-Riesz family at the nth angles theta, mean direction mu. */
void arcc_dens_fr(const double *theta, const int *nth, const double *pre, const double *pim,
                  const int *np, const double *mu, const int *nodes, double *out);

/* Simulation from the Fejer-Riesz family. As elsewhere the caller supplies the uniforms, so the
 * routine is deterministic and its output is byte-identical from R and from Python. */
void arcc_rand_fr(const double *unif, const int *n, const double *pre, const double *pim,
                  const int *np, const double *mu, const int *nodes, double *out);

/* ---- circular L-moments ------------------------------------------------------------------------ */

/* Circular L-moments of a sample x on [0,1): ell_m = int_0^1 Qhat(u) exp(-2 pi i m u) du for
 * m = 1..nm. With Qhat the step function through the order statistics this is exactly a linear
 * combination of order statistics with fixed complex weights, so it is a genuine L-statistic.
 * out has length 2*nm, holding (Re ell_1, Im ell_1, Re ell_2, ...). */
void arcc_lmom(const double *x, const int *n, const int *nm, double *out);

/* The exact linear inversion rho_m = 1 + 2 pi i m ell_m, the circular counterpart of the Legendre
 * L-moment map. out has length 2*nm. */
void arcc_rho_from_lmom(const double *ellre, const double *ellim, const int *nm, double *out);

/* Smallest value of the quantile density implied by a spectrum rho_1..rho_nm (with rho_0 = 1),
 * evaluated on `nodes` points. A negative value means the estimate has left the admissible set. */
void arcc_qmin_rho(const double *rhore, const double *rhoim, const int *nm, const int *nodes,
                   double *out);

/* Shrink a spectrum towards the circular uniform, rho_m -> t rho_m, by bisection on t in [0,1]
 * until the implied quantile density has minimum at least `margin`. Reports the factor used, which
 * is exactly one when the estimate was already admissible. out has length 2*nm + 1: the shrunken
 * spectrum followed by t. */
void arcc_admiss(const double *rhore, const double *rhoim, const int *nm, const double *margin,
                 const int *nodes, double *out);

/* Fejer-Riesz spectral factorisation: given rho_1..rho_nm with rho_0 = 1, recover the minimum-phase
 * coefficient vector p with |P(exp(2 pi i u))|^2 / sum|p_j|^2 = q(u). The degree-2*nm Laurent
 * polynomial is rooted by Durand-Kerner and the nm roots inside the unit disc are kept.
 * out has length 2*(nm+1). */
void arcc_factorise(const double *rhore, const double *rhoim, const int *nm, double *out);

/* Closed-form fit of the Fejer-Riesz family from n angles: sort, apply the L-moment weights, invert
 * the exact linear map, factorise. No numerical search enters anywhere.
 * The pair (p, mu) is NOT identified -- shifting the cut by delta is absorbed exactly by the phase
 * ramp p_j -> p_j exp(2 pi i j delta) -- so the gauge is fixed by cutting at the sample mean
 * direction. Two fits should therefore be compared by the circular density they imply, never by
 * their coefficient vectors, which need not agree.
 * out has length 2*(nm+1) + 3: the coefficients p, then mu, the shrink factor (one if the estimate
 * was admissible), and the admissibility flag. */
void arcc_fit_fr(const double *theta, const int *n, const int *nm, const double *margin,
                 const int *nodes, double *out);

/* ---- arc-length tempering of a circular base -------------------------------------------------- */

/* The generator of the linear theory, applied to a base density on the circle, is
 *      g(theta) = sqrt(1 + s^2 f(theta)^2) / int_0^{2pi} sqrt(1 + s^2 f(t)^2) dt .
 * The scale s is not cosmetic. Arc length adds a length to a density and so is not scale invariant;
 * on a circle of circumference 2 pi a density is of order 1/(2 pi), the constant swamps the base and
 * the transform returns the circular uniform. With s restored the family interpolates: s -> 0 gives
 * the uniform exactly and s -> infinity returns the base, monotonically in concentration, and the
 * mode is preserved at every s because t -> sqrt(1 + s^2 t^2) is strictly increasing on t >= 0.
 * The base here is von Mises with concentration kappa and mean direction mu. */
void arcc_temper_vm(const double *theta, const int *nth, const double *kappa, const double *mu,
                    const double *s, const int *nodes, double *out);

/* Trigonometric moment of order p of the tempered law; out has length 2 (real, imaginary). */
void arcc_temper_vm_trigmom(const int *p, const double *kappa, const double *mu, const double *s,
                            const int *nodes, double *out);

#endif
void arcc_exact_ci(const double *theta, const int *n, const int *B, const int *seed,
                   const double *cgrid, const int *ng, const int *stat, const int *group,
                   const double *lev, double *out);
void arcc_gof(const double *theta, const int *n, const double *c3, const double *mu,
              const int *B, const int *seed, double *out);
