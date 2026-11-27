# arcstat 0.4.1

## Bug fixes

* `k4eq_boot()` refitted every replicate WITH a linear drift term, whatever the point fit was: the
  setting of the point fit was not passed to the replicates (the back end held the flag at "drift" for
  every bootstrap). For a fit made with `k4eq_fit(drift = FALSE)` the replicates were therefore fitted
  by a different model from the one estimated, and the bootstrap spread of the reading was inflated
  (on two induction curves of 600 points, standard errors of 0.262 and 0.151 hours where the
  replicates refitted without drift give 0.138 and 0.103). Fits made with drift, the default, are not
  affected: their replicates are unchanged to the last bit.
* `k4_boot()` had the same fault: its replicates were always refitted with drift, so the bootstrap of a
  `k4_fit(drift = FALSE)` fit was wrong in the same way. Fits with drift are not affected.
* Both bootstraps now take a `drift` argument, by default read off the fit (a fit made without drift has
  a slope `m` of exactly zero), and refuse `drift = FALSE` for a fit that carries a slope. Their output
  gains a last column `m`, the fitted slope of each replicate, so that the setting can be checked.
* `k4_boot()` gains `gate`: `"best"` keeps a replicate only if its refit is physical, the rule of a point
  fit made with `gate = "best"`; and `gev` is now read from the `gev` column of a `k4_fit()` result when
  it is not given.
* `eqfit_bc_boot()` refuses a point fit that does not lie on the equivalence curve it is given: its
  replicates are always refitted on the curve, so a free fit would have been bootstrapped as a
  constrained one.
* `bb_evidence()` refuses a reference from `bb_ref_disc()` made for another sample size (the reference now
  carries its `n`). In Python its threshold is now R's type-7 quantile, so that the two fronts return the
  same evidence (they differed when the level fell between order statistics of the reference).
* The parallel loop of the exact interval of the circular family now honours `arc_threads()`.
* New tests check, for every bootstrap, that a replicate is fitted under the settings of its point fit
  (drift, submodel, box or bounds, locus piece and gate): a block as long as the series makes the
  replicate the data themselves, whose refit must return the point fit.

# arcstat 0.4.0

## New features

* Kappa induction curves with drift, many runs in one call: `k4_fit()` (free fits, with options for
  holding h, the h = 0 submodel and acceptance of physical fits only), `k4eq_fit()` (fits on the
  equivalence locus, whole or tabulated), `k4_locus_k()` (the locus), `k4_boot()` and `k4eq_boot()`
  (residual bootstraps), `k4_shared_fit()` (runs sharing one shape) and `k4_shared_perm()` (the
  permutation test of a shared shape against a shape per group).
* `wald_fit()`, `wald_cdf()` and `wald_readings()` for inverse Gaussian curves.
* `bc_host_fit()` for beta-companion curves with drift, free and on the equivalence curve; the dense
  search is the default.
* `k4_lmom_window()` for L-moment matching over the observed window; `k4_mv_share()` for the
  variance components of a multiverse of fitting pipelines.
* `arc_threads()` sets the number of threads of every parallel loop (option `arcstat.threads` or
  environment variable `ARCSTAT_THREADS`; two by default). Results do not depend on it.

## Changes

* `k4_readings()`, `eq_readings()`, `eq_readings_vsl()` and `k4_tau34()` are computed in closed form
  or by quadrature to about 1e-15; readings on heavy tails and near the ends of the unit interval
  change.
* `bc_q()`, `bc_cdf()` and `bc_pdf()` keep their accuracy as beta approaches -1; `eq_bstar()` is
  defined up to the end of the curve.
* The companion curve of `eqfit_bc()`, `eqfit_bc_boot()`, `eqfit_score_at()`, `eqfit_msim()` and
  `eqfit_nullT()` is solved at every point instead of interpolated, and the minimiser of the last two
  is a global search; their values change and `maxit` has no effect there.
* `k4_mv_boot()` searches h on both sides of zero in the box it is given, refits every replicate from
  its pipeline's own point fit, and has `drift`.
* `k4_fit_aleq()` integrates bands at the ends of the unit interval exactly; `qmin_rho()` and
  `admiss_rho()` refine the grid minimum.
* The Python package is built with OpenMP where the compiler supports it.

# arcstat 0.3.1

* `al_pvalue()` could return a **negative** probability far in the right tail. The Lugannani-Rice
  tail is an asymptotic expansion and its correction term overshoots there: at n = 150 on badly
  non-uniform samples the returned value was -6.5e-29 and -2.1e-89. Numerically zero, and no
  decision at any conventional level was affected, but a probability must lie in [0, 1] and callers
  comparing against a level must never see a negative number. The value is now clamped to [0, 1] in
  the C back-end, so both the R and the Python bindings inherit the fix. Found 17 September 2026
  while using the test as a calibration instrument for a time-series language model. A regression
  test now checks the tail directly.

# arcstat 0.3.0

* **Breaking, and it corrects a wrong answer.** `al_scale()` now matches the sample band arc length
  to the spacing-corrected model band length `al_band_model_star()` rather than to the population
  band length `al_band_model()`. The old behaviour was inconsistent: order statistics inside the
  band are spaced like Q'(u)E/n with E standard exponential, so the polygonal sum converges to

      S*(sigma) = int E sqrt(f0(z)^2 + sigma^2 E^2) dz = E S(sigma E),

  and not to S(sigma). The two differ by a constant factor, because the square root is nonlinear
  and E[E^2] = 2 rather than 1, and the gap does NOT close as the sample grows. At the standard
  normal on the default band the old estimator returned about 1.05 when the scale was 1, at every
  sample size. `al_scale()` now returns 1.001 at n = 20000 and its bias shrinks with n.

  Estimates of SCALE change. Regression slopes computed from this estimator do not, because a
  multiplicative bias in the level cancels in a log-linear slope; interval coverage does change,
  and becomes correct.

* New `al_band_model_star()`, the spacing-corrected model band arc length, on both fronts. Its
  exponential expectation uses an 80-node Gauss-Legendre rule under the map t = s/(1 - s), which
  agrees with adaptive double integration to 3e-13 across sigma from 0.1 to 10 and across four
  bands. A Gauss-Laguerre rule was tried first and rejected: it reaches only 2e-5 at sigma = 10,
  because the integrand has a kink at t of order f/sigma that falls below its first node.

* New `al_scale_raw()`, the uncorrected estimator, retained only so that the size of the
  correction can be measured. It is inconsistent and is documented as such.

# arcstat 0.2.0

* **Breaking.** `k4_fit_varpro()` now searches the fourth shape parameter `h`
  DIRECTLY rather than on the log scale, and the sentinel asking for the free
  four-parameter fit changes from a negative `hfix` to `NA` (`NaN` in C and in
  Python). Column four of `starts` now carries `h` itself, not `log h`.

  This is not a refactor. On the log scale `h` was necessarily positive and the
  boundary at `h = 0` lay at minus infinity, so the free fit could not reach the
  negative half-line at all and crept toward the boundary reporting a spurious
  small `h` instead of naming the submodel it had effectively chosen. The
  negative half-line holds real members -- the log-logistic at `h = -1`, the
  Burr III and five-parameter logistic slice at `h = -1/m` -- and an audit of the
  edible-oil traces found curves that genuinely prefer it. `h = 0`, the
  generalised extreme value member, is now an interior point of the free search.

  **It fails silently if you do not act.** Code that passed `hfix = -1` to mean
  "fit h freely" now gets `h` HELD at -1, which fits the log-logistic and returns
  a plausible number. Replace any negative `hfix` meant as a free-fit request
  with `NA_real_` in R, or `hfix=None` in Python. Passing `hfix = 0` for the
  generalised extreme value submodel is unaffected, as is `NA`.

  Both fronts move together: the shared C back end, the R wrapper and the Python
  wrapper all take the new sentinel, and the Python front had been sending -1.0
  internally for a free fit.

* The quantile-family routines `arck4_q()`, `arck4_tau34()`, `arck4_fit_lmom()`,
  `arceq2_bc_q()` and `arceq2_bc_pdf()` are now registered as C callables and
  declared in `inst/include/arcstatAPI.h`, so a dependent package can call them
  from its own C rather than through R. The registration is additive: nothing
  already exported changes, no numerics are affected, and dynamic symbol lookup
  stays disabled, which is why each entry point is registered explicitly.

* `arck4_readings()` now requires the density's grid maximum to be INTERIOR before
  reporting a mode. For `k >= 1` with `h < 1` both terms of
  `(log q)'(u) = (h-1)/u - (k-1) h u^(h-1)/(1-u^h)` are negative, so the quantile
  density is strictly decreasing, the density is strictly increasing, and the
  argmax lands on the support endpoint. The routine previously returned that
  endpoint as a mode together with two readings computed from it -- finite,
  plausible-looking, and meaningless. Both readings are now `NaN` there, which is
  correct: neither is defined when the curve has no interior steepest point.
  Found while re-deriving the equivalence locus of the four-parameter kappa
  family, where the artefact had been read as a second branch of the locus.

# arcstat 0.1.0

* First release. The computation runs on a shared pure-C back-end that is
  also bound from Python, and the two fronts are checked against each other
  value by value.
