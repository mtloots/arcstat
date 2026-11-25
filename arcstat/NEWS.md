# arcstat (development version 0.3.1.9004)

* New `wald_fit()`, `wald_cdf()` and `wald_readings()`: inverse Gaussian curves with drift, many runs
  in one call, and their two readings without a grid (exact for every phi, the Levy limit included).

* `k4_fit()` and `k4eq_fit()` have `drift`; `FALSE` holds the slope at zero.

* New `k4_lmom_window()`: L-moment matching of kappa shapes over the observed window, both signs of h,
  many runs in one call. Two windowed ratios can be matched exactly by more than one shape; the number
  of matches is returned.

# arcstat (development version 0.3.1.9003)

* New `k4_fit()`: free kappa fits with drift, many runs in one call, by Levenberg-Marquardt from a
  built-in start grid wider than the one used with `k4_fit_varpro()`, with options for holding h, for
  the h = 0 submodel selection and for the rule on physical fits. `gate = "best"` is the lowest
  physical minimum. A ladder of floors on h (`gate = TRUE`) depends on how well the search does: a
  better search finds lower minima that are not physical and steps past floors a weaker one accepts.

* `k4_mv_boot()` has `drift`: `med9-LS`, `med21-LS` and `arc` can carry the linear drift the other
  pipelines always carried (`FALSE`, as before, is the default). Its point fits use the search of
  `k4_fit()`.

* `bc_host_fit()`: `curve = NULL` is the whole equivalence curve, solved exactly; `eps` and `dense`.

# arcstat (development version 0.3.1.9002)

* **Number of threads.** Every parallel loop now uses one setting: the option `arcstat.threads`,
  otherwise the environment variable `ARCSTAT_THREADS`, otherwise at most two as before.
  `arc_threads()` reads and sets it. Results are bit-identical whatever the number (tested on 1, 2
  and 8). In Python, `set_threads()`, `get_threads()` and the environment variable; the Python back
  end is now built with OpenMP where the compiler has it.

* New `k4eq_fit()`: kappa fits with drift constrained to the equivalence locus, many runs in one
  call, with the linear part solved exactly and the starts of all runs sharing the threads. By default
  on the whole locus, computed inside; alternatively on tables of it. No fit is refused on physical
  grounds; the best fit on every piece of the locus is returned as well, so that a caller can refuse.
  A linearly interpolated table puts kinks in the objective, at which a descent can stop at a local
  minimum; every node and the cells beside the best nodes are searched.

* New `k4_locus_k()`: the equivalence locus itself, to about 1e-15, over h from -infinity to 1/2.

* New `bc_host_fit()`: beta-companion fits with drift, free and constrained to the equivalence curve,
  many runs in one call. The curve is exactly 0 and 1 off the support, the exponents may lie in
  [-1 + 1e-6, 20], and the free fit also starts from the constrained one.

* New `k4_boot()` and `k4eq_boot()`: moving-block (or single) residual bootstraps of the free and the
  constrained kappa fit, many runs in one call, each run in its own box and with h held at zero where
  its point fit holds it. Replicate b of a run depends only on the run's seed and b.

* **`k4_mv_boot()` could not reach a curve with negative h.** Its start grid had positive h only, and
  the presmoothed and banded pipelines searched log h from 0.02 up, whatever the box. The start grid now
  covers both signs, all eight pipelines search h directly in the box `bounds`, the presmoothed
  pipelines solve their linear part exactly, and a replicate starts from its pipeline's own point fit
  rather than from the base fit. Point estimates and replicates change; by how much depends on the
  share of curves with negative h.

# arcstat (development version 0.3.1.9001)

* **`k4_readings()` returned tangent readings that were too low on heavy right tails**, and with it
  every routine that reads a kappa curve: `eqfit_k4eq()`, `eqfit_k4eq_boot()`, `k4_mv_boot()`,
  `k4_b_sweep()` and `k4_ab_sweep()`. The steepest point was taken from a grid between the 0.1 and
  99.9 per cent quantiles; with `k` near -1 that range is so long that the steepest point falls
  between two nodes. On 97 fitted oil shapes the tangent reading moves up by more than 0.1 per cent
  on 45 and by more than 1 per cent on 26, at most by 15.8 per cent, and never down. Both readings
  are now closed forms: in r = u^h / w, w = (1 - u^h)/h, the steepest point is
  r = (1 - h)/(1 - k) and the third-derivative reading is a root of a quadratic. They agree with
  50-digit values to about 1e-15 relative over k from -2.5 to 1.2 and h from -60 to 2, including the
  limits. The `grid` argument is kept and ignored. The existence guards are now the exact conditions
  (`a` and `mode`: k < 1, h < 1, hk < 1; `b`: also h < 1/2, hk < 1/2); readings that used to be `NaN`
  only because they lay outside the 0.1 to 99.9 per cent range, or below the rounding threshold of
  the differenced density, are now returned (`b` on 6 of the 97 oil shapes).

* **`k4_mv_boot()` did not refit the shape on its replicates.** The replicate start carried log h
  where the variable-projection fit takes h, so for h < 1 it lay outside the box, the free fit never
  moved, and the selection fell back on the h = 0 submodel: the replicate columns `LS`,
  `transient-excised`, `grid-odd` and `grid-even` were those of `GEV-LS`. Present since h has been
  searched directly (27 August 2026). Point estimates `a_pt` were not affected by this.

* `eqfit_msim()`, `eqfit_nullT()`, `eqfit_score_at()`, `eqfit_bc()` and `eqfit_bc_boot()` evaluated
  the companion curve by linear interpolation of its inverse on 2049 levels: up to 1.7e-4 wrong in
  F, and, differentiated numerically, a step function, so that the score statistic was rough on
  every scale and its apparent minima were often artefacts of the interpolation (0.057 at a point
  where the exact statistic is 0.424). The curve is now solved to machine precision at each point,
  the Jacobian uses exact derivatives in location and scale, and the projection is an
  orthogonalisation without a ridge. One evaluation of the statistic is about four times faster.
  **Every value of the score statistic changes**; see the next item for the minimum.

* The minimiser behind `eqfit_msim()` and `eqfit_nullT()` is a global search: a dense grid over
  location, right end of the support and alpha, least-squares starts, descents from every grid-local
  minimum, tight polishing, and an explicit search of the faces on which a design point leaves the
  support, where the statistic is discontinuous. It costs about 7 seconds of processor time per
  draw at n = 400, against about 17 before. On 50 null draws at each of alpha = -0.52 and -0.60
  it was never above a separate, about ten times more expensive search by more than 1e-9, and on 40
  of each never above a second one. Values are mostly LARGER than in 0.3.1 and 0.3.1.9000, whose low
  values were interpolation artefacts: higher at 30 of 40 draws at -0.52 and 40 of 50 at -0.60. `maxit` no longer has an effect in these two
  functions.

* `k4_tau34()` (and `k4_fit_lmom()`, which uses it) used a Gauss rule on an integrand singular at
  both ends: errors in tau_3 of 2e-6 at (0.2, 0.4), 3e-3 at (-0.5, 0.3) and up to 0.1 near k = -1. It
  now uses a transformed trapezoid that is exact to 1e-15 and returns `NaN` where the moments do not
  exist.

* `k4_fit_aleq()` integrated a band that starts at 0 with a 40-node Gauss rule across the singularity
  of the quantile density: model band lengths wrong by 5 per cent at h = 0.3 and 25 per cent at
  h = 0.1. Bands within 0.05 of either end now use a double-exponential rule. Interior bands are
  unchanged.

* `bc_q()`, `bc_cdf()` and `bc_pdf()` lost accuracy as beta approaches -1 (relative errors 2e-3 in
  `bc_q()` at beta = -1 + 1e-12): the incomplete beta is now computed unscaled on the logit scale.
  `eq_bstar()` no longer returns `NaN` when the root lies within 0.001 of beta = -1. `eq_readings()`
  without a tilt is in closed form (1e-15 instead of 2e-9), and `eq_readings_vsl()` finds derivative
  readings below u = 0.001.

* `qmin_rho()` refines the grid minimum (it was up to 2e-6 too high), and with it `admiss_rho()`.

# arcstat (development version 0.3.1.9000)

* **`eq_readings()` returned readings that were too low, and so did `eqfit_bc()` and
  `eqfit_bc_boot()`, which use it.** The routine integrates the quantile density from a grid that
  starts at 1e-5 and used to set the integral to zero there. For a quantile density
  u^alpha (1-u)^beta with alpha between -1 and -1/2 the mass below the first node is not small
  (0.025 at alpha = -0.6), and the trapezoid rule overshoots on the first cells of a steep
  singular curve. Both the tangent reading `a` and the derivative reading `b` were therefore too
  low, by 1 to 5 per cent for alpha between -0.6435 and -0.505 and beta between -1 and the
  equivalence curve; the gap `D = b - a` was almost unaffected. The integral now starts at zero:
  the mass below the first node is added, the singular power factor is integrated exactly on
  every cell, and the roots `uc` and `ub` are located between neighbouring nodes instead of being
  snapped to the grid. Against 50-digit values the readings are now right to about 2e-9 relative
  on the whole range, including beta = -1 + 1e-9. The observation-domain readings that
  `eqfit_bc()` returns (`a`, `b`) and the bootstrap vector of `eqfit_bc_boot()` move with them.
  `eq_readings_vsl()` had the same omission for lambda below 1 and now uses the closed-form integral.
  The kappa readings (`k4_readings()`) use a different routine; see 0.3.1.9001 for theirs.

* The minimiser behind `eqfit_msim()` and `eqfit_nullT()` often stopped above the minimum of its
  own objective, so the statistic was too large at those draws: on 40 null draws at alpha = -0.52
  a restarted search found a lower value in 19, by up to 0.36. The five starts used to differ in
  one coordinate only and nearly always ended at the same point. After them the search now
  restarts from the best point along all three coordinates, at three scales, so the returned
  minimum can only be lower than before. Statistics from these two functions can therefore be
  smaller than in 0.3.1 at the same seed. The cost is about twice.

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
