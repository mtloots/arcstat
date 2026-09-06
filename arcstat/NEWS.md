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
