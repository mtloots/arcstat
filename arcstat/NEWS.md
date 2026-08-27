# arcstat (development)

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
