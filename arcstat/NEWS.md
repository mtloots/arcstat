# arcstat (development)

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
