# ellstat 0.1.0

First release. Implements the elliptic moment system for circular data:

* `ell_basis()`, `ell_moments()` --- the modulus-indexed basis and sample moments.
* `ell_moment_var()` --- the exact finite-sample identity
  `tr Var(hat E_p) = (1 - |E_p|^2)/n`, which holds for every sample size.
* `ell_transfer()` --- the triangular transfer to trigonometric moments that
  establishes characterisation.
* `ell_exact_uniform()`, `ell_exact_ci()` --- exact inference by rotation
  invariance and test inversion, valid for every sample size and every number
  of replicates.
* `dellvm()`, `rellvm()`, `fit_ellvm()` --- the elliptic von Mises family.
* `ell_dn2_moments()` --- closed-form coefficients of the arc-length law on an
  ellipse, the diagonal member of the system.
* `ell_efficiency()` --- the matched-modulus efficiency study.

All computation is in C. The same source backs the Python package `ellstat`;
`parity_ellstat.sh` checks that the two fronts agree to the last bit.
