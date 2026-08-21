# ellstat

The elliptic moment system for circular data.

Trigonometric moments are a single fixed system. Replacing the circular functions
by the Jacobi elliptic functions of modulus `k` gives a one-parameter family of
moment systems that reduces to the trigonometric system exactly at `k = 0`,
characterises every circular distribution at every modulus, and carries exact
finite-sample first- and second-order theory.

```r
library(ellstat)
th <- rellvm(200, mu = 1, kappa = 3, k = 0.9, seed = 42)
ell_moment_var(th, P = 3L, k = 0.9)   # realised and exact variance agree identically
ell_exact_uniform(th, k = 0.9)        # exact for every n and every B
fit_ellvm(th, k = 0.9)
```

Computation is in C. The same source backs the Python package `ellstat`, and
`parity_ellstat.sh` checks the two fronts agree to the last bit.
