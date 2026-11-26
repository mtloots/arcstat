## Submission

arcstat 0.4.0 adds batched fits, bootstraps and a permutation test for induction curves
(`k4_fit()`, `k4eq_fit()`, `k4_shared_fit()`, `k4_shared_perm()`, `k4_boot()`, `k4eq_boot()`,
`wald_fit()`, `bc_host_fit()`), parallel by OpenMP with a thread count set by `arc_threads()`
(default two threads; results do not depend on the number). Every exported function now has an
example, and the slowest runs in about two seconds.

## Test environments

* local macOS 26 (arm64), R 4.6.1, R CMD check --as-cran
* win-builder, R-devel
* GitHub Actions: macOS, Windows and Ubuntu (R release), Ubuntu (R devel)

## R CMD check results

0 errors | 0 warnings | 1 note

* The note is the CRAN incoming feasibility check.

## Reverse dependencies

There are no reverse dependencies.
