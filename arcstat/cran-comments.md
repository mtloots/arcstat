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

* The note is the CRAN incoming feasibility check, which lists "Hosking" as a possibly misspelled
  word in DESCRIPTION. It is the surname of the author of the four-parameter kappa distribution,
  cited there with its DOI.
* Locally the HTML manual check also notes that math rendering was skipped because 'V8' is not
  installed; win-builder checks the HTML manual without a note.

## Reverse dependencies

There are no reverse dependencies.
