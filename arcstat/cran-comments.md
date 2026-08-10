## Test environments

* local macOS (arm64), R 4.6.1
* GitHub Actions, ubuntu-latest, R release
* `R CMD check --as-cran` on the built tarball

## R CMD check results

0 errors | 0 warnings | 1 note

The note is "New submission", which is expected for a first submission.

## Notes for the reviewer

* This is a first submission.
* The package contains compiled C. It uses only libm and R's own headers, has no system
  requirements beyond a C compiler, and registers its native routines with
  `R_registerRoutines` and `R_useDynamicSymbols(dll, FALSE)`.
* The same C sources back a Python package of the same name. A parity harness in the
  repository compares every exported quantity from both front ends and requires the
  printed values to be identical; the tests in this package encode the mathematical
  identities the routines must satisfy, including two closed forms
  (`cf_arclength_family("normal")` is exactly 2 and `("exponential")` exactly pi).
* No package uses more than two cores at any point in examples, tests or vignettes.
