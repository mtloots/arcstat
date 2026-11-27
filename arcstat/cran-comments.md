## Submission

arcstat 0.4.1 is a bug-fix release that follows 0.4.0. Two bootstrap functions, `k4eq_boot()` and
`k4_boot()`, refitted every replicate with a linear drift term even when the point fit had been made
without one, so their standard errors for such fits were too large. Fits with drift (the default) are
unchanged. Both functions now take the drift setting of the point fit; the details are in NEWS.md.

## Test environments

* local macOS 26 (arm64), R 4.6.1, R CMD check --as-cran
* win-builder, R-devel
* GitHub Actions: macOS, Windows and Ubuntu (R release), Ubuntu (R devel)

## R CMD check results

0 errors | 0 warnings | 1 note

* The note is the CRAN incoming feasibility check, which lists "Hosking" as a possibly misspelled
  word in DESCRIPTION. It is the surname of the author of the four-parameter kappa distribution,
  cited there with its DOI. It may also report the number of days since the last update.
* Locally the HTML manual check also notes that math rendering was skipped because 'V8' is not
  installed; win-builder checks the HTML manual without a note.

## Reverse dependencies

There are no reverse dependencies.
