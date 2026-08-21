## Test environments

* local macOS (arm64), R release
* R CMD check --as-cran on the built tarball

## R CMD check results

There were no ERRORs or WARNINGs.

## Notes

This is a new submission.

The package contains no compiled dependencies beyond the C standard library and
libm. OpenMP is used where available and the package builds and runs correctly
without it. All random number generation is by a self-contained splitmix64
stream inside the C back end, seeded explicitly, so results are reproducible
across platforms and independent of the R RNG state.


## On URLs

The package URL and BugReports fields are omitted from this submission because
the public repository is not yet created; they will be added once it is, rather
than submitted pointing at a location that does not resolve.
