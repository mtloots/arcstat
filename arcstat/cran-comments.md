## Submission

arcstat 0.3.1 is a bug-fix release, four weeks after 0.3.0.

`al_pvalue()` evaluates a Lugannani-Rice saddlepoint tail. Far in the right tail the expansion's
correction term overshoots, and on badly non-uniform samples (n = 150) it returned -6.5e-29 and
-2.1e-89 rather than a value in [0, 1]. No decision at any conventional level changes, but a
probability must not be negative, so the value is now clamped to [0, 1] in the C back-end, with a
regression test. Nothing else changes; the release is made promptly so that users comparing
p-values against a level never see a negative number.

## Test environments

* local macOS (arm64), R 4.5, R CMD check --as-cran

## R CMD check results

See below; no ERRORs or WARNINGs.

## Reverse dependencies

There are no reverse dependencies.
