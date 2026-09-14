## Submission of 0.3.0, a correction

This update follows 0.2.0 within a day, and I am sorry to come back so soon. The reason is that
0.2.0 returns a wrong number, and I would rather say so than let it stand.

**What was wrong.** `al_scale()` matched the sample band arc length to the population band arc
length. Those are not the same quantity. Order statistics inside the band are spaced like
`Q'(u)E/n` with `E` standard exponential, so the polygonal sum converges to

    S*(sigma) = int E sqrt(f0(z)^2 + sigma^2 E^2) dz,

not to `S(sigma)`. The square root is not linear and `E[E^2] = 2`, so the two differ by a constant
factor that does not vanish as the sample grows. For a standard normal response on the default band
the old estimator returned about 1.05 when the scale was 1, at every sample size. It is an
inconsistent estimator, not a small-sample effect.

**What changed.** `al_scale()` now matches to the spacing-corrected `al_band_model_star()`, which is
new and exported. At n = 20000 on the same example it returns 1.001, and its bias shrinks with n.
`al_scale_raw()` is also new, and keeps the old behaviour under a name that says what it is, so that
the size of the correction can be measured rather than merely asserted.

Estimates of scale from this package change. Regression slopes computed from it do not, because a
multiplicative bias in the level cancels in a log-linear slope. Interval coverage changes, and
becomes correct. NEWS.md states this in the same terms.

**Reverse dependencies.** None; the package has been on CRAN less than a day.

## Test environments

* local: macOS 26.5, R 4.6.1, Apple clang 21.0.0
* `R CMD check --as-cran` on the tarball: Status 2 NOTEs, both described below

## R CMD check results

Two NOTEs, no ERRORs and no WARNINGs.

* *Days since last update: 0.* This is the update described above. I judged a correct answer to be
  worth the short interval, and I will not submit again soon.
* *Skipping checking HTML validation / math rendering.* HTML Tidy on this machine is too old and the
  'V8' package is unavailable locally. Neither reflects package content.

## On the DESCRIPTION reference

Still none to give, and the position has changed since 0.2.0. The manuscript describing these
methods was declined by the Journal of Statistical Planning and Inference; the work has since been
merged into a single paper, which corrects exactly the defect this release fixes, and is being
prepared for the Electronic Journal of Statistics. There is no DOI to cite yet. I will add the
citation as soon as there is one.
