## Resubmission

This is a resubmission of a first submission, following Leonore Hochhauser's review of 0.1.0 on
31 August 2026. Both points raised have been addressed.

**`\value` tags.** `\value` was missing from `arcc_gof.Rd` and `eqfit_bc.Rd`. Both now document
the class and structure of what is returned and what it means. `eqfit_bc.Rd` documents thirteen
exported functions that share one help page, so its `\value` names each of them and the object it
returns. A sweep of the package confirms no other `.Rd` for an exported function is missing the tag.

**References in DESCRIPTION.** There are none to give. The methods in this package are described in
a manuscript currently under review at the Journal of Statistical Planning and Inference
(JSPI-D-26-00392), and in a University of Pretoria PhD thesis from which nothing was published as an
article. Neither carries a DOI I can point to, and I would rather cite nothing than cite something a
reader cannot reach. The CRAN cookbook notes that a reference in the Description field is optional,
so the field is unchanged in that respect. I will add the citation as soon as the paper has a DOI.

I also took the opportunity to single-quote the software name 'Python' in the Description, which the
cookbook asks for and which the earlier version did not do.

**Version.** This resubmission is 0.2.0 rather than 0.1.0. The additional change is a breaking
reparameterisation of the fourth shape parameter in `k4_fit_varpro()`, searched directly rather than
on the log scale, which is documented at the head of NEWS.md. The log scale confined that parameter
to the positive half-line, so the free fit could not reach members that genuinely lie on the negative
one. Submitting 0.2.0 avoids publishing semantics I would have to deprecate immediately.

## Test environments

* local macOS 26 (arm64), R 4.6.1, `R CMD check --as-cran` on the built tarball
* GitHub Actions, R release on ubuntu-latest
* win-builder, R-devel, on this exact tarball

NOTE TO SELF, 31 Aug 2026: the earlier version of this file claimed GitHub Actions coverage of
R-devel, R oldrel-1, macOS and Windows. THAT WAS NOT TRUE FOR 0.2.0. The workflow's CRAN-like
matrix job carries `if: github.event_name == 'workflow_dispatch'`, so on an ordinary push only the
quick ubuntu-release job runs, and the matrix has to be started by hand. Before claiming multi-
platform coverage again, either dispatch that workflow and wait for it, or cite win-builder.

## R CMD check results

0 errors | 0 warnings | 1 note

The note is "New submission", which is expected.

## Notes for the reviewer

* The package contains compiled C. It uses only libm and R's own headers, has no system
  requirements beyond a C compiler, and registers its native routines with
  `R_registerRoutines` and `R_useDynamicSymbols(dll, FALSE)`.
* The same C sources back a 'Python' package of the same name. A parity harness in the
  repository compares every exported quantity from both front ends and requires the printed
  values to be identical; the tests here encode the mathematical identities the routines must
  satisfy, including two closed forms (`cf_arclength_family("normal")` is exactly 2 and
  `("exponential")` exactly pi).
* No example, test or vignette uses more than two cores.
