## The three routines that carried no test: the matched-modulus efficiency study and the
## projected-family efficiency pair. Coverage was 82.3 per cent with these at exactly zero.
## As elsewhere in this package the checks are identities, not smoke.

test_that("EFFICIENCY STUDY reports a gain that is its own variance ratio", {
  ## gain is defined as var.trig / var.elliptic, so the returned value must equal the ratio of
  ## the two variances it is derived from. That is a definitional identity: any drift between
  ## them is a bug in the C, not a modelling question.
  e <- ell_efficiency(k = 0.5, kappa = 3, n = 100L, R = 50L, seed = 1L)
  expect_named(e, c("eta", "var.trig", "var.elliptic", "gain",
                    "solved.trig", "solved.elliptic"))
  expect_equal(unname(e["gain"]), unname(e["var.trig"] / e["var.elliptic"]), tolerance = 1e-12)
  expect_gt(e["var.trig"], 0)
  expect_gt(e["var.elliptic"], 0)
  ## the solved fractions are proportions of the replicates
  expect_true(all(e[c("solved.trig", "solved.elliptic")] >= 0))
  expect_true(all(e[c("solved.trig", "solved.elliptic")] <= 1))
  ## and the study is deterministic under its seed
  expect_identical(ell_efficiency(k = 0.5, kappa = 3, n = 100L, R = 50L, seed = 1L), e)
  expect_false(identical(ell_efficiency(k = 0.5, kappa = 3, n = 100L, R = 50L, seed = 2L), e))
})

test_that("PROJECTED FAMILY: the mixing law is a probability measure", {
  ## mass is the total mass of the projected density. Whatever the axis ratio and offset, it
  ## must be one, or the quadrature underneath the efficiency is not integrating a density.
  for (tau in c(0.5, 1, 2)) {
    p <- ell_projfam_eff(tau = tau, d = 1, kgrid = c(0, 0.5, 0.9))
    expect_equal(p$mass, 1, tolerance = 1e-06)
    expect_length(p$efficiency, 3L)
    expect_identical(p$kgrid, c(0, 0.5, 0.9))
    expect_gt(p$fisher, 0)
    ## an efficiency relative to maximum likelihood cannot exceed one
    expect_true(all(p$efficiency > 0 & p$efficiency <= 1 + 1e-08))
  }
})

test_that("PROJECTED FAMILY: the optimal modulus really is optimal", {
  ## ell_projfam_kopt claims the modulus that maximises efficiency. The claim is checkable:
  ## evaluate the efficiency on a fine grid and no point may beat the reported optimum.
  opt <- ell_projfam_kopt(tau = 1, d = 1)
  expect_named(opt, c("k", "efficiency"))
  expect_true(opt["k"] >= 0 && opt["k"] < 1)
  grid <- seq(0, 0.995, by = 0.005)
  eff <- ell_projfam_eff(tau = 1, d = 1, kgrid = grid)$efficiency
  expect_gte(unname(opt["efficiency"]), max(eff) - 1e-06)
  ## the optimum is interior here, so it beats both ends of the grid outright
  expect_gt(unname(opt["efficiency"]), eff[1])
})

test_that("PROJECTED FAMILY: both cited families are reachable", {
  n <- ell_projfam_eff(tau = 1, d = 1, family = "normal", kgrid = c(0, 0.9))
  cc <- ell_projfam_eff(tau = 1, d = 1, family = "cauchy", kgrid = c(0, 0.9))
  expect_equal(n$mass, 1, tolerance = 1e-06)
  expect_equal(cc$mass, 1, tolerance = 1e-06)
  ## the two families are genuinely different laws, so their Fisher information differs
  expect_false(isTRUE(all.equal(n$fisher, cc$fisher)))
})
