## The circular family, the arc-length generator, the fitting front ends and the equivalence
## simulation layer. These carried the lowest coverage in the package after test-backends.R:
## arccirc.c at 15.6 per cent, arceqfit.c at 14.5, arcdistc.c at 38.3. As elsewhere the checks
## are identities and invariances, not smoke.

test_that("CIRCULAR family: the uniform member has no resultant and concentration rises with c3", {
  ## c3 = 0 is the uniform circular law, whose mean resultant length is exactly zero. Increasing
  ## c3 concentrates the law, so rho must rise monotonically, and it cannot leave [0, 1].
  expect_lt(abs(rho_arccirc(arccirc(0, 0))), 1e-12)
  rhos <- vapply(c(0, 0.3, 0.6, 0.9, 1.03), function(c) rho_arccirc(arccirc(c, 0)), numeric(1))
  expect_true(all(diff(rhos) > 0))
  expect_true(all(rhos >= 0 & rhos <= 1))
  ## c3 is bounded by the family's own admissible maximum
  expect_gt(arcc_c3max(), 1)
})

test_that("CIRCULAR family: draws lie on the circle and the fit recovers a concentrated member", {
  set.seed(9)
  o <- arccirc(c3 = 1.0, mu = 1.0)
  s <- rarccirc(4000L, o)
  expect_length(s, 4000L)
  expect_true(all(s >= 0 & s < 2 * pi))
  ## the quantile density of a proper circular quantile function is strictly positive
  expect_true(all(qd_arccirc(seq(0.01, 0.99, by = 0.01), o) > 0))
  f <- fit_arccirc(s)
  expect_named(f, c("mu", "c3", "rho", "interior"))
  expect_lt(abs(f$mu - 1.0), 0.1)
  expect_lt(abs(f$c3 - 1.0), 0.1)
})

test_that("CIRCULAR Fourier form: a null first moment is the uniform law", {
  ## arccirc_fr is parameterised by Fourier coefficients. With a single coefficient carrying no
  ## first trigonometric moment the density is flat at 1/(2 pi) and the quantile function is the
  ## identity on [0, 1] scaled to the circle.
  fr <- arccirc_fr(c(0.3 + 0.1i))
  expect_lt(Mod(trigmom_arccirc_fr(1L, fr)), 1e-08)
  expect_equal(darccirc_fr(1.0, fr), 1/(2*pi), tolerance = 1e-06)
  expect_equal(Q_arccirc_fr(0.5, fr), 0.5, tolerance = 1e-06)
})

test_that("CIRCULAR moments: rho factorises and the tempered von Mises is a density", {
  ## factorise_rho returns Fejer-Riesz coefficients for a requested mean resultant length. The
  ## coefficients are normalised, which is exact and is asserted here; the achieved resultant is
  ## checked only in the small-rho regime.
  ## NOTE, 21 Aug 2026: the achieved resultant SATURATES at 0.440051 for every request of 0.5 or
  ## more -- identical for 0.5, 0.7 and 0.9. That is either the true admissible ceiling of the
  ## degree-one family or a silent clamp, and it is deliberately NOT asserted either way here
  ## until the intended behaviour is confirmed.
  for (r in c(0.2, 0.5, 0.9)) expect_equal(sum(Mod(factorise_rho(r))^2), 1, tolerance = 1e-10)
  expect_lt(abs(Mod(trigmom_arccirc_fr(1L, arccirc_fr(factorise_rho(0.2)))) - 0.2), 0.01)
  ## and the tempered von Mises must integrate to one over the circle
  expect_equal(integrate(function(t) dtemper_vm(t, kappa = 1, mu = 0, s = 1), 0, 2*pi)$value,
               1, tolerance = 1e-06)
  ## circular L-moments of a sample are complex and finite
  set.seed(21)
  l <- lmom_circ(rarccirc(500L, arccirc(0.8, 0.3)), 3L)
  expect_length(l, 3L)
  expect_true(all(is.finite(Mod(unlist(l)))))
})

test_that("ARC-LENGTH GENERATOR builds a proper distribution", {
  ## The generator turns the arc length of a density curve into a distribution. Whatever it
  ## produces must still be a distribution: the cumulative function runs from 0 to 1, never
  ## decreases, and the density it carries is non-negative.
  g <- arc_generator(function(z) dnorm(z), -6, 6)
  expect_true(all(c("x", "G", "g", "total_arclength") %in% names(g)))
  expect_equal(pgen(-6, g), 0, tolerance = 1e-06)
  expect_equal(pgen( 6, g), 1, tolerance = 1e-06)
  grid <- seq(-5, 5, by = 0.25)
  expect_true(all(diff(pgen(grid, g)) > 0))
  expect_true(all(dgen(grid, g) >= 0))
  expect_gt(g$total_arclength, 12)          # at least the width of the window
  set.seed(6)
  r <- rgen(3000L, g)
  expect_true(all(r >= -6 & r <= 6))
  ## the sample must reproduce its own generating cumulative function
  expect_lt(abs(mean(r <= 0) - pgen(0, g)), 0.03)
})

test_that("THE TEST front end returns a proper htest and separates uniform from clustered", {
  set.seed(31)
  t_unif <- al_test(runif(300))
  expect_s3_class(t_unif, "htest")
  expect_true(t_unif$p.value >= 0 && t_unif$p.value <= 1)
  expect_gt(t_unif$p.value, 0.01)
  ## a strongly bimodal sample departs from uniformity and must be detected
  t_bad <- al_test(c(runif(150, 0, 0.05), runif(150, 0.95, 1)))
  expect_lt(t_bad$p.value, t_unif$p.value)
  ## al_scale is a positive scale summary
  expect_gt(al_scale(rnorm(200)), 0)
})

test_that("QUANTILE-FAMILY fits hit the moments they target", {
  ## fit_arcq matches L-moment ratios by construction, so the fitted taus must equal the targets
  ## to optimiser tolerance -- an identity the routine either satisfies or has failed.
  set.seed(41)
  f <- fit_arcq(rnorm(400))
  expect_s3_class(f, "arcq")
  expect_equal(f$fit$tau_fitted, f$fit$tau_target, tolerance = 1e-04)
  expect_identical(f$fit$convergence, 0L)
  fc <- fit_arcq_cf(rnorm(400))
  expect_s3_class(fc, "arcq")
  expect_true(fc$fit$admissible)
  ## the four L-moments of a fitted object are finite and the scale is positive
  L <- lmoments(f, 4L)
  expect_length(L, 4L)
  expect_true(all(is.finite(L)))
  expect_gt(L[2], 0)
})

test_that("CLOSED FORMS of the characteristic-function arc length", {
  ## Two members have exact values, and they are the sharpest checks in the package.
  expect_equal(cf_arclength_family("normal"), 2, tolerance = 1e-08)
  expect_equal(cf_arclength_family("exponential"), pi, tolerance = 1e-08)
})

test_that("EQUIVALENCE fitting recovers a member that lies on the curve", {
  ## The equivalence curve is traced by eq_bstar, which returns a root only where the
  ## discrepancy changes sign over the scanned window -- two bands, one of which contains the
  ## symmetric member the construction turns on. Fitting noiseless data generated ON the curve
  ## must return the generating alpha and drive the sum of squares to numerical zero.
  al <- -0.60
  bstar <- eq_bstar(al)
  ag <- seq(-0.64, -0.51, length.out = 14)
  curve <- cbind(ag, vapply(ag, eq_bstar, numeric(1)))
  expect_true(all(is.finite(curve)))

  x <- seq(0.05, 0.95, length.out = 40)
  y <- bc_cdf(x, al, bstar)
  st <- matrix(c(0, 1, al,    bstar, 0, 0,
                 0, 1, -0.58, -0.35, 0, 0), nrow = 2, byrow = TRUE)
  f <- eqfit_bc(x, y, st, curve)
  expect_lt(f$sse, 1e-06)
  expect_lt(abs(f$th[5] - al), 0.05)
  ## the bootstrap of the same fit is shaped correctly and reproducible under its seed
  b <- eqfit_bc_boot(x, y, f$th, curve, B = 4L, seed = 4207L)
  expect_identical(eqfit_bc_boot(x, y, f$th, curve, B = 4L, seed = 4207L), b)
})
