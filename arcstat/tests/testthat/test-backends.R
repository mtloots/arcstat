## The C back-ends that carried NO test at all. Measured coverage on 21 Aug 2026 was 21.8 per
## cent overall, with arclen.c, bayesarc.c, arcmv.c, arcvc.c, arceq.c, arceqfit.c and
## arck4fit.c at exactly ZERO -- including the arc-length core the package is named for and the
## Bayesian test named in its own title. The checks below are identities and invariances that
## any correct implementation must satisfy, not smoke tests.

test_that("ARC LENGTH of the probability plot is minimised by the diagonal", {
  ## The plot of a perfectly uniform sample against its own ranks IS the diagonal of the unit
  ## square, whose arc length is exactly sqrt(2). No configuration can be shorter, so the
  ## uniform grid attains the support of the statistic and its p-value is one.
  u <- (1:200 - 0.5)/200
  expect_lt(abs(al_statistic(u) - sqrt(2)), 1e-03)
  ## al_moments returns the null mean and variance together with the SUPPORT INTERVAL of the
  ## statistic. Its lower end is the diagonal, exactly sqrt(2); the mean must lie strictly
  ## inside the interval, which is the only arrangement a proper distribution allows.
  m <- al_moments(200L)
  expect_named(m, c("mean", "var", "support"))
  expect_length(m[["support"]], 2L)
  expect_equal(m[["support"]][1], sqrt(2), tolerance = 1e-12)
  expect_gt(m[["mean"]], m[["support"]][1])
  expect_lt(m[["mean"]], m[["support"]][2])
  expect_gt(m[["var"]], 0)
  expect_equal(al_pvalue(al_statistic(u), 200L), 1, tolerance = 1e-08)

  ## clustering the sample lengthens the path, and the p-value must fall monotonically as the
  ## statistic rises
  s <- sort(c(al_statistic(u), al_statistic(sort(rbeta(200, 0.4, 0.4)))))
  p <- vapply(s, al_pvalue, numeric(1), n = 200L)
  expect_true(all(p >= 0 & p <= 1))
  expect_true(p[1] >= p[2])
})

test_that("VARIABLE-PROJECTION fit recovers a noiseless kappa response exactly", {
  ## With no noise the design is interpolated by the true parameters, so a correct optimiser
  ## must return them and drive the residual sum of squares to numerical zero. This is the
  ## strongest available check on arck4fit.c.
  x  <- seq(0, 10, length.out = 60)
  th <- c(0, 1, 5, 1.5, 0.2, 0.8)
  y  <- th[1] + th[2] * k4_cdf(x, th[3], th[4], th[5], th[6])
  st <- matrix(c(5, 1.5, 0.2, 0.8,
                 4, 2.0, 0.0, 1.0), nrow = 2, byrow = TRUE)
  f <- k4_fit_varpro(x, y, st)
  expect_named(f, c("theta", "drift", "rss"))
  expect_lt(f$rss, 1e-12)
  expect_equal(f$theta[3:6], th[3:6], tolerance = 1e-04)
})

test_that("THE MULTIVERSE returns eight named pipelines and is reproducible", {
  x  <- seq(0, 10, length.out = 60)
  y  <- k4_cdf(x, 5, 1.5, 0.2, 0.8)
  m0 <- k4_mv_boot(x, y, B = 0L, seed = 4207L)
  expect_length(m0$a_pt, 8L)
  expect_named(m0$a_pt, c("LS", "med9-LS", "med21-LS", "arc", "GEV-LS",
                          "transient-excised", "grid-odd", "grid-even"))
  expect_null(m0$A)                            # no replicate matrix when none are asked for
  m1 <- k4_mv_boot(x, y, B = 3L, seed = 4207L)
  expect_equal(dim(m1$A), c(3L, 8L))
  ## the seed drives per-replicate splitmix64 streams, so a rerun must be bit-for-bit identical
  expect_identical(k4_mv_boot(x, y, B = 3L, seed = 4207L)$A, m1$A)
  expect_false(identical(k4_mv_boot(x, y, B = 3L, seed = 99L)$A, m1$A))
})

test_that("VARIANCE COMPONENTS: the intraclass correlation hits its two endpoints", {
  ## Groups whose members are identical carry all their variation between groups, so the
  ## intraclass correlation is exactly one; independent noise with no group structure sits at
  ## the other end and cannot be distinguished from zero.
  icc_perfect <- k4_icc(rep(c(1, 2, 3), each = 4), rep(1:3, each = 4))
  expect_equal(unname(icc_perfect[1]), 1, tolerance = 1e-10)
  set.seed(11)
  icc_noise <- k4_icc(rnorm(300), rep(1:30, each = 10))
  expect_lt(abs(icc_noise[1]), 0.35)
})

test_that("LOCUS DISTANCE vanishes on the locus and grows away from it", {
  lh <- seq(0, 1, length.out = 21); lk <- lh^2
  expect_equal(k4_locus_dist(lh[7], lk[7], lh, lk), 0, tolerance = 1e-12)
  d_near <- k4_locus_dist(lh[7], lk[7] + 0.01, lh, lk)
  d_far  <- k4_locus_dist(lh[7], lk[7] + 0.50, lh, lk)
  expect_gt(d_near, 0)
  expect_gt(d_far, d_near)
})

test_that("EQUIVALENCE readings and the curve that defines them", {
  ## eq_bstar returns the beta at which the discrepancy vanishes for a given alpha; the
  ## discrepancy must therefore be zero there and change sign across it. The admissible domain
  ## is NEGATIVE alpha, which is what the documented example uses.
  bs <- eq_bstar(-0.60)
  expect_true(is.finite(bs))
  expect_lt(abs(eq_E(-0.60, bs)), 1e-10)
  expect_lt(eq_E(-0.60, bs - 0.05) * eq_E(-0.60, bs + 0.05), 0)
  ## the sweep is the vectorised form of the same quantity and must agree with it pointwise
  bb <- c(bs - 0.05, bs, bs + 0.05)
  expect_equal(eq_E_sweep(-0.60, bb), vapply(bb, function(b) eq_E(-0.60, b), numeric(1)),
               tolerance = 1e-10)
  ## eq_readings is defined ON the curve: away from it the construction has no solution and
  ## every reading is NaN, which is the honest signal rather than a silent number. Evaluated at
  ## (alpha, bstar(alpha)) the readings are finite and the discrepancy D vanishes.
  r <- eq_readings(-0.60, bs)
  expect_named(r, c("D", "a", "b", "uc", "ub"))
  expect_true(all(is.finite(unlist(r))))
  ## D is computed on the object's own 20001-point grid, so it vanishes to grid accuracy
  ## rather than to machine precision
  expect_lt(abs(unname(r[["D"]])), 1e-03)
  expect_true(all(is.nan(unlist(eq_readings(-0.60, bs + 0.5)))))
})

test_that("EQUIVALENCE simulation helpers are shaped and seeded correctly", {
  t1 <- eqfit_taus(2, 3, R = 5L, n = 50L, seed = 1L)
  expect_equal(dim(t1), c(5L, 4L))
  expect_named(as.data.frame(t1), c("t3", "t4", "al", "be"))
  expect_identical(eqfit_taus(2, 3, R = 5L, n = 50L, seed = 1L), t1)
  set.seed(5)
  bl <- eqfit_blocklen(rnorm(100))
  expect_true(is.finite(bl) && bl >= 1)
})

test_that("THE BAYESIAN TEST separates a uniform sample from a departure", {
  ## The evidence is the posterior probability that the arc-length functional exceeds the
  ## reference quantile. A sample that IS uniform must produce none of it; a beta(2,5) sample
  ## is visibly not uniform and must produce nearly all of it. The reference is shared so the
  ## two are compared against one threshold.
  set.seed(3)
  ref <- bb_ref_disc(60L)
  expect_equal(dim(ref)[2], 2L)
  expect_named(as.data.frame(ref), c("arc", "ks"))

  e_unif <- bb_evidence((1:60 - 0.5)/60, ref = ref, M = 400L, seed = 1L)
  e_beta <- bb_evidence(sort(rbeta(60, 2, 5)), ref = ref, M = 400L, seed = 1L)
  expect_true(e_unif >= 0 && e_unif <= 1)
  expect_lt(e_unif, 0.10)
  expect_gt(e_beta, 0.90)

  ## the posterior draws are a matrix over both functionals, and the seed reproduces them
  p <- bb_post_disc((1:60 - 0.5)/60, M = 200L, seed = 1L)
  expect_equal(dim(p), c(200L, 2L))
  expect_identical(bb_post_disc((1:60 - 0.5)/60, M = 200L, seed = 1L), p)
})

test_that("eq_readings matches high-precision values, including the singular left end", {
  ## Constants: mpmath at 50 digits (incomplete beta, closed-form roots of g' = 0 and
  ## 2 g'^2 = g''). Columns: alpha, beta, a, b, uc, ub. The integral of q must start at ZERO:
  ## for alpha below -1/2 the mass on (0, 1e-5) and the trapezoid overshoot on the first cells
  ## were once dropped, which lowered a and b by up to 5 per cent. beta = -1 + 1e-9 is the
  ## near-degenerate end of the equivalence curve.
  ref <- rbind(
    c(-0.6434, -0.966, 1.11192920277, 1.11193319294, 0.399776314154, 0.0709201279567),
    c(-0.6, -0.378118, 1.06789574393, 1.06789588842, 0.613422920343, 0.115438787487),
    c(-0.55, -0.9, 0.640186539725, 0.442993243486, 0.379310344828, 0.027298078942),
    c(-0.52, -1 + 1e-9, 0.512394550159, 0.228607978552, 0.342105263383, 0.00994857612993),
    c(-0.505, -1 + 1e-9, 0.465652925522, 0.104103098037, 0.33554817298, 0.00249685186998),
    c(-0.6, -1 + 1e-9, 0.843319480661, 0.756439847235, 0.375000000234, 0.0486043951283),
    c(-0.64, -0.99, 1.08223342321, 1.07374270604, 0.39263803681, 0.0678008353378)
  )
  for (i in seq_len(nrow(ref))) {
    r <- eq_readings(ref[i, 1], ref[i, 2])
    expect_equal(unname(r[["a"]]),  ref[i, 3], tolerance = 1e-7)
    expect_equal(unname(r[["b"]]),  ref[i, 4], tolerance = 1e-7)
    expect_equal(unname(r[["uc"]]), ref[i, 5], tolerance = 1e-7)
    expect_equal(unname(r[["ub"]]), ref[i, 6], tolerance = 1e-7)
  }
  ## tilted kernel (theta = 0.05, 0.02, 0.01), mpmath quadrature: a, uc
  rt <- eq_readings(-0.6, -0.9, c(0.05, 0.02, 0.01))
  expect_equal(unname(rt[["a"]]), 0.826378472447, tolerance = 1e-7)
  expect_equal(unname(rt[["uc"]]), 0.391886932994, tolerance = 1e-7)
  ## van Staden-Loots host, closed form (mpmath): lambda = 0.3, delta = 0.6
  rv <- eq_readings_vsl(0.3, 0.6)
  expect_equal(unname(rv[["a"]]), 0.652733314402, tolerance = 1e-7)
  expect_equal(unname(rv[["uc"]]), 0.440653844715, tolerance = 1e-7)
  qv <- function(u) 0.4 * u^(-0.7) + 0.6 * (1 - u)^(-0.7)
  Qv <- function(u) (0.4 * u^0.3 + 0.6 * (1 - (1 - u)^0.3)) / 0.3
  expect_equal(unname(rv[["b"]]), Qv(unname(rv[["ub"]])), tolerance = 1e-10)
  expect_equal(unname(rv[["D"]]), unname(rv[["b"]] - rv[["a"]]), tolerance = 1e-12)
  expect_equal(unname(rv[["a"]]), Qv(unname(rv[["uc"]])) - unname(rv[["uc"]]) * qv(unname(rv[["uc"]])),
               tolerance = 1e-10)
})

test_that("the profiled minimiser does not stop above the minimum it can reach", {
  skip_on_cran()
  ## Null draw 0 of the manifold test at alpha = -0.52, seed 912, n = 400, noise 2, design and
  ## curve as in the equivalence study. Versions up to 0.3.1 returned 0.292312001128; a restarted
  ## search reaches 0.1354 on the same data, and the statistic at the truth (1.865) bounds both.
  alpha1 <- -0.6449
  als <- seq(alpha1 + 1e-3, -0.502, length.out = 40)
  bs <- vapply(als, eq_bstar, numeric(1)); ok <- is.finite(bs)
  curve <- cbind(als[ok], bs[ok])
  th0 <- c(1, 400, 0, 20, -0.52, approx(curve[, 1], curve[, 2], xout = -0.52)$y)
  Tn <- eqfit_nullT(400L, 2, 1L, th0, c(0, 0, 0), curve, c(-0.641, -0.506), seed = 912L)
  expect_true(is.finite(Tn))
  expect_lt(Tn, 0.292312001128 - 1e-3)
  expect_lt(Tn, 1.86)
})
