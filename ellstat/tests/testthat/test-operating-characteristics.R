# Operating characteristics are the claims most easily broken and least easily
# noticed: a test can keep its level while losing all its power, and an interval
# can be derived from a correct proof and still fail to cover if the proof's
# hypotheses do not hold in the implementation. Both happened during development
# and neither was caught by a proof, a review or a parity check. Every claim of
# level, power or coverage in the manuscript is therefore pinned here.

test_that("the uniformity test holds its level AND has power", {
  # Level alone is not enough: a test that never rejects also never exceeds its
  # level. An earlier version calibrated by random rotations passed every level
  # check and was nearly blind, because |hat E_1| is almost rotation invariant.
  R <- 400L; n <- 60L
  lev <- mean(vapply(seq_len(R), function(r)
    ell_exact_uniform(rellvm(n, 0, 0, 0.9, seed = 7919 * r), k = 0.9,
                      B = 99L, seed = 104729 * r + 1)$p.value, 0) <= 0.05)
  expect_lt(abs(lev - 0.05), 4 * sqrt(0.05 * 0.95 / R))
  pow <- mean(vapply(seq_len(R), function(r)
    ell_exact_uniform(rellvm(n, 1, 1.5, 0.9, seed = 3313 * r), k = 0.9,
                      B = 99L, seed = 5171 * r + 7)$p.value, 0) <= 0.05)
  expect_gt(pow, 0.5)
})

test_that("the joint test holds its level and sees axial departures", {
  R <- 300L; n <- 60L
  lev <- mean(vapply(seq_len(R), function(r)
    ell_joint_uniform(rellvm(n, 0, 0, 0.9, seed = 7919 * r), P = 2L, k = 0.9,
                      B = 99L, seed = 971 * r + 3)$p.value, 0) <= 0.05)
  expect_lt(abs(lev - 0.05), 4 * sqrt(0.05 * 0.95 / R))

  mkax <- function(nn, e, seed) {
    g <- seq(0, 2 * pi, length.out = 1025)[-1025]
    f <- sqrt(1 - e^2 * cos(g)^2); f <- f / sum(f)
    set.seed(seed); sample(g, nn, replace = TRUE, prob = f)
  }
  p1 <- mean(vapply(seq_len(200L), function(r)
    ell_exact_uniform(mkax(200L, 0.9, 55 + r), k = 0.9, B = 99L,
                      seed = 811 * r + 1)$p.value, 0) <= 0.05)
  p2 <- mean(vapply(seq_len(200L), function(r)
    ell_joint_uniform(mkax(200L, 0.9, 55 + r), P = 2L, k = 0.9, B = 99L,
                      seed = 811 * r + 1)$p.value, 0) <= 0.05)
  expect_lt(p1, 0.15)     # blind to an antipodal departure, as it must be
  expect_gt(p2, 0.50)     # the joint test recovers it
})

test_that("the exact interval covers at every mean direction", {
  # This is the check that exposed two real faults. The null is COMPOSITE: the
  # statistic depends on the mean direction as well as the concentration, so a
  # single-mu calibration under-covers (0.587 was observed at mu = 1). And the
  # endpoints must be interpolated, not taken as extreme retained grid points,
  # which truncates inward by up to one grid step.
  gr <- seq(0.4, 8, length.out = 20L)
  for (mu in c(0, 1, 3)) {
    R <- 120L
    hit <- vapply(seq_len(R), function(r) {
      th <- rellvm(40L, mu = mu, kappa = 2, k = 0.9, seed = 613 * r + 5)
      ci <- ell_exact_ci(th, k = 0.9, kgrid = gr, B = 49L, level = 0.10,
                         seed = 977 * r + 3, ng = 128L, nmu = 6L)
      as.numeric(ci$lower <= 2 && 2 <= ci$upper)
    }, 0)
    # the construction is conservative by design, so coverage must be at least
    # nominal, allowing for Monte Carlo error
    expect_gt(mean(hit), 0.90 - 3 * sqrt(0.9 * 0.1 / R))
  }
})

test_that("the interval endpoints do not depend on the grid resolution", {
  # the interpolation exists precisely so that they do not
  th <- rellvm(60L, mu = 0.5, kappa = 3, k = 0.9, seed = 4242)
  coarse <- ell_exact_ci(th, 0.9, seq(0.5, 9, length.out = 15L), 99L, 0.10,
                         7L, 128L, 4L)
  fine   <- ell_exact_ci(th, 0.9, seq(0.5, 9, length.out = 60L), 99L, 0.10,
                         7L, 128L, 4L)
  expect_lt(abs(coarse$lower - fine$lower), 1.0)
  expect_lt(abs(coarse$upper - fine$upper), 1.5)
})
