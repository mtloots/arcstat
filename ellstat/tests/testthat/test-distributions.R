# Densities and their samplers: a density must integrate to one and a sampler
# must reproduce the density it claims. Compare against the BIN AVERAGE, never
# the density at the bin midpoint, which measures the bin width instead.

test_that("the projected family is a proper density and its sampler matches it", {
  th <- seq(0, 2 * pi, length.out = 20001)[-20001]
  f <- dprojell(th, tau = 0.3, d = 0.6, beta = 2, ng = 400L)
  expect_equal(sum(f) * (2 * pi / 20000), 1, tolerance = 1e-6)

  # A histogram estimates the BIN AVERAGE of the density, not its value at the
  # bin midpoint. For a peaked density on wide bins those differ by far more
  # than sampling error, so the comparison must be against the bin average or
  # the test measures the bin width rather than the sampler.
  nb <- 40L; sub <- 64L; N <- 400000L
  br <- seq(0, 2 * pi, length.out = nb + 1); w <- diff(br)[1]
  x <- rprojell(N, tau = 0.3, d = 0.6, beta = 2, seed = 5)
  emp <- as.numeric(table(cut(x, br))) / (N * w)
  fine <- seq(0, 2 * pi, length.out = nb * sub + 1)[-(nb * sub + 1)] +
          pi / (nb * sub)
  avg <- colMeans(matrix(dprojell(fine, 0.3, 0.6, 2, 400L), nrow = sub))
  se <- sqrt(max(avg) / (N * w))
  expect_lt(max(abs(emp - avg)), 4 * se)
})

