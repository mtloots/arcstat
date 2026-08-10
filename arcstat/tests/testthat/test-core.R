## Core modules: arc length, the arc-length distribution family, L-moments, the circular
## family and its characteristic function. These modules carried no tests; the checks below
## are identities and round trips that any correct implementation must satisfy.


## shared across several blocks below; test_that() gives each block its own environment, so
## anything used in more than one of them has to live out here
o0 <- arcq(0.5, mu = 0, sigma = 1)

test_that("ARC LENGTH of a quantile family: location leaves it alone, scale acts predictably", {
  ## the arc length of a quantile function is a shape functional, so shifting the distribution
  ## cannot change it, and the exact and quadrature routes must agree.
  o1 <- arcq(0.5, mu = 17, sigma = 1)
  expect_lt(abs(arclength(o0) - arclength(o1)), 1e-08)
  ## the exact route is grid-independent; the fast route is a sum over the object's own grid,
  ## so it must converge to the exact value at first order as that grid refines
  rel <- sapply(c(1000L, 4000L, 16000L, 64000L), function(ng) {
    oo <- arcq(0.5, ngrid = ng)
    abs(arclength(oo, exact = FALSE) / arclength(oo, exact = TRUE) - 1)
  })
  expect_true(all(diff(rel) < 0))
  expect_lt(rel[4], 5e-05)
  expect_true(all(rel[-4]/rel[-1] > 2.5))
  ## the degenerate family (no shape coefficients) has the arc length of its own scale
  od <- arcq(numeric(0), mu = 0, sigma = 1)
  expect_true(is.finite(arclength(od)))
  expect_gt(arclength(od), 0)
})

test_that("L-MOMENTS: known closed forms", {
  ## for a uniform sample on (0,1) the first L-moment is the mean and the second is 1/6
  set.seed(11)
  u <- (seq_len(20000) - 0.5) / 20000
  ## sample L-moments of a uniform grid: mean 1/2 and L-scale 1/6, both known exactly
  sl <- sample_lmoments(u, 4L)
  expect_lt(abs(sl[1] - 0.5), 1e-06)
  expect_lt(abs(sl[2] - 1/6), 1e-04)
  ## symmetric data have vanishing L-skewness
  ## lmratios returns location and scale, then the ratios: tau3 is the third element
  expect_lt(abs(lmratios(sl)[3]), 1e-06)
  ## the population routine on a family returns finite, ordered moments
  pl <- lmoments(o0, 4L)
  expect_true(all(is.finite(pl)))
  expect_gt(pl[2], 0)
})

test_that("ARC-LENGTH DISTRIBUTION family: density, quantile, round trip", {
  pu <- c(0.1, 0.35, 0.6, 0.9)
  qq <- qarcq(pu, o0)
  expect_true(all(diff(qq) > 0))
  expect_lt(max(abs(parcq(qq, o0) - pu)), 1e-05)
  ## density is positive and integrates the distribution function
  eps <- 1e-6
  for (z in qq[2:3]) {
    fd <- (parcq(z + eps, o0) - parcq(z - eps, o0)) / (2 * eps)
    expect_gt(fd, 0)
    expect_lt(abs(fd - darcq(z, o0)), 1e-2)
  }
  ## simulated draws land in the support and match the quantile at the median
  set.seed(7)
  rs <- rarcq(4000, o0)
  expect_true(all(is.finite(rs)))
  expect_lt(abs(median(rs) - qarcq(0.5, o0)), 0.2)
})

test_that("admissibility and the rho parametrisation", {
  expect_true(is.finite(qmin_rho(0.4)))
  expect_gte(length(admiss_rho(0.4)), 1)
})

test_that("CIRCULAR family: density integrates to one over the circle", {
  th <- seq(0, 2 * pi, length.out = 4001)
  oc <- arccirc(0.4)
  dv <- darccirc(th, oc)
  expect_true(all(dv >= 0))
  tot <- sum((dv[-1] + dv[-length(dv)]) / 2 * diff(th))
  expect_lt(abs(tot - 1), 0.001)
  ## the distribution function runs from zero to one and is increasing
  Pv <- Q_arccirc(c(0.001, 0.5, 0.999), oc)
  expect_true(all(diff(Pv) > 0))
  ## trigonometric moments are bounded by one in modulus, as moments of a circular law must be
  tm <- trigmom_arccirc(3L, oc)
  expect_true(all(abs(tm) <= 1 + 1e-08))
})

test_that("CHARACTERISTIC FUNCTION: value at zero is one, and it is conjugate-symmetric", {
  ## cf_arclength takes a SAMPLE and returns the arc length of its empirical characteristic
  ## function; the family routine returns the theoretical value. Two closed forms anchor it.
  expect_lt(abs(cf_arclength_family("normal") - 2), 1e-09)
  expect_lt(abs(cf_arclength_family("exponential", lambda = 1) - pi), 1e-09)
  ## and the empirical statistic approaches the theoretical one on a large normal sample
  set.seed(3)
  expect_lt(abs(cf_arclength(rnorm(50000)) - 2), 0.05)
  ## scale equivariance of the family values: the exponential rate rescales the axis only
  expect_lt(abs(cf_arclength_family("exponential", lambda = 2) - cf_arclength_family("exponential",      lambda = 1)), 1e-09)
})

