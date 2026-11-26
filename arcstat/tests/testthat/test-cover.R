## correctness checks on routines the other files reach only through examples

test_that("circular family: exact interval covers the truth, goodness of fit is a valid p-value", {
  set.seed(11)
  th <- rarccirc(300, arccirc(0.4))
  ci <- arcc_exact_ci(th, B = 199L)
  expect_true(ci$lower <= 0.4 && ci$upper >= 0.4)
  expect_identical(arcc_exact_ci(th, B = 199L, stat = 1L, group = 1L)$lower,
                   arcc_exact_ci(th, B = 199L, stat = 1L, group = 1L)$lower)
  g <- arcc_gof(th, fit_arccirc(th), B = 199L)
  expect_true(g$p > 0 && g$p <= 1 && g$stat > 0)
  ## a sample from another law is rejected
  expect_lt(arcc_gof(c(rnorm(150, 0, 0.2), rnorm(150, pi, 0.2)) %% (2 * pi), fit_arccirc(th), B = 199L)$p, 0.05)
})

test_that("Fejer-Riesz family: the density integrates to one and the fit recovers the coefficients", {
  fr <- arccirc_fr(c(0.3 + 0.1i), mu = 0.2)
  expect_equal(stats::integrate(function(t) darccirc_fr(t, fr), 0, 2 * pi)$value, 1, tolerance = 1e-6)
  expect_equal(diff(Q_arccirc_fr(c(0, 1), fr)), 1, tolerance = 1e-6)
  expect_true(all(qd_arccirc_fr(c(0.2, 0.7), fr) > 0))
  set.seed(5)
  th <- rarccirc_fr(4000, fr)
  f <- fit_arccirc_fr(th, nm = 1L)
  expect_true(f$admissible)
  expect_equal(Mod(trigmom_arccirc_fr(1L, arccirc_fr(f$p[-1], f$mu)) -
                   mean(exp(1i * th))), 0, tolerance = 0.05)
  l <- lmom_circ(th, 3L)
  r <- rho_from_lmom(l)
  expect_equal(length(r), length(l))
})

test_that("tempered von Mises: density and moments agree", {
  m <- trigmom_temper_vm(1, kappa = 2, mu = 0.5)
  num <- stats::integrate(function(t) cos(t) * dtemper_vm(t, 2, 0.5), 0, 2 * pi)$value +
    1i * stats::integrate(function(t) sin(t) * dtemper_vm(t, 2, 0.5), 0, 2 * pi)$value
  expect_equal(Mod(m - num), 0, tolerance = 1e-6)
})

test_that("closed forms of the characteristic-function arc length", {
  expect_equal(cf_arclength_family("cauchy", delta = 0), 2)
  expect_equal(cf_arclength_family("gamma", k = 1), 2 * sqrt(pi) * gamma(1.5) / gamma(1))
  expect_equal(cf_arclength_family("skewstable", alpha = 1.5, beta = 0), 2)
  expect_gt(cf_arclength_family("normal", delta = 1), 2)
  expect_equal(cf_arclength_family("poisson", lambda = 1), 2 * pi * exp(-1) * besselI(1, 0))
})

test_that("goodness-of-fit test by bootstrap and in the extreme tail", {
  set.seed(2)
  b <- al_test(runif(80), nboot = 200L, rnull = runif)
  expect_true(b$p.value > 0 && b$p.value <= 1)
  expect_lt(al_test(rep(c(0.1, 0.9), 100) + runif(200, 0, 1e-3))$p.value, 1e-3)
})

test_that("sweeps equal the readings they sweep", {
  k <- c(-0.5, 0, 0.3); h <- 0.2
  s <- k4_ab_sweep(k, h)
  r <- t(sapply(k, function(kk) k4_readings(c(0, 1, 0, 1, kk, h))))
  expect_equal(s$a, unname(r[, "a"]), tolerance = 1e-12)
  expect_equal(s$b, unname(r[, "b"]), tolerance = 1e-12)
  expect_equal(length(k4_b_sweep(k, h)), 3L)
})

test_that("kappa curve fits recover a noiseless curve", {
  x <- seq(0, 20, length.out = 150)
  th <- c(1, 300, 8, 2, 0.15, 0.4)
  y <- th[1] + th[2] * k4_cdf(x, th[3], th[4], th[5], th[6])
  f <- k4_fit_nls(x, y, c(1, log(290), 7.8, log(2.1), 0.1, log(0.45)))
  expect_lt(f[["obj"]], 1e-4)
  fd <- k4_fit_nls(x, y + 0.5 * x, c(1, 0.5, log(290), 7.8, log(2.1), 0.1, log(0.45)), drift = TRUE)
  expect_lt(fd[["obj"]], 1e-2)
  fa <- k4_fit_nalr(x, y, c(1, log(290), 7.8, log(2.1), 0.1, log(0.45)))
  expect_lt(abs(fa$theta[3] - 8), 0.5)
})

test_that("tropical band arc length: its mean is the Monte Carlo mean", {
  x <- seq(0, 20, length.out = 100); th <- c(0, 50, 8, 2, -0.3, 0.2); br <- c(0, 10, 20)
  m <- k4_band_trop_mean(x, th, sigma = 1, breaks = br)
  set.seed(3)
  mc <- rowMeans(replicate(2000, {
    y <- th[1] + th[2] * k4_cdf(x, 8, 2, -0.3, 0.2) + rnorm(100)
    k4_band_sample(x, y, br, p = Inf) }))
  expect_equal(m, mc, tolerance = 0.01)
})

test_that("kappa fits on a locus table and their bootstrap", {
  h <- seq(0.05, 0.45, by = 0.05); loc <- cbind(h, k4_locus_k(h))
  x <- seq(0, 20, length.out = 120)
  y <- 1 + 300 * k4_cdf(x, 9, 2, k4_locus_k(0.2), 0.2) + sin(4321 * seq_along(x))
  st <- as.matrix(expand.grid(1, log(300), c(8, 10), log(c(1.5, 2.5)), c(-1, 0, 1)))
  f <- eqfit_k4eq(x, y, st, loc, maxit = 600L)
  expect_equal(f$th[6], 0.2, tolerance = 0.05)
  w <- (f$th[6] - min(h)) / diff(range(h))
  p0 <- c(f$th[1], log(f$th[2]), f$th[3], log(f$th[4]), stats::qlogis(w))
  b <- eqfit_k4eq_boot(x, y, f$th, p0, loc, B = 6L, maxit = 300L)
  expect_equal(length(b), 6L)
  expect_lt(max(abs(b - f$a), na.rm = TRUE), 1)
})

test_that("manifold-test simulation is deterministic and its audit holds", {
  al <- seq(-0.64, -0.52, length.out = 13); curve <- cbind(al, vapply(al, eq_bstar, 0))
  th0 <- c(1, 300, 0, 20, -0.60, eq_bstar(-0.60))
  thref <- rbind(c(0, 20, -0.58, eq_bstar(-0.58)), c(0, 20, -0.60, eq_bstar(-0.60)), c(0, 20, -0.62, eq_bstar(-0.62)))
  m <- eqfit_msim(60L, 1, 2L, th0, c(0, 0, 0), curve, c(-0.64, -0.52), thref, seed = 3L, maxit = 100L)
  expect_true(all(unlist(m) >= 0))
  expect_true(all(m$Tlev <= m$Taud + 1e-9))
  expect_identical(m, eqfit_msim(60L, 1, 2L, th0, c(0, 0, 0), curve, c(-0.64, -0.52), thref,
                                 seed = 3L, maxit = 100L))
})
