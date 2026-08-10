test_that("test-arck4", {
  ## Exercises for the four-parameter kappa module: round trips, known identities, fit recovery.
  expect_lt(abs(k4_cdf(k4_q(0.35, 2, 1.5, 0.2, 0.4), 2, 1.5, 0.2, 0.4) -      0.35), 1e-10)
  ## density integrates the cdf: finite-difference check
  x <- k4_q(0.5, 0, 1, 0.2, 0.4); eps <- 1e-6
  x <- k4_q(0.5, 0, 1, 0.2, 0.4); eps <- 1e-6
  expect_lt(abs((k4_cdf(x + eps, 0, 1, 0.2, 0.4) - k4_cdf(x - eps, 0, 1,      0.2, 0.4))/(2 * eps) - k4_pdf(x, 0, 1, 0.2, 0.4)), 1e-05)
  ## theoretical tau ratios invert through the L-moment fit
  tt <- k4_tau34(0.2, 0.4)
  ft <- k4_fit_lmom(tt[1], tt[2])
  expect_lt(abs(ft["k"] - 0.2), 0.02)
  expect_lt(abs(ft["h"] - 0.4), 0.05)
  ## running median removes an isolated spike
  y <- c(rep(1, 10), 50, rep(1, 10))
  expect_lt(max(abs(k4_runmed(y, 9L) - 1)), 1e-12)
  ## band arcs: model on a flat segment equals the band width
  th <- c(0, 1, 100, 1, 0.2, 0.4)   # curve flat over [0,1]
  expect_lt(abs(k4_band_model(th, c(0, 1))[1] - 1), 1e-06)
  ## readings on a known configuration are finite and ordered sensibly
  r <- k4_readings(c(1, 400, 8, 2.5, 0.2, 0.4))
  expect_true(is.finite(r["a"]))
  expect_true(is.finite(r["b"]))
  expect_lt(r["a"], r["mode"])
  expect_lt(r["b"], r["mode"])
  ## NLS recovers a clean curve
  xg <- 20*(0:399)/399
  yg <- 1 + 400*k4_cdf(xg, 8, 2.5, 0.2, 0.4)
  f  <- k4_fit_nls(xg, yg, c(1, log(380), 7.5, log(2.2), 0.15, log(0.5)))
  expect_lt(abs(f$theta[3] - 8), 0.2)
  expect_lt(abs(f$theta[5] - 0.2), 0.05)
})
