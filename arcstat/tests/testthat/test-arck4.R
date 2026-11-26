test_that("test-arck4", {
  ## Exercises for the four-parameter kappa module: round trips, known identities, fit recovery.
  expect_lt(abs(k4_cdf(k4_q(0.35, 2, 1.5, 0.2, 0.4), 2, 1.5, 0.2, 0.4) -      0.35), 1e-10)
  ## density integrates the cdf: finite-difference check
  x <- k4_q(0.5, 0, 1, 0.2, 0.4); eps <- 1e-6
  x <- k4_q(0.5, 0, 1, 0.2, 0.4); eps <- 1e-6
  expect_lt(abs((k4_cdf(x + eps, 0, 1, 0.2, 0.4) - k4_cdf(x - eps, 0, 1,      0.2, 0.4))/(2 * eps) - k4_pdf(x, 0, 1,
      0.2, 0.4)), 1e-05)
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

test_that("kappa readings match 50-digit values, whatever the grid argument", {
  ## mpmath, 50 digits; the reference locates the steepest point and the third-derivative reading
  ## by bisection on numerically differentiated log quantile densities. Columns mu, sg, k, h, a,
  ## b, mode. The last four rows are fitted oil shapes at the fitting bound k = -0.98, where the
  ## grid search of versions up to 0.3.1 read the tangent up to 16 per cent too early.
  ref <- rbind(
    c(0, 1, -0.98, 0.2418, -0.618391710679162, -0.664133279975884, -0.376903116941754),
    c(0, 1, -0.98, 0.8088, -0.182328163643372, NaN, -0.094726927684711),
    c(6, 1.3, -0.98, 0.02, 5.02992976054898, 4.9785373175359, 5.36569333649175),
    c(6, 1.3, -0.999999, 0.1, 5.09324965009938, 5.0416030219296, 5.41500013495375),
    c(0, 1, 0.2, 0.4, -0.591499322579312, -0.776714806468733, 0.13782667601799),
    c(0, 1, 0, 0, -1, -0.962423650119207, 0),
    c(0, 1, -0.4155, -1.5959, -2.09123825802291, NaN, -1.08116330098459),
    c(0, 1, -0.9, -30, NaN, NaN, NaN),
    c(0, 1, 1e-9, 1e-9, -0.999999999, -0.962423649305943, 9.999999995e-10),
    c(0, 1, 0.95, 0.3, -0.149752443420154, 0.476837548835723, 0.968548780180409),
    c(2, 0.5, -1.5, 0.45, 1.78914682864312, 1.77195187103037, 1.84947287857182),
    c(4.92632094446957, 6.95102192620741, -0.98, 0.2418, 0.627866604553808, 0.309915953433158, 2.30645911455152),
    c(8.19125403708667, 10.1479488622584, -0.98, 0.2723, 2.11334668357306, 1.6327669378209, 4.52106459703162),
    c(13.9330859822076, 28.8494200479927, -0.98, 0.8088, 8.6730242026808, NaN, 11.2002690555755),
    c(11.4738105118771, 10.7081893223202, -0.98, 0.3776, 5.81731877858587, 5.23762306767254, 8.16365159941232))
  for (i in seq_len(nrow(ref))) {
    th <- c(1, 400, ref[i, 1:4])
    r <- unname(k4_readings(th))
    expect_identical(is.nan(r), is.nan(ref[i, 5:7]))
    ok <- !is.nan(r)
    expect_lt(max(c(0, abs(r - ref[i, 5:7])[ok])), 1e-12)
    expect_identical(unname(k4_readings(th, grid = 400000L)), r)
    expect_identical(unname(k4_readings(th, grid = 50L)), r)
  }
  ## the first guard: no interior steepest point, no reading
  expect_true(all(is.nan(k4_readings(c(0, 1, 0, 1, 1, 0.5)))))
  expect_true(all(is.nan(k4_readings(c(0, 1, 0, 1, 1.2, 0.3)))))
  expect_true(all(is.nan(k4_readings(c(0, 1, 0, 1, 0.1, 1)))))
  expect_true(all(is.nan(k4_readings(c(0, -1, 0, 1, 0.1, 0.2)))))
  ## the second: b exists exactly for h < 1/2 and h k < 1/2
  expect_true(is.nan(k4_readings(c(0, 1, 0, 1, 0.1, 0.5))[["b"]]))
  expect_true(is.finite(k4_readings(c(0, 1, 0, 1, 0.1, 0.4999))[["b"]]))
})

test_that("theoretical L-moment ratios are exact to the edge of existence", {
  ## mpmath, 60 digits, Hosking's gamma-function forms: tau3, tau4, l1, l2
  ref <- rbind(
    c(0.2, 0.4, 0.1373718513157561, 0.1023301585433455, 0.6060735086987798, 0.4879222286971861),
    c(-0.5, 0.3, 0.553162257372821, 0.4059396936298152, 1.67986130941115, 1.428191375575225),
    c(0.3, 0.01, -0.006343883233203597, 0.105492004532896, 0.3475767330701642, 0.5582123759495253),
    c(-0.98, 0.3, 0.9793008357201542, 0.966575435477793, 49.58640438135162, 49.05639225332278),
    c(-0.3, -0.9, 0.3051467962097121, 0.2442369232751757, 0.5921629194784419, 1.150269280221254))
  for (i in seq_len(nrow(ref)))
    expect_equal(k4_tau34(ref[i, 1], ref[i, 2]), ref[i, 3:6], tolerance = 1e-12)
  expect_true(all(is.nan(k4_tau34(-1, 0.3))))
})

test_that("the arc-length shape fit integrates a band that starts at zero", {
  ## quantiles of a known member on a fine grid: the sample band lengths are then the model's, and
  ## the shape comes back. The old rule had a relative error of 5 per cent on the first band.
  n <- 20000; u <- (1:n) / (n + 1)
  y <- k4_q(u, 0, 1, 0.1, 0.6)
  f <- k4_fit_aleq(y, rbind(c(0, 0.2), c(0.2, 0.4), c(0.4, 0.6), c(0.6, 0.8), c(0.8, 0.9)))
  expect_lt(abs(f[["k"]] - 0.1), 5e-3)
  expect_lt(abs(f[["h"]] - 0.6), 5e-3)
})

test_that("multiverse replicates restart from the fitted shape, h included", {
  skip_on_cran()
  ## the replicate start once carried log h where the fit takes h: for h < 1 it was outside the box,
  ## the free fit never moved and the selection always fell back on the h = 0 submodel
  set.seed(1)
  x <- seq(0, 20, length.out = 200)
  y <- 1 + 300 * k4_cdf(x, 8, 2, -0.3, 0.6) + rnorm(200, sd = 0.5)
  m <- k4_mv_boot(x, y, B = 4L)
  expect_gt(max(abs(m$A[, "LS"] - m$A[, "GEV-LS"])), 1e-5)      # identical before the fix
  expect_lt(max(abs(m$A[, "LS"] - m$a_pt[["LS"]])), 0.05)
})
