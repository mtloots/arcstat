## threads, the batched fits, and the bootstraps of arcbatch.c

kappa_run <- function(seed, k, h, n = 120L, drift = 0.2) {
  x <- seq(0, 20, length.out = n)
  list(x = x, y = 1 + drift * x + 300 * k4_cdf(x, 9, 2, k, h) + 2 * sin(seed * seq_len(n)))
}
R1 <- kappa_run(4321, -0.3, 0.2); R2 <- kappa_run(977, 0.1, -0.4, 90L)
X <- list(R1$x, R2$x); Y <- list(R1$y, R2$y)
box <- c(-0.98, 0.95, -60, 4)
curve <- cbind(seq(-0.62, -0.51, length.out = 8), vapply(seq(-0.62, -0.51, length.out = 8), eq_bstar, 0))
ob <- function(x, y, Fv) sum(stats::lm.fit(cbind(1, x, Fv), y)$residuals^2)
free <- function(r) {
  f <- k4_fit_varpro(r$x, r$y, rbind(c(9, log(2), -0.1, 0.1), c(9, log(2), 0, -0.5)), bounds = box)
  c(f$theta[1], f$drift, f$theta[2:6])
}

test_that("the thread setting is honoured and cleared", {
  old <- options(arcstat.threads = NULL); on.exit(options(old))
  expect_lte(arc_threads(), max(2L, suppressWarnings(as.integer(Sys.getenv("ARCSTAT_THREADS", "2"))), na.rm = TRUE))
  n <- arc_threads(3L)
  expect_true(n %in% c(1L, 3L))                  # 1 where OpenMP is absent
  expect_identical(getOption("arcstat.threads"), 3L)
  arc_threads(NULL)
  expect_null(getOption("arcstat.threads"))
  expect_error(arc_threads(0))
})

test_that("results are bit-identical on 1, 2 and 8 threads", {
  skip_on_cran()                                 # eight threads
  old <- options(arcstat.threads = NULL); on.exit(options(old))
  fit <- rbind(free(R1), free(R2))
  run <- function(nt) {
    options(arcstat.threads = nt)
    e <- k4eq_fit(X, Y)
    list(e, k4eq_fit(X, Y, cbind(seq(-1, 0.4, 0.2), k4_locus_k(seq(-1, 0.4, 0.2)))),
         bc_host_fit(X, Y, curve, keep = 3L), k4_fit(X, Y, gate = TRUE, gev = TRUE),
         k4eq_fit(X, Y, gate = "best"), k4_fit(X, Y, gate = "best", infimum = TRUE),
         k4eq_boot(X, Y, e$fit, B = 3L, gate = "best"), wald_fit(X, Y), k4_fit(X, Y, drift = FALSE),
         k4_shared_fit(c(X, X), c(Y, lapply(Y, rev)), sets = list(1:4, 1:2)),
         k4_lmom_window(X[1], Y[1], c(1, 0.2, 300), ngrid = 400L, nmodel = 400L),
         k4_boot(X, Y, fit, box, gev = c(FALSE, TRUE), B = 6L), k4eq_boot(X, Y, e$fit, B = 6L),
         k4_mv_boot(R2$x, R2$y, B = 2L, bounds = box, maxit = 300L),
         k4_fit_varpro(R1$x, R1$y, rbind(c(9, log(2), 0, 0.1), c(8, 0, 0, -1)), bounds = box),
         eqfit_estsim(0.4, 0.65, 3L, 60L, seed = 11L, maxit = 100L), eqfit_taus(0.4, 0.65, 4L, 50L))
  }
  one <- run(1L)
  expect_identical(run(2L), one)
  expect_identical(run(8L), one)
})

test_that("the locus has equal readings, and its two ends", {
  h <- c(-12, -3, -0.4, -1e-9, 0, 1e-7, 0.2, 0.45, 0.499)
  k <- k4_locus_k(h)
  d <- vapply(seq_along(h), function(i) { r <- k4_readings(c(0, 1, 0, 1, k[i], h[i])); r[["a"]] - r[["b"]] }, 0)
  expect_lt(max(abs(d[-1])), 1e-12)
  expect_lt(abs(d[1]), 1e-6)                     # the readings lose the frontier distance, not k
  expect_lt(max(abs(k[c(4, 6)] - k[5])), 1e-6)   # one curve through h = 0
  expect_equal(k[8], 0.932174548080419, tolerance = 1e-14)      # mpmath, 40 digits
  expect_equal(k[3] < 0 && h[3] * k[3] < 0.5, TRUE)
  expect_equal(k4_locus_k(-60) * -60, 0.5, tolerance = 1e-14)
  expect_lt(1 - k4_locus_k(0.4999), 2e-4)
  expect_true(all(is.nan(k4_locus_k(c(0.5, 2)))))
})

test_that("constrained kappa fit: objective, constraint, and no better point nearby", {
  f <- k4eq_fit(X, Y)
  for (r in 1:2) {
    p <- f$fit[r, ]
    expect_equal(p[["rss"]], ob(X[[r]], Y[[r]], k4_cdf(X[[r]], p[["mu"]], p[["sigma"]], p[["k"]], p[["h"]])),
                 tolerance = 1e-11)
    expect_equal(p[["a"]], p[["b"]], tolerance = 1e-10)
    g <- function(q) { k <- k4_locus_k(q[3]); if (!is.finite(k)) return(1e12)
      ob(X[[r]], Y[[r]], k4_cdf(X[[r]], q[1], exp(q[2]), k, q[3])) }
    o <- stats::optim(c(p[["mu"]], log(p[["sigma"]]), p[["h"]]), g, control = list(maxit = 2000, reltol = 1e-14))
    expect_gt(o$value, p[["rss"]] * (1 - 1e-7))
    expect_equal(p[["rss"]], min(f$pieces[r, "rss", ]))
  }
  ## a member of the locus, without noise, is recovered
  x <- seq(0, 20, length.out = 150)
  y <- 1 - 0.3 * x + 200 * k4_cdf(x, 9, 2, k4_locus_k(-0.7), -0.7)
  q <- k4eq_fit(x, y)$fit[1, ]
  expect_lt(q[["rss"]], 1e-8)
  expect_equal(q[["h"]], -0.7, tolerance = 1e-4)
})

test_that("a table is searched past its kinks", {
  h <- seq(-1.2, 0.45, 0.15); tab <- cbind(h, k4_locus_k(h))
  f <- k4eq_fit(X, Y, tab)$fit
  for (r in 1:2) {
    cell <- vapply(seq_len(nrow(tab) - 1L), function(i) {
      fit <- k4eq_fit(X[[r]], Y[[r]], tab[i + 0:1, ])$fit[1, "rss"]; unname(fit) }, 0)
    expect_lt(f[r, "rss"], min(cell) * (1 + 1e-7))
  }
  expect_error(k4eq_fit(X, Y, rbind(c(0, 0), c(0, 1))), "distinct")
})

test_that("companion host: objective, exact ends, constrained never below free", {
  h <- bc_host_fit(X, Y, curve)
  expect_true(all(h$con[, "rss"] >= h$free[, "rss"]))
  for (r in 1:2) for (m in c("free", "con")) {
    p <- h[[m]][r, ]
    z <- (X[[r]] - p[["s"]]) / p[["w"]]
    Fv <- stats::qbeta(pmin(pmax(z, 0), 1), p[["alpha"]] + 1, p[["beta"]] + 1)
    expect_equal(p[["rss"]], ob(X[[r]], Y[[r]], Fv), tolerance = 1e-7)
  }
  expect_equal(unname(h$con[, "beta"]), stats::approx(curve[, 1], curve[, 2], h$con[, "alpha"])$y, tolerance = 1e-12)
})

test_that("bootstraps: streams, block length, boxes, held submodel", {
  fit <- rbind(free(R1), free(R2))
  b <- k4_boot(X, Y, fit, rbind(box, c(-0.98, 0.95, -3, 4)), gev = c(TRUE, FALSE), B = 8L, seed = c(5L, 9L))
  expect_true(all(b$boot[, "h", 1] == 0) && all(is.nan(b$boot[, "sel", 1])))
  expect_true(all(b$boot[, "h", 2] >= -3) && all(b$boot[, "sel", 2] %in% c(0, 1)))
  alone <- k4_boot(R2$x, R2$y, fit[2, , drop = FALSE], c(-0.98, 0.95, -3, 4), B = 8L, seed = 9L)
  expect_identical(alone$boot[, , 1], b$boot[, , 2])     # a run does not see the others
  for (r in 1:2) {
    p <- fit[r, ]
    res <- Y[[r]] - (p[1] + p[2] * X[[r]] + p[3] * k4_cdf(X[[r]], p[4], p[5], p[6], p[7]))
    expect_identical(b$block[r], eqfit_blocklen(res))
  }
  i <- k4_boot(X, Y, fit, box, B = 3L, block = 1L)
  expect_identical(i$block, c(1L, 1L))
  e <- k4eq_fit(X, Y); q <- k4eq_boot(X, Y, e$fit, B = 8L)
  expect_equal(q$boot[, "a", ], q$boot[, "b", ], tolerance = 1e-9)
  expect_lt(max(abs(q$boot[, "k", ] - k4_locus_k(q$boot[, "h", ]))), 1e-14)
  expect_lt(max(abs(apply(q$boot[, "a", ], 2, stats::median) - e$fit[, "a"])), 0.5)
})

test_that("the multiverse reaches a curve with h < 0", {
  x <- seq(0, 20, length.out = 150)
  y <- 1 + 300 * k4_cdf(x, 9, 2, 0.1, -1.5) + 0.5 * sin(4321 * seq_along(x))
  a0 <- k4_readings(c(1, 300, 9, 2, 0.1, -1.5))[["a"]]
  a <- k4_mv_boot(x, y, B = 0L, bounds = c(-0.98, 0.95, -60, 4))$a_pt
  expect_lt(max(abs(a[c("LS", "med9-LS", "med21-LS", "grid-odd", "grid-even")] - a0)), 0.15)
  expect_lt(abs(a[["arc"]] - a0), 0.6)           # twelve bands of a smoothed trace: coarser
})

test_that("free kappa fit: objective, submodel, gate, and not above Nelder-Mead", {
  f <- k4_fit(X, Y); g <- k4_fit(X, Y, hfix = 0); s <- k4_fit(X, Y, gate = "best", gev = TRUE)
  expect_true(all(g[, "h"] == 0) && all(g[, "rss"] >= f[, "rss"]))
  for (r in 1:2) {
    p <- f[r, ]
    expect_equal(p[["rss"]], ob(X[[r]], Y[[r]], k4_cdf(X[[r]], p[["mu"]], p[["sigma"]], p[["k"]], p[["h"]])),
                 tolerance = 1e-11)
    v <- k4_fit_varpro(X[[r]], Y[[r]], as.matrix(expand.grid(c(7, 10), log(c(1, 3)), c(-0.4, 0.2), c(-1, 0.1))),
                       bounds = box)
    expect_lt(p[["rss"]], v$rss * (1 + 1e-9))
    d <- length(X[[r]]) * log(g[r, "rss"] / f[r, "rss"]) - 2
    expect_equal(s[[r, "gev"]], as.numeric(d < 0))
    expect_equal(s[[r, "physical"]], 1)
  }
  ## a falling curve: the lowest minimum is not physical, gate = "best" does not return it
  y <- 300 - 250 * k4_cdf(R1$x, 9, 2, 0.1, -0.5) + sin(4321 * seq_along(R1$x))
  expect_lt(k4_fit(R1$x, y)[1, "g1"], 0)
  b <- k4_fit(R1$x, y, gate = "best")[1, ]
  expect_true(is.nan(b[["rss"]]) || b[["g1"]] > 0)
})

test_that("the whole companion curve, and drift in every pipeline", {
  skip_on_cran()                                 # about 11 s: the dense host search, twice
  h <- bc_host_fit(X, Y); t <- bc_host_fit(X, Y, curve)
  expect_true(all(h$con[, "rss"] <= t$con[, "rss"] * (1 + 1e-4)) && all(h$con[, "rss"] >= h$free[, "rss"]))
  expect_equal(h$con[[1, "beta"]], eq_bstar(h$con[[1, "alpha"]]), tolerance = 1e-12)
  a <- k4_mv_boot(R1$x, R1$y, B = 0L, bounds = box, drift = TRUE)$a_pt
  n <- k4_mv_boot(R1$x, R1$y, B = 0L, bounds = box)$a_pt
  expect_identical(a[c("LS", "GEV-LS", "grid-odd")], n[c("LS", "GEV-LS", "grid-odd")])
  expect_lt(abs(a[["med9-LS"]] - a[["LS"]]), abs(n[["med9-LS"]] - n[["LS"]]))
})

test_that("drift-free fits are the two-column least squares at their shape, never below the drifted fit", {
  f <- k4_fit(X, Y); n <- k4_fit(X, Y, drift = FALSE); q <- k4eq_fit(X, Y, drift = FALSE)$fit
  expect_true(all(n[, "m"] == 0) && all(q[, "m"] == 0) && all(n[, "rss"] >= f[, "rss"]) && all(q[, "rss"] >= n[,
      "rss"] * (1 - 1e-9)))
  for (r in 1:2) for (p in list(n[r, ], q[r, ])) {
    Fv <- k4_cdf(X[[r]], p[["mu"]], p[["sigma"]], p[["k"]], p[["h"]])
    expect_equal(p[["rss"]], sum(stats::lm.fit(cbind(1, Fv), Y[[r]])$residuals^2), tolerance = 1e-11)
  }
})

test_that("Wald: distribution function, readings against the definition, fit", {
  pw <- function(x, mu, l) { s <- sqrt(l / x); pnorm(s * (x / mu - 1)) + exp(2 * l / mu + pnorm(-s * (x / mu + 1),
                                                                                                log.p = TRUE)) }
  x <- c(0.01, 0.5, 3, 9, 40)
  expect_equal(wald_cdf(x, 9, 4), pw(x, 9, 4), tolerance = 1e-12)
  expect_identical(wald_cdf(c(-1, 0), 9, 4), c(0, 0))
  lf <- function(x, mu, l) 0.5 * log(l) - 1.5 * log(x) - l * (x - mu)^2 / (2 * mu^2 * x)
  for (phi in c(1e-7, 0.01, 0.159, 1, 30, 1e5)) {
    r <- wald_readings(1, phi)[1, ]; m <- r[["mode"]]; b <- r[["b"]]
    d1 <- function(x) -1.5 / x - phi / 2 + phi / (2 * x^2)
    d2 <- function(x) 1.5 / x^2 - phi / x^3
    expect_lt(abs(d1(m)) * m, 1e-9 * (1 + phi))
    expect_lt(abs(d1(b)^2 + d2(b)) / (abs(d1(b)^2) + abs(d2(b))), 1e-9)
    expect_true(b < m && r[["a"]] < m && r[["a"]] > 0)
  }
  ## the Levy limit: the readings no longer depend on mu
  expect_equal(wald_readings(1e8, 15)[1, ], wald_readings(1e6, 15)[1, ], tolerance = 1e-4)
  y <- 1 + 0.2 * R1$x + 300 * wald_cdf(R1$x, 9, 150)
  f <- wald_fit(R1$x, y)[1, ]
  expect_equal(unname(f[c("mu", "lambda", "m")]), c(9, 150, 0.2), tolerance = 1e-5)
  g <- wald_fit(X, Y)
  for (r in 1:2) expect_equal(g[[r, "rss"]], ob(X[[r]], Y[[r]], pw(pmax(X[[r]], 1e-12), g[r, "mu"], g[r,
                                                                                                      "lambda"])),
      tolerance = 1e-9)
})

test_that("windowed L-moment matching is the R construction and recovers a shape of either sign", {
  x <- seq(0, 20, length.out = 400)
  for (h in c(0.3, -0.5)) {
    y <- 1 + 300 * k4_cdf(x, 9, 2, -0.2, h)
    L <- k4_lmom_window(x, y, c(1, 0, 300), bounds = c(-0.98, 0.95, -8, 4))[1, ]
    Fh <- cummax(pmin(pmax((as.numeric(stats::runmed(y, 21)) - 1) / 300, 0), 1)); kp <- Fh > 1e-6 & Fh < 1 - 1e-6
    a <- max(0.005, min(Fh[kp])); b <- min(0.995, max(Fh[kp])); u <- seq(a, b, length.out = 3000)
    tau <- function(u, Q) { w <- c(diff(u), utils::tail(diff(u), 1)); Q <- Q - sum(w * Q) / sum(w)
      c2 <- sum(w * Q * (2 * u - 1)); c(sum(w * Q * (6 * u^2 - 6 * u + 1)),
                                        sum(w * Q * (20 * u^3 - 30 * u^2 + 12 * u - 1))) / c2 }
    e <- tau(u, stats::approx(Fh[kp], x[kp], xout = u, ties = "ordered")$y)
    expect_equal(unname(L[c("t3", "t4")]), e, tolerance = 1e-9)
    um <- seq(a, b, length.out = 4000)
    expect_lt(sum((tau(um, k4_q(um, 0, 1, L[["k"]], L[["h"]])) - e)^2), 1e-12)
    expect_true(L[["nsol"]] >= 1 && L[["h_min"]] <= h + 0.1 && L[["h_max"]] >= h - 0.1)   # not always unique
  }
})

test_that("gated constrained fit: the ungated one where that is physical, none where there is none", {
  u <- k4eq_fit(X, Y); g <- k4eq_fit(X, Y, gate = "best")
  expect_true(all(g$physical[, "status"] %in% 1:2))
  expect_equal(g$fit[, "rss"], u$fit[, "rss"], tolerance = 1e-9)
  expect_true(all(g$physical[, "inf_rss"] <= g$fit[, "rss"]))
  ## a falling curve: nothing rises to a fit, and what is returned is not a fit
  y <- 300 - 250 * k4_cdf(R1$x, 9, 2, 0.1, -0.5) + sin(4321 * seq_along(R1$x))
  expect_lt(k4eq_fit(R1$x, y)$fit[1, "g1"], 0)
  n <- k4eq_fit(R1$x, y, gate = "best")
  expect_true(n$physical[1, "status"] %in% c(0, 3) == is.nan(n$fit[1, "rss"]))
  if (is.finite(n$fit[1, "rss"])) expect_true(n$fit[1, "g1"] > 0 && n$fit[1, "a"] > 0)
  b <- k4eq_boot(X, Y, g$fit, B = 6L, gate = "best"); b0 <- k4eq_boot(X, Y, g$fit, B = 6L)
  expect_identical(unname(b$lost), unname(colSums(is.nan(b$boot[, "a", ]))))
  expect_lt(max(abs(b$boot - b0$boot), na.rm = TRUE), 1e-4)
  f <- k4_fit(X, Y, gate = "best", infimum = TRUE)
  expect_identical(f[, 1:14], k4_fit(X, Y, gate = "best"))
  expect_true(all(f[, "status"] == 1) && all(f[, "inf_rss"] <= f[, "rss"]))
  expect_error(k4eq_fit(X, Y, gate = TRUE))
})

test_that("multiverse components are the covariance identity", {
  set.seed(3); A <- matrix(stats::rnorm(240), 40) + rep(1:6, each = 40) / 10 + stats::rnorm(40)
  C <- stats::cov(A); v <- mean(diag(C)); w <- v - mean(C[upper.tri(C)]); m <- stats::var(colMeans(A))
  expect_equal(k4_mv_share(A)[c("v", "w", "delta2", "tau2")], c(v = v, w = w, delta2 = max(0, m - w / 40),
                                                                tau2 = max(0, m - w * 41 / 40)), tolerance = 1e-10)
  B <- A; B[1:30, 2] <- NaN; B[1, 1] <- NaN
  expect_equal(k4_mv_share(B), k4_mv_share(A[-1, -2]), tolerance = 1e-12)
})

test_that("shared shape: the least squares of the model, never above a pinned refit, and the permutation", {
  skip_on_cran()                                 # about 8 s
  f <- function(mu, s, g, h, e) 1 + 0.2 * R1$x + g * k4_cdf(R1$x, mu, s, -0.4, h) + sin(e * seq_along(R1$x))
  y <- list(f(8, 2, 300, 0.3, 11), f(10, 2.5, 280, 0.3, 13), f(5, 1.2, 310, -0.4, 17), f(6, 1.5, 290, -0.4, 19))
  x <- rep(list(R1$x), 4)
  r <- k4_shared_fit(x, y, sets = list(1:4, 1:2, 3:4))
  for (s in 1:3) {
    m <- r$runs[[s]]; ix <- list(1:4, 1:2, 3:4)[[s]]
    v <- sum(sapply(seq_along(ix), function(j) ob(x[[ix[j]]], y[[ix[j]]], k4_cdf(x[[ix[j]]], m[j, "mu"], m[j,
      "sigma"], r$fit[s, "k"], r$fit[s, "h"]))))
    expect_equal(v, unname(r$fit[s, "rss"]), tolerance = 1e-10)
    st <- cbind(as.matrix(expand.grid(c(4, 6, 8, 10), log(c(1, 2, 4)))), r$fit[s, "k"], r$fit[s, "h"])
    pin <- sum(sapply(ix, function(j) k4_fit_varpro(x[[j]], y[[j]], st, 1500L, c(r$fit[s, "k"] - 1e-9, r$fit[s,
      "k"] + 1e-9, -3, 4), hfix = r$fit[s, "h"])$rss))
    expect_lt(r$fit[s, "rss"], pin * (1 + 1e-7))
  }
  expect_lt(max(abs(unname(r$fit[2:3, "h"]) - c(0.3, -0.4))), 0.25)    # shape is weakly determined by two runs
  p <- k4_shared_perm(x, y, group = c(1, 1, 2, 2), rss_free = 1, slack = c(1e-9, 0.5))
  expect_equal(unname(p$stat[1, c("rss_tr", "rss_w")]), unname(c(r$fit[1, "rss"], sum(r$fit[2:3, "rss"]))),
               tolerance = 1e-6)
  expect_equal(unname(p$stat[1, "n_alloc"]), 6)
  expect_equal(unname(p$p[1, 1]), 2 / 6)                  # the observed split and its mirror
  expect_gte(p$p[1, 2], p$p[1, 1])
  ## the permutation sums are those of the fits of the sets they are built from
  q <- k4_shared_fit(x, y, sets = list(c(1, 3), c(2, 4), c(1, 4), c(2, 3)))
  expect_equal(unname(p$stat[1, "min_rss_w"]), min(sum(r$fit[2:3, "rss"]), sum(q$fit[1:2, "rss"]), sum(q$fit[3:4,
      "rss"])), tolerance = 1e-6)
})

test_that("free fit reaches minima with k far below -1", {
  y <- 1 + 0.1 * R1$x + 400 * k4_cdf(R1$x, 9, 6, -4, 0.2) + sin(4321 * seq_along(R1$x))
  a <- k4_fit(R1$x, y, bounds = c(-30, 0.95, -60, 4)); d <- k4_fit(R1$x, y, bounds = c(-30, 0.95, -60, 4), dense = TRUE)
  expect_lt(a[1, "rss"], d[1, "rss"] * (1 + 1e-6))
  expect_lt(abs(a[1, "k"] + 4), 1)
})
