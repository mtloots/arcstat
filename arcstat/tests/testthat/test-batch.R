## threads, the batched fits, and the bootstraps of arcbatch.c

kappa_run <- function(seed, k, h, n = 120L, drift = 0.2) {
  x <- seq(0, 20, length.out = n)
  list(x = x, y = 1 + drift * x + 300 * k4_cdf(x, 9, 2, k, h) + 2 * sin(seed * seq_len(n)))
}
R1 <- kappa_run(4321, -0.3, 0.2); R2 <- kappa_run(977, 0.1, -0.4, 90L)
X <- list(R1$x, R2$x); Y <- list(R1$y, R2$y)
box <- c(-0.98, 0.95, -60, 4)
curve <- cbind(seq(-0.62, -0.51, length.out = 8), vapply(seq(-0.62, -0.51, length.out = 8), eq_bstar, 0))
ob <- function(x, y, F) sum(stats::lm.fit(cbind(1, x, F), y)$residuals^2)
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
    F <- stats::qbeta(pmin(pmax(z, 0), 1), p[["alpha"]] + 1, p[["beta"]] + 1)
    expect_equal(p[["rss"]], ob(X[[r]], Y[[r]], F), tolerance = 1e-7)
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
    expect_equal(p[["rss"]], ob(X[[r]], Y[[r]], k4_cdf(X[[r]], p[["mu"]], p[["sigma"]], p[["k"]], p[["h"]])), tolerance = 1e-11)
    v <- k4_fit_varpro(X[[r]], Y[[r]], as.matrix(expand.grid(c(7, 10), log(c(1, 3)), c(-0.4, 0.2), c(-1, 0.1))), bounds = box)
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
  h <- bc_host_fit(X, Y); t <- bc_host_fit(X, Y, curve)
  expect_true(all(h$con[, "rss"] <= t$con[, "rss"] * (1 + 1e-4)) && all(h$con[, "rss"] >= h$free[, "rss"]))
  expect_equal(h$con[[1, "beta"]], eq_bstar(h$con[[1, "alpha"]]), tolerance = 1e-12)
  a <- k4_mv_boot(R1$x, R1$y, B = 0L, bounds = box, drift = TRUE)$a_pt
  n <- k4_mv_boot(R1$x, R1$y, B = 0L, bounds = box)$a_pt
  expect_identical(a[c("LS", "GEV-LS", "grid-odd")], n[c("LS", "GEV-LS", "grid-odd")])
  expect_lt(abs(a[["med9-LS"]] - a[["LS"]]), abs(n[["med9-LS"]] - n[["LS"]]))
})
