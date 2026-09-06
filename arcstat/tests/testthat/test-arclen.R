
## ---- the spacing correction -------------------------------------------------------------------
## The sample band arc length converges to S*, not to S. These tests pin the correction down: the
## corrected model against an independent double integration, the consistency the correction buys,
## and the size of the bias the uncorrected estimator carries.

test_that("al_band_model_star matches independent double integration", {
  ind <- function(sigma, a, b) {
    z <- qnorm(c(a, b))
    integrate(function(zz) vapply(zz, function(z1) {
      f2 <- dnorm(z1)^2
      integrate(function(t) sqrt(f2 + sigma^2 * t^2) * exp(-t), 0, Inf, rel.tol = 1e-11)$value
    }, numeric(1)), z[1], z[2], rel.tol = 1e-11)$value
  }
  for (bd in list(c(0.05, 0.95), c(0.10, 0.90), c(0.25, 0.75))) {
    for (s in c(0.25, 1, 5)) {
      expect_equal(al_band_model_star(s, bd[1], bd[2], 4000), ind(s, bd[1], bd[2]),
                   tolerance = 1e-9)
    }
  }
})

test_that("S* exceeds S, and the two agree in neither limit", {
  for (s in c(0.25, 1, 5)) {
    expect_gt(al_band_model_star(s, 0.1, 0.9), al_band_model(s, 0.1, 0.9))
  }
})

test_that("al_scale is consistent where al_scale_raw is not", {
  set.seed(11)
  for (n in c(1000, 8000)) {
    v  <- vapply(seq_len(60), function(i) al_scale(rnorm(n)),     numeric(1))
    vr <- vapply(seq_len(60), function(i) al_scale_raw(rnorm(n)), numeric(1))
    expect_equal(mean(v), 1, tolerance = 0.02)      # corrected: on target
    expect_gt(mean(vr), 1.02)                       # uncorrected: biased high, and it does not
  }                                                 # shrink with n, which is the point
})

test_that("al_scale is exactly scale equivariant", {
  set.seed(12); x <- rnorm(2000)
  expect_equal(al_scale(7 * x), 7 * al_scale(x), tolerance = 1e-10)
})
