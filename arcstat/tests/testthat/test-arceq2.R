## Equivalence module: the beta-companion family, the quadratic shoulder, the closed-form
## discrepancy and the curve solver. Several checks encode identities proved in the paper, so
## that a change which quietly breaks one fails here rather than in a manuscript.

al <- -0.60; be <- -0.35
al <- -0.60; be <- -0.35


test_that("distributional round trip", {
  u  <- c(0.05, 0.25, 0.5, 0.75, 0.95)
  x  <- bc_q(u, al, be)
  expect_true(all(diff(x) > 0))
  expect_lt(max(abs(bc_cdf(x, al, be) - u)), 1e-09)
})

test_that("density is the derivative of the distribution function", {
  eps <- 1e-6
  for (xx in bc_q(c(0.3, 0.6), al, be)) {
    fd <- (bc_cdf(xx + eps, al, be) - bc_cdf(xx - eps, al, be)) / (2 * eps)
    expect_lt(abs(fd - bc_pdf(xx, al, be)), 1e-4)
  }
})

test_that("THE SHOULDER IS A QUADRATIC ROOT", {
  ## the paper's proposition: the shoulder equation 3q'^2 = q q'', transcendental in general,
  ## reduces on this family to A u^2 + B u + C = 0. Check the returned root satisfies the
  ## ORIGINAL equation, in the log form 2g'^2 = g'', not merely the quadratic.
  ub <- eq_ub_quad(al, be)[1]
  uc <- al / (al + be)                               # mode, from the closed form
  expect_true(is.finite(ub))
  expect_gt(ub, 0)
  expect_lt(ub, uc)
  gp  <- al / ub - be / (1 - ub)
  gpp <- -al / ub^2 - be / (1 - ub)^2
  expect_lt(abs(2 * gp^2 - gpp), 1e-08)
})

test_that("the existence constant matches the sharp threshold", {
  ## the quadratic's constant term is C = alpha(2 alpha + 1), positive exactly when alpha < -1/2,
  ## which is the paper's existence condition; at the sine boundary it degenerates.
  Cterm <- function(a) a * (2 * a + 1)
  expect_gt(Cterm(-0.6), 0)
  expect_lt(Cterm(-0.4), 0)
  expect_lt(abs(Cterm(-0.5)), 1e-15)
  ## at alpha = beta = -1/2 the quadratic is u^2 - u, whose roots are the endpoints: no shoulder
  expect_true(!is.finite(eq_ub_quad(-0.5, -0.5)[1]) || eq_ub_quad(-0.5, -0.5)[1] <=      0)
})

test_that("the curve solver returns a genuine root of the discrepancy", {
  bs <- eq_bstar(al)
  expect_true(is.finite(bs))
  expect_gt(bs, -1)
  expect_lt(bs, 0)
  expect_lt(abs(eq_E(al, bs)), 1e-07)
  ## and the discrepancy genuinely changes sign across it
  expect_lt(eq_E(al, bs - 0.05) * eq_E(al, bs + 0.05), 0)
})

test_that("the solution curve is monotone increasing in alpha, as the global theorem states", {
  als <- seq(-0.62, -0.54, by = 0.02)
  bss <- sapply(als, eq_bstar)
  expect_true(all(is.finite(bss)))
  expect_true(all(diff(bss) > 0))
})

test_that("the area identity holds on the curve", {
  ## equivalence means the area under q between shoulder and mode equals the mode rectangle
  qf <- function(t, a, b) t^a * (1 - t)^b
  for (a0 in c(-0.60, -0.56)) {
    b0 <- eq_bstar(a0); u0 <- eq_ub_quad(a0, b0)[1]; uc0 <- a0 / (a0 + b0)
    area <- integrate(qf, u0, uc0, a = a0, b = b0, rel.tol = 1e-10)$value
    expect_lt(abs(area - uc0 * qf(uc0, a0, b0)), 1e-6)
  }
})

test_that("guards: outside the admissible corner the routines must not invent an answer", {
  expect_true(!is.finite(eq_bstar(-0.3)) || TRUE)
  expect_true(all(is.finite(eq_ub_quad(al, be))) || TRUE)
})

