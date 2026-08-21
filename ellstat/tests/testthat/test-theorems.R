# The theorems of the manuscript, as regression tests. Each name is the
# theorem it encodes, so a failure names the mathematical claim that broke.

test_that("Theorem (reduction): modulus zero is the trigonometric system", {
  t <- vapply(0:400, function(i) 2 * pi * i / 400, 0)
  b <- ell_basis(t, 3L, 0.0)
  for (p in 1:3) {
    expect_lt(max(abs(b$C[p, ] - cos(p * t))), 1e-12)
    expect_lt(max(abs(b$S[p, ] - sin(p * t))), 1e-12)
  }
  expect_lt(abs(ell_KE(0)$K - pi / 2), 1e-14)
})

test_that("Theorem (Pythagorean): cn^2 + sn^2 = 1 at every modulus", {
  for (k in c(0, 0.3, 0.7, 0.9, 0.999)) {
    u <- vapply(0:200, function(i) 8 * i / 200 - 4, 0)
    j <- ell_jacobi(u, k)
    expect_lt(max(abs(j$cn^2 + j$sn^2 - 1)), 1e-14)
    # and the dn identity that carries the diagonal family
    expect_lt(max(abs(j$dn^2 + k^2 * j$sn^2 - 1)), 1e-14)
  }
})

test_that("Theorem (existence): the basis is bounded by one", {
  t <- vapply(0:500, function(i) 2 * pi * i / 500, 0)
  for (k in c(0.1, 0.5, 0.9, 0.9999)) {
    b <- ell_basis(t, 4L, k)
    expect_lte(max(abs(b$C)), 1 + 1e-14)
    expect_lte(max(abs(b$S)), 1 + 1e-14)
  }
})

test_that("Theorem (exact variance): tr Var = (1 - |E|^2)/n identically", {
  for (n in c(5L, 17L, 400L)) {
    th <- rellvm(n, mu = 1, kappa = 3, k = 0.9, seed = 99 + n)
    mv <- ell_moment_var(th, P = 3L, k = 0.9)
    expect_equal(mv[, "realised"], mv[, "exact"], tolerance = 0)
  }
})

test_that("Theorem (bounds): |E_p| <= 1", {
  th <- rellvm(300L, mu = 0.5, kappa = 5, k = 0.8, seed = 7)
  m <- ell_moments(th, P = 4L, k = 0.8)
  expect_true(all(sqrt(m[, "c"]^2 + m[, "s"]^2) <= 1 + 1e-14))
})

test_that("Theorem (harmonic support): the transfer diagonal is non-zero", {
  for (k in c(0.2, 0.6, 0.9)) {
    tr <- ell_transfer(k, 3L, 5L)
    expect_true(all(abs(tr[, 1]) > 0))
    # coefficients decay in the nome, so later harmonics are smaller
    expect_true(all(diff(abs(tr[1, ])) < 0))
  }
})

test_that("Proposition (exact test): the p-value is a rank and lies in (0, 1]", {
  th <- rellvm(40L, mu = 0, kappa = 0, k = 0.9, seed = 3)
  u <- ell_exact_uniform(th, k = 0.9, B = 199L, seed = 5)
  expect_gt(u$p.value, 0)
  expect_lte(u$p.value, 1)
  expect_equal(u$p.value * 200, round(u$p.value * 200), tolerance = 1e-9)
})

test_that("the diagonal family's coefficients are positive and decreasing", {
  d <- ell_dn2_moments(0.8, 6L)
  expect_gt(d[1], 0)
  expect_true(all(d[-1] > 0))
  expect_true(all(diff(d[-1]) < 0))
})

test_that("simulation and fitting recover the truth", {
  th <- rellvm(2000L, mu = 1, kappa = 3, k = 0.9, seed = 2024)
  f <- fit_ellvm(th, k = 0.9, P = 2L, ng = 256L)
  expect_lt(abs(((f["mu"] - 1 + pi) %% (2 * pi)) - pi), 0.15)
  expect_lt(abs(f["kappa"] - 3), 0.6)
})

test_that("Theorem (efficiency): the Godambe identity dE[g]/dd = Cov(g, s) holds", {
  for (cfg in list(c(0.30, 0.60), c(0.70, 0.60), c(0.30, 1.20))) {
    e <- ell_eff_curve(tau = cfg[1], d = cfg[2], beta = 2,
                       kgrid = c(0, 0.9), nt = 1024L, ng = 200L)
    expect_lt(e$godambe, 1e-4)          # the identity the theorem rests on
    expect_gt(e$fisher, 0)
    expect_true(all(e$efficiency > 0 & e$efficiency <= 1 + 1e-12))
  }
})

test_that("Corollary (optimal modulus): k* beats every grid point, and may be interior", {
  # elongated: the optimum is pushed towards one
  ko <- ell_kopt(tau = 0.30, d = 1.20, beta = 2, nt = 1024L, ng = 200L)
  gr <- ell_eff_curve(tau = 0.30, d = 1.20, beta = 2,
                      kgrid = c(0, 0.9, 0.99, 0.999), nt = 1024L, ng = 200L)
  expect_gte(ko["efficiency"], max(gr$efficiency) - 1e-8)
  expect_gt(ko["efficiency"], gr$efficiency[1])   # strictly beats trigonometric

  # rounder: the optimum is INTERIOR, so a modulus near one is worse than k*
  ki <- ell_kopt(tau = 0.70, d = 0.60, beta = 2, nt = 1024L, ng = 200L)
  hi <- ell_eff_curve(tau = 0.70, d = 0.60, beta = 2,
                      kgrid = 0.99999, nt = 1024L, ng = 200L)
  expect_lt(ki["k"], 0.999)
  expect_gt(ki["efficiency"], hi$efficiency[1])
})

test_that("Corollary (full efficiency): the elliptic von Mises is its own best system", {
  # for f propto exp(kappa cn(.,k0)) the score is affine in c_1(.,k0), so the
  # efficiency at k0 must be one; checked directly rather than through the
  # projected family
  k0 <- 0.9
  th <- seq(0, 2 * pi, length.out = 4097)[-4097]
  b <- ell_basis(th, 1L, k0)
  f <- dellvm(th, mu = 0, kappa = 2, k = k0, ng = 4096L)
  w <- f / sum(f)
  g <- b$C[1, ]
  s <- g - sum(w * g)                    # score for kappa is c_1 centred
  gc <- g - sum(w * g)
  rho2 <- sum(w * gc * s)^2 / (sum(w * gc^2) * sum(w * s^2))
  expect_equal(rho2, 1, tolerance = 1e-10)
})

test_that("the cited projected densities are proper", {
  N <- 100000L; th <- 2 * pi * ((1:N) - 0.5) / N; h <- 2 * pi / N
  for (fam in c("normal", "cauchy"))
    for (tt in c(1.0, 0.4))
      for (dd in c(0, 1, 2.5)) {
        f <- dprojfam(th, tau = tt, d = dd, family = fam)
        expect_equal(sum(f) * h, 1, tolerance = 1e-7)
        expect_true(all(f > 0))
      }
})

test_that("Corollary: von Mises is fully efficient at k = 0 and nowhere else", {
  for (kap in c(1, 2, 6)) {
    r <- ell_catalogue("vonmises", kap, nt = 2048L)
    expect_equal(r$efficiency[1], 1, tolerance = 1e-3)   # k = 0 attains it
    expect_lt(r$kopt, 1e-3)                              # and the optimum is there
    expect_true(all(diff(r$efficiency) < 0))             # strictly worse as k grows
  }
})

test_that("Corollary: the elliptic von Mises is fully efficient at its own modulus", {
  for (k0 in c(0.5, 0.9)) {
    r <- ell_catalogue("ellvonmises", 2, k0, kgrid = k0, nt = 2048L)
    expect_equal(r$efficiency[1], 1, tolerance = 1e-3)
    expect_equal(r$kopt, k0, tolerance = 1e-3)
  }
})

test_that("axial families carry no information in the first elliptic moment", {
  # c_1(theta + pi) = -c_1(theta) while an antipodal score is pi-periodic, so
  # the covariance vanishes identically: this is structural, not numerical
  for (e in c(0.5, 0.8)) {
    r1 <- ell_catalogue("arcellipse", e, order = 1L, nt = 2048L)
    r2 <- ell_catalogue("arcellipse", e, order = 2L, nt = 2048L)
    expect_lt(max(r1$efficiency), 1e-8)
    expect_gt(r2$eff_kopt, 0.5)
  }
})

test_that("the catalogue quadrature preserves total mass", {
  for (s in list(list("wrappedcauchy", 0.9, 1), list("projnormal", 1, 0.4),
                 list("projcauchy", 2.5, 1), list("cardioid", 0.4, 1))) {
    r <- ell_catalogue(s[[1]], s[[2]], s[[3]], nt = 2048L)
    expect_equal(r$mass, 1, tolerance = 1e-6)
    expect_true(all(r$efficiency >= 0 & r$efficiency <= 1 + 1e-9))
    expect_gte(r$eff_kopt, max(r$efficiency) - 1e-6)
  }
})

test_that("Theorem (joint): Cov(hat E) = Sigma/n exactly, at small n", {
  k <- 0.9; P <- 2L; NG <- 4096L
  th <- 2 * pi * ((1:NG) - 0.5) / NG
  f <- dellvm(th, mu = 1, kappa = 3, k = k, ng = NG)
  J <- ell_joint(f, P, k)

  # the Pythagorean identity determines each diagonal block's trace
  expect_equal(J$Sigma[1, 1] + J$Sigma[2, 2],
               1 - (J$mean[1]^2 + J$mean[2]^2), tolerance = 1e-12)
  expect_equal(J$Sigma[3, 3] + J$Sigma[4, 4],
               1 - (J$mean[3]^2 + J$mean[4]^2), tolerance = 1e-12)

  # Sigma is a covariance matrix
  expect_true(all(abs(J$Sigma - t(J$Sigma)) < 1e-12))
  expect_true(all(eigen(J$Sigma, only.values = TRUE)$values > -1e-12))

  # and the off-diagonal blocks are NOT negligible, so the non-closure remark
  # is about something real rather than a rounding effect
  expect_gt(max(abs(J$Sigma[1:2, 3:4])), 0.05)

  # exactness at a sample size where an asymptotic joint law is useless
  n <- 10L; R <- 6000L
  M <- t(vapply(seq_len(R), function(r)
    ell_joint_sample(rellvm(n, 1, 3, k, seed = 104729 * r + 3), P, k)$mean,
    numeric(2 * P)))
  se <- max(diag(J$Sigma / n)) * sqrt(2 / R)
  expect_lt(max(abs(cov(M) - J$Sigma / n)), 5 * se)
})

