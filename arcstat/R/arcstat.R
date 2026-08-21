## arcstat: arc-length statistics -- goodness of fit, distributions, and a Bayesian test.
## One shared pure-C back-end (compiled at install), also bound from Python (arcstat on PyPI).
#' @useDynLib arcstat, .registration = TRUE, .fixes = "C_"
#' @importFrom stats punif approx runif optim quantile rgamma
#' @keywords internal
"_PACKAGE"

## ================= 1. Goodness-of-fit test ============================================
## alGOF: arc-length goodness-of-fit test. R binding to the bundled C back-end (compiled at install).
## Successor to the withdrawn 'alR' package (whose arc-length regression objective was invalid).

#' Arc-length goodness-of-fit statistic (C back-end)
#' @param u numeric probability-integral transforms in [0,1].
#' @return the scalar arc length of the probability-plot ogive.
#' @export
al_statistic <- function(u) {
  if (any(u < 0 | u > 1)) stop("`u` must be probability-integral transforms in [0,1]")
  .C(C_al_statistic, n = as.integer(length(u)), u = as.double(u), out = double(1))$out
}

#' Saddlepoint right-tail probability of the arc-length statistic (C back-end)
#' @param s observed statistic value;
#' @param n sample size.
#' @return the right-tail p-value.
#' @export
al_pvalue <- function(s, n)
  .C(C_al_pvalue, s = as.double(s), n = as.integer(n), out = double(1))$out

#' Exact mean, variance and support of the arc-length statistic (C back-end)
#' @param n sample size.
#' @return list with mean, var, support.
#' @export
al_moments <- function(n) {
  r <- .C(C_al_moments, n = as.integer(n), out = double(4))$out
  list(mean = r[1], var = r[2], support = r[3:4])
}

#' Arc-length goodness-of-fit test
#'
#' Tests H0: F = F0 via the arc length of the probability plot, with the analytic saddlepoint null.
#' Powerful against local density structure (multimodality, clustering, oscillation, heaping) that the
#' empirical-distribution tests miss; weak against smooth location/scale departures.
#' @param x data;
#' @param null CDF F0 giving the PIT (default punif);
#' @param nboot optional bootstrap reps.
#' @param rnull generator matching null (needed with nboot).
#' @return an object of class "htest".
#' @export
al_test <- function(x, null = stats::punif, nboot = NULL, rnull = NULL) {
  x <- x[is.finite(x)]; n <- length(x)
  u <- pmin(pmax(null(sort(x)), 1e-12), 1 - 1e-12); S <- al_statistic(u)
  if (is.null(nboot)) {
    p <- al_pvalue(S, n)
    method <- "Arc-length goodness-of-fit test (saddlepoint null)"
    if (is.na(p) || p <= 1e-12 || p > 1) {   # extreme tail: saddlepoint out of range; brief exact simulation
      Snull <- replicate(20000L, al_statistic(sort(stats::runif(n))))
      p <- (sum(Snull >= S) + 1) / 20001
      method <- "Arc-length goodness-of-fit test (saddlepoint null; extreme tail by simulation)"
    }
  }
  else {
    if (is.null(rnull)) stop("supply `rnull` when using `nboot`")
    Snull <- replicate(nboot, al_statistic(pmin(pmax(null(sort(rnull(n))), 1e-12), 1 - 1e-12)))
    p <- (sum(Snull >= S) + 1) / (nboot + 1)
    method <- sprintf("Arc-length goodness-of-fit test (parametric bootstrap, %d reps)", nboot)
  }
  structure(list(statistic = c(S = S), p.value = p, method = method,
                 data.name = deparse(substitute(x)), alternative = "local density structure present"),
            class = "htest")
}

## ================= 2. Arc-length distributions ========================================
## arcdist: arc-length transforms in distribution theory
## The numerical primitives (arc-length quadrature, sample L-moments) are computed by a shared C
## back-end compiled at install and bound through .C; the same sources back the Python package
## arcdist-py. Grid-based object numerics stay in R; support is compact by construction.
NULL

## ---- shifted Legendre polynomials on [0,1] (P_0..P_k) --------------------------------
.shifted_legendre <- function(u, k) {
  # returns a matrix length(u) x (k+1): columns P_0(u),...,P_k(u)
  P <- matrix(1, length(u), k + 1L)
  if (k >= 1L) P[, 2L] <- 2 * u - 1
  if (k >= 2L) for (r in 2L:k) {
    # recurrence for shifted Legendre: (r) P_r = (2r-1)(2u-1) P_{r-1} - (r-1) P_{r-2}
    P[, r + 1L] <- ((2 * r - 1) * (2 * u - 1) * P[, r] - (r - 1) * P[, r - 1L]) / r
  }
  P
}

## ================= Quantile arc-length family =========================================

#' Construct a quantile arc-length distribution
#'
#' The quantile density is Q'(u) = sigma * [1 + sum_k coef[k] P_k(u)]_+, with P_k the shifted
#' Legendre polynomials; the quantile function is Q(u) = mu + integral_0^u Q'(v) dv. Support is
#' compact. coef = numeric(0) gives the uniform on [mu, mu+sigma].
#'
#' @param coef shape coefficients (c_1, c_2, ...); may be empty for the uniform.
#' @param mu location (left endpoint of support).
#' @param sigma positive scale.
#' @param ngrid grid resolution for the internal numerics.
#' @return an object of class "arcq".
#' @examples
#' ## an empty coefficient vector gives the uniform on [mu, mu + sigma]
#' arclength(arcq(numeric(0), mu = 0, sigma = 1))
#'
#' ## the quantile function is increasing and the distribution function inverts it
#' o <- arcq(0.5)
#' pu <- c(0.1, 0.35, 0.6, 0.9)
#' max(abs(parcq(qarcq(pu, o), o) - pu))
#' @export
arcq <- function(coef = numeric(0), mu = 0, sigma = 1, ngrid = 4000L) {
  stopifnot(sigma > 0, ngrid > 10)
  u <- seq(0, 1, length.out = ngrid); du <- u[2] - u[1]
  k <- length(coef)
  shape <- rep(1, ngrid)
  if (k > 0) { P <- .shifted_legendre(u, k); shape <- 1 + as.numeric(P[, -1L, drop = FALSE] %*% coef) }
  qd <- sigma * pmax(0, shape)                 # quantile density Q'(u) >= 0
  Q  <- mu + cumsum(qd) * du - qd[1] * du       # Q(0) = mu
  structure(list(coef = coef, mu = mu, sigma = sigma, u = u, du = du,
                 qd = qd, Q = Q, a = sqrt(1 + qd^2)), class = "arcq")
}

#' @describeIn arcq quantile function
#' @param p probabilities.
#' @param obj an "arcq" object.
#' @export
qarcq <- function(p, obj) approx(obj$u, obj$Q, xout = p, rule = 2)$y

#' @describeIn arcq cumulative distribution function
#' @param x quantiles.
#' @export
parcq <- function(x, obj) approx(obj$Q, obj$u, xout = x, rule = 2)$y

#' @describeIn arcq density
#' @export
darcq <- function(x, obj) {
  uu <- parcq(x, obj); qd <- approx(obj$u, obj$qd, xout = uu, rule = 2)$y
  ifelse(qd > 0, 1 / qd, 0)
}

#' @describeIn arcq random generation
#' @param n number of draws.
#' @export
rarcq <- function(n, obj) qarcq(stats::runif(n), obj)

#' @describeIn arcq total arc length of the quantile curve (curve complexity)
#'
#' The integrand sqrt(1+Q'(u)^2) is analytic on [0,1], so Gauss-Legendre converges geometrically
#' where the equally spaced grid converges linearly: sixteen nodes reach machine precision on a
#' quadratic quantile density, against a relative error of 3e-04 for four thousand grid points.
#' Set exact = FALSE to recover the old grid value.
#'
#' The integral is generally not elementary. With Q' of degree one the integrand is the square root
#' of a quadratic and the antiderivative is an inverse hyperbolic sine; with Q' of degree two it is
#' the square root of a quartic, an elliptic integral, which by Liouville's theorem has no
#' elementary antiderivative; beyond that it is hyperelliptic.
#'
#' @param exact use Gauss-Legendre quadrature (default) rather than the internal grid.
#' @param nodes number of Gauss-Legendre nodes.
#' @return A single number, the arc length of the quantile function.
#' @examples
#' ## Arc length is a SHAPE functional, so shifting the distribution cannot change it.
#' arclength(arcq(0.5, mu = 0, sigma = 1))
#' arclength(arcq(0.5, mu = 17, sigma = 1))
#'
#' ## The quadrature route is spectral: the value is settled at a handful of nodes, where
#' ## the equally spaced grid it replaced converges only linearly.
#' o <- arcq(c(0.5, -0.3))
#' sapply(c(8L, 16L, 32L), function(k) arclength(o, nodes = k))
#' @export
arclength <- function(obj, exact = TRUE, nodes = 24L) {
  if (!exact) return(sum(obj$a * obj$du))
  .C(C_arcq_arclength, coef = as.double(obj$coef), k = as.integer(length(obj$coef)),
     sigma = as.double(obj$sigma), nodes = as.integer(nodes), out = double(1))$out    # shared C back-end
}

## Gauss-Legendre nodes and weights on [0,1] (Newton on the Legendre polynomial)
.gauss_legendre <- function(n) {
  i <- seq_len(n); x <- cos(pi * (i - 0.25) / (n + 0.5))
  for (it in 1:100) {
    p0 <- rep(1, n); p1 <- x
    for (k in 2:n) { p2 <- ((2 * k - 1) * x * p1 - (k - 1) * p0) / k; p0 <- p1; p1 <- p2 }
    dp <- n * (x * p1 - p0) / (x^2 - 1); dx <- p1 / dp; x <- x - dx
    if (max(abs(dx)) < 1e-15) break
  }
  p0 <- rep(1, n); p1 <- x
  for (k in 2:n) { p2 <- ((2 * k - 1) * x * p1 - (k - 1) * p0) / k; p0 <- p1; p1 <- p2 }
  dp <- n * (x * p1 - p0) / (x^2 - 1); w <- 2 / ((1 - x^2) * dp^2)
  list(x = (x + 1) / 2, w = w / 2)
}

#' Theoretical L-moments of an arcq distribution
#'
#' L_r is the shifted-Legendre projection of Q (Hosking's convention): L_r =
#' integral_0^1 Q(u) P*_(r-1)(u) du, r>=1, with P*_r the shifted Legendre polynomial (P*_r(1)=1).
#' @param obj an "arcq" object.
#' @param nmom number of L-moments.
#' @return the first nmom L-moments; ratios tau_r = L_r/L_2 for r>=3 via lmratios().
#' @export
lmoments <- function(obj, nmom = 4L) {
  P <- .shifted_legendre(obj$u, nmom - 1L)
  n <- length(obj$u)
  w <- rep(obj$du, n); w[1] <- w[1]/2; w[n] <- w[n]/2   # trapezoid: drops the O(1/n) endpoint bias
  m <- sum(obj$Q * w)                                   # L_1 is the mean
  if (nmom < 2L) return(m)
  ## L_r for r >= 2 is location-invariant, so it is taken from the CENTRED quantile function.
  ## Projecting the uncentred Q instead multiplies the quadrature error of \int P*_{r-1} = 0 by the
  ## location, and that residual does not vanish on a finite grid: for a distribution with mean 576
  ## and L-scale 0.75 it reverses the sign of L-skewness. Centring is exact, not an approximation.
  Qc <- obj$Q - m
  c(m, vapply(2:nmom, function(r) sum(Qc * P[, r] * w), numeric(1)))
}

#' L-moment ratios (tau_3, tau_4, ...) from a vector of L-moments
#' @param L L-moments.
#' @return A numeric vector: the first two entries are the location and scale L-moments
#'   themselves, followed by the ratios \eqn{\tau_3, \tau_4, \ldots}.
#' @examples
#' ## a symmetric sample has vanishing L-skewness
#' u <- (seq_len(20000) - 0.5) / 20000
#' lmratios(sample_lmoments(u, 4L))[3]
#' @export
lmratios <- function(L) c(L[1], L[2], L[-(1:2)] / L[2])

## ---- sample L-moments (Hosking) -----------------------------------------------------
#' Sample L-moments of data
#' @param x data.
#' @param nmom number of L-moments.
#' @return A numeric vector of \code{nmom} sample L-moments.
#' @examples
#' ## for a uniform grid on (0,1) the first L-moment is 1/2 and the second is 1/6
#' u <- (seq_len(20000) - 0.5) / 20000
#' sample_lmoments(u, 4L)[1:2]
#' @export
sample_lmoments <- function(x, nmom = 4L) {
  .C(C_sample_lmoments_c, x = as.double(x), n = as.integer(length(x)),
     nmom = as.integer(nmom), out = double(nmom))$out                                 # shared C back-end
}

#' Closed-form L-moment fit of the order-two arcq family
#'
#' The family's L-moments are an exact linear function of its shape coefficients, so matching them
#' inverts explicitly: c2 = 35 t4 / (3 + 7 t4) and c1 = t3 (5 - c2), with sigma and mu then absorbing
#' the sample L-mean and L-scale. No optimisation, starting value or tolerance is involved, and the
#' estimator is consistent and asymptotically normal by the delta method. This is the estimator of
#' record for the order-two family; `fit_arcq` remains available for higher orders, where the
#' inversion is done numerically.
#' @param x data.
#' @return an `arcq` object, with a `fit` component holding the coefficients, `mu`, `sigma` and an
#'   `admissible` flag; a warning is issued when the estimate leaves the admissible set.
#' @export
fit_arcq_cf <- function(x) {
  o <- .C(C_arcq_fit_cf_c, x = as.double(x), n = as.integer(length(x)),
          out = double(5))$out                                                        # shared C back-end
  if (anyNA(o)) stop("closed-form fit undefined: the sample L-scale is zero")
  adm <- o[5] > 0.5
  if (!adm)                       # the inversion is derived on the admissible set only
    warning("closed-form estimate is outside the admissible set: the positive part is active, ",
            "so the fitted law's L-moments are not the values matched")
  fit <- arcq(o[1:2], o[3], o[4])
  fit$fit <- list(coef = o[1:2], mu = o[3], sigma = o[4], admissible = adm)
  fit
}

#' Fit an arcq distribution to data by matching L-moments
#'
#' Matches the first (order+2) L-moments: mu and sigma absorb L_1 and L_2, and the `order` shape
#' coefficients are chosen to match tau_3, tau_4, ... by least squares.
#' @param x data.
#' @param order number of shape coefficients (>=1).
#' @return An \code{arcq} object with the fitted location, scale and shape coefficients.
#' @examples
#' set.seed(1)
#' fit <- fit_arcq(rnorm(500), order = 2L)
#' arclength(fit)
#' @export
fit_arcq <- function(x, order = 2L) {
  Ls <- sample_lmoments(x, order + 2L); tau_target <- Ls[-(1:2)] / Ls[2]
  obj_ratios <- function(coef) {
    o <- arcq(coef, 0, 1); Lt <- lmoments(o, order + 2L); Lt[-(1:2)] / Lt[2]
  }
  fn <- function(coef) sum((obj_ratios(coef) - tau_target)^2)
  opt <- stats::optim(rep(0, order), fn, method = "BFGS", control = list(reltol = 1e-10))
  o0 <- arcq(opt$par, 0, 1); L0 <- lmoments(o0, 2L)
  sigma <- Ls[2] / L0[2]; mu <- Ls[1] - sigma * L0[1]
  fit <- arcq(opt$par, mu, sigma)
  fit$fit <- list(coef = opt$par, mu = mu, sigma = sigma, tau_target = tau_target,
                  tau_fitted = obj_ratios(opt$par), convergence = opt$convergence)
  fit
}

## ================= Arc-length generator ==============================================

#' Arc-length generator: transform a bounded-support base into a new distribution
#'
#' G(x) = S_F^[a,x] / S_F^[a,b], the normalised cumulative arc length of the base cumulative
#' distribution function, whose density is proportional to sqrt(1 + f(x)^2). Requires bounded support.
#' @param dens base density function (vectorised).
#' @param lower,upper bounded support.
#' @param ngrid grid resolution.
#' @return an object of class "arcgen".
#' @export
arc_generator <- function(dens, lower, upper, ngrid = 4000L) {
  stopifnot(is.finite(lower), is.finite(upper), upper > lower)
  x <- seq(lower, upper, length.out = ngrid); dx <- x[2] - x[1]
  f <- dens(x); speed <- sqrt(1 + f^2)
  L <- cumsum(speed) * dx - speed[1] * dx; S <- L[ngrid]
  G <- L / S; g <- speed / S
  structure(list(x = x, dx = dx, G = G, g = g, base_f = f, total_arclength = S,
                 lower = lower, upper = upper), class = "arcgen")
}

#' @describeIn arc_generator cumulative distribution function of the generated law
#' @param q quantiles.
#' @param obj an "arcgen" object.
#' @export
pgen <- function(q, obj) approx(obj$x, obj$G, xout = q, rule = 2)$y

#' @describeIn arc_generator density of the generated law
#' @export
dgen <- function(q, obj) approx(obj$x, obj$g, xout = q, rule = 2)$y

#' @describeIn arc_generator random generation from the generated law
#' @param n draws.
#' @export
rgen <- function(n, obj) approx(obj$G, obj$x, xout = stats::runif(n), rule = 2)$y

## ================= 3. Bayesian arc-length test ========================================
## R (.C) binding to the shared Bayesian arc-length C back-end (libbayesarc).

#' Bayesian-bootstrap posterior arc-length discrepancies
#'
#' Draws \code{M} Bayesian-bootstrap (Dirichlet-process, concentration to zero) posterior samples of the
#' arc-length and Kolmogorov--Smirnov discrepancies of the probability-plot ogive for the probability
#' integral transforms \code{u}.
#' @param u numeric probability-integral transforms in [0,1].
#' @param M number of posterior draws.
#' @param seed integer seed for the C back-end.
#' @return an \code{M} by 2 matrix with columns \code{arc} and \code{ks}.
#' @export
bb_post_disc <- function(u, M, seed = 1L) {
  u <- sort(u)
  r <- .C(C_bb_post, n = as.integer(length(u)), u = as.double(u), M = as.integer(M),
          seed = as.integer(seed), out = double(2L * M))$out
  matrix(r, ncol = 2, byrow = TRUE, dimnames = list(NULL, c("arc", "ks")))
}

#' Null reference distribution of the arc-length discrepancy
#'
#' Pools the posterior discrepancies of \code{D} uniform data sets (\code{m} draws each) to give the
#' distribution of the discrepancy under the null hypothesis at sample size \code{n}.
#' @param n sample size.
#' @param D number of uniform reference data sets.
#' @param m posterior draws per data set.
#' @param seed integer seed for the C back-end.
#' @return a matrix with columns \code{arc} and \code{ks}.
#' @export
bb_ref_disc <- function(n, D = 300L, m = 16L, seed = 7L) {
  r <- .C(C_bb_ref, n = as.integer(n), D = as.integer(D), m = as.integer(m),
          seed = as.integer(seed), out = double(2L * D * m))$out
  matrix(r, ncol = 2, byrow = TRUE, dimnames = list(NULL, c("arc", "ks")))
}

#' Bayesian arc-length goodness-of-fit evidence
#'
#' The posterior probability that the discrepancy exceeds its \code{(1 - level)} quantile under the null;
#' values near one are strong evidence against the hypothesised distribution. Sensitive to local density
#' structure (multimodality, clustering, heaping) that the vertical discrepancies miss.
#' @param u numeric probability-integral transforms in [0,1].
#' @param which discrepancy to use, \code{"arc"} (default) or \code{"ks"}.
#' @param level nominal level defining the reference quantile.
#' @param M number of posterior draws.
#' @param ref optional precomputed reference from \code{bb_ref_disc}; computed if \code{NULL}.
#' @param seed integer seed for the C back-end.
#' @return the posterior evidence against the null, in [0,1].
#' @export
bb_evidence <- function(u, which = "arc", level = 0.95, M = 1500L, ref = NULL, seed = 1L) {
  n <- length(u); if (is.null(ref)) ref <- bb_ref_disc(n)
  thr <- as.numeric(quantile(ref[, which], level))
  mean(bb_post_disc(u, M, seed)[, which] > thr)
}

#' Arc length of the empirical characteristic function (C back-end)
#'
#' The empirical characteristic function of a sample traces a curve in the complex plane; this returns
#' the arc length of that curve over a window \eqn{[0, T]}, doubled for the symmetric half. The sample is
#' standardised internally, so the value is scale free. It is a descriptive, asymmetry-sensitive summary:
#' a symmetric law tends to the value two, and asymmetry adds length.
#' @param x numeric data vector.
#' @param T upper end of the frequency window.
#' @param ngrid number of grid points on \eqn{[0, T]}.
#' @return the windowed arc length of the empirical characteristic function.
#' @export
cf_arclength <- function(x, T = 6, ngrid = 1200L) {
  x <- as.double(x[is.finite(x)])
  .C(C_cf_arclength_emp, x = x, n = length(x), T = as.double(T),
     ng = as.integer(ngrid), out = double(1))$out
}

#' Closed-form characteristic-function arc length for named families
#'
#' Total arc length of the characteristic-function curve for the families with a closed form:
#' \code{"exponential"} (\eqn{\pi}), \code{"gamma"} (shape \code{k}), \code{"normal"} (drift ratio
#' \code{delta}), \code{"cauchy"} (drift ratio \code{delta}), \code{"skewstable"} (index \code{alpha},
#' skewness \code{beta}), and \code{"poisson"} (rate \code{lambda}, per period). Total arc length is scale
#' free; a symmetric monotone (Polya) law carries the value two.
#' @param family one of the names above.
#' @param k gamma shape.
#' @param delta location-to-scale ratio for the normal and Cauchy.
#' @param alpha,beta stability index and skewness for the skew-stable family.
#' @param lambda Poisson rate (per-period arc length).
#' @return the closed-form arc length.
#' @examples
#' ## two closed forms that anchor the whole construction
#' cf_arclength_family("normal")                    # exactly 2
#' cf_arclength_family("exponential", lambda = 1)   # exactly pi
#'
#' ## total arc length is scale free, so the exponential rate cannot matter
#' cf_arclength_family("exponential", lambda = 2)
#' @export
cf_arclength_family <- function(family = c("exponential", "gamma", "normal", "cauchy",
                                           "skewstable", "poisson"),
                                k = 1, delta = 0, alpha = 1.5, beta = 0, lambda = 1) {
  family <- match.arg(family)
  switch(family,
    exponential = pi,
    gamma = 2 * sqrt(pi) * gamma(k / 2 + 1) / gamma((k + 1) / 2),
    normal = { x <- delta^2 / 4
      if (x == 0) 2 else 2 * x * exp(x) * (besselK(x, 0) + besselK(x, 1)) },
    cauchy = 2 * sqrt(1 + delta^2),
    skewstable = 2 * sqrt(1 + (beta * tan(pi * alpha / 2))^2),
    poisson = 2 * pi * lambda * exp(-lambda) * besselI(lambda, 0))
}

#' Band arc length of the normal reference curve (C back-end)
#'
#' @param sigma scale of the normal reference;
#' @param a,b the probability band;
#' @param nodes quadrature panels.
#' @return the model band arc length.
#' @export
al_band_model <- function(sigma, a = 0.05, b = 0.95, nodes = 400L)
  .C(C_al_band_model, sigma = as.double(sigma), a = as.double(a), b = as.double(b),
     nodes = as.integer(nodes), out = double(1))$out

#' Sample band arc length of the empirical distribution curve (C back-end)
#'
#' @param x data;
#' @param a,b the probability band.
#' @return the sample band arc length, NA when fewer than two points fall in the band.
#' @export
al_band_sample <- function(x, a = 0.05, b = 0.95) {
  x <- as.double(x[is.finite(x)])
  .C(C_al_band_sample, x = x, n = length(x), a = as.double(a), b = as.double(b), out = double(1))$out
}

#' Scale by arc-length band matching, in the scale-equivariant standardised form (C back-end)
#'
#' The raw matching equation fixes an aspect ratio between the probability and response axes and is
#' therefore unit dependent; this estimator divides by the MAD, matches on the standardised scale and
#' rescales, which is exactly scale equivariant.
#' @param x data;
#' @param a,b the probability band, whose tail mass sets the breakdown point.
#' @return the estimated scale, NA when the MAD vanishes or the matching equation has no root.
#' @export
al_scale <- function(x, a = 0.05, b = 0.95) {
  x <- as.double(x[is.finite(x)])
  .C(C_al_scale, x = x, n = length(x), a = as.double(a), b = as.double(b), out = double(1))$out
}

## ---------------------------------------------------------------------------------------------
## arccirc: the circular (wrapped) arc-length family
## ---------------------------------------------------------------------------------------------
## The linear quantile arc-length family has Q(0) = 0 and Q(1) = 1 identically, so theta = 2 pi Q(U)
## covers the circle exactly once and the wrapping series has a single term. Smoothness across the
## join costs two linear constraints on the shifted-Legendre coefficients, which at order two leave
## only the uniform; order three is the first smooth non-uniform member, with one shape parameter.
## All numerical work is in the shared C back-end, which also serves the Python front end.

#' Largest admissible shape parameter of the smooth circular family
#'
#' @return the exact bound 3 sqrt(3) / 5 on |c3|, beyond which the quantile density goes negative.
#' @export
arcc_c3max <- function() .C(C_arcc_c3max, out = double(1))$out

#' The smooth circular arc-length family
#'
#' @param c3 shape parameter, |c3| <= arcc_c3max(); c3 = 0 is the circular uniform.
#' @param mu mean direction in radians.
#' @return an object of class "arccirc".
#' @export
arccirc <- function(c3 = 0, mu = 0) {
  stopifnot(length(c3) == 1L, length(mu) == 1L)
  if (abs(c3) > arcc_c3max() + 1e-12)
    stop("c3 outside the admissible range: |c3| must not exceed 3*sqrt(3)/5")
  structure(list(c3 = c3, mu = mu), class = "arccirc")
}

#' @describeIn arccirc quantile density q(u) of the underlying linear family
#' @param u probabilities in [0,1].
#' @param obj an "arccirc" object.
#' @export
qd_arccirc <- function(u, obj)
  .C(C_arcc_qd3, u = as.double(u), nu = as.integer(length(u)),
     c3 = as.double(obj$c3), out = double(length(u)))$out

#' @describeIn arccirc quantile function Q(u) of the underlying linear family
#' @export
Q_arccirc <- function(u, obj)
  .C(C_arcc_Q3, u = as.double(u), nu = as.integer(length(u)),
     c3 = as.double(obj$c3), out = double(length(u)))$out

#' @describeIn arccirc circular density
#' @param theta angles in radians.
#' @export
darccirc <- function(theta, obj)
  .C(C_arcc_dens3, theta = as.double(theta), nth = as.integer(length(theta)),
     c3 = as.double(obj$c3), mu = as.double(obj$mu), out = double(length(theta)))$out

#' @describeIn arccirc trigonometric moment of order p
#' @param p integer order.
#' @param nodes number of rectangle-rule nodes; the integrand is C^1 periodic, so convergence is
#'   geometric and the default is ample.
#' @export
trigmom_arccirc <- function(p, obj, nodes = 4096L) {
  z <- .C(C_arcc_trigmom3, p = as.integer(p), c3 = as.double(obj$c3),
          mu = as.double(obj$mu), nodes = as.integer(nodes), out = double(2))$out
  complex(real = z[1], imaginary = z[2])
}

#' @describeIn arccirc mean resultant length
#' @export
rho_arccirc <- function(obj, nodes = 4096L) Mod(trigmom_arccirc(1L, obj, nodes))

#' @describeIn arccirc random generation
#' @param n sample size.
#' @export
rarccirc <- function(n, obj)
  .C(C_arcc_rand3, unif = as.double(stats::runif(n)), n = as.integer(n),
     c3 = as.double(obj$c3), mu = as.double(obj$mu), out = double(n))$out

#' Method-of-moments fit of the smooth circular arc-length family
#'
#' The modulus of the first trigonometric moment is even in c3 and strictly increasing in |c3|, so it
#' identifies the magnitude; the rotation-invariant psi = arg(phi_2) - 2 arg(phi_1), reduced to the
#' principal branch, identifies the sign. Reliable sign recovery needs a few thousand observations.
#'
#' @param theta angles in radians.
#' @param nodes rectangle-rule nodes used inside the inversion.
#' @return a list with mu, c3, the observed |phi_1|, and whether the inversion stayed interior.
#' @export
fit_arccirc <- function(theta, nodes = 4096L) {
  o <- .C(C_arcc_fit3, theta = as.double(theta), n = as.integer(length(theta)),
          nodes = as.integer(nodes), out = double(4))$out
  list(mu = o[1], c3 = o[2], rho = o[3], interior = o[4] == 1)
}

#' Fejer-Riesz form of the circular arc-length family
#'
#' q(u) = |P(exp(2 pi i u))|^2 normalised to integrate to one. Every such q is non-negative and
#' infinitely differentiable around the circle by construction, so neither an admissibility condition
#' nor a smoothness constraint arises, and the coefficients are unrestricted.
#'
#' @param p complex vector of polynomial coefficients.
#' @param mu mean direction.
#' @return an object of class "arccirc_fr".
#' @export
arccirc_fr <- function(p, mu = 0) structure(list(p = as.complex(p), mu = mu), class = "arccirc_fr")

#' @describeIn arccirc_fr quantile density
#' @param u probabilities in [0,1].
#' @param obj an "arccirc_fr" object.
#' @export
qd_arccirc_fr <- function(u, obj)
  .C(C_arcc_qd_fr, u = as.double(u), nu = as.integer(length(u)),
     pre = as.double(Re(obj$p)), pim = as.double(Im(obj$p)),
     np = as.integer(length(obj$p)), out = double(length(u)))$out

#' @describeIn arccirc_fr trigonometric moment of order p
#' @param pord integer order.
#' @param nodes rectangle-rule nodes.
#' @export
trigmom_arccirc_fr <- function(pord, obj, nodes = 4096L) {
  z <- .C(C_arcc_trigmom_fr, p = as.integer(pord),
          pre = as.double(Re(obj$p)), pim = as.double(Im(obj$p)),
          np = as.integer(length(obj$p)), mu = as.double(obj$mu),
          nodes = as.integer(nodes), out = double(2))$out
  complex(real = z[1], imaginary = z[2])
}

## ---- Fejer-Riesz family: quantile function, density, simulation, and circular L-moments --------

#' @describeIn arccirc_fr quantile function of the underlying linear family
#' @param nodes Simpson cells used to accumulate Q.
#' @export
Q_arccirc_fr <- function(u, obj, nodes = 4096L)
  .C(C_arcc_Q_fr, u = as.double(u), nu = as.integer(length(u)),
     pre = as.double(Re(obj$p)), pim = as.double(Im(obj$p)), np = as.integer(length(obj$p)),
     nodes = as.integer(nodes), out = double(length(u)))$out

#' @describeIn arccirc_fr circular density
#' @param theta angles in radians.
#' @export
darccirc_fr <- function(theta, obj, nodes = 4096L)
  .C(C_arcc_dens_fr, theta = as.double(theta), nth = as.integer(length(theta)),
     pre = as.double(Re(obj$p)), pim = as.double(Im(obj$p)), np = as.integer(length(obj$p)),
     mu = as.double(obj$mu), nodes = as.integer(nodes), out = double(length(theta)))$out

#' @describeIn arccirc_fr random generation
#' @param n sample size.
#' @export
rarccirc_fr <- function(n, obj, nodes = 4096L)
  .C(C_arcc_rand_fr, unif = as.double(stats::runif(n)), n = as.integer(n),
     pre = as.double(Re(obj$p)), pim = as.double(Im(obj$p)), np = as.integer(length(obj$p)),
     mu = as.double(obj$mu), nodes = as.integer(nodes), out = double(n))$out

#' Circular L-moments
#'
#' On the line the L-moments are Legendre projections of the quantile function and are linear in the
#' Legendre coefficients of the quantile density. On the circle the natural basis is Fourier, and the
#' m-th circular L-moment is the Fourier projection
#' \eqn{\ell_m = \int_0^1 Q(u) e^{-2\pi i m u}\,du}. With the empirical quantile function this is a
#' linear combination of order statistics with fixed complex weights, hence a genuine L-statistic.
#' The map to the quantile-density spectrum is exact and inverts,
#' \eqn{\ell_m = i(1-\rho_m)/(2\pi m)}, equivalently \eqn{\rho_m = 1 + 2\pi i m \ell_m}.
#'
#' The general device of projecting a quantile function onto an orthogonal basis is not new: it is
#' Sillitto (1969) in one dimension and Decurninge (2014) for multivariate quantile maps. What is
#' specific here is the Fourier form and its exact inversion to the Fejer-Riesz parameters.
#'
#' @param x observations on the unit interval, obtained by cutting the circle at the mean direction.
#' @param nm number of L-moments.
#' @return complex vector of length nm.
#' @export
lmom_circ <- function(x, nm = 3L) {
  o <- .C(C_arcc_lmom, x = as.double(x), n = as.integer(length(x)),
          nm = as.integer(nm), out = double(2 * nm))$out
  complex(real = o[seq(1, 2 * nm, by = 2)], imaginary = o[seq(2, 2 * nm, by = 2)])
}

#' @describeIn lmom_circ the exact linear inversion to the quantile-density spectrum
#' @param ell complex vector of circular L-moments.
#' @export
rho_from_lmom <- function(ell) {
  nm <- length(ell)
  o <- .C(C_arcc_rho_from_lmom, ellre = as.double(Re(ell)), ellim = as.double(Im(ell)),
          nm = as.integer(nm), out = double(2 * nm))$out
  complex(real = o[seq(1, 2 * nm, by = 2)], imaginary = o[seq(2, 2 * nm, by = 2)])
}

#' @describeIn lmom_circ smallest implied quantile density; negative means the estimate has left the
#'   admissible set
#' @param rho complex vector of spectrum values.
#' @param nodes grid points.
#' @export
qmin_rho <- function(rho, nodes = 4096L)
  .C(C_arcc_qmin_rho, rhore = as.double(Re(rho)), rhoim = as.double(Im(rho)),
     nm = as.integer(length(rho)), nodes = as.integer(nodes), out = double(1))$out

#' @describeIn lmom_circ shrink a spectrum towards the circular uniform until admissible
#' @param margin required minimum of the implied quantile density.
#' @export
admiss_rho <- function(rho, margin = 0, nodes = 4096L) {
  nm <- length(rho)
  o <- .C(C_arcc_admiss, rhore = as.double(Re(rho)), rhoim = as.double(Im(rho)),
          nm = as.integer(nm), margin = as.double(margin), nodes = as.integer(nodes),
          out = double(2 * nm + 1))$out
  list(rho = complex(real = o[seq(1, 2 * nm, by = 2)], imaginary = o[seq(2, 2 * nm, by = 2)]),
       shrink = o[2 * nm + 1])
}

#' @describeIn lmom_circ Fejer-Riesz spectral factorisation, recovering the coefficient vector
#' @export
factorise_rho <- function(rho) {
  nm <- length(rho)
  o <- .C(C_arcc_factorise, rhore = as.double(Re(rho)), rhoim = as.double(Im(rho)),
          nm = as.integer(nm), out = double(2 * (nm + 1)))$out
  complex(real = o[seq(1, 2 * (nm + 1), by = 2)], imaginary = o[seq(2, 2 * (nm + 1), by = 2)])
}

#' Closed-form fit of the Fejer-Riesz circular arc-length family
#'
#' Sort, apply the circular L-moment weights, invert the exact linear map, factorise: no numerical
#' search enters anywhere. The pair (p, mu) is not identified, because shifting the cut is absorbed
#' exactly by a phase ramp on the coefficients, so the gauge is fixed by cutting at the sample mean
#' direction; compare two fits by the density they imply, not by their coefficients. Following the practice of the linear
#' family, an estimate that leaves the admissible set is reported rather than silently replaced; the
#' shrunken spectrum is returned as a labelled fallback because, unlike the linear case, the
#' factorisation itself fails without one.
#'
#' @param theta angles in radians.
#' @param nm number of circular L-moments, equal to the degree of the fitted family.
#' @param margin required minimum of the implied quantile density.
#' @param nodes quadrature nodes.
#' @return a list with the coefficients, mean direction, shrink factor and admissibility.
#' @export
fit_arccirc_fr <- function(theta, nm = 3L, margin = 0, nodes = 4096L) {
  o <- .C(C_arcc_fit_fr, theta = as.double(theta), n = as.integer(length(theta)),
          nm = as.integer(nm), margin = as.double(margin), nodes = as.integer(nodes),
          out = double(2 * (nm + 1) + 3))$out
  b <- 2 * (nm + 1)
  list(p = complex(real = o[seq(1, b, by = 2)], imaginary = o[seq(2, b, by = 2)]),
       mu = o[b + 1], shrink = o[b + 2], admissible = o[b + 3] == 1)
}

#' Arc-length tempering of a von Mises base
#'
#' The generator of the linear theory applied to a circular base, with the scale the construction
#' implicitly carries. Arc length adds a length to a density and is therefore not scale invariant;
#' on a circle of circumference \eqn{2\pi} a density is of order \eqn{1/(2\pi)}, so without the scale
#' the constant dominates and the transform returns the circular uniform. With it, the family
#' interpolates from the uniform at \eqn{s\to0} to the base itself as \eqn{s\to\infty}, monotonically
#' in concentration, preserving the mode throughout.
#'
#' @param theta angles in radians.
#' @param kappa von Mises concentration of the base.
#' @param mu mean direction.
#' @param s tempering scale.
#' @param nodes quadrature nodes for the normaliser.
#' @return density values at theta.
#' @export
dtemper_vm <- function(theta, kappa, mu = 0, s = 1, nodes = 4096L)
  .C(C_arcc_temper_vm, theta = as.double(theta), nth = as.integer(length(theta)),
     kappa = as.double(kappa), mu = as.double(mu), s = as.double(s),
     nodes = as.integer(nodes), out = double(length(theta)))$out

#' @describeIn dtemper_vm trigonometric moment of the tempered law
#' @param p integer order.
#' @export
trigmom_temper_vm <- function(p, kappa, mu = 0, s = 1, nodes = 4096L) {
  z <- .C(C_arcc_temper_vm_trigmom, p = as.integer(p), kappa = as.double(kappa),
          mu = as.double(mu), s = as.double(s), nodes = as.integer(nodes), out = double(2))$out
  complex(real = z[1], imaginary = z[2])
}

## ---- four-parameter kappa module (induction-curve paper) ----------------------------------------

#' Four-parameter kappa quantile, distribution and density functions
#' @param u,x numeric vectors of probabilities or quantiles
#' @param mu,sg,k,h kappa parameters (location, scale, two shapes)
#' @return numeric vector
#' @export
k4_q <- function(u, mu=0, sg=1, k, h)
  .C("arck4_q", as.double(u), as.integer(length(u)), as.double(mu), as.double(sg),
     as.double(k), as.double(h), out=double(length(u)))$out

#' @rdname k4_q
#' @export
k4_cdf <- function(x, mu=0, sg=1, k, h)
  .C("arck4_cdf", as.double(x), as.integer(length(x)), as.double(mu), as.double(sg),
     as.double(k), as.double(h), out=double(length(x)))$out

#' @rdname k4_q
#' @export
k4_pdf <- function(x, mu=0, sg=1, k, h)
  .C("arck4_pdf", as.double(x), as.integer(length(x)), as.double(mu), as.double(sg),
     as.double(k), as.double(h), out=double(length(x)))$out

#' Theoretical L-moment ratios of the standard kappa distribution
#' @param k,h shape parameters
#' @param nodes number of Gauss-Legendre nodes
#' @return c(tau3, tau4, l1, l2)
#' @export
k4_tau34 <- function(k, h, nodes=200L)
  .C("arck4_tau34", as.double(k), as.double(h), as.integer(nodes), out=double(4))$out

#' Running median with shrinking symmetric windows at the edges
#' @param y numeric vector
#' @param w odd window width
#' @return A numeric vector the same length as \code{y}.
#' @examples
#' ## the window shrinks symmetrically at the ends rather than padding, so the smoothed
#' ## series keeps the length of the original
#' y <- c(1, 8, 2, 9, 3, 10, 4)
#' length(k4_runmed(y, 3L)) == length(y)
#' @export
k4_runmed <- function(y, w=9L)
  .C("arck4_runmed", as.double(y), as.integer(length(y)), as.integer(w),
     out=double(length(y)))$out

## R's .C refuses a non-finite argument, so p = Inf (the tropical norm) is carried to C as -1.
#' @noRd
.p_code <- function(p) if (length(p) == 1L && is.infinite(p) && p > 0) -1 else as.double(p)

#' Band arc lengths of a scaled kappa curve and of a data polyline
#'
#' The arc-length element is a norm of the vector \eqn{(dx, dy)}. Which norm is a modelling
#' choice: \code{p = 1} gives the ordinary sum, \code{p = 2} the Euclidean arc length, and
#' \code{p = Inf} the max-plus (tropical) form \eqn{\max(|dx|, |dy|)}. The three agree to within
#' a factor of \eqn{\sqrt{2}}.
#' @param theta c(g0, g1, mu, sg, k, h)
#' @param breaks band break points (length J+1)
#' @param nodes number of Gauss-Legendre nodes per band
#' @param p the norm used for the arc-length element; 2 is the Euclidean default, Inf the tropical.
#' @return A numeric vector of band arc lengths, one per band.
#' @export
k4_band_model <- function(theta, breaks, nodes=60L, p=2)
  .C("arck4_band_model", as.double(theta), as.double(breaks),
     as.integer(length(breaks)-1L), as.integer(nodes), .p_code(p),
     out=double(length(breaks)-1L))$out

#' @rdname k4_band_model
#' @param x,y data ordered in x
#' @export
k4_band_sample <- function(x, y, breaks, p=2)
  .C("arck4_band_sample", as.double(x), as.double(y), as.integer(length(x)),
     as.double(breaks), as.integer(length(breaks)-1L), .p_code(p),
     out=double(length(breaks)-1L))$out

#' The two standard induction-period readings of a fitted kappa curve
#' @param theta c(g0, g1, mu, sg, k, h)
#' @param grid number of grid points for the dense evaluation
#' @return c(a = tangent reading, b = third-derivative reading, mode)
#' @export
k4_readings <- function(theta, grid=4000L){
  r <- .C("arck4_readings", as.double(theta), as.integer(grid), out=double(3))$out
  names(r) <- c("a","b","mode"); r }

#' Deterministic kappa fits: L-moment shape inversion, quantile-domain arc-length shape fit,
#' curve-domain NLS and NALR (banded arc lengths of the running-median presmoothed polyline)
#' @param t3,t4 target sample L-moment ratios
#' @param nodes number of Gauss-Legendre nodes for the theoretical ratios
#' @param y response vector (sorted internally where required)
#' @param bands two-column matrix of quantile bands
#' @param x data vector
#' @param J number of curve bands
#' @param lambda anchor weight
#' @param w running-median window
#' @param p the norm used for the arc-length element; 2 is the Euclidean default, Inf the tropical
#' @return A list of fitted parameters, whose components depend on which fitting routine
#'   is called; all return the six-vector \code{theta} of the scaled kappa curve.
#' @export
k4_fit_lmom <- function(t3, t4, nodes=200L){
  r <- .C("arck4_fit_lmom", as.double(t3), as.double(t4), as.integer(nodes), out=double(3))$out
  c(k=r[1], h=r[2], obj=r[3]) }

#' @rdname k4_fit_lmom
#' @export
k4_fit_aleq <- function(y, bands, start=c(0,0,log(0.3))){
  r <- .C("arck4_fit_aleq", as.double(sort(y)), as.integer(length(y)),
          as.double(t(bands)), as.integer(nrow(bands)), as.double(start), out=double(4))$out
  c(k=r[1], h=r[2], alpha=r[3], obj=r[4]) }

#' @rdname k4_fit_lmom
#' @param start transformed start vector
#' @param drift model the slow linear rise of the Rancimat water-trap conductivity,
#'   \code{y = g0 + m x + g1 F(x)}. Both induction-period readings are invariant to it, so the
#'   equivalence theory is unaffected and only the fit improves.
#' @export
k4_fit_nls <- function(x, y, start, drift = FALSE){
  if (drift) {
    if (length(start) == 6L) start <- c(start[1], 0, start[-1])
    r <- .C("arck4_fit_nls_drift", as.double(x), as.double(y), as.integer(length(x)),
            as.double(start), out=double(8))$out
    return(list(theta = r[c(1,3:7)], drift = r[2], obj = r[8]))
  }
  r <- .C("arck4_fit_nls", as.double(x), as.double(y), as.integer(length(x)),
          as.double(start), out=double(7))$out
  list(theta=r[1:6], obj=r[7]) }

#' @rdname k4_fit_lmom
#' @export
k4_fit_nalr <- function(x, y, start, J=12L, lambda=1, w=9L, p=2){
  r <- .C("arck4_fit_nalr", as.double(x), as.double(y), as.integer(length(x)),
          as.integer(J), as.double(lambda), as.integer(w), .p_code(p),
          as.double(start), out=double(7))$out
  list(theta=r[1:6], obj=r[7]) }

#' Quantile-domain induction readings and equivalence discrepancy for tilted beta-kernel
#' quantile densities q(u) = u^alpha (1-u)^beta exp(sum theta_j P_j(u))
#' @param alpha,beta kernel exponents
#' @param theta tilt coefficients on shifted Legendre polynomials (orders 1 to 3)
#' @param ngrid evaluation grid size
#' @return c(D, a, b, u_c, u_b); NA when the geometry is invalid
#' @export
eq_readings <- function(alpha, beta, theta=numeric(0), ngrid=20001L){
  r <- .C("arceq_readings", as.double(alpha), as.double(beta), as.integer(length(theta)),
          as.double(theta), as.integer(ngrid), out=double(5))$out
  names(r) <- c("D","a","b","uc","ub"); r }

#' @rdname eq_readings
#' @param lambda,delta van Staden-Loots parameters (kurtosis and skew weights)
#' @export
eq_readings_vsl <- function(lambda, delta, ngrid=20001L){
  r <- .C("arceq_readings_vsl", as.double(lambda), as.double(delta),
          as.integer(ngrid), out=double(5))$out
  names(r) <- c("D","a","b","uc","ub"); r }

## ---- closed-form machinery of the incomplete-beta equivalence family (C back-end) ----

#' Beta-companion quantile function (standardised support)
#'
#' The two-exponent quantile density \eqn{q(u) = u^{\alpha}(1-u)^{\beta}} integrated from the
#' origin, so \code{bc_q} is the quantile function of the beta-companion family on its
#' standardised support. Both exponents are negative on the family of interest.
#' @param u probabilities in \eqn{(0,1)}.
#' @param alpha,beta the two exponents of the quantile density.
#' @return A numeric vector of quantiles, the same length as \code{u}.
#' @examples
#' ## the quantile function is increasing, and the distribution function inverts it
#' u <- c(0.05, 0.25, 0.5, 0.75, 0.95)
#' x <- bc_q(u, alpha = -0.60, beta = -0.35)
#' all(diff(x) > 0)
#' max(abs(bc_cdf(x, -0.60, -0.35) - u))
#' @export
bc_q <- function(u, alpha, beta)
  .C("arceq2_bc_q", as.integer(length(u)), as.double(u), as.double(alpha),
     as.double(beta), out = double(length(u)))$out

#' Beta-companion distribution function
#'
#' The inverse of \code{\link{bc_q}}, obtained in the back-end by a continued-fraction
#' incomplete beta and a safeguarded root find.
#' @param x quantiles at which to evaluate.
#' @inheritParams bc_q
#' @return A numeric vector of probabilities, the same length as \code{x}.
#' @examples
#' bc_cdf(bc_q(c(0.3, 0.7), -0.60, -0.35), -0.60, -0.35)
#' @export
bc_cdf <- function(x, alpha, beta)
  .C("arceq2_bc_cdf", as.integer(length(x)), as.double(x), as.double(alpha),
     as.double(beta), out = double(length(x)))$out

#' Beta-companion density
#'
#' The derivative of \code{\link{bc_cdf}}, equal to the reciprocal of the quantile density
#' evaluated at the corresponding probability.
#' @param x points at which to evaluate the density.
#' @inheritParams bc_q
#' @return A numeric vector of density values, the same length as \code{x}.
#' @examples
#' ## the density is the derivative of the distribution function
#' xx <- bc_q(0.3, -0.60, -0.35); eps <- 1e-6
#' (bc_cdf(xx + eps, -0.60, -0.35) - bc_cdf(xx - eps, -0.60, -0.35)) / (2 * eps)
#' bc_pdf(xx, -0.60, -0.35)
#' @export
bc_pdf <- function(x, alpha, beta)
  .C("arceq2_bc_pdf", as.integer(length(x)), as.double(x), as.double(alpha),
     as.double(beta), out = double(length(x)))$out

#' Quadratic shoulder and mode of the kernel family
#'
#' The shoulder is the root of \eqn{3q'^2 = q q''}, transcendental on a general family. On this
#' one it reduces to a quadratic \eqn{Au^2 + Bu + C = 0}, whose constant term \eqn{C =
#' \alpha(2\alpha+1)} is positive exactly when \eqn{\alpha < -1/2}: that inequality is the
#' family's existence condition for a shoulder.
#' @inheritParams bc_q
#' @return A length-two vector, the shoulder \eqn{u_b} and the mode \eqn{u_c}.
#' @examples
#' ## the returned root satisfies the ORIGINAL transcendental equation, in the log form
#' ## 2g'^2 = g'', not merely the quadratic that replaced it
#' al <- -0.60; be <- -0.35
#' ub <- eq_ub_quad(al, be)[1]
#' gp <- al / ub - be / (1 - ub); gpp <- -al / ub^2 - be / (1 - ub)^2
#' abs(2 * gp^2 - gpp)
#' @export
eq_ub_quad <- function(alpha, beta)
  .C("arceq2_ub_quad", as.double(alpha), as.double(beta), out = double(2))$out

#' Closed-form equivalence discrepancy
#'
#' The signed amount by which the area identity \eqn{\int q = u_c q(u_c)} fails at
#' \eqn{(\alpha, \beta)}. Equivalence holds where it vanishes.
#' @inheritParams bc_q
#' @return A single number, the discrepancy.
#' @examples
#' ## the discrepancy vanishes on the equivalence curve and changes sign across it
#' bs <- eq_bstar(-0.60)
#' eq_E(-0.60, bs)
#' eq_E(-0.60, bs - 0.05) * eq_E(-0.60, bs + 0.05) < 0
#' @export
eq_E <- function(alpha, beta)
  .C("arceq2_E", as.double(alpha), as.double(beta), out = double(1))$out

#' Equivalence curve
#'
#' Solves \code{\link{eq_E}} for \eqn{\beta} at a given \eqn{\alpha} by a deterministic
#' scan followed by bisection, so the returned curve is reproducible rather than dependent on a
#' starting value.
#' @param alpha the first exponent of the quantile density. Documented here rather than
#'   inherited, because a combined two-name parameter tag cannot be inherited one name at a time.
#' @return A single number, \eqn{\beta^*(\alpha)}.
#' @examples
#' ## the solution curve is increasing in alpha, as the global theorem states
#' als <- seq(-0.62, -0.54, by = 0.02)
#' bss <- sapply(als, eq_bstar)
#' all(diff(bss) > 0)
#' @export
eq_bstar <- function(alpha)
  .C("arceq2_bstar", as.double(alpha), out = double(1))$out

#' The equivalence-system paper's fitting and simulation loops, in the back end
#'
#' \code{eqfit_bc} fits the beta-companion sigmoid by deterministic Nelder--Mead from a matrix
#' of starts (free: six columns; constrained to the equivalence curve: five), returning the
#' fitted six-parameter theta, the SSE and both readings. \code{eqfit_bc_boot} runs the iid
#' residual bootstrap of the constrained fit, warm-started, OpenMP over replicates with
#' per-replicate splitmix64 streams. \code{eqfit_k4eq} and \code{eqfit_k4eq_boot} are the
#' kappa-family analogues constrained to a locus table h -> k. \code{eqfit_msim} runs the
#' complete manifold-test simulation (truth and two displacement arms, with the optimiser
#' audit certificate). \code{eqfit_estsim} runs the three-estimator study. All outputs are
#' byte-identical across thread counts and between the R and Python fronts.
#' @param x,y the curve being fitted
#' @param curve two-column matrix (alpha, beta) of the traced equivalence curve
#' @param locus two-column matrix (h, k) of the kappa equivalence locus
#' @param th a fitted theta, used to warm-start the bootstrap
#' @param starts matrix of optimiser starts, one per row
#' @param B bootstrap or reference replicates
#' @param R,R2 replicate counts for the estimator simulation and for the null simulation
#' @param n,nsizes the sample size, and the vector of sample sizes swept over
#' @param alpha,betas a single alpha, and the vector of betas swept across it
#' @param st5 the constrained five-parameter start, held on the equivalence curve
#' @param awin two-element admissible window in alpha, passed as lower and upper bounds
#' @param thref three by four matrix of reference thetas, one row per reference member
#' @param seed seed for the per-replicate splitmix64 streams, so a run reproduces exactly
#' @param maxit Nelder--Mead iteration cap per pass
#' @export
eqfit_bc <- function(x, y, starts, curve, maxit = 4000L){
  starts <- as.matrix(starts); np <- ncol(starts)
  stopifnot(np %in% c(5L, 6L))
  r <- .C("arceqfit_bc", as.double(x), as.double(y), as.integer(length(x)),
          as.double(t(starts)), as.integer(nrow(starts)), as.integer(np), as.integer(maxit),
          as.double(curve[,1]), as.double(curve[,2]), as.integer(nrow(curve)),
          as.double(min(curve[,1])), as.double(max(curve[,1])), out = double(9))$out
  list(th = r[1:6], sse = r[7], a = r[8], b = r[9])
}
#' @rdname eqfit_bc
#' @export
eqfit_bc_boot <- function(x, y, th, curve, B = 30L, seed = 4207L, maxit = 2000L){
  .C("arceqfit_bc_boot", as.double(x), as.double(y), as.integer(length(x)), as.double(th),
     as.integer(B), as.integer(seed), as.integer(maxit),
     as.double(curve[,1]), as.double(curve[,2]), as.integer(nrow(curve)),
     as.double(min(curve[,1])), as.double(max(curve[,1])), out = double(B))$out
}
#' @rdname eqfit_bc
#' @export
eqfit_k4eq <- function(x, y, starts, locus, maxit = 4000L){
  starts <- as.matrix(starts); stopifnot(ncol(starts) == 5L)
  r <- .C("arceqfit_k4eq", as.double(x), as.double(y), as.integer(length(x)),
          as.double(t(starts)), as.integer(nrow(starts)), as.integer(maxit),
          as.double(locus[,1]), as.double(locus[,2]), as.integer(nrow(locus)),
          as.double(min(locus[,1])), as.double(max(locus[,1])), out = double(8))$out
  list(th = r[1:6], sse = r[7], a = r[8])
}
#' @rdname eqfit_bc
#' @param p0 warm-start parameter vector for the bootstrap refits
#' @export
eqfit_k4eq_boot <- function(x, y, th, p0, locus, B = 25L, seed = 4207L, maxit = 2000L){
  .C("arceqfit_k4eq_boot", as.double(x), as.double(y), as.integer(length(x)), as.double(th),
     as.double(p0), as.integer(B), as.integer(seed), as.integer(maxit),
     as.double(locus[,1]), as.double(locus[,2]), as.integer(nrow(locus)),
     as.double(min(locus[,1])), as.double(max(locus[,1])), out = double(B))$out
}
#' @rdname eqfit_bc
#' @param nM,sde,Rm design size, noise standard deviation and replicates of the manifold study
#' @param th0 the true manifold member; \code{st5} the three profiling start transforms
#'   (mu, sigma, alpha-centre); \code{thref} a 3 x 4 matrix of reference manifold members
#'   (mu, sigma, alpha, beta) for the three displaced arms, used for the power certificate;
#'   \code{awin} the
#'   profiling alpha window \code{c(alo, ahi)}
#' @export
eqfit_msim <- function(nM, sde, Rm, th0, st5, curve, awin, thref, seed = 4207L, maxit = 2000L){
  stopifnot(is.matrix(thref), nrow(thref) == 3L, ncol(thref) == 4L)
  r <- .C("arceqfit_msim", as.integer(nM), as.double(sde), as.integer(Rm), as.integer(seed),
          as.double(th0), as.double(st5), as.integer(maxit),
          as.double(curve[,1]), as.double(curve[,2]), as.integer(nrow(curve)),
          as.double(awin[1]), as.double(awin[2]), as.double(t(thref)),
          out = double(8 * Rm))$out
  list(Tlev = r[1:Rm], Taud = r[Rm + 1:Rm], Tp2 = r[2*Rm + 1:Rm],
       Tp15 = r[3*Rm + 1:Rm], Tp4 = r[4*Rm + 1:Rm],
       Tref2 = r[5*Rm + 1:Rm], Tref15 = r[6*Rm + 1:Rm], Tref4 = r[7*Rm + 1:Rm])
}
#' @rdname eqfit_bc
#' @param mu,sg,al,be a single parameter point for \code{eqfit_score_at}
#' @export
eqfit_score_at <- function(x, y, mu, sg, al, be){
  .C("arceqfit_score_at", as.double(x), as.double(y), as.integer(length(x)),
     as.double(mu), as.double(sg), as.double(al), as.double(be), out = double(1))$out
}
#' @rdname eqfit_bc
#' @param res a residual vector
#' @export
eqfit_blocklen <- function(res){
  .C("arceqfit_blocklen", as.double(res), as.integer(length(res)), out = integer(1))$out
}
#' @rdname eqfit_bc
#' @export
eqfit_nullT <- function(nM, sde, R2, th0, st5, curve, awin, seed = 4207L, maxit = 2000L){
  .C("arceqfit_nullT", as.integer(nM), as.double(sde), as.integer(R2), as.integer(seed),
     as.double(th0), as.double(st5), as.integer(maxit),
     as.double(curve[,1]), as.double(curve[,2]), as.integer(nrow(curve)),
     as.double(awin[1]), as.double(awin[2]), out = double(R2))$out
}
#' @rdname eqfit_bc
#' @param a0,b0 beta parameters of the sampled member
#' @export
eqfit_estsim <- function(a0, b0, R, nsizes, seed = 4207L, maxit = 500L){
  nn <- length(nsizes)
  r <- .C("arceqfit_estsim", as.double(a0), as.double(b0), as.integer(R),
          as.integer(nsizes), as.integer(nn), as.integer(seed), as.integer(maxit),
          out = double(R * 3 * 2 * nn))$out
  array(r, dim = c(R, 3, 2, nn))
}
#' @rdname eqfit_bc
#' @param k,h shape grid for the wall sweep; \code{betas} the beta grid of an E section
#' @export
k4_b_sweep <- function(k, h){
  .C("arceqfit_bsweep", as.double(k), as.integer(length(k)), as.double(h),
     out = double(length(k)))$out
}
#' @rdname eqfit_bc
#' @export
eqfit_taus <- function(a0, b0, R, n, seed = 4207L){
  r <- .C("arceqfit_taus", as.double(a0), as.double(b0), as.integer(R), as.integer(n),
          as.integer(seed), out = double(4 * R))$out
  cbind(t3 = r[1:R], t4 = r[R + 1:R], al = r[2*R + 1:R], be = r[3*R + 1:R])
}
#' @rdname eqfit_bc
#' @export
k4_ab_sweep <- function(k, h){
  r <- .C("arceqfit_absweep", as.double(k), as.integer(length(k)), as.double(h),
          a = double(length(k)), b = double(length(k)))
  list(a = r$a, b = r$b)
}
#' @rdname eqfit_bc
#' @export
eq_E_sweep <- function(alpha, betas){
  .C("arceqfit_Esweep", as.double(alpha), as.double(betas), as.integer(length(betas)),
     out = double(length(betas)))$out
}

#' The fitting-method multiverse of one induction run: eight admissible pipelines under one
#' moving-block residual bootstrap
#'
#' All eight pipelines report the standard tangent reading from the same trace: the
#' variable-projection least-squares fit with its GEV submodel selection (\code{LS}), running
#' medians of window nine and twenty-one followed by transformed-parameter NLS
#' (\code{med9-LS}, \code{med21-LS}), the banded arc-length estimator (\code{arc}), the GEV
#' submodel in its own right (\code{GEV-LS}), the selection fit past the run's measured
#' transient cutoff (\code{transient-excised}), and the selection fit on each half-density
#' index grid (\code{grid-odd}, \code{grid-even}). Replicates start from the base fit;
#' point estimates use the full start grid. Resampling indices come from per-replicate
#' splitmix64 streams, so the result is identical whatever the OpenMP thread count and
#' identical between the R and Python fronts.
#'
#' @param x,y the logged curve
#' @param cut transient cutoff in the units of \code{x}; \code{-Inf} excises nothing
#' @param B bootstrap replicates; \code{0} returns point estimates only
#' @param seed integer seed for the per-replicate streams
#' @param bounds the admissible shape box \code{(k_lo, k_hi, h_lo, h_hi)}
#' @param maxit Nelder--Mead iteration cap per start
#' @return \code{list(a_pt, A)}: the eight named point estimates, and the \code{B} by eight
#'   matrix of replicate readings (\code{NULL} when \code{B = 0})
#' @export
k4_mv_boot <- function(x, y, cut = -Inf, B = 50L, seed = 4207L,
                       bounds = c(-0.98, 0.95, 0, 4), maxit = 1500L){
  n <- length(x); stopifnot(length(y) == n, B >= 0)
  r <- .C("arck4_mv_boot", as.double(x), as.double(y), as.integer(n),
          as.double(if (is.finite(cut)) cut else -1e300), as.integer(B),
          as.integer(seed), as.double(bounds), as.integer(maxit),
          a_pt = double(8), A = double(max(1L, B) * 8))
  meth <- c("LS", "med9-LS", "med21-LS", "arc", "GEV-LS", "transient-excised",
            "grid-odd", "grid-even")
  a_pt <- r$a_pt; names(a_pt) <- meth
  A <- if (B > 0) matrix(r$A, nrow = B, ncol = 8, byrow = TRUE,
                         dimnames = list(NULL, meth)) else NULL
  list(a_pt = a_pt, A = A)
}

#' Variable-projection fit of the drifted kappa response
#'
#' Fits \eqn{y = g_0 + m x + g_1 F(x;\mu,\sigma,k,h)}. The mean curve is linear in the three
#' coefficients, so they are solved exactly for any shape and only the four shape parameters are
#' searched; the multi-start is run in the C back end, in parallel over starts, and the best is
#' chosen serially so the answer does not depend on the thread count.
#' @param x,y the curve.
#' @param starts a matrix with four columns, \code{(mu, log sigma, k, log h)}, one row per start.
#' @param maxit iterations per Nelder-Mead descent.
#' @param bounds admissible shape box \code{c(k_lo, k_hi, h_lo, h_hi)}. The default lower bound on
#'   \eqn{k} is the one the edible-oil analysis settled on: a wider box lowers the residual sum of
#'   squares but moves the fitted induction periods away from the laboratory values, so the choice is
#'   the analyst's and is stated here rather than fixed in the back end. A lower bound of zero
#'   on \eqn{h} is admissible and meaningful: it admits the generalised extreme value
#'   distribution, the \eqn{h \to 0} member of the family, and with \eqn{k = 0} the Gumbel.
#' @param hfix hold \eqn{h} at this value and fit the remaining three parameters, which is how a
#'   submodel is fitted in its own right; \code{NA} (the default) fits all four. Because \eqn{h} is
#'   searched on the log scale, \eqn{h = 0} is unreachable by the free search, which creeps toward
#'   it and reports a spurious small value instead of naming the submodel. Choosing between the two
#'   fits is model selection between a nested pair, and is left to the caller.
#' @return a list with \code{theta} in the six-vector convention \code{(g0, g1, mu, sg, k, h)},
#'   the fitted \code{drift}, and the residual sum of squares \code{rss}.
#' @export
k4_fit_varpro <- function(x, y, starts, maxit = 1500L,
                          bounds = c(-0.98, 0.95, 0, 4), hfix = NA_real_) {
  starts <- as.matrix(starts)
  stopifnot(ncol(starts) == 4L, length(bounds) == 4L,
            bounds[1] < bounds[2], bounds[3] < bounds[4], bounds[3] >= 0)
  r <- .C("arck4_fit_varpro", as.double(x), as.double(y), as.integer(length(x)),
          as.double(t(starts)), as.integer(nrow(starts)), as.integer(maxit),
          as.double(bounds), as.double(if (is.na(hfix)) -1 else hfix),
          out = double(8))$out
  list(theta = c(r[1], r[3], r[4], r[5], r[6], r[7]), drift = r[2], rss = r[8])
}

#' Exact mean of the tropical band arc length under Gaussian error
#'
#' The tropical (max-plus) arc-length element \eqn{\max(dx, |dy|)} has an expectation elementary in
#' \eqn{\Phi} and \eqn{\varphi} when the curve is observed with independent Gaussian error, whereas
#' the Euclidean element's expectation is a confluent hypergeometric function. The observed band arc
#' length can therefore be compared with its own mean rather than with a clean-curve quantity it does
#' not estimate.
#'
#' Only the mean is returned. The variance of a band sum is not the sum of the elements' variances:
#' consecutive increments share an observation, and although the noise increments correlate at
#' \eqn{-1/2} the absolute value discards the sign and leaves the elements positively correlated, so
#' summing as if independent understates a band's variance by roughly a fifth. That correction is not
#' elementary and is deliberately not supplied.
#' @param x ordered design points.
#' @param theta the six-vector \code{(g0, g1, mu, sg, k, h)}.
#' @param sigma error standard deviation of a single observation.
#' @param breaks the \code{J+1} band edges.
#' @return the \code{J} expected band arc lengths.
#' @export
k4_band_trop_mean <- function(x, theta, sigma, breaks) {
  J <- length(breaks) - 1L
  stopifnot(length(theta) == 6L, sigma >= 0, J >= 1L)
  .C("arck4_band_trop_mean", as.double(x), as.integer(length(x)), as.double(theta),
     as.double(sigma), as.double(breaks), as.integer(J), out = double(J))$out
}

#' One-way variance components and the intraclass correlation
#'
#' Decomposes the variance of \code{y} into within-group and between-group parts under a one-way
#' random-effects model, and returns the intraclass correlation. Groups may be of unequal size;
#' the between-group mean square is divided by the unbalanced constant
#' \eqn{k_0 = (N - \sum n_i^2 / N)/(G-1)} rather than by the mean group size, which is the
#' balanced-design shortcut and biases the ratio when group sizes differ.
#'
#' A negative between-group variance estimate is truncated at zero, the usual convention.
#'
#' @param y numeric observations.
#' @param g grouping vector; coerced with \code{factor}, so any labels will do.
#' @return named numeric vector: \code{icc}, \code{sd_within}, \code{sd_between}.
#' @examples
#' ## replicate runs on the same specimen agree more closely than runs on different specimens
#' y <- c(1.1, 1.2, 1.0,  5.1, 5.3, 4.9,  9.0, 9.2, 8.8)
#' g <- rep(1:3, each = 3)
#' k4_icc(y, g)["icc"] > 0.9
#' @export
k4_icc <- function(y, g) {
  ok <- is.finite(y)
  y  <- as.double(y[ok]); g <- factor(g[ok])
  out <- .C("arcvc_icc_oneway", y, as.integer(as.integer(g) - 1L),
            as.integer(length(y)), as.integer(nlevels(g)), out = double(3))$out
  stats::setNames(out, c("icc", "sd_within", "sd_between"))
}

#' Distance from fitted shapes to the equivalence locus
#'
#' Shortest Euclidean distance in the \eqn{(h,k)} shape plane from each fitted shape to the
#' equivalence locus, supplied as a polyline. Distance is measured to the segments of the polyline
#' rather than to its vertices: a vertex-only search overstates the distance by up to half the
#' vertex spacing, which matters when the distance is compared against the price of the constraint.
#'
#' @param h,k numeric vectors of fitted shape parameters, of equal length.
#' @param locus_h,locus_k the locus polyline, of equal length.
#' @return numeric vector of distances, \code{NA} where the shape is not finite.
#' @examples
#' ## a shape sitting on the locus is at distance zero
#' lh <- seq(0.05, 0.39, length.out = 20); lk <- seq(-0.04, 0.81, length.out = 20)
#' k4_locus_dist(lh[5], lk[5], lh, lk) < 1e-12
#' @export
k4_locus_dist <- function(h, k, locus_h, locus_k) {
  stopifnot(length(h) == length(k), length(locus_h) == length(locus_k))
  .C("arcvc_locus_dist", as.double(h), as.double(k), as.integer(length(h)),
     as.double(locus_h), as.double(locus_k), as.integer(length(locus_h)),
     out = double(length(h)))$out
}

#' Exact confidence interval for the shape parameter of the wrapped order-three family
#'
#' Inverts a Monte Carlo test on a rotation-invariant statistic. The rotation is eliminated exactly
#' by invariance, so no nuisance parameter is profiled; the level is exact for every sample size and
#' every number of reference draws rather than asymptotically, which matters because the admissible
#' boundary is attained and because the likelihood of this family is unbounded. Both statistics are
#' reflection invariants, so the interval brackets \code{|c3|}; the sign requires a reflection-odd
#' statistic. Recorded resolution is handled by grouping the reference draws as the data are grouped.
#'
#' @param theta observed angles in radians
#' @param cgrid grid over the admissible \code{|c3|} interval; the grid is the parameter space, so
#'   it is exhaustive rather than a search
#' @param B reference draws per grid point
#' @param stat 0 for the first trigonometric moment (about twice as efficient) or 1 for the
#'   arc-length spacings functional. Both are reflection invariants.
#' @param seed seed for the splitmix64 stream, so a run reproduces exactly
#' @param group 0 for continuous data, or the number of equal cells the data are grouped into
#'   (36 for ten-degree grouping)
#' @param level test level, so the interval has coverage \code{1 - level}
#' @return list with the interval, the observed statistic and the p-value curve that was inverted
#' @export
arcc_exact_ci <- function(theta, cgrid = seq(0, arcc_c3max(), length.out = 81), B = 999L,
                          stat = 0L, group = 0L, level = 0.10, seed = 4207L){
  ng <- length(cgrid)
  o <- .C(C_arcc_exact_ci, as.double(theta), as.integer(length(theta)), as.integer(B),
          as.integer(seed), as.double(cgrid), as.integer(ng), as.integer(stat),
          as.integer(group), as.double(level), out = double(3L + ng))$out
  list(lower = o[1], upper = o[2], stat = o[3], pcurve = o[3 + seq_len(ng)], cgrid = cgrid,
       excludes_uniform = !is.na(o[1]) && o[1] > 0)
}

#' Goodness of fit of a fitted circular arc-length member
#'
#' The probability-integral transform through the fitted member sends the sample to uniform under the
#' hypothesis, so the reference law of the arc-length spacings statistic is distribution-free and the
#' Monte Carlo p-value is exact in level for any number of reference draws.
#'
#' @param theta angles in radians
#' @param fit a fitted object from \code{fit_arccirc}, supplying the member through which the
#'   probability-integral transform is taken
#' @param B reference draws
#' @param seed seed for the splitmix64 stream, so a run reproduces exactly
#' @export
arcc_gof <- function(theta, fit, B = 999L, seed = 4207L){
  o <- .C(C_arcc_gof, as.double(theta), as.integer(length(theta)), as.double(fit$c3),
          as.double(fit$mu), as.integer(B), as.integer(seed), out = double(2))$out
  list(stat = o[1], p = o[2])
}
