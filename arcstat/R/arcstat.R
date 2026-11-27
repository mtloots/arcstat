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
#' @examples
#' u <- sort(runif(100))
#' al_statistic(u)
#' @export
al_statistic <- function(u) {
  if (any(u < 0 | u > 1)) stop("`u` must be probability-integral transforms in [0,1]")
  .C(C_al_statistic, n = as.integer(length(u)), u = as.double(u), out = double(1))$out
}

#' Saddlepoint right-tail probability of the arc-length statistic (C back-end)
#' @param s observed statistic value;
#' @param n sample size.
#' @return the right-tail p-value.
#' @examples
#' u <- sort(stats::runif(100))
#' al_pvalue(al_statistic(u), 100)
#' @export
al_pvalue <- function(s, n)
  .C(C_al_pvalue, s = as.double(s), n = as.integer(n), out = double(1))$out

#' Exact mean, variance and support of the arc-length statistic (C back-end)
#' @param n sample size.
#' @return list with mean, var, support.
#' @examples
#' m <- al_moments(100L)
#' m$mean; m$support
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
#' @examples
#' set.seed(1)
#' al_test(runif(150))$p.value
#' al_test(rbeta(150, 0.5, 0.5))$p.value
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

#' Theoretical L-moments of an arcq distribution
#'
#' L_r is the shifted-Legendre projection of Q (Hosking's convention): L_r =
#' integral_0^1 Q(u) P*_(r-1)(u) du, r>=1, with P*_r the shifted Legendre polynomial (P*_r(1)=1).
#' @param obj an "arcq" object.
#' @param nmom number of L-moments.
#' @return the first nmom L-moments; ratios tau_r = L_r/L_2 for r>=3 via lmratios().
#' @examples
#' lmoments(arcq(0.5))
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
#' @examples
#' set.seed(1)
#' fit_arcq_cf(rnorm(400))
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
#' @examples
#' g <- arc_generator(function(z) dnorm(z), -6, 6)
#' pgen(0, g)
#' dgen(c(-1, 0, 1), g)
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
#' @examples
#' p <- bb_post_disc((1:60 - 0.5) / 60, M = 200L)
#' summary(p)
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
#' @return a matrix with columns \code{arc} and \code{ks}, and the sample size as attribute \code{n}.
#' @examples
#' ref <- bb_ref_disc(60L)
#' str(ref)
#' @export
bb_ref_disc <- function(n, D = 300L, m = 16L, seed = 7L) {
  r <- .C(C_bb_ref, n = as.integer(n), D = as.integer(D), m = as.integer(m),
          seed = as.integer(seed), out = double(2L * D * m))$out
  structure(matrix(r, ncol = 2, byrow = TRUE, dimnames = list(NULL, c("arc", "ks"))), n = as.integer(n))
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
#' @param ref optional precomputed reference from \code{bb_ref_disc}, for the sample size of \code{u}
#'   (a reference for another size is refused); computed if \code{NULL}.
#' @param seed integer seed for the C back-end.
#' @return the posterior evidence against the null, in [0,1].
#' @examples
#' ref <- bb_ref_disc(60L)
#' bb_evidence((1:60 - 0.5) / 60, ref = ref, M = 300L)
#' @export
bb_evidence <- function(u, which = "arc", level = 0.95, M = 1500L, ref = NULL, seed = 1L) {
  n <- length(u); if (is.null(ref)) ref <- bb_ref_disc(n)
  if (!is.null(attr(ref, "n")) && attr(ref, "n") != n)
    stop("`ref` was computed for samples of size ", attr(ref, "n"), ", not ", n)
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
#' @examples
#' set.seed(1)
#' ## the normal carries arc length 2 (see cf_arclength_family)
#' cf_arclength(rnorm(5000))
#' @export
cf_arclength <- function(x, T = 6, ngrid = 1200L) {
  x <- as.double(x[is.finite(x)])
  .C(C_cf_arclength_emp, x = x, n = length(x), tmax = as.double(T), # nolint: T_and_F_symbol_linter.
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
#' @examples
#' ## the normal reference band arc length grows with the scale
#' al_band_model(c(0.5, 1, 2))
#' @export
al_band_model <- function(sigma, a = 0.05, b = 0.95, nodes = 400L)
  .C(C_al_band_model, sigma = as.double(sigma), a = as.double(a), b = as.double(b),
     nodes = as.integer(nodes), out = double(1))$out

#' Sample band arc length of the empirical distribution curve (C back-end)
#'
#' @param x data;
#' @param a,b the probability band.
#' @return the sample band arc length, NA when fewer than two points fall in the band.
#' @examples
#' set.seed(1)
#' al_band_sample(rnorm(400, sd = 2))
#' @export
al_band_sample <- function(x, a = 0.05, b = 0.95) {
  x <- as.double(x[is.finite(x)])
  .C(C_al_band_sample, x = x, n = length(x), a = as.double(a), b = as.double(b), out = double(1))$out
}

#' Spacing-corrected model band arc length (C back-end)
#'
#' The sample band arc length does not converge to the population band arc length
#' \code{\link{al_band_model}}. Consecutive order statistics inside the band are spaced like
#' \eqn{Q'(u)E/n} with \eqn{E} standard exponential, so the polygonal sum converges instead to
#' \deqn{S^*(\sigma) = \int_{z_a}^{z_b} E\sqrt{f_0(z)^2 + \sigma^2 E^2}\,dz = E\,S(\sigma E),}
#' which differs from \eqn{S(\sigma)} by a constant factor because the square root is nonlinear and
#' \eqn{E[E^2] = 2}. The gap does not close as the sample grows. This is the functional
#' \code{\link{al_scale}} matches, and matching to it is what makes that estimator consistent.
#' @param sigma scale at which to evaluate the model band length.
#' @param a,b the probability band.
#' @param nodes number of Simpson nodes across the band; the exponential expectation uses a
#'   64-point Gauss-Laguerre rule, which is exact to machine precision for this integrand.
#' @return the spacing-corrected model band arc length at \code{sigma}.
#' @examples
#' al_band_model_star(1, 0.1, 0.9)
#' @export
al_band_model_star <- function(sigma, a = 0.05, b = 0.95, nodes = 400) {
  .C(C_al_band_model_star, sigma = as.double(sigma), a = as.double(a), b = as.double(b),
     nodes = as.integer(nodes), out = double(1))$out
}

#' Scale by arc-length band matching, in the scale-equivariant standardised form (C back-end)
#'
#' The raw matching equation fixes an aspect ratio between the probability and response axes and is
#' therefore unit dependent; this estimator divides by the MAD, matches on the standardised scale and
#' rescales, which is exactly scale equivariant.
#'
#' The sample band length is matched to \code{\link{al_band_model_star}}, not to
#' \code{\link{al_band_model}}. Matching to the population band length is inconsistent: the sample
#' quantity converges to \eqn{S^*} rather than to \eqn{S}, so that estimator is biased by a constant
#' factor that does not vanish with the sample size, about five per cent high at the standard normal
#' on the default band. \code{\link{al_scale_raw}} retains the uncorrected form for comparison.
#' @param x data;
#' @param a,b the probability band, whose tail mass sets the breakdown point.
#' @return the estimated scale, NA when the MAD vanishes or the matching equation has no root.
#' @examples
#' set.seed(1)
#' x <- rnorm(400, sd = 3)
#' al_scale(x)
#' ## scale equivariance
#' al_scale(7 * x) / al_scale(x)
#' @export
al_scale <- function(x, a = 0.05, b = 0.95) {
  x <- as.double(x[is.finite(x)])
  .C(C_al_scale, x = x, n = length(x), a = as.double(a), b = as.double(b), out = double(1))$out
}

#' The uncorrected band-matching scale estimator (C back-end)
#'
#' Matches the sample band arc length to the population band length \code{\link{al_band_model}}
#' rather than to \code{\link{al_band_model_star}}. INCONSISTENT: it estimates the wrong quantity,
#' by a constant factor that does not vanish as the sample grows. Retained only so that the size of
#' the correction can be measured. Use \code{\link{al_scale}}.
#' @param x data;
#' @param a,b the probability band.
#' @return the uncorrected scale estimate, NA when the MAD vanishes or there is no root.
#' @examples
#' set.seed(1)
#' al_scale_raw(rnorm(400, sd = 3))
#' @export
al_scale_raw <- function(x, a = 0.05, b = 0.95) {
  x <- as.double(x[is.finite(x)])
  .C(C_al_scale_raw, x = x, n = length(x), a = as.double(a), b = as.double(b), out = double(1))$out
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
#' @examples
#' arcc_c3max()
#' @export
arcc_c3max <- function() .C(C_arcc_c3max, out = double(1))$out

#' The smooth circular arc-length family
#'
#' @param c3 shape parameter, |c3| <= arcc_c3max(); c3 = 0 is the circular uniform.
#' @param mu mean direction in radians.
#' @return an object of class "arccirc".
#' @examples
#' o <- arccirc(0.4, mu = 0.3)
#' ## the density integrates to one
#' stats::integrate(function(t) darccirc(t, o), 0, 2 * pi)$value
#' rho_arccirc(o)
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
#' @examples
#' th <- rarccirc(300, arccirc(0.5, 0.2))
#' fit_arccirc(th)[c("c3", "mu")]
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
#' @examples
#' fr <- arccirc_fr(c(0.3 + 0.1i))
#' darccirc_fr(c(0, 1, 2), fr)
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
#' @examples
#' l <- lmom_circ(rarccirc(500L, arccirc(0.8, 0.3)), 3L)
#' l
#' qmin_rho(0.4)
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
#'
#' A factorisation exists only where the implied quantile density is non-negative, which is what
#' \code{qmin_rho} measures; for a single moment with \eqn{\rho_0 = 1} the condition reduces to
#' \eqn{|\rho| \le 1/2}. Asked for an inadmissible spectrum this routine used to run the root
#' finder anyway and return the coefficients of a factorisation that does not exist, silently and
#' at any modulus. It now does what the fitting path in the C back-end has always done: shrink
#' towards the circular uniform until the spectrum is admissible, and REPORT that it did, so a
#' caller can tell an answer to the question asked from an answer to a nearby one.
#' @return a list with the coefficient vector \code{p}, the \code{shrink} factor applied to the
#'   spectrum (one when nothing was needed), and the \code{admissible} flag for the spectrum AS
#'   SUPPLIED.
#' @export
factorise_rho <- function(rho, margin = 0, nodes = 4096L) {
  adm    <- qmin_rho(rho, nodes) >= 0
  shrink <- 1
  if (!adm) {
    a <- admiss_rho(rho, margin, nodes)
    rho <- a$rho; shrink <- a$shrink
    warning("spectrum is not admissible: no Fejer-Riesz factorisation exists for it, so it was ",
            "shrunk towards the circular uniform by a factor of ", format(shrink, digits = 6),
            ". The returned coefficients describe the SHRUNKEN spectrum, not the one supplied.",
            call. = FALSE)
  }
  nm <- length(rho)
  o <- .C(C_arcc_factorise, rhore = as.double(Re(rho)), rhoim = as.double(Im(rho)),
          nm = as.integer(nm), out = double(2 * (nm + 1)))$out
  list(p = complex(real = o[seq(1, 2 * (nm + 1), by = 2)],
                   imaginary = o[seq(2, 2 * (nm + 1), by = 2)]),
       shrink = shrink, admissible = adm)
}

#' Closed-form fit of the Fejer-Riesz circular arc-length family
#'
#' Sort, apply the circular L-moment weights, invert the exact linear map, factorise: no numerical
#' search enters anywhere. The pair (p, mu) is not identified, because shifting the cut is absorbed
#' exactly by a phase ramp on the coefficients, so the gauge is fixed by cutting at the sample mean
#' direction; compare two fits by the density they imply, not by their coefficients. Following the practice of the
#'   linear
#' family, an estimate that leaves the admissible set is reported rather than silently replaced; the
#' shrunken spectrum is returned as a labelled fallback because, unlike the linear case, the
#' factorisation itself fails without one.
#'
#' @param theta angles in radians.
#' @param nm number of circular L-moments, equal to the degree of the fitted family.
#' @param margin required minimum of the implied quantile density.
#' @param nodes quadrature nodes.
#' @return a list with the coefficients, mean direction, shrink factor and admissibility.
#' @examples
#' th <- rarccirc(300, arccirc(0.5, 0.2))
#' str(fit_arccirc_fr(th))
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
#' @examples
#' stats::integrate(function(t) dtemper_vm(t, kappa = 1), 0, 2 * pi)$value
#' trigmom_temper_vm(1, kappa = 1)
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
#' @examples
#' u <- c(0.1, 0.5, 0.9)
#' x <- k4_q(u, 8, 2, -0.3, 0.2)
#' k4_cdf(x, 8, 2, -0.3, 0.2)
#' k4_pdf(x, 8, 2, -0.3, 0.2)
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
#' @param nodes number of quadrature nodes, at least 200. The rule is the trapezoid after a map that
#'   removes both end singularities, and 200 nodes give full double precision.
#' @return c(tau3, tau4, l1, l2); NaN where the moments do not exist (\eqn{k \le -1} or
#'   \eqn{hk \le -1} with \eqn{h < 0})
#' @examples
#' k4_tau34(0.1, 0.4)
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
#' @examples
#' th <- c(0, 300, 8, 2, -0.3, 0.2)
#' br <- seq(0, 20, by = 5)
#' k4_band_model(th, br)
#' x <- seq(0, 20, length.out = 200)
#' k4_band_sample(x, th[1] + th[2] * k4_cdf(x, 8, 2, -0.3, 0.2), br)
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
#'
#' Both readings are closed forms. The steepest point \code{mode} exists exactly when
#' \eqn{k < 1}, \eqn{h < 1} and \eqn{hk < 1}; the third-derivative reading \code{b} exists when in
#' addition \eqn{h < 1/2} and \eqn{hk < 1/2}. Where a reading does not exist it is \code{NaN}. Versions
#' up to 0.3.1 searched a grid between the 0.1 and 99.9 per cent quantiles; with a heavy right tail
#' (\eqn{k} near \eqn{-1}) the grid missed the steepest point and the tangent reading came out too low.
#' @param theta c(g0, g1, mu, sg, k, h)
#' @param grid ignored, kept so that existing calls work; the result does not depend on a grid
#' @return c(a = tangent reading, b = third-derivative reading, mode)
#' @examples
#' ## the tangent at the steepest point meets the baseline at a
#' th <- c(1, 400, 8, 2.5, -0.9, 0.2)
#' r <- unname(k4_readings(th)); cdf <- k4_cdf(r[3], 8, 2.5, -0.9, 0.2)
#' r[3] - cdf / k4_pdf(r[3], 8, 2.5, -0.9, 0.2) - r[1]
#' @export
k4_readings <- function(theta, grid=4000L){
  r <- .C("arck4_readings", as.double(theta), as.integer(grid), out=double(3))$out
  names(r) <- c("a","b","mode"); r }

#' Deterministic kappa fits: L-moment shape inversion, quantile-domain arc-length shape fit,
#' curve-domain NLS and NALR (banded arc lengths of the running-median presmoothed polyline)
#' @param t3,t4 target sample L-moment ratios
#' @param nodes number of quadrature nodes for the theoretical ratios (at least 200 are used)
#' @param y response vector (sorted internally where required)
#' @param bands two-column matrix of quantile bands. A band that ends at one has infinite model
#'   length unless \eqn{k > 0}, so such a band excludes \eqn{k \le 0}.
#' @param x data vector
#' @param J number of curve bands
#' @param lambda anchor weight
#' @param w running-median window
#' @param p the norm used for the arc-length element; 2 is the Euclidean default, Inf the tropical
#' @return A list of fitted parameters, whose components depend on which fitting routine
#'   is called; all return the six-vector \code{theta} of the scaled kappa curve.
#' @examples
#' ## invert the L-moment ratios of a member
#' tt <- k4_tau34(0.1, 0.4)
#' k4_fit_lmom(tt[1], tt[2])
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
#' @details The quantile function \eqn{Q(u)} is the integral of \eqn{q} from zero, so \code{a}
#'   (tangent reading) and \code{b} (derivative reading) are quantile values on the same footing
#'   as \code{\link{bc_q}}. Without a tilt (and for \eqn{\alpha > -1}) everything is in closed form
#'   and \code{ngrid} has no effect. With a tilt the readings are located on a grid from
#'   \eqn{10^{-5}} to \eqn{1 - 10^{-5}} and refined; for \eqn{\alpha > -1} the mass below the first
#'   node is added and the singular power factor is integrated exactly on every cell, with an error
#'   of about \eqn{10^{-8}}; a derivative reading below \eqn{u = 3\times 10^{-4}} is not found there.
#'   Versions up to 0.3.1 started the integral at the first node with a trapezoid, which lowered
#'   \code{a} and \code{b} by up to 5 per cent for \eqn{\alpha} between \eqn{-0.6435} and
#'   \eqn{-0.505}; the gap \code{D} was almost unaffected.
#' @return c(D, a, b, u_c, u_b); NA when the geometry is invalid
#' @examples
#' ## on the equivalence curve the two readings coincide
#' al <- -0.60
#' r <- eq_readings(al, eq_bstar(al))
#' r
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
#' fitted six-parameter theta, the SSE and both readings. \code{eqfit_bc_boot} runs the
#' moving-block residual bootstrap (block length from the residual autocorrelation, as
#' \code{eqfit_blocklen}) of the constrained fit, warm-started, OpenMP over replicates with
#' per-replicate splitmix64 streams; its \code{th} must be a fit on \code{curve}, since every
#' replicate is refitted on the curve, and anything else is an error. \code{eqfit_k4eq} and
#' \code{eqfit_k4eq_boot} are the kappa-family analogues constrained to a locus table h -> k (both
#' without drift); \code{eqfit_k4eq_boot} must be given the locus table of its point fit. \code{eqfit_msim} runs the
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
#' @param maxit Nelder--Mead iteration cap per pass. Kept but without effect in \code{eqfit_msim}
#'   and \code{eqfit_nullT}, whose profiled minimum is a global search with fixed internal limits.
#' @return All of these return plain R objects rather than a class, so they can be compared
#'   element by element against the Python front end.
#'   \code{eqfit_bc} returns a list with \code{th}, the fitted six-parameter vector,
#'   \code{sse}, the residual sum of squares at that fit, and \code{a} and \code{b}, the two
#'   curve readings; \code{eqfit_k4eq} returns the same list without \code{b}.
#'   \code{eqfit_bc_boot}, \code{eqfit_k4eq_boot} and \code{eqfit_nullT} return a numeric vector
#'   holding one statistic per replicate, of length \code{B} for the first two and \code{R2} for
#'   the third.
#'   \code{eqfit_msim} returns a list of eight numeric vectors, each of length \code{Rm}: the test
#'   statistic under the truth (\code{Tlev}), under the displacement arms (\code{Taud},
#'   \code{Tp2}, \code{Tp15}, \code{Tp4}) and at the three reference members (\code{Tref2},
#'   \code{Tref15}, \code{Tref4}).
#'   \code{eqfit_estsim} returns a numeric array of dimension
#'   \code{c(R, 3, 2, length(nsizes))}, indexed by replicate, estimator, parameter and sample size.
#'   \code{eqfit_taus} returns a numeric matrix of \code{R} rows with columns \code{t3},
#'   \code{t4}, \code{al} and \code{be}, the two L-moment ratios and the parameter pair they map to.
#'   \code{k4_ab_sweep} returns a list of the two numeric vectors \code{a} and \code{b}, one entry
#'   per element of \code{k}; \code{k4_b_sweep} and \code{eq_E_sweep} return a single numeric
#'   vector of the same length as \code{k} and \code{betas} respectively.
#'   \code{eqfit_score_at} returns a single number, the score at the given parameter point, and
#'   \code{eqfit_blocklen} a single integer, the selected moving-block length.
#' @examples
#' curve <- cbind(a <- seq(-0.62, -0.51, length.out = 8), sapply(a, eq_bstar))
#' x <- 20 * (0:119) / 119
#' y <- 1 + 300 * (1:120) / 121 + 3 * sin(4321 * (1:120))
#' st <- rbind(c(1, log(300), 0, log(20), 0.3, 0.5))
#' f <- eqfit_bc(x, y, st, curve, maxit = 400L)
#' f$sse
#' eqfit_blocklen(sin(0.05 * (1:300)))
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
  ## the replicates are refitted on the curve, so the point fit must be on it too
  bc <- stats::approx(curve[, 1], curve[, 2], th[5])$y
  if (length(th) != 6L || !is.finite(bc) || abs(th[6] - bc) > 1e-9 * max(1, abs(bc)))
    stop("`th` must be a fit on `curve` (the constrained fit of eqfit_bc): its replicates are refitted on the curve")
  .arc_sync_threads()
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
  .arc_sync_threads()
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
  .arc_sync_threads()
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
  .arc_sync_threads()
  .C("arceqfit_nullT", as.integer(nM), as.double(sde), as.integer(R2), as.integer(seed),
     as.double(th0), as.double(st5), as.integer(maxit),
     as.double(curve[,1]), as.double(curve[,2]), as.integer(nrow(curve)),
     as.double(awin[1]), as.double(awin[2]), out = double(R2))$out
}
#' @rdname eqfit_bc
#' @param a0,b0 beta parameters of the sampled member
#' @export
eqfit_estsim <- function(a0, b0, R, nsizes, seed = 4207L, maxit = 500L){
  .arc_sync_threads()
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
  .arc_sync_threads()
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
#' medians of window nine and twenty-one followed by least squares of the curve without drift
#' (\code{med9-LS}, \code{med21-LS}), the banded arc-length estimator (\code{arc}), the GEV
#' submodel in its own right (\code{GEV-LS}), the selection fit past the run's measured
#' transient cutoff (\code{transient-excised}), and the selection fit on each half-density
#' index grid (\code{grid-odd}, \code{grid-even}). Every pipeline searches the box
#' \code{bounds}, with \eqn{h} itself as the coordinate and starts on both sides of zero, so pass
#' the box of the run's own point fit. Point estimates use the full start grid; a replicate starts
#' from its pipeline's own point fit. Resampling indices come from per-replicate
#' splitmix64 streams, so the result is identical whatever the OpenMP thread count and
#' identical between the R and Python fronts.
#'
#' @param x,y the logged curve
#' @param cut transient cutoff in the units of \code{x}; \code{-Inf} excises nothing
#' @param B bootstrap replicates; \code{0} returns point estimates only
#' @param seed integer seed for the per-replicate streams
#' @param bounds the admissible shape box \code{(k_lo, k_hi, h_lo, h_hi)}
#' @param maxit Nelder--Mead iteration cap of a replicate refit
#' @param drift \code{TRUE} lets \code{med9-LS}, \code{med21-LS} and \code{arc} carry the linear drift
#'   that the other five pipelines always carry. With \code{FALSE} the spread between pipelines is
#'   mostly the effect of leaving the drift out.
#' @return \code{list(a_pt, A)}: the eight named point estimates, and the \code{B} by eight
#'   matrix of replicate readings (\code{NULL} when \code{B = 0})
#' @examples
#' x <- seq(0, 20, length.out = 120)
#' y <- 1 + 300 * k4_cdf(x, 8, 2.5, 0.2, 0.4) + 2 * sin(9871 * seq_along(x))
#' k4_mv_boot(x, y, B = 0L, maxit = 300L)$a_pt
#' @export
k4_mv_boot <- function(x, y, cut = -Inf, B = 50L, seed = 4207L,
                       bounds = c(-0.98, 0.95, 0, 4), maxit = 1500L, drift = FALSE){
  .arc_sync_threads()
  n <- length(x); stopifnot(length(y) == n, B >= 0)
  r <- .C("arck4_mv_boot", as.double(x), as.double(y), as.integer(n),
          as.double(if (is.finite(cut)) cut else -1e300), as.integer(B),
          as.integer(seed), as.double(bounds), as.integer(maxit), as.integer(drift),
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
#' chosen serially so the answer does not depend on the thread count (see \code{\link{arc_threads}}).
#' @param x,y the curve.
#' @param starts a matrix with four columns, \code{(mu, log sigma, k, h)}, one row per start; h is carried directly so
#'   negative values are admissible.
#' @param maxit iterations per Nelder-Mead descent.
#' @param bounds admissible shape box \code{c(k_lo, k_hi, h_lo, h_hi)}. The default lower bound on
#'   \eqn{k} is the one the edible-oil analysis settled on: a wider box lowers the residual sum of
#'   squares but moves the fitted induction periods away from the laboratory values, so the choice is
#'   the analyst's and is stated here rather than fixed in the back end. A lower bound of zero
#'   on \eqn{h} is admissible and meaningful: it admits the generalised extreme value
#'   distribution, the \eqn{h \to 0} member of the family, and with \eqn{k = 0} the Gumbel.
#' @param hfix hold \eqn{h} at this value and fit the remaining three parameters, which is how a
#'   submodel is fitted in its own right; \code{NA} (the default) fits all four. \eqn{h} is searched
#'   directly, so \eqn{h = 0} is an ordinary point of the free search. Choosing between the two
#'   fits is model selection between a nested pair, and is left to the caller.
#' @return a list with \code{theta} in the six-vector convention \code{(g0, g1, mu, sg, k, h)},
#'   the fitted \code{drift}, and the residual sum of squares \code{rss}.
#' @examples
#' x <- seq(0, 20, length.out = 120)
#' y <- 1 + 0.2 * x + 300 * k4_cdf(x, 9, 2, -0.3, -0.5) + sin(4321 * seq_along(x))
#' st <- as.matrix(expand.grid(c(8, 10), log(c(1, 3)), c(-0.5, 0), c(-1, 0.2)))
#' f <- k4_fit_varpro(x, y, st, bounds = c(-0.98, 0.95, -60, 4))
#' f$theta; f$rss
#' @export
k4_fit_varpro <- function(x, y, starts, maxit = 1500L,
                          bounds = c(-0.98, 0.95, 0, 4), hfix = NA_real_) {
  .arc_sync_threads()
  starts <- as.matrix(starts)
  stopifnot(ncol(starts) == 4L, length(bounds) == 4L,
            bounds[1] < bounds[2], bounds[3] < bounds[4])
  ## h is searched directly (column four of starts holds h itself, not log h): the negative
  ## half-line carries real members -- the log-logistic at h = -1, the Burr III / 5PL slice at
  ## h = -1/m -- and h = 0 (GEV) is interior. hfix = NA asks for the free fit; any finite value,
  ## including negative ones, holds h there.
  r <- .C("arck4_fit_varpro", as.double(x), as.double(y), as.integer(length(x)),
          as.double(t(starts)), as.integer(nrow(starts)), as.integer(maxit),
          as.double(bounds), as.double(hfix),
          out = double(8), NAOK = TRUE)$out
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
#' @examples
#' x <- seq(0, 20, length.out = 200)
#' k4_band_trop_mean(x, c(0, 300, 8, 2, -0.3, 0.2), sigma = 1, breaks = c(0, 10, 20))
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
#' @examples
#' set.seed(1)
#' th <- rarccirc(200, arccirc(0.4))
#' arcc_exact_ci(th, B = 99L)[c("lower", "upper")]
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
#' @return A list of two numbers. \code{stat} is the observed arc-length spacings statistic of the
#'   angles after the probability-integral transform through the fitted member, and \code{p} is its
#'   Monte Carlo p-value against the distribution-free reference law. Because the transform sends
#'   the sample to uniform under the hypothesis, the p-value is exact in level for any \code{B}.
#' @examples
#' th <- rarccirc(200, arccirc(0.4, 0.3))
#' arcc_gof(th, fit_arccirc(th), B = 99L)
#' @export
arcc_gof <- function(theta, fit, B = 999L, seed = 4207L){
  o <- .C(C_arcc_gof, as.double(theta), as.integer(length(theta)), as.double(fit$c3),
          as.double(fit$mu), as.integer(B), as.integer(seed), out = double(2))$out
  list(stat = o[1], p = o[2])
}

## ================= threads ============================================================

.arc_sync_threads <- function() {
  n <- getOption("arcstat.threads", 0L)
  n <- if (is.numeric(n) && length(n) == 1L && is.finite(n) && n >= 1) as.integer(n) else 0L
  .C(C_arc_set_threads, n, now = integer(1))$now
}

#' Number of threads used by the back end
#'
#' Every parallel loop of the package uses the same number of threads: the option
#' \code{arcstat.threads} if it is set, otherwise the environment variable \code{ARCSTAT_THREADS},
#' otherwise at most two. Results do not depend on the number.
#' @param n a positive integer to set the option, \code{NULL} to clear it; missing only asks.
#' @return the number of threads in force (1 where the package was built without OpenMP), invisibly
#'   when it was set.
#' @examples
#' arc_threads()
#' @export
arc_threads <- function(n) {
  if (missing(n)) return(.arc_sync_threads())
  stopifnot(is.null(n) || (is.numeric(n) && length(n) == 1L && is.finite(n) && n >= 1))
  options(arcstat.threads = if (is.null(n)) NULL else as.integer(n))
  invisible(.arc_sync_threads())
}

## ================= fits and bootstraps over many runs ==================================

## one run as two vectors, or many as two lists: concatenated, with offsets
.arc_runs <- function(x, y) {
  if (!is.list(x)) x <- list(x)
  if (!is.list(y)) y <- list(y)
  n <- lengths(x)
  stopifnot(length(x) == length(y), length(x) >= 1L, all(n == lengths(y)), all(n >= 8L))
  X <- as.double(unlist(x, use.names = FALSE)); Y <- as.double(unlist(y, use.names = FALSE))
  if (!all(is.finite(X)) || !all(is.finite(Y))) stop("x and y must be finite")
  list(X = X, Y = Y, off = as.integer(c(0, cumsum(n))), R = length(n), names = names(x))
}
.arc_locus <- function(locus) {
  if (is.null(locus)) return(list(np = 0L, h = 0, k = 0, off = 0L))
  if (!is.list(locus) || is.data.frame(locus)) locus <- list(locus)
  locus <- lapply(locus, function(m) { m <- as.matrix(m)
    stopifnot(ncol(m) == 2L, nrow(m) >= 2L, all(is.finite(m)))
    m <- m[order(m[, 1]), , drop = FALSE]
    if (any(diff(m[, 1]) <= 0)) stop("the h of a locus table must be distinct")
    m })
  list(np = length(locus), h = as.double(unlist(lapply(locus, function(m) m[, 1]))),
       k = as.double(unlist(lapply(locus, function(m) m[, 2]))),
       off = as.integer(c(0, cumsum(vapply(locus, nrow, 1L)))))
}

#' The equivalence locus of the four-parameter kappa family
#'
#' For every \eqn{h<1/2} there is one \eqn{k} at which the tangent reading and the third-derivative
#' reading coincide. The locus is a single curve through \eqn{h=0}; it tends to the corner
#' \eqn{(h,k)=(1/2,1)}, and for \eqn{h<0} it approaches the frontier \eqn{hk=1/2} from inside,
#' exponentially fast in \eqn{h}: below \eqn{h} of about \eqn{-15} the value returned is the double
#' nearest to \eqn{1/(2h)}.
#' @param h shape values.
#' @return the \eqn{k} of each, \code{NaN} where \eqn{h\ge 1/2}.
#' @examples
#' k <- k4_locus_k(c(-1, 0, 0.3))
#' r <- k4_readings(c(0, 1, 0, 1, k[3], 0.3))
#' r[["a"]] - r[["b"]]
#' @export
k4_locus_k <- function(h)
  .C(C_arcbt_locus_k, as.double(h), as.integer(length(h)), k = double(length(h)), NAOK = TRUE)$k

.arc_gate <- function(gate) {
  if (identical(gate, "best")) return(1L)
  if (identical(gate, FALSE)) return(0L)
  stop("`gate` must be FALSE or \"best\"")
}
.arc_status <- function(v, d)
  matrix(v, d$R, 5, byrow = TRUE, dimnames = list(d$names, c("status", "inf_rss", "inf_k", "inf_h", "inf_a")))

#' Equivalence-constrained kappa fits with drift, many runs in one call
#'
#' Fits \eqn{y = g_0 + m x + g_1 F(x;\mu,\sigma,k,h)} with \eqn{(h,k)} held on the equivalence locus,
#' so that the two standard readings coincide. The linear part \eqn{(g_0,m,g_1)} is solved exactly at
#' every trial shape; \eqn{(\mu,\log\sigma,h)} is searched by Levenberg--Marquardt from a grid of
#' starts, the best of which are carried to convergence and polished by a simplex. The starts of all
#' runs share one pool of threads (\code{\link{arc_threads}}); results do not depend on the number.
#'
#' Unless \code{gate = "best"}, no fit is refused on physical grounds. The sign of \eqn{g_1} and the position of the
#'   reading are
#' the caller's to judge; a caller who gates the free fit in that way must gate both fits, or a ratio
#' of the two sums of squares can fall below one.
#' @param x,y one run as two vectors, or many runs as two lists of vectors.
#' @param locus \code{NULL} for the whole locus, computed inside (\code{\link{k4_locus_k}}); or a
#'   two-column matrix \code{(h, k)}, or a list of such matrices (pieces), interpolated linearly and
#'   not extrapolated. An interpolated table puts a kink in the objective at each of its nodes; all
#'   nodes are scanned, so the best one is found, but the fit is then a property of the table.
#' @param bounds the box \code{(k_lo, k_hi, h_lo, h_hi)} the shape must lie in.
#' @param drift \code{FALSE} holds the slope \eqn{m} at zero.
#' @param gate \code{FALSE} takes the lowest minimum on the locus; \code{"best"} the lowest interior
#'   minimum that is physical: \code{g1 > 0} and the tangent reading between 0 and \code{max(x)}. Where
#'   there is none the fit is \code{NaN}: the edge of the physical set, or the end of the searched part
#'   of the locus, is never returned as a fit. The
#'   result then has \code{physical}, one row per run: \code{status}, 1 for an interior minimum that is
#'   the lowest physical value found, 2 for an interior minimum with lower physical values towards the
#'   edge of the set, 3 for physical points without an interior minimum, 0 for no physical point; and
#'   \code{inf_rss, inf_k, inf_h, inf_a}, the lowest physical sum of squares found and the shape and
#'   reading there, which in states 2 and 3 is a limit and not a fit.
#' @param it1,it2 Levenberg--Marquardt steps per start, and for the best starts.
#' @param keep how many distinct best starts of each run are carried to convergence.
#' @return a list: \code{fit}, a matrix with one row per run and columns \code{g0, m, g1, mu, sigma,
#'   k, h, rss, a, b} and \code{piece} (the table piece chosen; for the whole locus, which is
#'   searched from both sides of zero, 1 for \eqn{h\le 0} and 2 for \eqn{h\ge 0}); \code{a} and
#'   \code{b} are the two readings of the fitted shape, whatever the sign of \code{g1}; and
#'   \code{pieces}, the best fit on every piece in the same columns (runs by columns by pieces), from
#'   which a caller can take the best one that is acceptable on other grounds.
#' @examples
#' x <- seq(0, 20, length.out = 150)
#' y <- 1 + 0.2 * x + 300 * k4_cdf(x, 9, 2, k4_locus_k(0.2), 0.2) + sin(4321 * seq_along(x))
#' f <- k4eq_fit(x, y)
#' f$fit[1, c("k", "h", "rss", "a", "b")]
#' @export
k4eq_fit <- function(x, y, locus = NULL, bounds = c(-0.98, 0.95, -60, 4), drift = TRUE,
                     gate = FALSE, it1 = 12L, it2 = 150L, keep = 6L) {
  d <- .arc_runs(x, y); L <- .arc_locus(locus); np <- if (L$np > 0L) L$np else 2L
  stopifnot(length(bounds) == 4L, bounds[1] < bounds[2], bounds[3] < bounds[4], keep >= 1L)
  .arc_sync_threads(); g <- .arc_gate(gate)
  r <- .C(C_arcbt_k4eq, d$X, d$Y, d$off, d$R, L$np, L$h, L$k, L$off, as.double(bounds),
          as.integer(drift), as.integer(it1), as.integer(it2), as.integer(keep),
          g, out = double(11 * d$R), p = double(11 * np * d$R),
          s = double(5 * d$R), NAOK = TRUE)
  nm <- c("g0", "m", "g1", "mu", "sigma", "k", "h", "rss", "a", "b", "piece")
  z <- list(fit = matrix(r$out, d$R, 11, byrow = TRUE, dimnames = list(d$names, nm)),
            pieces = aperm(array(r$p, c(11L, np, d$R), list(nm, NULL, d$names)), c(3L, 1L, 2L)))
  if (g) z$physical <- .arc_status(r$s, d)
  z
}

#' Beta-companion fits with drift, free and equivalence-constrained, many runs in one call
#'
#' Fits \eqn{y = g_0 + m x + g_1 F((x-s)/w;\alpha,\beta)}, \eqn{F} the inverse of the regularised
#' incomplete beta function with parameters \eqn{\alpha+1,\beta+1}, once with \eqn{\alpha,\beta} free
#' in \eqn{[-1+}\code{eps}\eqn{, 20]} and once with \eqn{\beta} on the equivalence curve. \eqn{F} is exactly 0
#' and 1 off the support. The method is that of \code{\link{k4eq_fit}}. The constrained optimum is a
#' point of the free family and is also a start of the free fit, so the constrained sum of squares is
#' never below the free one.
#' @inheritParams k4eq_fit
#' @param curve \code{NULL} for the whole equivalence curve with beta of at least \code{-1 + eps},
#'   computed inside (\code{\link{eq_bstar}}); or a two-column matrix \code{(alpha, beta)} of it,
#'   increasing in alpha, interpolated linearly, in which alpha is held.
#' @param eps how close to -1 the exponents may come.
#' @param dense \code{TRUE} (the default) uses 16 (free) and 8 (constrained) times as many starts as
#'   \code{FALSE}, which is about three times faster. On the 98 induction runs of the archive the quick
#'   search stopped above the dense one on 84 constrained fits (by up to a factor of 14 on one run);
#'   the dense search left no gap of more than 7.4e-6 relative to a heavier search (16 kept starts,
#'   600 iterations) except where the fit runs along a boundary ridge: an exponent tending to -1 while
#'   the scale \code{w} grows, on which the sum of squares falls slowly and the minimum is an infimum.
#'   There the heavier search went lower by at most 0.16 per cent on 27 constrained fits (24 with beta
#'   at the boundary) and 0.65 per cent on 3 free fits. The readings moved by less than 4e-4 except on
#'   the degenerate Marula runs, whose readings lie far before the data (below -37) and moved by up to
#'   240; on such runs a reading is not determined by the data.
#' @return a list of two matrices, \code{free} and \code{con}, with one row per run and columns
#'   \code{g0, m, g1, s, w, alpha, beta, rss, a, b}; \code{a} and \code{b} are the two readings.
#' @examples
#' curve <- cbind(a <- seq(-0.62, -0.51, length.out = 6), sapply(a, eq_bstar))
#' x <- seq(0, 20, length.out = 50)
#' y <- 2 + 0.1 * x + 250 * bc_cdf((x - 1) / 17 * beta(0.4, 0.6), -0.6, -0.4) + cos(999 * x)
#' h <- bc_host_fit(x, y, curve)
#' h$con[, "rss"] / h$free[, "rss"]
#' @export
bc_host_fit <- function(x, y, curve = NULL, it1 = 15L, it2 = 200L, keep = 6L, eps = 1e-6,
                        dense = TRUE) {
  d <- .arc_runs(x, y)
  curve <- if (is.null(curve)) matrix(0, 0L, 2L) else as.matrix(curve)
  stopifnot(nrow(curve) != 1L, eps >= 1e-12, eps < 0.01, ncol(curve) == 2L, all(is.finite(curve)),
            all(diff(curve[, 1]) > 0))
  .arc_sync_threads()
  r <- .C(C_arcbt_host, d$X, d$Y, d$off, d$R, as.double(curve[, 1]), as.double(curve[, 2]),
          nrow(curve), as.double(eps), as.integer(dense), as.integer(it1), as.integer(it2), as.integer(keep),
          f = double(10 * d$R), c = double(10 * d$R), NAOK = TRUE)
  nm <- list(d$names, c("g0", "m", "g1", "s", "w", "alpha", "beta", "rss", "a", "b"))
  list(free = matrix(r$f, d$R, 10, byrow = TRUE, dimnames = nm),
       con  = matrix(r$c, d$R, 10, byrow = TRUE, dimnames = nm))
}

.arc_seeds <- function(seed, R) {
  if (length(seed) == 1L) seed <- seed + seq_len(R) - 1L
  stopifnot(length(seed) == R, all(is.finite(seed)))
  as.integer(seed)
}
## the drift setting of every run of a point fit: by default read off the fit, whose slope m is exactly 0
## when the fit was made with drift = FALSE (the back end holds it there); a run given drift = FALSE whose
## fit carries a slope is refused, because its replicates would be drawn from one mean and fitted by another
.arc_drift <- function(drift, m, R) {
  m <- as.double(m)
  if (is.null(drift)) return(as.integer(!(is.finite(m) & m == 0)))
  drift <- rep_len(as.logical(drift), R)
  stopifnot(!anyNA(drift))
  if (any(!drift & is.finite(m) & m != 0))
    stop("`drift = FALSE` for a fit with a nonzero slope m: the bootstrap must use the point fit's own setting")
  as.integer(drift)
}
.arc_bootarr <- function(v, B, R, cols, runs)
  aperm(array(v, c(length(cols), B, R), list(cols, NULL, runs)), c(2L, 1L, 3L))

#' Residual bootstraps of kappa fits, many runs in one call
#'
#' \code{k4_boot} resamples the residuals of a fitted free kappa curve in moving blocks and refits it
#' by the variable-projection search of \code{\link{k4_fit_varpro}} from the fitted shape, each run in
#' its own box and, where the point fit is the generalised extreme value submodel, with \eqn{h} held at
#' zero. \code{k4eq_boot} does the same for a fit of \code{\link{k4eq_fit}}, refitted on the piece of
#' the locus it lies on. Replicate \eqn{b} of a run draws from its own stream, determined by the run's
#' seed and \eqn{b} alone, so the result does not depend on the number of threads or on which other runs
#' are in the call.
#'
#' Every replicate is fitted under the settings of the point fit: the same box or bounds and locus, the
#' same submodel, the same rule on physical fits (\code{gate}) and the same drift. The drift setting is
#' read off the fit by default: a fit made with \code{drift = FALSE} has a slope \code{m} of exactly zero,
#' and its replicates are then drawn from the mean without drift and refitted without it. (Up to version
#' 0.4.0 the replicates were always refitted with a drift term, which inflated the spread of a fit made
#' without one.) The box, bounds, locus and gate cannot be read off the fit and must be passed as they
#' were given to the point fit.
#'
#' The block length is the first lag at which the residual autocorrelation is inside
#' \eqn{\pm 2/\sqrt n} (at least 2; \eqn{n^{1/3}} if none of the first 200 is); \code{block = 1}
#' resamples singly and any larger value fixes the length.
#' @inheritParams k4eq_fit
#' @param fit a matrix with one row per run whose columns \code{g0, m, g1, mu, sigma, k, h} are the
#'   fit (for \code{k4eq_boot} the \code{fit} of \code{\link{k4eq_fit}}, with its \code{piece}).
#' @param box a matrix with one row per run, or one row for all: the box \code{(k_lo, k_hi, h_lo,
#'   h_hi)} of the run's own point fit.
#' @param gev logical, one per run: the point fit holds \eqn{h = 0}. \code{NULL} (the default) reads it from
#'   a column \code{gev} of \code{fit} (as returned by \code{\link{k4_fit}}), and is \code{FALSE} without one.
#' @param drift logical, one per run or one for all: the point fit carries the drift \eqn{m x}.
#'   \code{NULL} (the default) takes it from the fit, \code{FALSE} where \code{m} is exactly zero. A run
#'   given \code{FALSE} whose fit has a nonzero slope is an error.
#' @param B replicates.
#' @param seed one seed per run, or one for all, in which case run \eqn{r} uses \code{seed + r - 1}.
#' @param block 0 for the block length of the residual autocorrelation, 1 for single resampling, or
#'   a block length.
#' @param maxit,it2 Nelder--Mead iterations per pass of a free refit; Levenberg--Marquardt steps of a
#'   constrained one.
#' @param gate the rule of the point fit. For \code{k4eq_boot}, \code{"best"} keeps a constrained
#'   replicate only if its refit, a local search from the fitted point, ends at an interior physical
#'   minimum; for \code{k4_boot}, \code{"best"} keeps a replicate only if its refit is physical
#'   (\code{g1 > 0} and the tangent reading between 0 and \code{max(x)}). The others are \code{NaN}.
#'   Pass \code{"best"} when the point fit was made with \code{gate = "best"}.
#' @return a list: \code{boot}, an array of \code{B} rows by columns by runs. The columns are the two
#'   readings \code{a} and \code{b}, the fitted \code{k} and \code{h}, and for \code{k4_boot} also
#'   \code{sel}, 1 if \eqn{n\log(rss_0/rss)-2<0} would choose \eqn{h=0} on that replicate and
#'   \code{NaN} where \eqn{h} is held; for both, last, \code{m}, the fitted slope of the replicate, exactly
#'   zero when the point fit has no drift; and \code{block}, the block length of every run. \code{k4eq_boot} adds
#'   \code{lost}, the number of replicates of every run without a refit.
#' @examples
#' x <- seq(0, 20, length.out = 120)
#' y <- 1 + 300 * k4_cdf(x, 9, 2, -0.3, 0.2) + sin(4321 * seq_along(x))
#' f <- k4_fit_varpro(x, y, rbind(c(9, log(2), -0.3, 0.2)), bounds = c(-0.98, 0.95, -60, 4))
#' fit <- rbind(c(f$theta[1], f$drift, f$theta[2:6]))
#' b <- k4_boot(x, y, fit, box = c(-0.98, 0.95, -60, 4), B = 5)
#' b$boot[, c("a", "b"), 1]
#' @export
k4_boot <- function(x, y, fit, box, gev = NULL, B = 400L, seed = 4207L, block = 0L, maxit = 800L,
                    drift = NULL, gate = FALSE) {
  d <- .arc_runs(x, y)
  fm <- as.matrix(fit)
  stopifnot(nrow(fm) == d$R, ncol(fm) >= 7L)
  if (is.null(gev)) gev <- if ("gev" %in% colnames(fm)) fm[, "gev"] == 1 else FALSE
  dr <- .arc_drift(drift, fm[, 2L], d$R)
  fit <- matrix(as.double(fm[, 1:7, drop = FALSE]), nrow = d$R)
  box <- if (is.matrix(box)) box else matrix(box, 1L)
  if (nrow(box) == 1L) box <- box[rep(1L, d$R), , drop = FALSE]
  stopifnot(ncol(box) == 4L, nrow(box) == d$R, B >= 1L, block >= 0L)
  .arc_sync_threads()
  r <- .C(C_arcbt_k4_boot, d$X, d$Y, d$off, d$R, as.double(t(fit)), as.double(t(box)),
          as.integer(rep_len(as.logical(gev), d$R)), dr, .arc_gate(gate), as.integer(B), .arc_seeds(seed, d$R),
          as.integer(block), as.integer(maxit), bl = integer(d$R),
          out = double(6 * B * d$R), NAOK = TRUE)
  list(boot = .arc_bootarr(r$out, B, d$R, c("a", "b", "k", "h", "sel", "m"), d$names), block = r$bl)
}
#' @rdname k4_boot
#' @export
k4eq_boot <- function(x, y, fit, locus = NULL, bounds = c(-0.98, 0.95, -60, 4), B = 400L,
                      seed = 4207L, block = 0L, it2 = 100L, gate = FALSE, drift = NULL) {
  d <- .arc_runs(x, y); L <- .arc_locus(locus)
  fit <- as.matrix(fit)
  stopifnot(nrow(fit) == d$R, ncol(fit) >= 7L, length(bounds) == 4L, B >= 1L, block >= 0L)
  piece <- if (L$np == 0L) rep(1L, d$R) else if ("piece" %in% colnames(fit)) fit[, "piece"] else
    stop("with a locus table `fit` needs its `piece` column")
  stopifnot(L$np == 0L || all(piece >= 1 & piece <= L$np))
  dr <- .arc_drift(drift, fit[, 2L], d$R)
  .arc_sync_threads()
  r <- .C(C_arcbt_k4eq_boot, d$X, d$Y, d$off, d$R, as.double(t(fit[, 1:7, drop = FALSE])),
          as.integer(piece), L$np, L$h, L$k, L$off, as.double(bounds), as.integer(B),
          .arc_seeds(seed, d$R), as.integer(block), as.integer(it2), .arc_gate(gate), dr, bl = integer(d$R),
          out = double(5 * B * d$R), NAOK = TRUE)
  b <- .arc_bootarr(r$out, B, d$R, c("a", "b", "k", "h", "m"), d$names)
  list(boot = b, block = r$bl,
       lost = colSums(!is.finite(b[, "h", , drop = FALSE]), dims = 2L))
}

#' Free kappa fits with drift, many runs in one call
#'
#' The fit of \code{\link{k4_fit_varpro}} from a built-in start grid, by the method of
#' \code{\link{k4eq_fit}} (Levenberg--Marquardt on the projected residual, best starts carried to
#' convergence and polished), with the submodel selection and the rule on acceptable fits as options.
#' @inheritParams k4eq_fit
#' @param bounds the box \code{(k_lo, k_hi, h_lo, h_hi)}; the submodel is always fitted in it.
#' @param floors lower bounds on \eqn{h}, deepest first, tried in turn when \code{gate = TRUE}; otherwise
#'   \code{bounds[3]} is used.
#' @param gate \code{FALSE} takes the lowest minimum whatever it is; \code{TRUE} takes the first of \code{floors}
#'   whose best fit is physical, meaning that
#'   \code{g1 > 0} and the tangent reading lies between 0 and \code{max(x)}, and the fit at the last
#'   floor if none is; \code{"best"} takes the lowest minimum that is physical in the box
#'   \code{bounds}, and \code{NaN} where none is found. \code{TRUE} depends on the search: the better
#'   it is, the more often an unphysical minimum is the lowest at every floor.
#' @param gev fit also the generalised extreme value submodel \eqn{h=0} and choose it when
#'   \eqn{n\log(rss_0/rss)-2<0} (with \code{gate}, only if it is physical).
#' @param hfix hold \eqn{h} there; \code{NA} searches it.
#' @param infimum with \code{gate = "best"}, also report the status and the infimum over the physical
#'   set of the four-parameter family, as columns \code{status, inf_rss, inf_k, inf_h, inf_a} (see
#'   \code{\link{k4eq_fit}}). They describe the search with h free; the fit does not change.
#' @param dense use about ten times as many starts at negative \eqn{k} (in the coordinates of the lower end
#'   of the support, which the search also uses for \eqn{k < -0.15}).
#' @return a matrix with one row per run and columns \code{g0, m, g1, mu, sigma, k, h, rss, a, b, gev,
#'   dAIC, h_floor, physical}. \code{gev} is 1 where the submodel was chosen.
#' @examples
#' x <- seq(0, 20, length.out = 120)
#' y <- 1 + 0.2 * x + 300 * k4_cdf(x, 9, 2, -0.3, -0.5) + sin(4321 * seq_along(x))
#' k4_fit(x, y)[1, c("k", "h", "rss", "a", "b")]
#' @export
k4_fit <- function(x, y, bounds = c(-0.98, 0.95, -60, 4), floors = c(-60, -8, -3, -1.2, -0.5, 0),
                   gate = FALSE, gev = FALSE, hfix = NA_real_, drift = TRUE, it1 = 12L, it2 = 200L, keep = 8L,
                   infimum = FALSE, dense = FALSE) {
  d <- .arc_runs(x, y)
  stopifnot(length(bounds) == 4L, bounds[1] < bounds[2], bounds[3] < bounds[4], length(hfix) == 1L)
  gate <- if (identical(gate, "best")) 2L else as.integer(isTRUE(gate))
  if (infimum && (gate != 2L || is.finite(hfix))) stop("`infimum` goes with gate = \"best\" and free h")
  fl <- if (gate == 1L) floors else bounds[3]
  stopifnot(length(fl) >= 1L, all(is.finite(fl)), all(fl < bounds[4]))
  .arc_sync_threads()
  r <- .C(C_arcbt_k4fit, d$X, d$Y, d$off, d$R, as.double(bounds), as.double(fl), length(fl),
          as.integer(gate), as.integer(gev), as.double(hfix), as.integer(drift), as.integer(it1),
          as.integer(it2), as.integer(keep), as.integer(infimum), as.integer(dense), out = double(14 * d$R),
          s = double(5 * d$R), NAOK = TRUE)
  z <- matrix(r$out, d$R, 14, byrow = TRUE, dimnames = list(d$names, c("g0", "m", "g1", "mu", "sigma", "k",
         "h", "rss", "a", "b", "gev", "dAIC", "h_floor", "physical")))
  if (infimum) cbind(z, .arc_status(r$s, d)) else z
}

#' Wald (inverse Gaussian) curve fits, many runs in one call
#'
#' Fits \eqn{y = g_0 + m x + g_1 F(x;\mu,\lambda)} with \eqn{F} the inverse Gaussian distribution
#' function, by the method of \code{\link{k4eq_fit}}, and reads the fitted curve exactly. Time is
#' measured from zero, so \code{x} must be in the units and from the origin the curve is meant in.
#' @inheritParams k4eq_fit
#' @return \code{wald_fit}: a matrix with one row per run and columns \code{g0, m, g1, mu, lambda, phi,
#'   rss, a, b, mode}; \code{phi} is \eqn{\lambda/\mu}; \code{a} and \code{b} are the tangent and the
#'   third-derivative reading of the fitted shape, whatever the sign of \code{g1} and wherever they
#'   lie, and can be of any size when the fit runs to the end of its range (\eqn{\mu} up to
#'   \eqn{10^8}). \code{wald_cdf}: a vector. \code{wald_readings}: a matrix with columns \code{a},
#'   \code{b}, \code{mode}; both readings exist for every \eqn{\phi>0}.
#' @examples
#' x <- seq(0, 20, length.out = 120)
#' y <- 1 + 0.2 * x + 300 * wald_cdf(x, 9, 150) + sin(4321 * seq_along(x))
#' wald_fit(x, y)[1, c("mu", "lambda", "a", "b")]
#' wald_readings(9, 150)
#' @export
wald_fit <- function(x, y, drift = TRUE, it1 = 25L, it2 = 200L, keep = 6L) {
  d <- .arc_runs(x, y); .arc_sync_threads()
  r <- .C(C_arcbt_wald, d$X, d$Y, d$off, d$R, as.integer(drift), as.integer(it1), as.integer(it2),
          as.integer(keep), out = double(10 * d$R), NAOK = TRUE)$out
  matrix(r, d$R, 10, byrow = TRUE, dimnames = list(d$names, c("g0", "m", "g1", "mu", "lambda", "phi",
         "rss", "a", "b", "mode")))
}
#' @rdname wald_fit
#' @param mu,lambda positive parameters, one pair for \code{wald_cdf}, any number for
#'   \code{wald_readings}.
#' @export
wald_cdf <- function(x, mu, lambda) {
  stopifnot(length(mu) == 1L, length(lambda) == 1L, mu > 0, lambda > 0)
  .C(C_arcbt_wald_cdf, as.double(x), as.integer(length(x)), as.double(mu), as.double(lambda),
     out = double(length(x)), NAOK = TRUE)$out
}
#' @rdname wald_fit
#' @export
wald_readings <- function(mu, lambda) {
  n <- max(length(mu), length(lambda)); mu <- rep_len(mu, n); lambda <- rep_len(lambda, n)
  r <- .C(C_arcbt_wald_readings, as.double(mu), as.double(lambda), as.integer(n),
          out = double(3 * n), NAOK = TRUE)$out
  matrix(r, n, 3, byrow = TRUE, dimnames = list(NULL, c("a", "b", "mode")))
}

#' Windowed L-moment matching of kappa shapes, many runs in one call
#'
#' Puts each run on the probability scale by the linear part of a fit only, takes the third and fourth
#' L-moment ratios of its quantile function over the observed window of levels, and finds the kappa
#' shape whose quantile function has the same ratios over the same window. The fitted shape is not
#' used. \eqn{h} is searched on both sides of zero.
#' @inheritParams k4eq_fit
#' @param lin a matrix with one row per run whose columns \code{g0, m, g1} (or the first three) are the
#'   linear part.
#' @param window running-median window (odd; the ends as in \code{\link[stats]{runmed}}).
#' @param ngrid,nmodel number of equally spaced levels for the data and for the model.
#' @return a matrix with one row per run and columns \code{k, h, obj, t3, t4, lo, hi, nsol,
#'   h_min, h_max}: the shape, the squared distance left, the two ratios of the data, the window, the
#'   number of distinct exact matches found and their smallest and largest \code{h}. Two ratios over a
#'   window do not always determine the shape: where \code{nsol} exceeds one the shape returned is only
#'   the one with the smallest \eqn{|h|}. \code{NaN} where fewer than 100
#'   points lie inside the window.
#' @examples
#' x <- seq(0, 20, length.out = 300)
#' y <- 1 + 300 * k4_cdf(x, 9, 2, -0.2, 0.3)
#' k4_lmom_window(x, y, c(1, 0, 300))
#' @export
k4_lmom_window <- function(x, y, lin, bounds = c(-0.98, 0.95, -8, 4), window = 21L, ngrid = 3000L,
                           nmodel = 4000L) {
  d <- .arc_runs(x, y)
  lin <- if (is.matrix(lin)) lin else matrix(lin, 1L)
  if (!is.null(colnames(lin)) && all(c("g0", "m", "g1") %in% colnames(lin)))
    lin <- lin[, c("g0", "m", "g1"), drop = FALSE]
  stopifnot(nrow(lin) == d$R, ncol(lin) >= 3L, length(bounds) == 4L, window %% 2 == 1, window >= 1,
            all(diff(d$off) > window), ngrid >= 50L, nmodel >= 50L)
  .arc_sync_threads()
  r <- .C(C_arcbt_lmomw, d$X, d$Y, d$off, d$R, as.double(t(lin[, 1:3, drop = FALSE])), as.double(bounds),
          as.integer(window), as.integer(ngrid), as.integer(nmodel), out = double(10 * d$R), NAOK = TRUE)$out
  matrix(r, d$R, 10, byrow = TRUE, dimnames = list(d$names, c("k", "h", "obj", "t3", "t4", "lo", "hi",
         "nsol", "h_min", "h_max")))
}

#' Variance components of a multiverse of pipelines
#'
#' From the replicate readings of \code{\link{k4_mv_boot}}. All pipelines are refitted on the same
#' replicates, so the column means share their resampling error. With \eqn{v} the mean variance of a
#' column over replicates and \eqn{w} that mean less the mean covariance of two columns (the
#' method-by-replicate variance), the variance of the \eqn{M} column means has expectation
#' \eqn{\delta^2 + w/B}, where \eqn{\delta^2} is the variance of the pipelines' differences \emph{on
#' this curve}. Subtracting \eqn{v/B}, as for independent methods, is too much when the columns are
#' positively correlated. \eqn{\delta^2} still contains the sampling error of the curve itself, which
#' affects the pipelines differently: \eqn{E\delta^2 = \tau^2 + w}. Use \code{delta2} to describe the
#' spread on the curve at hand and \code{tau2} for a method component that would remain with the
#' sampling error taken out. The step from \eqn{\delta^2} to \eqn{\tau^2} takes the bootstrap
#' covariance as that of the estimates.
#' @param A a matrix of readings, replicates by pipelines.
#' @param frac pipelines with a smaller share of finite readings are left out; then only replicates
#'   finite in all the others are used.
#' @return a named vector: \code{M} and \code{B} used, \code{v}, \code{w}, \code{delta2},
#'   \code{share_delta} (\eqn{\delta^2/(\delta^2+v)}), \code{tau2} and \code{share_tau}; \code{NaN} with
#'   fewer than three pipelines or replicates.
#' @examples
#' set.seed(1); e <- rnorm(200)
#' A <- cbind(e, e + 0.1, e + 0.2)      # differences without error: w = 0
#' k4_mv_share(A)
#' @export
k4_mv_share <- function(A, frac = 0.6) {
  A <- as.matrix(A)
  stats::setNames(.C(C_arcbt_mvshare, as.double(t(A)), nrow(A), ncol(A), as.double(frac), out = double(8),
                     NAOK = TRUE)$out,
                  c("M", "B", "v", "w", "delta2", "share_delta", "tau2", "share_tau"))
}

#' Kappa curves that share one shape, for many sets of runs in one call
#'
#' Every run has its own baseline, drift, plateau, location and scale, and the runs of a set share one
#' shape: \eqn{y = g_{0j} + m_j x + g_{1j} F((x-\mu_j)/\sigma_j; k, h)}. This is the model of time
#' rescaling, in which the curves differ in the clock and the amplitude only. The sum of squares is
#' minimised over all of it. At a fixed shape the runs separate, so the profile of every run over a
#' grid of shapes is computed once and shared by all the sets; the lowest minima of the summed profile
#' start a joint search, and the best results are checked run by run from cold starts and polished.
#'
#' The search is not certified. Its known weak case is a run that is fitted best by a curve with a very
#' long scale and a plateau far beyond the data, so that the curve is close to a power of a linear
#' function and its sum of squares falls slowly along a ridge (plateaux of \eqn{10^4} to \eqn{10^8}
#' times the rise observed). On sets containing such runs the result can lie above the minimum, by up to
#' 1.3 per cent on the archive sets where this was measured (0.9 per cent with \code{dense = TRUE});
#' sets of ordinary runs were not affected.
#' @inheritParams k4eq_fit
#' @param x,y many runs as two lists of vectors.
#' @param sets a list of integer vectors, each the indices of the runs of one set; by default one set
#'   of all the runs.
#' @param drift \code{FALSE} holds every \eqn{m_j} at zero.
#' @param dense a grid four times as fine, with more starts.
#' @return a list: \code{fit}, a matrix with one row per set and columns \code{rss, k, h}; and
#'   \code{runs}, a list with one matrix per set, its rows the runs of the set and its columns
#'   \code{g0, m, g1, mu, sigma, rss}. No run is refused on physical grounds.
#' @examples
#' x <- seq(0, 20, length.out = 60)
#' f <- function(mu, s, g) 1 + 0.2 * x + g * k4_cdf(x, mu, s, -0.4, 0.3) + sin(4321 * mu * x)
#' r <- k4_shared_fit(list(x, x), list(f(8, 2, 300), f(10, 2.5, 280)),
#'                    bounds = c(-0.9, 0.5, -1, 1))
#' r$fit
#' @export
k4_shared_fit <- function(x, y, sets = NULL, bounds = c(-0.98, 0.95, -3, 4), drift = TRUE, dense = FALSE) {
  d <- .arc_runs(x, y)
  if (is.null(sets)) sets <- list(seq_len(d$R))
  n <- lengths(sets); s <- as.integer(unlist(sets))
  stopifnot(is.list(sets), length(sets) >= 1L, all(n >= 1L), all(s >= 1L & s <= d$R), length(bounds) == 4L,
            bounds[1] < bounds[2], bounds[3] < bounds[4])
  .arc_sync_threads()
  r <- .C(C_arcsh_fit, d$X, d$Y, d$off, d$R, length(sets), s - 1L, as.integer(c(0, cumsum(n))),
          as.double(bounds), as.integer(drift), as.integer(dense), out = double(3 * length(sets)),
          run = double(6 * length(s)), NAOK = TRUE)
  m <- matrix(r$run, ncol = 6L, byrow = TRUE, dimnames = list(NULL, c("g0", "m", "g1", "mu", "sigma", "rss")))
  list(fit = matrix(r$out, ncol = 3L, byrow = TRUE, dimnames = list(names(sets), c("rss", "k", "h"))),
       runs = lapply(split(seq_along(s), rep(seq_along(n), n)), function(i) m[i, , drop = FALSE]))
}

#' Permutation test of one shared kappa shape against a shape per group
#'
#' For the runs of each problem, \eqn{RSS_{TR}} is the sum of squares of one shape for all of them and
#' \eqn{RSS_W} that of a shape for each group, \eqn{\Lambda_T = RSS_{TR}/RSS_W} and, given the sum of
#' squares of the free fits, \eqn{\Lambda_W = RSS_W/RSS_{free}}. Every set of runs of a group's size is
#' fitted once by the search of \code{\link{k4_shared_fit}} and every assignment of the group labels is
#' evaluated from those fits; the p-value is the share of assignments with an \eqn{RSS_W} not above the
#' observed one (\eqn{\Lambda_T} decreases in it).
#' @inheritParams k4_shared_fit
#' @param group a group label per run.
#' @param problem a problem label per run, when many tests are made in one call; the runs of a problem
#'   must be consecutive.
#' @param rss_free sum of squares of the free fits, one per problem, or \code{NULL}.
#' @param slack numerical allowances: for each, assignments with an \eqn{RSS_W} of at most
#'   \eqn{(1 + slack)} times the observed one are counted against rejection. The fits are the minima a
#'   search found, and an allowance shows whether the verdict depends on its last per cent.
#' @param maxalloc no p-value with more assignments than this, or with more than 14 runs.
#' @return a list: \code{stat}, a matrix with one row per problem and columns \code{rss_tr, rss_w,
#'   lambda_t, lambda_w, n_alloc, min_rss_w}; and \code{p}, a matrix with one column per \code{slack}.
#' @examples
#' x <- seq(0, 20, length.out = 50)
#' f <- function(mu, h) 1 + 300 * k4_cdf(x, mu, 2, -0.4, h) + sin(4321 * mu * x)
#' k4_shared_perm(rep(list(x), 4), list(f(8, 0.3), f(8.4, 0.3), f(6, -0.4), f(6.3, -0.4)),
#'                group = c(1, 1, 2, 2), bounds = c(-0.9, 0.5, -1, 1))
#' @export
k4_shared_perm <- function(x, y, group, problem = NULL, rss_free = NULL, bounds = c(-0.98, 0.95, -3, 4),
                           drift = TRUE, slack = c(1e-9, 0.01), dense = FALSE, maxalloc = 2e6) {
  d <- .arc_runs(x, y)
  pr <- if (is.null(problem)) rep(1L, d$R) else match(problem, unique(problem))
  stopifnot(length(group) == d$R, length(pr) == d$R, !is.unsorted(pr), all(slack >= 0), length(slack) >= 1L,
            length(bounds) == 4L, bounds[1] < bounds[2], bounds[3] < bounds[4])
  Q <- max(pr); g <- integer(d$R)
  for (q in seq_len(Q)) g[pr == q] <- match(group[pr == q], unique(group[pr == q])) - 1L
  stopifnot(is.null(rss_free) || length(rss_free) == Q, all(g < 32L))
  .arc_sync_threads()
  r <- .C(C_arcsh_perm, d$X, d$Y, d$off, Q, as.integer(c(0, cumsum(tabulate(pr, Q)))), g, as.double(bounds),
          as.integer(drift), as.integer(dense), as.double(slack), length(slack), as.integer(maxalloc),
          out = double(5 * Q), p = double(Q * length(slack)), NAOK = TRUE)
  o <- matrix(r$out, Q, 5L, byrow = TRUE)
  nm <- if (is.null(problem)) NULL else as.character(unique(problem))
  list(stat = matrix(c(o[, 1:3], if (is.null(rss_free)) rep(NA_real_, Q) else o[, 2] / rss_free, o[, 4:5]), Q, 6L,
         dimnames = list(nm, c("rss_tr", "rss_w", "lambda_t", "lambda_w", "n_alloc", "min_rss_w"))),
       p = matrix(r$p, Q, length(slack), byrow = TRUE, dimnames = list(nm, paste0("slack_", slack))))
}
