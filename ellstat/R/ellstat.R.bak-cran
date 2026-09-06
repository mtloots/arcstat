#' Complete elliptic integrals of the first and second kind
#' @param k modulus, in [0, 1)
#' @return list with components \code{K} and \code{E}
#' @examples
#' ell_KE(c(0, 0.5, 0.9))
#' @export
ell_KE <- function(k) {
  n <- length(k)
  z <- .C("C_ell_ke", as.double(k), as.integer(n),
          K = double(n), E = double(n))
  list(K = z$K, E = z$E)
}

#' Nome q = exp(-pi K'/K)
#' @param k modulus
#' @examples
#' ell_nome(c(0.5, 0.9))
#' @export
ell_nome <- function(k) {
  n <- length(k)
  .C("C_ell_nome", as.double(k), as.integer(n), out = double(n))$out
}

#' Jacobi elliptic functions
#' @param u argument
#' @param k modulus
#' @examples
#' ell_jacobi(seq(0, 4, length.out = 5), k = 0.7)
#' @export
ell_jacobi <- function(u, k) {
  n <- length(u)
  z <- .C("C_ell_jac", as.double(u), as.integer(n), as.double(k),
          sn = double(n), cn = double(n), dn = double(n))
  list(sn = z$sn, cn = z$cn, dn = z$dn)
}

#' The elliptic moment basis c_p, s_p
#' @param t angles in [0, 2*pi)
#' @param P highest order
#' @param k modulus
#' @return list of two P-by-length(t) matrices
#' @examples
#' b <- ell_basis(seq(0, 2 * pi, length.out = 9), P = 2L, k = 0.9)
#' # at modulus zero the basis is the trigonometric one
#' max(abs(ell_basis(seq(0, 2 * pi, length.out = 9), 1L, 0)$C[1, ] -
#'         cos(seq(0, 2 * pi, length.out = 9))))
#' @export
ell_basis <- function(t, P, k) {
  nt <- length(t)
  z <- .C("C_ell_basis", as.double(t), as.integer(nt), as.integer(P),
          as.double(k), C = double(P * nt), S = double(P * nt))
  list(C = matrix(z$C, nrow = P, byrow = TRUE),
       S = matrix(z$S, nrow = P, byrow = TRUE))
}

#' Sample elliptic moments
#' @param t angles
#' @param P highest order
#' @param k modulus
#' @return P-by-2 matrix of (c, s) sample moments
#' @examples
#' theta <- rellvm(200, mu = 1, kappa = 2, k = 0.9, seed = 1)
#' ell_moments(theta, P = 2L, k = 0.9)
#' @export
ell_moments <- function(t, P = 2L, k = 0.9) {
  z <- .C("C_ell_smom", as.double(t), as.integer(length(t)), as.integer(P),
          as.double(k), out = double(2 * P))$out
  matrix(z, ncol = 2, byrow = TRUE, dimnames = list(NULL, c("c", "s")))
}

#' Exact first and second order sample theory
#'
#' Returns the modulus of the sample elliptic moment, its realised trace
#' variance and the exact value (1 - |E|^2)/n. The two agree identically.
#' @param t angles
#' @param P highest order
#' @param k modulus
#' @examples
#' theta <- rellvm(200, mu = 1, kappa = 2, k = 0.9, seed = 1)
#' mv <- ell_moment_var(theta, P = 2L, k = 0.9)
#' # the realised and exact columns agree identically, not approximately
#' max(abs(mv[, "realised"] - mv[, "exact"]))
#' @export
ell_moment_var <- function(t, P = 3L, k = 0.9) {
  z <- .C("C_ell_moment_var", as.double(t), as.integer(length(t)),
          as.integer(P), as.double(k), out = double(3 * P))$out
  matrix(z, ncol = 3, byrow = TRUE,
         dimnames = list(NULL, c("modulus", "realised", "exact")))
}

#' Transfer coefficients from trigonometric to elliptic moments
#' @param k modulus
#' @param P highest elliptic order
#' @param M number of retained harmonics per order
#' @examples
#' ell_transfer(k = 0.9, P = 2L, M = 4L)
#' @export
ell_transfer <- function(k, P = 3L, M = 5L) {
  z <- .C("C_ell_transfer", as.double(k), as.integer(P), as.integer(M),
          out = double(P * M))$out
  matrix(z, nrow = P, byrow = TRUE)
}

#' Fourier coefficients of the arc-length family on an ellipse
#'
#' On its own clock the family is dn^2/(4E); element one is E/K and element
#' element \code{n + 1} is \eqn{(2\pi^2/K^2) n q^n/(1-q^{2n})}.
#' @param k modulus
#' @param P number of harmonics
#' @examples
#' ell_dn2_moments(k = 0.8, P = 4L)
#' @export
ell_dn2_moments <- function(k, P = 4L) {
  .C("C_ell_dn2mom", as.double(k), as.integer(P), out = double(P + 1))$out
}

#' Elliptic von Mises density
#' @param t angles
#' @param mu mean direction
#' @param kappa concentration
#' @param k modulus
#' @param ng quadrature points for the normalising constant
#' @examples
#' theta <- seq(0, 2 * pi, length.out = 9)
#' dellvm(theta, mu = 1, kappa = 2, k = 0.9)
#' @export
dellvm <- function(t, mu = 0, kappa = 1, k = 0.9, ng = 1024L) {
  .C("C_ell_evm_d", as.double(t), as.integer(length(t)), as.double(mu),
     as.double(kappa), as.double(k), as.integer(ng),
     out = double(length(t)))$out
}

#' Simulate from the elliptic von Mises
#' @param n sample size
#' @param mu mean direction
#' @param kappa concentration
#' @param k modulus
#' @param seed stream seed
#' @examples
#' theta <- rellvm(10, mu = 1, kappa = 2, k = 0.9, seed = 42)
#' round(theta, 3)
#' @export
rellvm <- function(n, mu = 0, kappa = 1, k = 0.9, seed = 1L) {
  .C("C_ell_evm_r", as.integer(n), as.double(mu), as.double(kappa),
     as.double(k), as.numeric(seed), out = double(n))$out
}

#' Population elliptic moments of the elliptic von Mises
#' @inheritParams dellvm
#' @param P highest order
#' @examples
#' ell_evm_moments(mu = 1, kappa = 2, k = 0.9, P = 2L)
#' @export
ell_evm_moments <- function(mu, kappa, k, P = 2L, ng = 1024L) {
  z <- .C("C_ell_evm_pop", as.double(mu), as.double(kappa), as.double(k),
          as.integer(P), as.integer(ng), out = double(2 * P))$out
  matrix(z, ncol = 2, byrow = TRUE, dimnames = list(NULL, c("c", "s")))
}

#' Minimum-distance fit of the elliptic von Mises
#' @param t angles
#' @param k modulus
#' @param P number of elliptic moments matched
#' @param ng quadrature points
#' @return named vector of mu and kappa
#' @examples
#' theta <- rellvm(500, mu = 1, kappa = 3, k = 0.9, seed = 7)
#' fit_ellvm(theta, k = 0.9, P = 2L)
#' @export
fit_ellvm <- function(t, k = 0.9, P = 2L, ng = 512L) {
  z <- .C("C_ell_evm_fit", as.double(t), as.integer(length(t)), as.double(k),
          as.integer(P), as.integer(ng), out = double(2))$out
  c(mu = z[1], kappa = z[2])
}

#' Exact test of circular uniformity by rotation invariance
#'
#' The rank p-value is exact for every sample size and every B.
#' @param t angles
#' @param k modulus
#' @param B number of rotations
#' @param seed stream seed
#' @examples
#' theta <- rellvm(60, mu = 0, kappa = 0, k = 0.9, seed = 3)  # uniform
#' ell_exact_uniform(theta, k = 0.9, B = 199L, seed = 1L)
#' @export
ell_exact_uniform <- function(t, k = 0.9, B = 999L, seed = 1L) {
  z <- .C("C_ell_exact_unif", as.double(t), as.integer(length(t)), as.double(k),
          as.integer(B), as.numeric(seed), out = double(2))$out
  list(statistic = z[1], p.value = z[2])
}

#' Exact confidence interval for the concentration by test inversion
#'
#' Inverts the Monte Carlo test of a fixed concentration. The null is composite,
#' the statistic depending on the mean direction as well as the concentration,
#' so the p-value is maximised over a grid of mean directions; provided that
#' grid contains the true direction the interval covers at no less than its
#' nominal level, at the price of conservatism. Endpoints are obtained by
#' interpolating the p-value curve across the level, so the answer does not
#' depend on the resolution of \code{kgrid} to first order.
#'
#' @param t angles
#' @param k modulus
#' @param kgrid grid of concentrations to invert over
#' @param B replicates per grid point
#' @param level nominal level
#' @param seed stream seed
#' @param ng quadrature points
#' @param nmu size of the grid of mean directions over which the p-value is
#'   maximised. The null is composite: the statistic depends on the mean
#'   direction as well as the concentration, so calibrating at a single value
#'   under-covers. Larger values give a more conservative and more expensive
#'   interval.
#' @examples
#' theta <- rellvm(80, mu = 1, kappa = 3, k = 0.9, seed = 11)
#' ell_exact_ci(theta, k = 0.9, kgrid = seq(0.5, 8, length.out = 20L),
#'              B = 99L, level = 0.05, nmu = 6L)[c("lower", "upper")]
#' @export
ell_exact_ci <- function(t, k = 0.9, kgrid = seq(0.05, 10, length.out = 60L),
                         B = 499L, level = 0.05, seed = 1L, ng = 512L,
                         nmu = 8L) {
  ngr <- length(kgrid)
  z <- .C("C_ell_exact_ci", as.double(t), as.integer(length(t)), as.double(k),
          as.double(kgrid), as.integer(ngr), as.integer(B), as.double(level),
          as.numeric(seed), as.integer(ng), as.integer(nmu),
          out = double(3 + ngr))$out
  list(lower = z[1], upper = z[2], statistic = z[3],
       curve = data.frame(kappa = kgrid, p.value = z[4:(3 + ngr)]))
}

#' Matched-modulus efficiency study
#'
#' Races the trigonometric estimator of the concentration against the elliptic
#' one, both by inverting a single monotone moment equation, on data from the
#' elliptic von Mises of modulus \code{k}.
#' @param k modulus
#' @param kappa true concentration
#' @param n sample size
#' @param R replicates
#' @param ng quadrature points
#' @param seed stream seed
#' @return named vector: eta, var.trig, var.elliptic, gain
#' @examples
#' ell_efficiency(k = 0.9, kappa = 3, n = 100L, R = 200L)
#' @export
ell_efficiency <- function(k, kappa = 3, n = 250L, R = 2000L, ng = 1024L,
                           seed = 1L) {
  z <- .C("C_ell_effstudy", as.double(k), as.double(kappa), as.integer(n),
          as.integer(R), as.integer(ng), as.numeric(seed),
          out = double(5))$out
  c(eta = z[1], var.trig = z[2], var.elliptic = z[3], gain = z[2] / z[3],
    solved.trig = z[4], solved.elliptic = z[5])
}

#' Density of the projected family
#'
#' The angular law of \code{arg(Y)} where \code{Y} given \code{U} is uniform on
#' the ellipse with semi-axes \code{U(1, tau)} offset by \code{d} along its
#' axis, and \code{U} has the generalised gamma mixing law of shape
#' \code{beta}. Projecting a uniform law is geometry, so the conditional
#' density is elementary and the mixture is a single quadrature.
#' @param t angles
#' @param tau axis ratio of the ellipse
#' @param d offset of the origin from the centre
#' @param beta shape of the mixing law
#' @param ng quadrature points for the mixture
#' @examples
#' dprojell(seq(0, 2 * pi, length.out = 9), tau = 0.3, d = 0.6, beta = 2)
#' @export
dprojell <- function(t, tau = 0.3, d = 0.6, beta = 2, ng = 400L) {
  .C("C_ell_pj_dens", as.double(t), as.integer(length(t)), as.double(tau),
     as.double(d), as.double(beta), as.integer(ng),
     out = double(length(t)))$out
}

#' Simulate from the projected family
#' @inheritParams dprojell
#' @param n sample size
#' @param seed stream seed
#' @examples
#' rprojell(10, tau = 0.3, d = 0.6, beta = 2, seed = 5)
#' @export
rprojell <- function(n, tau = 0.3, d = 0.6, beta = 2, seed = 1L) {
  .C("C_ell_pj_rand", as.integer(n), as.double(tau), as.double(d),
     as.double(beta), as.numeric(seed), out = double(n))$out
}

#' Efficiency of the elliptic moment estimator against the likelihood
#'
#' Returns the asymptotic efficiency \code{corr^2(c_1(., k), score)} of the
#' moment estimator of the offset, for each modulus on \code{kgrid}. By the
#' characterisation theorem this is exact and requires no simulation. The
#' returned \code{godambe} entry is the residual
#' \code{|dE[g]/dd - Cov(g, score)|}, which checks the theorem itself.
#' @inheritParams dprojell
#' @param kgrid moduli at which to evaluate the efficiency
#' @param nt quadrature points in the angle
#' @param eps finite-difference step for the score
#' @examples
#' e <- ell_eff_curve(tau = 0.3, d = 0.6, beta = 2,
#'                    kgrid = c(0, 0.9, 0.99), nt = 512L, ng = 100L)
#' round(e$efficiency, 3)
#' @export
ell_eff_curve <- function(tau = 0.3, d = 0.6, beta = 2,
                          kgrid = c(0, 0.9, 0.99, 0.999, 0.99999),
                          nt = 4096L, ng = 400L, eps = 1e-4) {
  ngr <- length(kgrid)
  z <- .C("C_ell_pj_effcurve", as.double(tau), as.double(d), as.double(beta),
          as.double(kgrid), as.integer(ngr), as.integer(nt), as.integer(ng),
          as.double(eps), out = double(ngr + 2))$out
  list(efficiency = z[seq_len(ngr)], fisher = z[ngr + 1],
       godambe = z[ngr + 2], kgrid = kgrid)
}

#' The optimal modulus
#'
#' Maximises \code{corr^2(c_1(., k), score)} over the modulus by golden
#' section. The optimum may be interior: pushing the modulus towards one is
#' not always right.
#' @inheritParams ell_eff_curve
#' @return named vector of the optimal modulus and its efficiency
#' @examples
#' ell_kopt(tau = 0.3, d = 0.6, beta = 2, nt = 512L, ng = 100L)
#' @export
ell_kopt <- function(tau = 0.3, d = 0.6, beta = 2, nt = 4096L, ng = 400L,
                     eps = 1e-4) {
  z <- .C("C_ell_pj_kopt", as.double(tau), as.double(d), as.double(beta),
          as.integer(nt), as.integer(ng), as.double(eps), out = double(2))$out
  c(k = z[1], efficiency = z[2])
}

#' Density of a cited projected family
#'
#' @param t angles
#' @param tau axis ratio of the scatter, \code{diag(1, tau^2)}
#' @param d offset of the mean from the origin
#' @param family "normal" for the projected (offset) normal, "cauchy" for the
#'   projected Cauchy
#' @examples
#' theta <- seq(0, 2 * pi, length.out = 9)
#' dprojfam(theta, tau = 1, d = 1, family = "normal")
#' dprojfam(theta, tau = 1, d = 1, family = "cauchy")
#' @export
dprojfam <- function(t, tau = 1, d = 1, family = c("normal", "cauchy")) {
  w <- if (match.arg(family) == "normal") 0L else 1L
  .C("C_ell_projfam_dens", as.double(t), as.integer(length(t)), as.double(tau),
     as.double(d), w, out = double(length(t)))$out
}

#' Simulate from a cited projected family
#' @inheritParams dprojfam
#' @param n sample size
#' @param seed stream seed
#' @examples
#' rprojfam(10, tau = 1, d = 1, family = "normal", seed = 2)
#' @export
rprojfam <- function(n, tau = 1, d = 1, family = c("normal", "cauchy"),
                     seed = 1L) {
  w <- if (match.arg(family) == "normal") 0L else 1L
  .C("C_ell_projfam_rand", as.integer(n), as.double(tau), as.double(d), w,
     as.numeric(seed), out = double(n))$out
}

#' Efficiency and optimal modulus for a cited projected family
#' @inheritParams dprojfam
#' @param kgrid moduli at which to evaluate the efficiency
#' @param nt quadrature points in the angle
#' @param eps finite-difference step for the score
#' @examples
#' ell_projfam_eff(tau = 0.4, d = 1, family = "normal",
#'                 kgrid = c(0, 0.9, 0.99), nt = 1024L)$efficiency
#' @export
ell_projfam_eff <- function(tau = 1, d = 1, family = c("normal", "cauchy"),
                            kgrid = c(0, 0.9, 0.99, 0.999, 0.99999),
                            nt = 4096L, eps = 1e-4) {
  w <- if (match.arg(family) == "normal") 0L else 1L
  ngr <- length(kgrid)
  z <- .C("C_ell_projfam_eff", as.double(tau), as.double(d), w,
          as.double(kgrid), as.integer(ngr), as.integer(nt), as.double(eps),
          out = double(ngr + 2))$out
  list(efficiency = z[seq_len(ngr)], fisher = z[ngr + 1], mass = z[ngr + 2],
       kgrid = kgrid)
}

#' @rdname ell_projfam_eff
#' @export
ell_projfam_kopt <- function(tau = 1, d = 1, family = c("normal", "cauchy"),
                             nt = 4096L, eps = 1e-4) {
  w <- if (match.arg(family) == "normal") 0L else 1L
  z <- .C("C_ell_projfam_kopt", as.double(tau), as.double(d), w,
          as.integer(nt), as.double(eps), out = double(2))$out
  c(k = z[1], efficiency = z[2])
}

#' The catalogue of circular distributions
#'
#' After Hosking's table of L-moments for standard distributions. For each
#' family the concentration parameter is the estimand, and the routine returns
#' the asymptotic efficiency of the first elliptic moment relative to maximum
#' likelihood at each modulus, together with the optimal modulus.
#'
#' @param family one of "vonmises", "wrappednormal", "wrappedcauchy",
#'   "cardioid", "projnormal", "projcauchy", "ellvonmises", "arcellipse"
#' @param p1 the concentration parameter (the estimand)
#' @param p2 auxiliary parameter: the scatter ratio for the projected families,
#'   the modulus for the elliptic von Mises, unused otherwise
#' @param order the elliptic moment order to use. Antipodally symmetric
#'   (axial) families carry no information in the first moment, because
#'   \eqn{c_1(\theta+\pi)=-c_1(\theta)} while their score is \eqn{\pi}-periodic,
#'   so the covariance vanishes identically; use \code{order = 2} for those,
#'   exactly as one uses the second trigonometric moment for axial data.
#' @param kgrid moduli at which to evaluate
#' @param nt quadrature points
#' @param eps finite-difference step for the score
#' @examples
#' # the von Mises law is fully efficient at modulus zero and nowhere else
#' r <- ell_catalogue("vonmises", p1 = 2, kgrid = c(0, 0.9), nt = 1024L)
#' round(r$efficiency, 4)
#' @export
ell_catalogue <- function(family, p1, p2 = 1, order = 1L,
                          kgrid = c(0, 0.9, 0.99, 0.999, 0.99999),
                          nt = 4096L, eps = 1e-4) {
  fams <- c("vonmises", "wrappednormal", "wrappedcauchy", "cardioid",
            "projnormal", "projcauchy", "ellvonmises", "arcellipse")
  f <- match(match.arg(family, fams), fams) - 1L
  ngr <- length(kgrid)
  z <- .C("C_ell_cat_eff", as.integer(f), as.double(p1), as.double(p2),
          as.integer(order), as.double(kgrid), as.integer(ngr), as.integer(nt),
          as.double(eps), out = double(ngr + 2))$out
  k <- .C("C_ell_cat_kopt", as.integer(f), as.double(p1), as.double(p2),
          as.integer(order), as.integer(nt), as.double(eps), out = double(2))$out
  list(efficiency = z[seq_len(ngr)], fisher = z[ngr + 1], mass = z[ngr + 2],
       kopt = k[1], eff_kopt = k[2], kgrid = kgrid)
}

#' Exact joint mean and covariance of the sample elliptic moments
#'
#' Stacking the basis as \eqn{(c_1, s_1, \dots, c_P, s_P)}, the sample
#' elliptic moment vector is a sample mean of a bounded vector, so its mean is
#' exact and its covariance is \eqn{\Sigma/n} for every sample size. This is
#' the exact counterpart of Hosking's asymptotic joint distribution of sample
#' L-moments.
#' @param f density on a uniform grid over \eqn{[0, 2\pi)}
#' @param P highest order
#' @param k modulus
#' @examples
#' theta <- seq(0, 2 * pi, length.out = 512)
#' f <- dellvm(theta, mu = 1, kappa = 3, k = 0.9, ng = 512L)
#' J <- ell_joint(f, P = 2L, k = 0.9)
#' round(J$Sigma, 4)
#' @export
ell_joint <- function(f, P = 2L, k = 0.9) {
  D <- 2L * P
  z <- .C("C_ell_joint_cov", as.double(f), as.integer(length(f)),
          as.integer(P), as.double(k), mean = double(D), Sigma = double(D * D))
  list(mean = z$mean, Sigma = matrix(z$Sigma, D, D, byrow = TRUE))
}

#' @rdname ell_joint
#' @param t angles
#' @examples
#' theta <- rellvm(200, mu = 1, kappa = 3, k = 0.9, seed = 9)
#' ell_joint_sample(theta, P = 2L, k = 0.9)$mean
#' @export
ell_joint_sample <- function(t, P = 2L, k = 0.9) {
  D <- 2L * P
  z <- .C("C_ell_joint_samp", as.double(t), as.integer(length(t)),
          as.integer(P), as.double(k), mean = double(D), Sigma = double(D * D))
  list(mean = z$mean, Sigma = matrix(z$Sigma, D, D, byrow = TRUE))
}

#' Simultaneous exact test of uniformity from the joint moment vector
#'
#' Under uniformity the joint covariance is fully determined, so the quadratic
#' form \eqn{n\hat E'\Sigma_0^{-1}\hat E} has a null distribution free of
#' unknowns. The null being simple, it is simulated directly and the rank
#' p-value is exact for every sample size and every \code{B}.
#' @param t angles
#' @param P number of elliptic moments used jointly
#' @param k modulus
#' @param B replicates
#' @param ng quadrature points for the null covariance
#' @param seed stream seed
#' @examples
#' theta <- rellvm(80, mu = 0, kappa = 0, k = 0.9, seed = 4)
#' ell_joint_uniform(theta, P = 2L, k = 0.9, B = 199L, seed = 1L)
#' @export
ell_joint_uniform <- function(t, P = 2L, k = 0.9, B = 999L, ng = 4096L,
                              seed = 1L) {
  z <- .C("C_ell_joint_unif", as.double(t), as.integer(length(t)),
          as.integer(P), as.double(k), as.integer(B), as.integer(ng),
          as.numeric(seed), out = double(2))$out
  list(statistic = z[1], p.value = z[2])
}
