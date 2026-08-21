# The foreign-function interface is unchecked by R: adding an argument to a C
# entry point touches the C signature, the extern prototype, the registration
# table and every front's call, and .C() complains about none of it. A mismatch
# reads garbage for the missing argument and writes later arguments through the
# wrong pointer, which segfaults the session. These tests pin the arity of every
# registered routine so that a change to the C cannot pass silently.

test_that("every registered .C routine has the expected number of arguments", {
  reg <- getDLLRegisteredRoutines("ellstat")$.C
  expect_gt(length(reg), 0)
  arity <- vapply(reg, function(z) z$numParameters, 0L)
  expected <- c(
    C_ell_ke = 4L, C_ell_nome = 3L, C_ell_jac = 6L, C_ell_basis = 6L,
    C_ell_smom = 5L, C_ell_transfer = 4L, C_ell_dn2mom = 3L,
    C_ell_evm_d = 7L, C_ell_evm_r = 6L, C_ell_moment_var = 5L,
    C_ell_popmom_grid = 5L, C_ell_evm_pop = 6L, C_ell_evm_fit = 6L,
    C_ell_exact_unif = 6L, C_ell_exact_ci = 11L, C_ell_effstudy = 7L,
    C_ell_pj_dens = 7L, C_ell_pj_rand = 6L, C_ell_pj_effcurve = 9L,
    C_ell_pj_kopt = 7L, C_ell_projfam_dens = 6L, C_ell_projfam_rand = 6L,
    C_ell_projfam_eff = 8L, C_ell_projfam_kopt = 6L, C_ell_cat_eff = 9L,
    C_ell_cat_kopt = 7L, C_ell_joint_cov = 6L, C_ell_joint_samp = 6L,
    C_ell_joint_unif = 8L)
  common <- intersect(names(expected), names(arity))
  expect_setequal(names(expected), names(arity))
  expect_equal(arity[common], expected[common])
})

test_that("every exported function runs and returns the documented shape", {
  # a crash here is the symptom of an arity mismatch; a wrong shape is the
  # symptom of a front that was edited without its back end
  th <- rellvm(30L, mu = 1, kappa = 2, k = 0.9, seed = 1)
  expect_length(th, 30L)
  expect_true(all(th >= 0 & th < 2 * pi))

  expect_named(ell_KE(0.5), c("K", "E"))
  expect_length(ell_nome(c(0.2, 0.8)), 2L)
  expect_named(ell_jacobi(c(0, 1), 0.5), c("sn", "cn", "dn"))
  expect_equal(dim(ell_basis(seq(0, 6, length.out = 11), 2L, 0.5)$C), c(2L, 11L))
  expect_equal(dim(ell_moments(th, 2L, 0.9)), c(2L, 2L))
  expect_equal(dim(ell_moment_var(th, 2L, 0.9)), c(2L, 3L))
  expect_equal(dim(ell_transfer(0.9, 2L, 4L)), c(2L, 4L))
  expect_length(ell_dn2_moments(0.7, 3L), 4L)
  expect_length(dellvm(c(0, 1, 2), 0, 1, 0.9, 128L), 3L)
  expect_equal(dim(ell_evm_moments(0, 1, 0.9, 2L, 128L)), c(2L, 2L))
  expect_named(fit_ellvm(th, 0.9, 1L, 128L), c("mu", "kappa"))
  expect_named(ell_exact_uniform(th, 0.9, 49L, 1L), c("statistic", "p.value"))
  expect_named(ell_joint(rep(1, 64), 1L, 0.9), c("mean", "Sigma"))
  expect_named(ell_joint_uniform(th, 1L, 0.9, 49L, 256L, 1L),
               c("statistic", "p.value"))
  expect_length(dprojell(c(0, 1), 0.4, 0.5, 2, 64L), 2L)
  expect_length(dprojfam(c(0, 1), 1, 1, "normal"), 2L)
  expect_length(dprojfam(c(0, 1), 1, 1, "cauchy"), 2L)
  expect_length(rprojfam(5L, 1, 1, "normal", 1L), 5L)
  expect_named(ell_kopt(0.4, 0.6, 2, 256L, 64L), c("k", "efficiency"))
  expect_true(is.list(ell_eff_curve(0.4, 0.6, 2, c(0, 0.5), 256L, 64L)))
  expect_true(is.list(ell_catalogue("vonmises", 2, 1, 1L, c(0, 0.5), 512L)))
  ci <- ell_exact_ci(th, 0.9, seq(0.5, 5, length.out = 8L), 19L, 0.1, 1L, 128L, 4L)
  expect_named(ci, c("lower", "upper", "statistic", "curve"))
})
