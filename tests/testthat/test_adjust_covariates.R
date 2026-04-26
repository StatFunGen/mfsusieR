# Unit tests for `mf_adjust_for_covariates()` and the
# helper functions in `R/adjust_covariates.R`.
#
# OLS path: closed-form sanity (orthogonality, idempotency).
# Wavelet-EB path: math-correct posterior shrinkage. The output
# is NOT bit-identical to fsusieR::EBmvFR(adjust = TRUE) because
# we corrected three documented upstream bugs (slot-name-vs-
# content for MLE_wc2 and the dimensionally-wrong get_ER2.EBmvFR
# formula); see `inst/notes/refactor-exceptions.md` for the full
# entries.

# --- test fixtures -----------------------------------------------

build_cov_sim <- function(seed = 1L, n = 100L, T_m = 64L, K = 3L) {
  set.seed(seed)
  Z    <- matrix(rnorm(n * K), n, K)
  beta <- matrix(rnorm(K * T_m), K, T_m)
  Y    <- Z %*% beta + matrix(rnorm(n * T_m, sd = 0.5), n, T_m)
  list(Y = Y, Z = Z, beta = beta, n = n, T_m = T_m, K = K)
}

# --- OLS path ----------------------------------------------------

test_that("mf_residualize_ols residuals are orthogonal to Z", {
  sim <- build_cov_sim()
  out <- mf_residualize_ols(sim$Y, sim$Z)
  ortho <- max(abs(crossprod(sim$Z, out$Y_adjusted)))
  expect_lt(ortho, 1e-10)
})

test_that("mf_residualize_ols is idempotent", {
  sim  <- build_cov_sim()
  out1 <- mf_residualize_ols(sim$Y, sim$Z)
  out2 <- mf_residualize_ols(out1$Y_adjusted, sim$Z)
  expect_equal(out2$Y_adjusted, out1$Y_adjusted, tolerance = 1e-10)
})

test_that("mf_residualize_ols X_adjusted is orthogonal to Z", {
  sim <- build_cov_sim()
  X   <- matrix(rnorm(sim$n * 5L), sim$n, 5L)
  out <- mf_residualize_ols(sim$Y, sim$Z, X = X)
  expect_lt(max(abs(crossprod(sim$Z, out$X_adjusted))), 1e-10)
})

test_that("dispatcher routes method = 'ols' to mf_residualize_ols", {
  sim <- build_cov_sim()
  a <- mf_adjust_for_covariates(sim$Y, sim$Z, method = "ols")
  b <- mf_residualize_ols(sim$Y, sim$Z)
  expect_equal(a$Y_adjusted, b$Y_adjusted, tolerance = 0)
  expect_equal(a$method, "ols")
})

# --- Wavelet-EB path ---------------------------------------------

test_that("mf_residualize_wavelet_eb returns sane structure", {
  sim <- build_cov_sim()
  out <- suppressWarnings(
    mf_adjust_for_covariates(sim$Y, sim$Z, method = "wavelet_eb"))
  expect_named(out, c("Y_adjusted", "X_adjusted", "fitted_func",
                      "sigma2", "niter", "converged", "method"))
  expect_equal(dim(out$Y_adjusted),  c(sim$n, sim$T_m))
  expect_equal(dim(out$fitted_func), c(sim$K, sim$T_m))
  expect_true(out$sigma2 > 0)              # math-correct ER2 is positive
  expect_true(out$niter <= 100L)
  expect_equal(out$method, "wavelet_eb")
})

test_that("wavelet_eb reduces to small residuals when Z explains Y", {
  # Construct Y as exactly Z*beta (no noise) and adjust away.
  sim <- build_cov_sim(seed = 5L)
  Y_clean <- sim$Z %*% sim$beta
  out <- suppressWarnings(
    mf_adjust_for_covariates(Y_clean, sim$Z, method = "wavelet_eb"))
  # Most of the variance is removed; not exact because wavelet-EB
  # adds shrinkage. Compare to OLS as an upper bound on RSS.
  ols <- mf_residualize_ols(Y_clean, sim$Z)
  expect_lt(sd(out$Y_adjusted), 5 * sd(ols$Y_adjusted) + 1e-3)
})

test_that("wavelet_eb honours the X argument with FWL residualization", {
  sim <- build_cov_sim()
  X   <- matrix(rnorm(sim$n * 8L), sim$n, 8L)
  out <- suppressWarnings(
    mf_adjust_for_covariates(sim$Y, sim$Z, X = X, method = "wavelet_eb"))
  expect_equal(dim(out$X_adjusted), c(sim$n, 8L))
  expect_lt(max(abs(crossprod(sim$Z, out$X_adjusted))), 1e-10)
})

test_that("wavelet_eb agrees with port source up to documented bug magnitude", {
  skip_if_not_installed("fsusieR")
  for (seed in c(1L, 42L)) {
    sim <- build_cov_sim(seed = seed)
    ours <- suppressWarnings(
      mf_adjust_for_covariates(sim$Y, sim$Z, method = "wavelet_eb"))
    ref  <- suppressWarnings(
      fsusieR::EBmvFR(sim$Y, X = sim$Z, adjust = TRUE, verbose = FALSE))
    # The two outputs differ by the magnitude of the documented
    # upstream bugs (Pattern A). Empirically the difference is
    # < 0.1 in absolute terms across the test sweep; tighten if
    # the algorithm is tightened.
    expect_lt(max(abs(ours$Y_adjusted  - ref$Y_adjusted)),  0.1)
    expect_lt(max(abs(ours$fitted_func - ref$fitted_func)), 0.05)
  }
})

# --- Error paths -------------------------------------------------

test_that("unsupported preprocessing flags error cleanly", {
  sim <- build_cov_sim()
  expect_error(
    mf_adjust_for_covariates(sim$Y, sim$Z, thresh_lowcount = 0.01),
    "not supported in v1")
  expect_error(
    mf_adjust_for_covariates(sim$Y, sim$Z, quantile_trans = TRUE),
    "not supported in v1")
})

test_that("non-power-of-two T is rejected by wavelet_eb", {
  sim <- build_cov_sim(T_m = 50L)
  expect_error(
    mf_adjust_for_covariates(sim$Y, sim$Z, method = "wavelet_eb"),
    "power of two")
})

test_that("shape mismatches error cleanly", {
  expect_error(
    mf_adjust_for_covariates(matrix(1, 10L, 8L),
                             matrix(1, 12L, 3L)),
    "same number of rows")
  expect_error(
    mf_adjust_for_covariates(matrix(1, 10L, 8L),
                             matrix(1, 10L, 3L),
                             X = matrix(1, 11L, 5L)),
    "same number of rows")
})
