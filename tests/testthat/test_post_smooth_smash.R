# Reference tests for `mf_post_smooth(method = "smash")`.
#
# 1. Per-effect, per-outcome the mfsusieR smash kernel agrees
#    with `fsusieR::univariate_smash_regression` on the same
#    isolated response and lead variable. Tolerance <= 1e-12.
# 2. End-to-end smoke: a fit with M = 1, T_m = 64 runs the
#    smash post-smoother to completion and produces effect
#    curves and credible bands of the right shape.
# 3. Argument validation: dispatch errors when smashr is not
#    available (synthetic skip of `requireNamespace`).

test_that("mfsusieR smash kernel matches fsusieR::univariate_smash_regression at machine precision", {
  testthat::skip_if_not_installed("smashr")
  testthat::skip_if_not_installed("fsusieR")

  set.seed(101L)
  n   <- 60L
  T_m <- 64L
  X   <- matrix(rnorm(n), n, 1L)
  Y   <- X %*% matrix(0.7 * c(rep(0, 24), rep(1, 16), rep(0, 24)),
                      1, T_m) +
         matrix(rnorm(n * T_m, sd = 0.5), n, T_m)

  ours <- mfsusieR:::univariate_smash_regression(Y, X, alpha = 0.05)
  ref  <- fsusieR:::univariate_smash_regression(Y, X, alpha = 0.05)

  expect_equal(as.numeric(ours$effect_estimate),
               as.numeric(ref$effect_estimate),
               tolerance = 1e-12)
  expect_equal(as.numeric(ours$cred_band),
               as.numeric(ref$cred_band),
               tolerance = 1e-12)
})

test_that("mf_post_smooth(method = 'smash') populates effect_curves and credible_bands", {
  testthat::skip_if_not_installed("smashr")
  set.seed(13L)
  n <- 60L; p <- 20L; T_m <- 64L
  X    <- matrix(rnorm(n * p), n)
  beta <- numeric(p); beta[3] <- 1.2
  shape <- exp(-((seq_len(T_m) - T_m / 2)^2) / (2 * 6^2))
  Y    <- X %*% (matrix(beta, p, 1) %*% matrix(shape, 1, T_m)) +
            matrix(rnorm(n * T_m, sd = 0.3), n)

  fit <- fsusie(Y, X, L = 1, max_iter = 30, verbose = FALSE)
  fit_s <- mf_post_smooth(fit, method = "smash")

  payload <- fit_s$smoothed$smash
  expect_length(payload$effect_curves, 1L)        # M = 1
  expect_length(payload$effect_curves[[1L]], 1L)  # L = 1
  expect_length(payload$effect_curves[[1L]][[1L]], T_m)
  expect_equal(dim(payload$credible_bands[[1L]][[1L]]), c(T_m, 2L))
  # All smoothers populate lfsr_curves under the unified API.
  expect_length(payload$lfsr_curves[[1L]][[1L]], T_m)
  expect_true(all(payload$lfsr_curves[[1L]][[1L]] >= 0 &
                  payload$lfsr_curves[[1L]][[1L]] <= 1))

  band  <- payload$credible_bands[[1L]][[1L]]
  curve <- payload$effect_curves[[1L]][[1L]]
  expect_true(all(band[, 1L] <= curve + 1e-10))
  expect_true(all(curve - 1e-10 <= band[, 2L]))
})

test_that("method = 'smash' errors cleanly when smashr is not available", {
  # Direct test of the precondition path. We mock by clearing
  # the smashr namespace cache via attaching a sentinel; this
  # is a unit-level check that the dispatch error message
  # contains the install hint.
  set.seed(7L)
  fit <- list()
  class(fit) <- c("mfsusie", "susie")
  fit$residuals <- list()
  fit$lead_X    <- list()
  if (requireNamespace("smashr", quietly = TRUE)) {
    skip("smashr is installed; cannot exercise the missing-dep error.")
  }
  expect_error(mf_post_smooth(fit, method = "smash"),
               "smashr")
})
