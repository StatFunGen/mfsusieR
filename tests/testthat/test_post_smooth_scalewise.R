# Smoke + shape tests for `mf_post_smooth(method = "scalewise")`.
# Scalewise soft-thresholds the lead variable's wavelet
# coefficients per scale at `threshold_factor * sqrt(2 log T) * sd`,
# inverts the DWT to position space, and computes pointwise
# credible bands from the per-position posterior sd.
#
# No upstream apple-to-apple fidelity assertion: fsusieR's
# `post_processing = "smash"` uses smashr's empirical-Bayes
# wavelet shrinkage, a different algorithm. Scalewise here is the
# in-package soft-threshold variant; a "smash" port would land
# alongside the smashr / ebnm dependencies (see PR group 6b.1).

test_that("mf_post_smooth(method = 'scalewise') populates effect_curves and credible_bands", {
  testthat::skip_if_not_installed("ashr")
  set.seed(13)
  n <- 60; p <- 20; T_m <- 64L
  X    <- matrix(rnorm(n * p), n)
  beta <- numeric(p); beta[3] <- 1.2
  shape <- exp(-((seq_len(T_m) - T_m / 2)^2) / (2 * 6^2))
  Y    <- X %*% (matrix(beta, p, 1) %*% matrix(shape, 1, T_m)) +
            matrix(rnorm(n * T_m, sd = 0.3), n)

  fit  <- fsusie(Y, X, L = 1, max_iter = 30, verbose = FALSE)
  fit_s <- mf_post_smooth(fit, method = "scalewise")

  expect_length(fit_s$effect_curves, 1L)        # M = 1
  expect_length(fit_s$effect_curves[[1L]], 1L)  # L = 1
  expect_length(fit_s$effect_curves[[1L]][[1L]], T_m)
  expect_equal(dim(fit_s$credible_bands[[1L]][[1L]]), c(T_m, 2L))
  expect_null(fit_s$lfsr_curves)

  # Credible band envelopes the effect curve.
  band <- fit_s$credible_bands[[1L]][[1L]]
  curve <- fit_s$effect_curves[[1L]][[1L]]
  expect_true(all(band[, 1L] <= curve + 1e-10))
  expect_true(all(curve - 1e-10 <= band[, 2L]))

  # Signal direction recovered on the high-signal block (rough,
  # not bit-fidelity): the smoothed curve agrees in sign with the
  # ground-truth Gaussian shape on positions in the [-2 sigma,
  # 2 sigma] window.
  sig_idx <- which(shape > exp(-2))
  expect_true(mean(sign(curve[sig_idx]) == sign(beta[3] * shape[sig_idx])) > 0.7)
})

test_that("mf_post_smooth(method = 'scalewise') gracefully handles T_m = 1 (scalar) outcomes", {
  # Scalar outcome: scalewise reduces to passthrough of the
  # per-effect mean / sd; the band is symmetric around the mean.
  testthat::skip_if_not_installed("ashr")
  set.seed(17)
  n <- 80; p <- 15; T_m <- 1L
  X    <- matrix(rnorm(n * p), n)
  beta <- numeric(p); beta[5] <- 1.5
  Y    <- X %*% matrix(beta, p, T_m) + matrix(rnorm(n * T_m, sd = 0.4), n)

  fit  <- fsusie(Y, X, L = 1, max_iter = 30, verbose = FALSE)
  fit_s <- mf_post_smooth(fit, method = "scalewise")

  expect_length(fit_s$effect_curves[[1L]][[1L]], T_m)
  expect_equal(dim(fit_s$credible_bands[[1L]][[1L]]), c(T_m, 2L))
  band <- fit_s$credible_bands[[1L]][[1L]]
  expect_true(band[1L, 1L] < band[1L, 2L])
})

test_that("mf_post_smooth(method = 'scalewise') honors threshold_factor", {
  # Larger threshold_factor produces a shrunker curve (smaller in
  # absolute value at every position); equality at threshold = 0
  # would be the unshrunk lead-variable mean.
  testthat::skip_if_not_installed("ashr")
  set.seed(19)
  n <- 60; p <- 20; T_m <- 64L
  X    <- matrix(rnorm(n * p), n)
  beta <- numeric(p); beta[2] <- 0.8
  shape <- exp(-((seq_len(T_m) - T_m / 2)^2) / (2 * 8^2))
  Y    <- X %*% (matrix(beta, p, 1) %*% matrix(shape, 1, T_m)) +
            matrix(rnorm(n * T_m, sd = 0.4), n)
  fit  <- fsusie(Y, X, L = 1, max_iter = 20, verbose = FALSE)
  c_low  <- mf_post_smooth(fit, method = "scalewise",
                           threshold_factor = 0.5)$effect_curves[[1L]][[1L]]
  c_high <- mf_post_smooth(fit, method = "scalewise",
                           threshold_factor = 2.0)$effect_curves[[1L]][[1L]]
  # Higher threshold shrinks at least as much in L1.
  expect_lte(sum(abs(c_high)), sum(abs(c_low)) + 1e-10)
})
