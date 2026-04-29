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

  fit  <- fsusie(wavelet_qnorm = FALSE, Y, X, L = 1, max_iter = 30, verbose = FALSE)
  fit_s <- mf_post_smooth(fit, method = "scalewise")
  payload <- fit_s$smoothed$scalewise

  expect_length(payload$effect_curves, 1L)        # M = 1
  expect_length(payload$effect_curves[[1L]], 1L)  # L = 1
  expect_length(payload$effect_curves[[1L]][[1L]], T_m)
  expect_equal(dim(payload$credible_bands[[1L]][[1L]]), c(T_m, 2L))
  # All smoothers populate lfsr_curves under the unified API.
  expect_length(payload$lfsr_curves[[1L]][[1L]], T_m)
  expect_true(all(payload$lfsr_curves[[1L]][[1L]] >= 0 &
                  payload$lfsr_curves[[1L]][[1L]] <= 1))

  # Credible band envelopes the effect curve.
  band  <- payload$credible_bands[[1L]][[1L]]
  curve <- payload$effect_curves[[1L]][[1L]]
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

  fit  <- fsusie(wavelet_qnorm = FALSE, Y, X, L = 1, max_iter = 30, verbose = FALSE)
  fit_s <- mf_post_smooth(fit, method = "scalewise")
  payload <- fit_s$smoothed$scalewise

  expect_length(payload$effect_curves[[1L]][[1L]], T_m)
  expect_equal(dim(payload$credible_bands[[1L]][[1L]]), c(T_m, 2L))
  band <- payload$credible_bands[[1L]][[1L]]
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
  fit  <- fsusie(wavelet_qnorm = FALSE, Y, X, L = 1, max_iter = 20, verbose = FALSE)
  c_low  <- mf_post_smooth(fit, method = "scalewise",
                           threshold_factor = 0.5)$smoothed$scalewise$effect_curves[[1L]][[1L]]
  c_high <- mf_post_smooth(fit, method = "scalewise",
                           threshold_factor = 2.0)$smoothed$scalewise$effect_curves[[1L]][[1L]]
  # Higher threshold shrinks at least as much in L1.
  expect_lte(sum(abs(c_high)), sum(abs(c_low)) + 1e-10)
})

# Closed-form variance-curve check on a small dyadic input.
# Verifies that `mf_invert_variance_curve` matches the squared-
# inverse-DWT-matrix reference at machine precision. This is
# the corrected formula now used by `.post_smooth_scalewise`.

test_that("mf_invert_variance_curve matches W^2 %*% var_w at machine precision", {
  set.seed(31L)
  T_basis <- 8L
  filt     <- 1L
  fam      <- "DaubExPhase"

  # Reference: build W^T explicitly column-by-column, square each
  # entry, multiply by `var_w`. This is the textbook variance
  # propagation for an orthonormal linear operator with diagonal
  # input covariance.
  W_T <- matrix(0, T_basis, T_basis)
  for (k in seq_len(T_basis)) {
    e_k <- numeric(T_basis); e_k[k] <- 1
    inv_k <- mfsusieR:::mf_invert_dwt(matrix(e_k, nrow = 1L),
                                      column_center = rep(0, T_basis),
                                      column_scale  = rep(1, T_basis),
                                      filter_number = filt,
                                      family        = fam)
    W_T[, k] <- as.numeric(inv_k)
  }
  var_w  <- runif(T_basis, 0.1, 2.0)
  ref    <- as.numeric((W_T^2) %*% var_w)
  ours   <- mfsusieR:::mf_invert_variance_curve(
              var_w, T_basis = T_basis,
              filter_number = filt, family = fam)
  expect_equal(ours, ref, tolerance = 1e-14)
})

test_that("mf_invert_variance_curve agrees with closed-form for a dirac input", {
  # If var_w is a dirac at index k, the position-space variance
  # is (column k of W^T)^2.
  T_basis <- 8L
  filt     <- 1L
  fam      <- "DaubExPhase"
  for (k in seq_len(T_basis)) {
    var_w     <- numeric(T_basis); var_w[k] <- 1
    e_k       <- numeric(T_basis); e_k[k]   <- 1
    inv_k     <- mfsusieR:::mf_invert_dwt(matrix(e_k, nrow = 1L),
                                          column_center = rep(0, T_basis),
                                          column_scale  = rep(1, T_basis),
                                          filter_number = filt,
                                          family        = fam)
    expected  <- as.numeric(inv_k)^2
    got       <- mfsusieR:::mf_invert_variance_curve(
                   var_w, T_basis = T_basis,
                   filter_number = filt, family = fam)
    expect_equal(got, expected, tolerance = 1e-14)
  }
})
