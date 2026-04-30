# Bit-identity tests for `mf_post_smooth(method = "TI")` against
# the upstream univariate_TI_regression. Tolerance = machine
# precision on the orthonormal wavelet path.

skip_if_no_fsusier_for_TI <- function() {
  testthat::skip_if_not_installed("fsusieR")
  testthat::skip_if_not_installed("ashr")
}

test_that("mf_post_smooth(method = 'TI') matches univariate_TI_regression bit-for-bit on a single-effect fit", {
  skip_if_no_fsusier_for_TI()
  set.seed(11)
  n <- 60; p <- 20; T_m <- 64L
  X    <- matrix(rnorm(n * p), n)
  beta <- numeric(p); beta[3] <- 1.2
  shape <- exp(-((seq_len(T_m) - T_m / 2)^2) / (2 * 6^2))
  Y    <- X %*% (matrix(beta, p, 1) %*% matrix(shape, 1, T_m)) +
            matrix(rnorm(n * T_m, sd = 0.3), n)

  # mfsusieR side: fit + smooth. Production mfsusieR uses
  # `X_eff = X %*% alpha` (alpha-weighted aggregate); upstream
  # `univariate_TI_regression` uses the lead variant. Swap
  # `fit$X_eff` to the lead column so the bit-identity comparison
  # against the upstream routine still holds. Test hack only;
  # production code keeps the alpha-weighted form.
  fit <- fsusie(wavelet_qnorm = FALSE, wavelet_standardize = FALSE, Y, X, L = 1, max_iter = 30, verbose = FALSE)
  L <- nrow(fit$alpha)
  fit$X_eff <- lapply(seq_len(L), function(l) {
    lead_l <- which.max(fit$alpha[l, ])
    X[, lead_l]
  })
  fit_s <- mf_post_smooth(fit, method = "TI",
                          wavelet_filter = 1L,
                          wavelet_family = "DaubExPhase")

  # Reconstruct the per-effect "isolated" position-space response
  # the same way the smoother does, then run the upstream
  # univariate routine on it for the lead variable.
  l <- 1L; m <- 1L
  Y_pos  <- fit$Y_grid[[m]]
  pos_m  <- fit$dwt_meta$pos[[m]]
  if (ncol(Y_pos) > length(pos_m)) {
    Y_pos <- Y_pos[, seq_along(pos_m), drop = FALSE]
  }
  x_lead <- fit$X_eff[[l]]

  ref <- fsusieR:::univariate_TI_regression(
    Y = Y_pos, X = matrix(x_lead, ncol = 1L),
    filter.number = 1, family = "DaubExPhase", alpha = 0.05)

  payload <- fit_s$smoothed$TI
  expect_equal(payload$effect_curves[[m]][[l]],
               as.numeric(ref$effect_estimate),
               tolerance = 0,
               info = "TI effect curve bit-identical")

  # Upstream cred_band rows are c("up", "low"); ours is
  # c(lower, upper) by column.
  expect_equal(payload$credible_bands[[m]][[l]][, 2L],
               as.numeric(ref$cred_band["up", ]),
               tolerance = 0,
               info = "TI upper band bit-identical")
  expect_equal(payload$credible_bands[[m]][[l]][, 1L],
               as.numeric(ref$cred_band["low", ]),
               tolerance = 0,
               info = "TI lower band bit-identical")
})
