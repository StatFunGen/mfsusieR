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

  # mfsusieR side: fit + smooth.
  fit  <- fsusie(Y, X, L = 1, max_iter = 30, verbose = FALSE)
  fit_s <- mf_post_smooth(fit, method = "TI",
                          wavelet_filter = 1L,
                          wavelet_family = "DaubExPhase")

  # Reconstruct the per-effect "isolated" position-space response
  # the same way the smoother does, then run the upstream
  # univariate routine on it for the lead variable.
  l <- 1L; m <- 1L
  meta <- fit$dwt_meta
  D_w  <- fit$residuals[[m]] + fit$fitted[[m]]
  Y_pos <- mfsusieR:::mf_invert_dwt(
    D_packed      = D_w,
    column_center = meta$column_center[[m]],
    column_scale  = meta$column_scale[[m]],
    filter_number = meta$wavelet_filter,
    family        = meta$wavelet_family
  )
  if (ncol(Y_pos) > length(meta$pos[[m]])) {
    Y_pos <- Y_pos[, seq_along(meta$pos[[m]]), drop = FALSE]
  }
  x_lead <- fit$lead_X[[l]]

  ref <- fsusieR:::univariate_TI_regression(
    Y = Y_pos, X = matrix(x_lead, ncol = 1L),
    filter.number = 1, family = "DaubExPhase", alpha = 0.05)

  expect_equal(fit_s$effect_curves[[m]][[l]],
               as.numeric(ref$effect_estimate),
               tolerance = 0,
               info = "TI effect curve bit-identical")

  # Upstream cred_band rows are c("up", "low"); ours is
  # c(lower, upper) by column.
  expect_equal(fit_s$credible_bands[[m]][[l]][, 2L],
               as.numeric(ref$cred_band["up", ]),
               tolerance = 0,
               info = "TI upper band bit-identical")
  expect_equal(fit_s$credible_bands[[m]][[l]][, 1L],
               as.numeric(ref$cred_band["low", ]),
               tolerance = 0,
               info = "TI lower band bit-identical")
})
