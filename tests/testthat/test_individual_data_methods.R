# Apple-to-apple comparison against:
#   mvf.susie.alpha/R/operation_on_multfsusie_obj.R#L89-L132 (cal_partial_resid.multfsusie)
#   mvf.susie.alpha/R/computational_routine.R#L21-L83 (cal_Bhat_Shat_multfsusie)
#   susieR/R/individual_data_methods.R#L104-L145 (compute_residuals/compute_ser_statistics.individual; T=1 paradigm)
#
# 5a contracts:
#   1. compute_residuals.mf_individual returns model with
#      per-modality residuals[[m]] = X^T (D[[m]] - fitted_without_l[[m]]).
#   2. compute_ser_statistics.mf_individual returns betahat[[m]],
#      shat2[[m]] per modality.
#   3. update_fitted_values.mf_individual reconstructs
#      model$fitted[[m]] from fitted_without_l[[m]] plus the
#      effect-l contribution.

# ---- Helpers -----------------------------------------------------------

make_data <- function(n = 30, J = 8, T_per_modality = c(64L, 128L)) {
  set.seed(mfsusier_test_seed())
  X <- matrix(rnorm(n * J), nrow = n)
  Y <- lapply(T_per_modality, function(T_m) matrix(rnorm(n * T_m), nrow = n))
  mfsusieR:::create_mf_individual(X = X, Y = Y, verbose = FALSE)
}

make_model <- function(data, L = 3, sigma2_scalar = 1) {
  params <- list(L = L, prior_weights = NULL, prior = NULL,
                 residual_variance = lapply(seq_len(data$M),
                                            function(m) sigma2_scalar))
  var_y <- lapply(seq_len(data$M),
                  function(m) rep(sigma2_scalar, data$T_padded[m]))
  mfsusieR:::initialize_susie_model.mf_individual(data, params, var_y)
}

# ---- predictor_weights cache -------------------------------------------

test_that("data$predictor_weights = colSums(X^2) cached on mf_individual", {
  data <- make_data()
  expect_equal(data$predictor_weights, colSums(data$X^2), tolerance = 1e-12)
})

# ---- compute_residuals.mf_individual ----------------------------------

test_that("compute_residuals returns model with per-modality residuals XtR", {
  data <- make_data()
  model <- make_model(data, L = 3)

  out <- mfsusieR:::compute_residuals.mf_individual(data, NULL, model, l = 1)

  expect_identical(length(out$residuals), data$M)
  for (m in seq_len(data$M)) {
    expect_identical(dim(out$residuals[[m]]),
                     c(data$J, data$T_padded[m]))
  }
})

test_that("compute_residuals XtR = X' (D - fitted_without_l) per modality", {
  data <- make_data()
  model <- make_model(data, L = 3)

  out <- mfsusieR:::compute_residuals.mf_individual(data, NULL, model, l = 1)

  for (m in seq_len(data$M)) {
    expected_R <- data$D[[m]] - out$fitted_without_l[[m]]
    expect_equal(out$raw_residuals[[m]], expected_R, tolerance = 1e-12)
    expected_XtR <- crossprod(data$X, expected_R)
    expect_equal(out$residuals[[m]], expected_XtR, tolerance = 1e-12)
  }
})

test_that("compute_residuals at zero-init effect l gives raw_residual = D[[m]]", {
  # When alpha[l, ] * mu[[l]][[m]] is all zero (initial state),
  # fitted_without_l[[m]] == fitted[[m]] == 0, so raw_residuals = D.
  data <- make_data()
  model <- make_model(data, L = 3)

  out <- mfsusieR:::compute_residuals.mf_individual(data, NULL, model, l = 1)

  for (m in seq_len(data$M)) {
    expect_equal(out$raw_residuals[[m]], data$D[[m]], tolerance = 1e-12)
    expect_equal(out$residuals[[m]],
                 crossprod(data$X, data$D[[m]]), tolerance = 1e-12)
  }
})

# ---- compute_ser_statistics.mf_individual -----------------------------

test_that("compute_ser_statistics returns per-modality betahat, shat2", {
  data <- make_data()
  model <- make_model(data, L = 3, sigma2_scalar = 1)
  model <- mfsusieR:::compute_residuals.mf_individual(data, NULL, model, 1)

  out <- mfsusieR:::compute_ser_statistics.mf_individual(data, NULL, model, 1)

  expect_named(out, c("betahat", "shat2", "optim_init", "optim_bounds",
                      "optim_scale"))
  expect_identical(length(out$betahat), data$M)
  expect_identical(length(out$shat2), data$M)
  for (m in seq_len(data$M)) {
    expect_identical(dim(out$betahat[[m]]),
                     c(data$J, data$T_padded[m]))
    expect_identical(dim(out$shat2[[m]]),
                     c(data$J, data$T_padded[m]))
  }
})

test_that("compute_ser_statistics betahat = XtR / predictor_weights", {
  data <- make_data()
  model <- make_model(data, L = 3, sigma2_scalar = 1)
  model <- mfsusieR:::compute_residuals.mf_individual(data, NULL, model, 1)

  out <- mfsusieR:::compute_ser_statistics.mf_individual(data, NULL, model, 1)

  for (m in seq_len(data$M)) {
    expect_equal(out$betahat[[m]],
                 model$residuals[[m]] / data$predictor_weights,
                 tolerance = 1e-12)
  }
})

test_that("compute_ser_statistics shat2 broadcasts scalar sigma2 per modality", {
  data <- make_data()
  sigma2 <- 0.7
  model <- make_model(data, L = 3, sigma2_scalar = sigma2)
  model <- mfsusieR:::compute_residuals.mf_individual(data, NULL, model, 1)

  out <- mfsusieR:::compute_ser_statistics.mf_individual(data, NULL, model, 1)

  for (m in seq_len(data$M)) {
    expected <- matrix(sigma2 / data$predictor_weights,
                       nrow = data$J, ncol = data$T_padded[m])
    expect_equal(out$shat2[[m]], expected, tolerance = 1e-12)
  }
})

test_that("compute_ser_statistics shat2 broadcasts per-scale sigma2", {
  data <- make_data(T_per_modality = 64L)
  S_m <- length(data$scale_index[[1]])
  sigma2_per_scale <- seq_len(S_m) * 0.1     # distinct per scale
  model <- make_model(data, L = 3)
  model$sigma2[[1]] <- sigma2_per_scale
  model <- mfsusieR:::compute_residuals.mf_individual(data, NULL, model, 1)

  out <- mfsusieR:::compute_ser_statistics.mf_individual(data, NULL, model, 1)

  # Each position gets the sigma2 of its scale per scale_index.
  expected_sigma2_per_pos <- numeric(ncol(data$D[[1]]))
  for (s in seq_along(data$scale_index[[1]])) {
    expected_sigma2_per_pos[data$scale_index[[1]][[s]]] <- sigma2_per_scale[s]
  }
  expected <- outer(1 / data$predictor_weights, expected_sigma2_per_pos)
  expect_equal(out$shat2[[1]], expected, tolerance = 1e-12)
})

# ---- update_fitted_values.mf_individual -------------------------------

test_that("update_fitted_values reconstructs fitted = fitted_without_l + Xb_l", {
  data <- make_data()
  model <- make_model(data, L = 3)
  model <- mfsusieR:::compute_residuals.mf_individual(data, NULL, model, 1)

  # Inject deterministic alpha and mu so we can verify the algebra.
  model$alpha[1, ] <- 1 / data$J
  for (m in seq_len(data$M)) {
    model$mu[[1]][[m]] <- matrix(0.5, nrow = data$J, ncol = data$T_padded[m])
  }

  out <- mfsusieR:::update_fitted_values.mf_individual(data, NULL, model, 1)

  for (m in seq_len(data$M)) {
    b_l_m <- (1 / data$J) * matrix(0.5, nrow = data$J, ncol = data$T_padded[m])
    Xb_l_m <- data$X %*% b_l_m
    expected <- out$fitted_without_l[[m]] + Xb_l_m
    expect_equal(out$fitted[[m]], expected, tolerance = 1e-12)
  }
})

# ---- Round-trip: residuals -> SER -> update_fitted_values -------------

test_that("zero-effect IBSS sweep keeps fitted at zero", {
  data <- make_data()
  model <- make_model(data, L = 3)

  model <- mfsusieR:::compute_residuals.mf_individual(data, NULL, model, 1)
  out <- mfsusieR:::update_fitted_values.mf_individual(data, NULL, model, 1)

  # alpha_l . mu_l == 0 at init, so update keeps fitted == 0 within ULP.
  for (m in seq_len(data$M)) {
    expect_equal(out$fitted[[m]],
                 matrix(0, nrow = data$n, ncol = data$T_padded[m]),
                 tolerance = 1e-12)
  }
})

# ---- C3 fidelity vs mvf.susie.alpha::cal_Bhat_Shat_multfsusie ---------

test_that("compute_ser_statistics betahat matches mvf.susie.alpha at 1e-12 (initial residual)", {
  skip_if_no_mvf_alpha()
  data <- make_data(n = 50, J = 10, T_per_modality = 64L)
  model <- make_model(data, L = 3, sigma2_scalar = 1)

  # Initial residual (l = 1, all-zero effects): R_m = D[[m]] entirely.
  model <- mfsusieR:::compute_residuals.mf_individual(data, NULL, model, 1)
  ours <- mfsusieR:::compute_ser_statistics.mf_individual(data, NULL, model, 1)

  # Reference: per modality, fsusieR::cal_Bhat_Shat on (D[[m]], X)
  # gives the same betahat (since cal_Bhat_Shat does X^T Y / d).
  ref <- fsusieR:::cal_Bhat_Shat(Y = data$D[[1]], X = data$X,
                                  v1 = rep(1, data$n), lowc_wc = NULL)
  expect_equal(ours$betahat[[1]], ref$Bhat, tolerance = 1e-12)
})
