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

make_data <- function(n = 30, J = 8, T_per_outcome = c(64L, 128L)) {
  set.seed(mfsusier_test_seed())
  X <- matrix(rnorm(n * J), nrow = n)
  Y <- lapply(T_per_outcome, function(T_m) matrix(rnorm(n * T_m), nrow = n))
  mfsusieR:::create_mf_individual(X = X, Y = Y, verbose = FALSE)
}

make_model <- function(data, L = 3, sigma2_scalar = 1) {
  params <- list(L = L, prior_weights = NULL, prior = NULL,
                 residual_variance = lapply(seq_len(data$M),
                                            function(m) sigma2_scalar))
  var_y <- lapply(seq_len(data$M),
                  function(m) rep(sigma2_scalar, data$T_basis[m]))
  mfsusieR:::initialize_susie_model.mf_individual(data, params, var_y)
}

# ---- xtx_diag cache -------------------------------------------

test_that("data$xtx_diag = colSums(X^2) cached on mf_individual", {
  data <- make_data()
  expect_equal(data$xtx_diag, colSums(data$X^2), tolerance = 1e-12)
})

# ---- compute_residuals.mf_individual ----------------------------------

test_that("compute_residuals returns model with per-modality residuals XtR", {
  data <- make_data()
  model <- make_model(data, L = 3)

  out <- mfsusieR:::compute_residuals.mf_individual(data, NULL, model, l = 1)

  expect_identical(length(out$residuals), data$M)
  for (m in seq_len(data$M)) {
    expect_identical(dim(out$residuals[[m]]),
                     c(data$p, data$T_basis[m]))
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
                     c(data$p, data$T_basis[m]))
    expect_identical(dim(out$shat2[[m]]),
                     c(data$p, data$T_basis[m]))
  }
})

test_that("compute_ser_statistics betahat = XtR / xtx_diag", {
  data <- make_data()
  model <- make_model(data, L = 3, sigma2_scalar = 1)
  model <- mfsusieR:::compute_residuals.mf_individual(data, NULL, model, 1)

  out <- mfsusieR:::compute_ser_statistics.mf_individual(data, NULL, model, 1)

  for (m in seq_len(data$M)) {
    expect_equal(out$betahat[[m]],
                 model$residuals[[m]] / data$xtx_diag,
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
    expected <- matrix(sigma2 / data$xtx_diag,
                       nrow = data$p, ncol = data$T_basis[m])
    expect_equal(out$shat2[[m]], expected, tolerance = 1e-12)
  }
})

test_that("compute_ser_statistics shat2 broadcasts per-scale sigma2", {
  data <- make_data(T_per_outcome = 64L)
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
  expected <- outer(1 / data$xtx_diag, expected_sigma2_per_pos)
  expect_equal(out$shat2[[1]], expected, tolerance = 1e-12)
})

# ---- update_fitted_values.mf_individual -------------------------------

test_that("update_fitted_values reconstructs fitted = fitted_without_l + Xb_l", {
  data <- make_data()
  model <- make_model(data, L = 3)
  model <- mfsusieR:::compute_residuals.mf_individual(data, NULL, model, 1)

  # Inject deterministic alpha and mu so we can verify the algebra.
  model$alpha[1, ] <- 1 / data$p
  for (m in seq_len(data$M)) {
    model$mu[[1]][[m]] <- matrix(0.5, nrow = data$p, ncol = data$T_basis[m])
  }

  out <- mfsusieR:::update_fitted_values.mf_individual(data, NULL, model, 1)

  for (m in seq_len(data$M)) {
    b_lm <- (1 / data$p) * matrix(0.5, nrow = data$p, ncol = data$T_basis[m])
    Xb_lm <- data$X %*% b_lm
    expected <- out$fitted_without_l[[m]] + Xb_lm
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
                 matrix(0, nrow = data$n, ncol = data$T_basis[m]),
                 tolerance = 1e-12)
  }
})

# ---- loglik.mf_individual ---------------------------------------------

# Helper: build a model WITH a scale-mixture prior + cross-modality combiner.
make_model_with_prior <- function(data, L = 3, sigma2_scalar = 1,
                                  prior_variance_grid = c(0.1, 0.5, 2)) {
  prior <- mfsusieR:::mf_prior_scale_mixture(
    data, prior_variance_grid = prior_variance_grid, null_prior_weight = 1
  )
  params <- list(L = L, prior_weights = NULL, prior = prior,
                 cross_outcome_prior = NULL,
                 residual_variance = lapply(seq_len(data$M),
                                            function(m) sigma2_scalar))
  var_y <- lapply(seq_len(data$M),
                  function(m) rep(sigma2_scalar, data$T_basis[m]))
  mfsusieR:::initialize_susie_model.mf_individual(data, params, var_y)
}

test_that("initialize stashes G_prior and cross_outcome_combiner on the model", {
  data <- make_data()
  model <- make_model_with_prior(data)
  expect_identical(length(model$G_prior), data$M)
  for (m in seq_len(data$M)) {
    expect_identical(length(model$G_prior[[m]]),
                     length(data$scale_index[[m]]))
  }
  expect_s3_class(model$cross_outcome_combiner,
                  "mf_prior_cross_outcome_independent")
})

test_that("loglik.mf_individual updates alpha (sums to 1), lbf finite, lbf_variable length J", {
  data <- make_data()
  model <- make_model_with_prior(data)
  model <- mfsusieR:::compute_residuals.mf_individual(data, NULL, model, 1)
  ser <- mfsusieR:::compute_ser_statistics.mf_individual(data, NULL, model, 1)

  out <- mfsusieR:::loglik.mf_individual(data, NULL, model, V = 1, ser, l = 1)

  expect_equal(sum(out$alpha[1, ]), 1, tolerance = 1e-12)
  expect_true(all(out$alpha[1, ] >= 0))
  expect_true(is.finite(out$lbf[1]))
  expect_identical(length(out$lbf_variable[1, ]), data$p)
  expect_true(all(is.finite(out$lbf_variable[1, ])))
})

test_that("loglik.mf_individual returns a scalar lbf_model when l = NULL", {
  data <- make_data()
  model <- make_model_with_prior(data)
  model <- mfsusieR:::compute_residuals.mf_individual(data, NULL, model, 1)
  ser <- mfsusieR:::compute_ser_statistics.mf_individual(data, NULL, model, 1)

  scalar <- mfsusieR:::loglik.mf_individual(data, NULL, model, V = 1, ser, l = NULL)

  expect_length(scalar, 1L)
  expect_true(is.finite(scalar))
  # And it equals the model$lbf[l] obtained via the l-supplied call.
  out <- mfsusieR:::loglik.mf_individual(data, NULL, model, V = 1, ser, l = 1)
  expect_equal(scalar, out$lbf[1], tolerance = 1e-12)
})

# ---- calculate_posterior_moments.mf_individual -----------------------

test_that("calculate_posterior_moments shapes mu, mu2 are J x T_basis[m] per modality", {
  data <- make_data()
  model <- make_model_with_prior(data)
  model <- mfsusieR:::compute_residuals.mf_individual(data, NULL, model, 1)

  out <- mfsusieR:::calculate_posterior_moments.mf_individual(data, NULL, model, V = 1, l = 1)

  for (m in seq_len(data$M)) {
    expect_identical(dim(out$mu[[1]][[m]]),
                     c(data$p, data$T_basis[m]))
    expect_identical(dim(out$mu2[[1]][[m]]),
                     c(data$p, data$T_basis[m]))
    # Second moment satisfies E[X^2] >= (E[X])^2 elementwise.
    expect_true(all(out$mu2[[1]][[m]] >= out$mu[[1]][[m]]^2 - 1e-12))
  }
})

test_that("calculate_posterior_moments with broadcast scalar sigma2 matches per-(modality, scale) ashr loop", {
  skip_if_not_installed("ashr")
  data <- make_data()
  model <- make_model_with_prior(data)
  model <- mfsusieR:::compute_residuals.mf_individual(data, NULL, model, 1)

  out <- mfsusieR:::calculate_posterior_moments.mf_individual(data, NULL, model, V = 1, l = 1)

  pw <- data$xtx_diag
  for (m in seq_len(data$M)) {
    Bhat_m <- model$residuals[[m]] / pw
    Shat_m <- sqrt(matrix(model$sigma2[[m]] / pw,
                          nrow = data$p, ncol = ncol(model$residuals[[m]])))
    indx_m <- data$scale_index[[m]]
    G_m <- model$G_prior[[m]]
    mu_expected  <- matrix(0, nrow = data$p, ncol = data$T_basis[m])
    mu2_expected <- matrix(0, nrow = data$p, ncol = data$T_basis[m])
    for (s in seq_along(indx_m)) {
      idx <- indx_m[[s]]
      g_s <- G_m[[s]]$fitted_g
      g_scaled <- ashr::normalmix(pi = g_s$pi, mean = g_s$mean, sd = g_s$sd)
      for (j in seq_len(data$p)) {
        d <- ashr::set_data(Bhat_m[j, idx], Shat_m[j, idx])
        mu_expected[j, idx]  <- ashr::postmean(g_scaled, d)
        post_sd_ji <- ashr::postsd(g_scaled, d)
        mu2_expected[j, idx] <- post_sd_ji^2 + mu_expected[j, idx]^2
      }
    }
    expect_equal(out$mu[[1]][[m]],  mu_expected,  tolerance = 1e-12)
    expect_equal(out$mu2[[1]][[m]], mu2_expected, tolerance = 1e-12)
  }
})

# ---- susieR-degenerate single-component path (C1 sanity) ------------

test_that("loglik / calculate_posterior_moments with K=1, no null reduce to susieR scalar formula", {
  # Single non-null component, sd_grid = 1, pi_grid = 1.
  data <- make_data()
  prior <- list(
    G_prior = lapply(seq_len(data$M), function(m) {
      g <- lapply(data$scale_index[[m]], function(idx) {
        list(fitted_g = list(pi = 1, sd = 1, mean = 0), idx = idx)
      })
      attr(g, "class") <- "mixture_normal_per_scale"
      g
    }),
    V_grid = lapply(seq_len(data$M), function(m) 1),
    pi     = lapply(seq_len(data$M), function(m) {
      matrix(1, nrow = length(data$scale_index[[m]]), ncol = 1)
    }),
    null_prior_weight    = 0,
    prior_variance_scope = "per_scale"
  )
  class(prior) <- "mf_prior_scale_mixture"
  params <- list(L = 1, prior_weights = NULL, prior = prior,
                 cross_outcome_prior = NULL,
                 residual_variance = lapply(seq_len(data$M), function(m) 1))
  var_y <- lapply(seq_len(data$M), function(m) rep(1, data$T_basis[m]))
  model <- mfsusieR:::initialize_susie_model.mf_individual(data, params, var_y)

  V <- 2.5
  model <- mfsusieR:::compute_residuals.mf_individual(data, NULL, model, 1)
  out <- mfsusieR:::calculate_posterior_moments.mf_individual(data, NULL, model, V = V, l = 1)

  pw <- data$xtx_diag
  for (m in seq_len(data$M)) {
    XtR <- model$residuals[[m]]
    Shat2 <- matrix(model$sigma2[[m]] / pw,
                    nrow = data$p, ncol = ncol(XtR))
    Bhat <- XtR / pw
    # susieR scalar formula: post_var = (1/V + pw/sigma2)^(-1)
    #                      = V*Shat^2 / (V + Shat^2)
    post_var  <- V * Shat2 / (V + Shat2)
    post_mean <- (V / (V + Shat2)) * Bhat
    expect_equal(out$mu[[1]][[m]],  post_mean,             tolerance = 1e-12)
    expect_equal(out$mu2[[1]][[m]], post_var + post_mean^2, tolerance = 1e-12)
  }
})

# ---- SER_posterior_e_loglik.mf_individual ---------------------------

test_that("SER_posterior_e_loglik returns scalar finite value", {
  data <- make_data()
  model <- make_model_with_prior(data)
  model <- mfsusieR:::compute_residuals.mf_individual(data, NULL, model, 1)
  model <- mfsusieR:::calculate_posterior_moments.mf_individual(data, NULL, model, V = 1, l = 1)

  out <- mfsusieR:::SER_posterior_e_loglik.mf_individual(data, NULL, model, 1)
  expect_length(out, 1L)
  expect_true(is.finite(out))
})

test_that("SER_posterior_e_loglik matches the per-modality, per-position formula", {
  data <- make_data()
  model <- make_model_with_prior(data)
  model <- mfsusieR:::compute_residuals.mf_individual(data, NULL, model, 1)
  model <- mfsusieR:::calculate_posterior_moments.mf_individual(data, NULL, model, V = 1, l = 1)

  # Hand re-derive the formula and compare.
  alpha_l <- model$alpha[1, ]
  pw <- data$xtx_diag
  expected <- 0
  for (m in seq_len(data$M)) {
    Eb_m  <- alpha_l * model$mu[[1]][[m]]
    Eb2_m <- alpha_l * model$mu2[[1]][[m]]
    R_m   <- model$raw_residuals[[m]]
    n     <- nrow(R_m)
    s2_t  <- mfsusieR:::mf_sigma2_per_position(data, model, m)

    XEb_m <- data$X %*% Eb_m
    quad_t <- colSums(R_m^2) - 2 * colSums(R_m * XEb_m) + colSums(pw * Eb2_m)
    expected <- expected + sum(-0.5 * n * log(2 * pi * s2_t) -
                               0.5 * quad_t / s2_t)
  }
  expect_equal(mfsusieR:::SER_posterior_e_loglik.mf_individual(data, NULL, model, 1),
               expected, tolerance = 1e-12)
})

test_that("SER_posterior_e_loglik handles per-scale sigma2 (model$sigma2 vector)", {
  data <- make_data(T_per_outcome = 64L)  # single modality so per-scale shape is unambiguous
  model <- make_model_with_prior(data)
  # Per-scale sigma2 vector (length S_m).
  S_m <- length(data$scale_index[[1]])
  model$sigma2[[1]] <- seq_len(S_m) * 0.1

  model <- mfsusieR:::compute_residuals.mf_individual(data, NULL, model, 1)
  model <- mfsusieR:::calculate_posterior_moments.mf_individual(data, NULL, model, V = 1, l = 1)

  out <- mfsusieR:::SER_posterior_e_loglik.mf_individual(data, NULL, model, 1)
  expect_length(out, 1L)
  expect_true(is.finite(out))

  # Cross-check the broadcast helper directly on the per-scale path.
  s2_t <- mfsusieR:::mf_sigma2_per_position(data, model, 1)
  expect_identical(length(s2_t), data$T_basis[1])
  for (s in seq_len(S_m)) {
    expect_true(all(s2_t[data$scale_index[[1]][[s]]] == model$sigma2[[1]][s]))
  }
})

# ---- compute_kl.mf_individual ----------------------------------------

test_that("compute_kl writes a finite scalar into model$KL[l]", {
  data <- make_data()
  model <- make_model_with_prior(data)
  model <- mfsusieR:::compute_residuals.mf_individual(data, NULL, model, 1)
  ser <- mfsusieR:::compute_ser_statistics.mf_individual(data, NULL, model, 1)
  model <- mfsusieR:::loglik.mf_individual(data, NULL, model, V = 1, ser, l = 1)
  model <- mfsusieR:::calculate_posterior_moments.mf_individual(data, NULL, model, V = 1, l = 1)

  out <- mfsusieR:::compute_kl.mf_individual(data, NULL, model, 1)

  expect_true(is.finite(out$KL[1]))
})

test_that("compute_kl identity: KL[l] == -lbf[l] - L_null + L_post", {
  data <- make_data()
  model <- make_model_with_prior(data)
  model <- mfsusieR:::compute_residuals.mf_individual(data, NULL, model, 1)
  ser <- mfsusieR:::compute_ser_statistics.mf_individual(data, NULL, model, 1)
  model <- mfsusieR:::loglik.mf_individual(data, NULL, model, V = 1, ser, l = 1)
  model <- mfsusieR:::calculate_posterior_moments.mf_individual(data, NULL, model, V = 1, l = 1)

  L_null <- 0
  for (m in seq_len(data$M)) {
    R_m <- model$raw_residuals[[m]]
    s2_t <- mfsusieR:::mf_sigma2_per_position(data, model, m)
    L_null <- L_null + sum(-0.5 * nrow(R_m) * log(2 * pi * s2_t) -
                           0.5 * colSums(R_m^2) / s2_t)
  }
  L_post <- mfsusieR:::SER_posterior_e_loglik.mf_individual(data, NULL, model, 1)
  expected <- -(model$lbf[1] + L_null) + L_post

  out <- mfsusieR:::compute_kl.mf_individual(data, NULL, model, 1)
  expect_equal(out$KL[1], expected, tolerance = 1e-12)
})

# ---- neg_loglik wrapper ----------------------------------------------

test_that("neg_loglik.mf_individual exponentiates V_param and negates loglik", {
  data <- make_data()
  model <- make_model_with_prior(data)
  model <- mfsusieR:::compute_residuals.mf_individual(data, NULL, model, 1)
  ser <- mfsusieR:::compute_ser_statistics.mf_individual(data, NULL, model, 1)

  V_param <- log(0.7)
  out <- mfsusieR:::neg_loglik.mf_individual(data, NULL, model, V_param, ser)
  expected_V <- exp(V_param)
  expected <- -mfsusieR:::loglik.mf_individual(data, NULL, model, expected_V, ser, l = NULL)
  expect_equal(out, expected, tolerance = 1e-12)
})

# ---- .onLoad: S3 methods are registered on susieR's namespace -------

test_that("S3 methods on mf_individual are registered on susieR's generics", {
  # `registerS3method` populates susieR's S3 dispatch table; verify
  # via `getS3method`, not top-level `exists`. Using
  # `optional = TRUE` so a missing entry returns NULL instead of erroring.
  for (g in c("compute_residuals", "compute_ser_statistics",
              "calculate_posterior_moments", "loglik", "neg_loglik",
              "compute_kl", "SER_posterior_e_loglik",
              "update_fitted_values", "initialize_susie_model")) {
    method_fn <- getS3method(g, "mf_individual",
                             optional = TRUE,
                             envir = asNamespace("susieR"))
    expect_true(is.function(method_fn),
                info = paste("Missing registration:", g, ".mf_individual"))
  }
})

test_that("susieR helpers are cached on the mfsusieR namespace", {
  pkg_ns <- asNamespace("mfsusieR")
  for (fn in c("warning_message", "SER_posterior_e_loglik")) {
    expect_true(is.function(get(fn, envir = pkg_ns)),
                info = paste("Missing cached binding:", fn))
  }
})

# ---- C3 fidelity vs mvf.susie.alpha::cal_Bhat_Shat_multfsusie ---------

test_that("compute_ser_statistics betahat matches fsusieR cal_Bhat_Shat at 1e-12 (initial residual)", {
  skip_if_not_installed("fsusieR")
  data <- make_data(n = 50, J = 10, T_per_outcome = 64L)
  model <- make_model(data, L = 3, sigma2_scalar = 1)

  # Initial residual (l = 1, all-zero effects): R_m = D[[m]] entirely.
  model <- mfsusieR:::compute_residuals.mf_individual(data, NULL, model, 1)
  ours <- mfsusieR:::compute_ser_statistics.mf_individual(data, NULL, model, 1)

  # Reference: per outcome, fsusieR's cal_Bhat_Shat on (D[[m]], X)
  # gives the same betahat (X^T Y / d).
  ref <- fsusieR:::cal_Bhat_Shat(Y = data$D[[1]], X = data$X,
                                  v1 = rep(1, data$n), lowc_wc = NULL)
  expect_equal(ours$betahat[[1]], ref$Bhat, tolerance = 1e-12)
})
