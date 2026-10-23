# Regression tests for the 2026-10-23 review bug-fix bundle
# (`optimization/results.md`, section 2). BUG-2 and BUG-3 tests live
# next to their kernels (test_per_scale_normal_degeneracy.R,
# test_variance_and_elbo.R, test_mvf_alpha_fidelity.R).

# ---- BUG-4: low-count mask is applied to the M-step cache ----------

test_that("cached and uncached M-step likelihoods agree on low-count columns", {
  set.seed(mfsusier_test_seed())
  n <- 60; p <- 15
  X <- matrix(rnorm(n * p), n)
  data <- mfsusieR:::create_mf_individual(X, list(matrix(rnorm(n * 32), n)))
  # The masking code only reads `data$lowc_idx`; inject a mask.
  data$lowc_idx[[1]] <- c(20L, 21L, 22L)
  prior  <- mfsusieR:::mf_prior_scale_mixture(data,
              prior_variance_scope = "per_outcome")
  params <- list(L = 2, prior = prior, cross_iter_prior = TRUE,
                 residual_variance = NULL)
  model <- mfsusieR:::ibss_initialize.mf_individual(data, params)
  model <- mfsusieR:::compute_residuals.mf_individual(data, params, model, 1)
  ss    <- mfsusieR:::compute_ser_statistics.mf_individual(data, params, model, 1)
  G     <- model$G_prior[[1]][[1]]
  b     <- ss$betahat[[1]][, G$idx, drop = FALSE]
  s     <- sqrt(ss$shat2[[1]][, G$idx, drop = FALSE])
  L_cached <- mfsusieR:::mf_em_likelihood_per_scale(
    b, s, G$fitted_g$sd,
    sdmat_cache     = model$iter_cache$sdmat[[1]][[1]],
    log_sdmat_cache = model$iter_cache$log_sdmat[[1]][[1]])
  L_plain <- mfsusieR:::mf_em_likelihood_per_scale(b, s, G$fitted_g$sd)
  expect_equal(L_cached, L_plain, tolerance = 1e-12)
})

# ---- BUG-5: wavelet_qnorm does not touch the caller's RNG ----------

test_that("mf_quantile_normalize restores the caller's RNG state", {
  Y <- matrix(c(1, 1, 2, 2, 3, 3), 3, 2)
  set.seed(123); before <- .Random.seed
  invisible(mfsusieR:::mf_quantile_normalize(Y))
  expect_identical(.Random.seed, before)
  # Output is still the seeded, reproducible transform.
  expect_identical(mfsusieR:::mf_quantile_normalize(Y),
                   mfsusieR:::mf_quantile_normalize(Y))
})

# ---- BUG-6: NA rows do not disable low-count detection ------------

test_that("mf_low_count_indices ignores NA rows", {
  set.seed(mfsusier_test_seed())
  Wd <- matrix(rnorm(40), 10, 4); Wd[, 2] <- 0
  Wd_na <- Wd; Wd_na[1, ] <- NA
  expect_identical(mfsusieR:::mf_low_count_indices(Wd), 2L)
  expect_identical(mfsusieR:::mf_low_count_indices(Wd_na), 2L)
})

# ---- BUG-7: OLS covariate adjustment with missing rows -------------

test_that("mf_residualize_ols fits on complete rows and keeps NA rows NA", {
  set.seed(mfsusier_test_seed())
  n <- 50; T_m <- 8
  Z <- cbind(1, matrix(rnorm(n * 2), n))
  Y <- Z %*% matrix(rnorm(3 * T_m), 3) + matrix(rnorm(n * T_m), n)
  X <- matrix(rnorm(n * 4), n)
  Y_na <- Y; Y_na[c(3, 7), ] <- NA
  out <- mfsusieR:::mf_residualize_ols(Y_na, Z, X)
  obs <- setdiff(seq_len(n), c(3, 7))
  ref <- mfsusieR:::mf_residualize_ols(Y[obs, ], Z[obs, ], X[obs, ])
  expect_true(all(is.na(out$Y_adjusted[c(3, 7), ])))
  expect_equal(out$Y_adjusted[obs, ], ref$Y_adjusted, tolerance = 1e-12)
  expect_equal(out$X_adjusted[obs, ], ref$X_adjusted, tolerance = 1e-12)
  expect_equal(out$fitted_func, ref$fitted_func, tolerance = 1e-12)
  # Complete data: unchanged closed form.
  full <- mfsusieR:::mf_residualize_ols(Y, Z)
  expect_equal(full$Y_adjusted,
               Y - Z %*% solve(crossprod(Z), crossprod(Z, Y)), tolerance = 1e-12)
})

test_that("mf_adjust_for_covariates(wavelet_eb) rejects missing values", {
  set.seed(mfsusier_test_seed())
  Y <- matrix(rnorm(20 * 16), 20); Y[2, ] <- NA
  expect_error(mf_adjust_for_covariates(Y, matrix(rnorm(20), 20)),
               "does not support missing values")
})

# ---- BUG-9: cs_min_lead_pip keeps the CS fields aligned ------------

test_that("cs_min_lead_pip subsets cs, purity, cs_index and coverage together", {
  set.seed(mfsusier_test_seed())
  n <- 150; p <- 40; T_m <- 16
  X <- matrix(rnorm(n * p), n)
  shape <- exp(-((seq_len(T_m) - 8)^2) / 8)
  Y <- outer(X[, 5], shape) * 1.2 + outer(X[, 30], shape) * 0.7 +
       matrix(rnorm(n * T_m), n)
  fit0 <- mfsusie(X, list(Y), L = 5, verbose = FALSE)
  skip_if(length(fit0$sets$cs) < 2L, "fixture produced fewer than 2 CS")
  lead <- vapply(fit0$sets$cs, function(s) max(fit0$pip[s]), numeric(1))
  skip_if(length(unique(round(lead, 6))) < 2L, "lead PIPs tie")
  thr <- mean(sort(lead, decreasing = TRUE)[1:2])
  fit <- mfsusie(X, list(Y), L = 5, verbose = FALSE, cs_min_lead_pip = thr)
  keep <- lead >= thr
  expect_identical(fit$sets$cs, fit0$sets$cs[keep])
  expect_identical(fit$sets$cs_index, fit0$sets$cs_index[keep])
  expect_identical(fit$sets$coverage, fit0$sets$coverage[keep])
  expect_identical(names(fit$sets$cs), paste0("L", fit$sets$cs_index))
  expect_equal(nrow(fit$sets$purity), length(fit$sets$cs))
})

