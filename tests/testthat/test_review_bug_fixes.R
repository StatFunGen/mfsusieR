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

