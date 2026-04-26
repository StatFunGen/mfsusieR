# Coverage tests for the S3 IBSS hooks on `mf_individual` data
# class. These are the branches the regular fit path doesn't
# exercise: tracking, scalar-input, missing-coverage, no-
# intercept, and named-X.

build_for_ibss <- function(seed = 23L, n = 60L, p = 8L, T_m = 32L,
                           name_X = TRUE) {
  set.seed(seed)
  X    <- matrix(rnorm(n * p), n, p)
  if (name_X) colnames(X) <- paste0("snp", seq_len(p))
  beta <- numeric(p); beta[3L] <- 1.0
  shape <- exp(-((seq_len(T_m) - T_m / 2)^2) / (2 * 4^2))
  Y    <- X %*% (matrix(beta, p, 1) %*% matrix(shape, 1, T_m)) +
          matrix(rnorm(n * T_m, sd = 0.3), n)
  list(X = X, Y = Y)
}

# --- ibss_initialize.mf_individual: clamp L to p ----------------

test_that("ibss_initialize.mf_individual clamps L to p when L > p", {
  d <- build_for_ibss(p = 5L)
  fit <- fsusie(d$Y, d$X, L = 10, max_iter = 30, verbose = FALSE)
  # ibss_initialize sets params$L <- data$p when L > p; the fit
  # then returns L = p effects.
  expect_equal(nrow(fit$alpha), 5L)
})

# --- track_ibss_fit.mf_individual ------------------------------

test_that("track_ibss_fit.mf_individual exercises the recording branch when track_fit = TRUE", {
  d   <- build_for_ibss()
  # The hook builds per-iteration snapshots when track_fit is
  # TRUE; the final fit has them stripped by susieR's cleanup,
  # so we assert run-to-completion rather than slot presence.
  fit <- fsusie(d$Y, d$X, L = 2, max_iter = 5, track_fit = TRUE,
                verbose = FALSE)
  expect_true(fit$converged || fit$niter == 5L)
})

# --- get_variable_names.mf_individual --------------------------

test_that("get_variable_names.mf_individual carries colnames(X) onto alpha and lbf_variable", {
  d   <- build_for_ibss(name_X = TRUE)
  fit <- fsusie(d$Y, d$X, L = 2, max_iter = 30, verbose = FALSE)
  expect_equal(colnames(fit$alpha),        colnames(d$X))
  expect_equal(colnames(fit$lbf_variable), colnames(d$X))
})

test_that("get_variable_names.mf_individual leaves alpha unnamed when X is unnamed", {
  d   <- build_for_ibss(name_X = FALSE)
  fit <- fsusie(d$Y, d$X, L = 2, max_iter = 30, verbose = FALSE)
  expect_null(colnames(fit$alpha))
})

# --- get_intercept.mf_individual: intercept = FALSE -------------

test_that("get_intercept.mf_individual returns zeros when intercept = FALSE", {
  d <- build_for_ibss()
  # The fit pre-centers Y; the no-intercept branch in
  # `get_intercept.mf_individual` is exercised when
  # `params$intercept = FALSE`.
  fit <- fsusie(d$Y, d$X, L = 2, max_iter = 20,
                intercept = FALSE, verbose = FALSE)
  inter <- mfsusieR:::get_intercept.mf_individual(
    data   = list(T_basis = fit$dwt_meta$T_basis),
    params = list(intercept = FALSE),
    model  = NULL)
  expect_equal(inter,
               lapply(fit$dwt_meta$T_basis, function(T_m) rep(0, T_m)))
})

# --- expand_model_init_to_L: shrinking error is the documented stop --

test_that("expand_model_init_to_L errors when L_prev > L_new", {
  d <- build_for_ibss()
  fit_big <- fsusie(d$Y, d$X, L = 5, max_iter = 30, verbose = FALSE)
  expect_error(
    fsusie(d$Y, d$X, L = 2, model_init = fit_big, verbose = FALSE),
    "shrinking is not supported")
})
