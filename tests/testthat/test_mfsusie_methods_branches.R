# Coverage tests for the error paths, scalar (T_m = 1) fallbacks,
# and edge branches of `mf_post_smooth`, `summary.mfsusie`, and
# friends. These are mfsusieR-internal contract tests; the
# upstream functional fine-mapping packages do not expose
# `predict`/`coef`/`summary` S3 methods on the same fit shape, so
# the assertions are internal-consistency rather than apple-to-
# apple bit identity.

build_small_fit <- function(seed = 41L, n = 60L, p = 12L, T_m = 32L) {
  set.seed(seed)
  X    <- matrix(rnorm(n * p), n, p)
  beta <- numeric(p); beta[3L] <- 1.0
  shape <- exp(-((seq_len(T_m) - T_m / 2)^2) / (2 * 4^2))
  Y    <- X %*% (matrix(beta, p, 1) %*% matrix(shape, 1, T_m)) +
          matrix(rnorm(n * T_m, sd = 0.3), n)
  fsusie(Y, X, L = 2, max_iter = 30, verbose = FALSE)
}

# --- mf_post_smooth error paths --------------------------------

test_that("mf_post_smooth errors when input is not an mfsusie fit", {
  expect_error(mf_post_smooth(list()), "must be an `mfsusie`")
})

test_that("mf_post_smooth errors when level is out of (0, 1)", {
  fit <- build_small_fit()
  expect_error(mf_post_smooth(fit, level = 0), "`level`")
  expect_error(mf_post_smooth(fit, level = 1), "`level`")
})

test_that("mf_post_smooth errors when residuals slot is missing for non-scalewise methods", {
  fit <- build_small_fit()
  fit_no_resid <- fit
  fit_no_resid$residuals <- NULL
  for (method in c("TI", "HMM", "smash")) {
    if (method == "smash" && !requireNamespace("smashr", quietly = TRUE)) next
    expect_error(mf_post_smooth(fit_no_resid, method = method),
                 "requires `fit\\$residuals`")
  }
})

test_that("mf_post_smooth(scalewise) succeeds without the residuals slot", {
  fit <- build_small_fit()
  fit_no_resid <- fit
  fit_no_resid$residuals <- NULL
  fit_no_resid$lead_X    <- NULL
  out <- mf_post_smooth(fit_no_resid, method = "scalewise")
  expect_true(!is.null(out$smoothed$scalewise$effect_curves))
  expect_true(!is.null(out$smoothed$scalewise$credible_bands))
})

# --- Summary / print branches ----------------------------------

test_that("summary.mfsusie carries cs table when credible sets exist; print formats it", {
  fit <- build_small_fit()
  s <- summary(fit)
  expect_true(!is.null(s$cs))
  expect_gte(nrow(s$cs), 1L)
  out <- capture.output(print(s))
  expect_true(any(grepl("Credible sets:", out)))
})

test_that("summary.mfsusie returns NULL cs and prints `No credible sets` when fit has none", {
  set.seed(101L)
  n <- 40L; p <- 8L; T_m <- 16L
  X <- matrix(rnorm(n * p), n, p)
  Y <- matrix(rnorm(n * T_m), n, T_m)   # pure noise
  fit_null <- fsusie(Y, X, L = 2, max_iter = 20, verbose = FALSE)
  s <- summary(fit_null)
  expect_null(s$cs)
  out <- capture.output(print(s))
  expect_true(any(grepl("No credible sets", out)))
})

# --- Scalar (T_m = 1) post-smooth fallbacks --------------------

build_scalar_mfsusie <- function(seed = 53L, n = 60L, p = 8L) {
  set.seed(seed)
  X    <- matrix(rnorm(n * p), n, p)
  beta <- numeric(p); beta[2L] <- 1.0
  Y_scalar     <- as.vector(X %*% beta) + rnorm(n, sd = 0.3)
  Y_functional <- X %*% (matrix(beta, p, 1) %*%
                         matrix(exp(-((1:32 - 16)^2 / 8)), 1, 32L)) +
                  matrix(rnorm(n * 32, sd = 0.3), n)
  mfsusie(X, list(matrix(Y_scalar, n, 1L), Y_functional),
          L = 2, max_iter = 30, verbose = FALSE)
}

test_that("post-smooth methods handle scalar (T_m = 1) outcomes via the scalewise fallback", {
  fit <- build_scalar_mfsusie()
  for (method in c("TI", "HMM", "scalewise")) {
    s  <- mf_post_smooth(fit, method = method)
    ec <- s$smoothed[[method]]$effect_curves
    # First outcome is scalar; effect_curves[[1]][[l]] is length 1.
    expect_length(ec[[1L]][[1L]], 1L)
    # Second outcome is functional; effect_curves[[2]][[l]] is length 32.
    expect_length(ec[[2L]][[1L]], 32L)
  }
  if (requireNamespace("smashr", quietly = TRUE)) {
    s  <- mf_post_smooth(fit, method = "smash")
    ec <- s$smoothed$smash$effect_curves
    expect_length(ec[[1L]][[1L]], 1L)
    expect_length(ec[[2L]][[1L]], 32L)
  }
})

# --- coef.mfsusie / predict.mfsusie -----------------------------

test_that("coef.mfsusie returns a per-outcome list of L x T_basis matrices", {
  fit <- build_small_fit()
  cf  <- coef(fit)
  expect_length(cf, 1L)              # M = 1
  expect_equal(nrow(cf[[1L]]), 2L)   # L = 2
  expect_equal(ncol(cf[[1L]]), 32L)  # T_basis
})

test_that("predict.mfsusie(newx) returns one matrix per outcome with the right shape", {
  fit <- build_small_fit()
  X_new <- matrix(rnorm(20L * 12L), 20L, 12L)
  pr <- predict(fit, newx = X_new)
  expect_length(pr, 1L)
  expect_equal(dim(pr[[1L]]), c(20L, 32L))
})
