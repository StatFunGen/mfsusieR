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
  fsusie(wavelet_qnorm = FALSE, wavelet_standardize = FALSE, Y, X, L = 2, max_iter = 30, verbose = FALSE)
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

test_that("mf_post_smooth errors when smoothing inputs are missing for non-scalewise methods", {
  fit <- build_small_fit()
  fit_no_inputs <- fit
  fit_no_inputs$Y_grid <- NULL
  fit_no_inputs$X_eff  <- NULL
  for (method in c("TI", "HMM", "smash")) {
    if (method == "smash" && !requireNamespace("smashr", quietly = TRUE)) next
    expect_error(mf_post_smooth(fit_no_inputs, method = method),
                 "smoothing inputs")
  }
})

test_that("mf_post_smooth(scalewise) succeeds without the smoothing inputs", {
  fit <- build_small_fit()
  fit_no_inputs <- fit
  fit_no_inputs$Y_grid <- NULL
  fit_no_inputs$X_eff  <- NULL
  out <- mf_post_smooth(fit_no_inputs, method = "scalewise")
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
  fit_null <- fsusie(wavelet_qnorm = FALSE, wavelet_standardize = FALSE, Y, X, L = 2, max_iter = 20, verbose = FALSE)
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
  mfsusie(wavelet_qnorm = FALSE, wavelet_standardize = FALSE, X, list(matrix(Y_scalar, n, 1L), Y_functional),
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

# --- Multi-method state ----------------------------------------

test_that("mf_post_smooth refuses to overwrite an existing method without overwrite_previous", {
  fit <- build_small_fit()
  fit_s <- mf_post_smooth(fit, method = "TI")
  expect_error(mf_post_smooth(fit_s, method = "TI"),
               "already on this fit")
  # With overwrite_previous = TRUE the call replaces it.
  fit_s2 <- mf_post_smooth(fit_s, method = "TI",
                           overwrite_previous = TRUE)
  expect_named(fit_s2$smoothed, "TI")
})

test_that("mf_post_smooth accumulates entries for distinct methods", {
  fit <- build_small_fit()
  fit_s <- mf_post_smooth(fit, method = "TI")
  fit_s <- mf_post_smooth(fit_s, method = "HMM")
  fit_s <- mf_post_smooth(fit_s, method = "scalewise")
  expect_setequal(names(fit_s$smoothed), c("TI", "HMM", "scalewise"))
})

test_that("coef.mfsusie(smooth_method) returns smoothed curves on the L x T_basis shape", {
  fit <- build_small_fit()
  fit_s <- mf_post_smooth(fit, method = "TI")
  cf  <- coef(fit_s, smooth_method = "TI")
  expect_equal(attr(cf, "smooth_method"), "TI")
  # Same outer shape as the raw coef path: list[M] of L x T_basis matrices.
  expect_length(cf, 1L)
  expect_equal(dim(cf[[1L]]), c(2L, 32L))
})

test_that("coef.mfsusie(smooth_method) errors when the named method is not applied", {
  fit <- build_small_fit()
  expect_error(coef(fit, smooth_method = "TI"),
               "No methods have been applied")
  fit_s <- mf_post_smooth(fit, method = "TI")
  expect_error(coef(fit_s, smooth_method = "HMM"),
               "Applied methods")
})

test_that(".pick_smooth_method honors `requested` and errors on missing method", {
  fit <- build_small_fit()
  fit_s <- mf_post_smooth(fit, method = "TI")
  expect_equal(mfsusieR:::.pick_smooth_method(fit_s, "TI"), "TI")
  expect_error(mfsusieR:::.pick_smooth_method(fit_s, "smash"),
               "is not on this fit")
})

test_that("mfsusie_plot emits a hint when multiple smoothings are stacked on the fit", {
  null_dev <- if (capabilities("png")) {
    f <- tempfile(fileext = ".png")
    grDevices::png(f, width = 600, height = 400)
    on.exit({ grDevices::dev.off(); unlink(f) }, add = TRUE)
  } else {
    grDevices::pdf(NULL); on.exit(grDevices::dev.off(), add = TRUE)
  }
  fit   <- build_small_fit()
  fit_s <- mf_post_smooth(fit, method = "TI")
  fit_s <- mf_post_smooth(fit_s, method = "HMM")
  # testthat edition 3 captures the message stream itself; using
  # `capture.output(..., type = "message")` returns an empty vector
  # because testthat's reporter intercepts first. `expect_message()`
  # is the idiomatic check.
  expect_message(mfsusie_plot(fit_s),
                 regexp = "Other smoothings on this fit")
})

test_that("mfsusie_plot_lfsr errors with a clean message when no HMM smoothing is on the fit", {
  null_dev <- if (capabilities("png")) {
    f <- tempfile(fileext = ".png")
    grDevices::png(f, width = 600, height = 400)
    on.exit({ grDevices::dev.off(); unlink(f) }, add = TRUE)
  } else {
    grDevices::pdf(NULL); on.exit(grDevices::dev.off(), add = TRUE)
  }
  fit   <- build_small_fit()
  fit_s <- mf_post_smooth(fit, method = "TI")
  expect_error(mfsusie_plot_lfsr(fit_s),
               "HMM-smoothed lfsr")
})
