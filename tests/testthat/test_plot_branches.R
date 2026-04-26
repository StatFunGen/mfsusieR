# Coverage tests for the multi-modality, single-outcome focus,
# and lfsr-overlay paths in `mfsusie_plot()` and
# `mfsusie_plot_lfsr()`. The plot functions have no upstream
# functional fine-mapping equivalent, so the assertions are
# layout / smoke only; numerical assertions live in the kernels
# the plot calls into.

null_device <- function() {
  if (capabilities("png")) {
    f <- tempfile(fileext = ".png")
    grDevices::png(f, width = 600, height = 400)
    on.exit({ grDevices::dev.off(); unlink(f) }, add = TRUE)
  } else {
    grDevices::pdf(NULL)
    on.exit(grDevices::dev.off(), add = TRUE)
  }
  invisible(NULL)
}

build_multi_outcome_fit <- function(seed = 5L) {
  set.seed(seed)
  n <- 80L; p <- 10L
  X <- matrix(rnorm(n * p), n, p)
  beta <- numeric(p); beta[3L] <- 1.0
  shape1 <- exp(-((seq_len(32L) - 16L)^2) / (2 * 3^2))
  shape2 <- exp(-((seq_len(64L) - 24L)^2) / (2 * 5^2))
  Y1 <- X %*% (matrix(beta, p, 1) %*% matrix(shape1, 1, 32L)) +
        matrix(rnorm(n * 32L, sd = 0.3), n)
  Y2 <- X %*% (matrix(beta, p, 1) %*% matrix(shape2, 1, 64L)) +
        matrix(rnorm(n * 64L, sd = 0.3), n)
  mfsusie(X, list(Y1, Y2), L = 3, max_iter = 30, verbose = FALSE)
}

build_three_outcome_fit <- function(seed = 7L) {
  set.seed(seed)
  n <- 80L; p <- 10L
  X <- matrix(rnorm(n * p), n, p)
  beta <- numeric(p); beta[3L] <- 1.0
  shape <- exp(-((seq_len(32L) - 16L)^2) / (2 * 3^2))
  Y1 <- X %*% (matrix(beta, p, 1) %*% matrix(shape, 1, 32L)) +
        matrix(rnorm(n * 32L, sd = 0.3), n)
  Y2 <- X %*% (matrix(beta, p, 1) %*% matrix(shape, 1, 32L)) +
        matrix(rnorm(n * 32L, sd = 0.4), n)
  Y3 <- X %*% (matrix(beta, p, 1) %*% matrix(shape, 1, 32L)) +
        matrix(rnorm(n * 32L, sd = 0.5), n)
  mfsusie(X, list(Y1, Y2, Y3), L = 3, max_iter = 30, verbose = FALSE)
}

# --- Multi-modality grid layout ---------------------------------

test_that("mfsusie_plot() lays out the M > 1 panel grid", {
  null_device()
  fit <- build_multi_outcome_fit()
  expect_silent(mfsusie_plot(fit))
})

test_that("mfsusie_plot() M = 3 triggers the rows*cols-n_panels remainder branch", {
  # M + 1 = 4 panels, cols = ceiling(sqrt(4)) = 2, rows = 2.
  # remainder = 4 - 4 = 0 (no plot.new() filler). M = 3 with
  # M + 1 = 4 hits the rectangular path.
  null_device()
  fit <- build_three_outcome_fit()
  expect_silent(mfsusie_plot(fit))
})

test_that("mfsusie_plot(m = ...) draws a single outcome panel and accepts m bounds", {
  null_device()
  fit <- build_multi_outcome_fit()
  expect_silent(mfsusie_plot(fit, m = 2L))
  expect_error(mfsusie_plot(fit, m = 5L), "must be in 1..2")
  expect_error(mfsusie_plot(fit, m = 0L), "must be in 1..2")
})

test_that("mfsusie_plot() rejects non-mfsusie input", {
  null_device()
  expect_error(mfsusie_plot(list()), "must be an `mfsusie`")
})

# --- TI / HMM post-smoothed multi-outcome plotting --------------

test_that("mfsusie_plot() draws bands and lfsr overlays after multi-outcome HMM smoothing", {
  null_device()
  fit <- build_multi_outcome_fit()
  fit_h <- mf_post_smooth(fit, method = "HMM")
  expect_silent(mfsusie_plot(fit_h))
  expect_silent(mfsusie_plot(fit_h, show_lfsr_curve = FALSE))
  expect_silent(mfsusie_plot(fit_h, show_affected_region = FALSE))
})

test_that("mfsusie_plot() in errorbar style on a multi-outcome smoothed fit", {
  null_device()
  fit <- build_multi_outcome_fit()
  fit_s <- mf_post_smooth(fit, method = "TI",
                          wavelet_filter = 1L,
                          wavelet_family = "DaubExPhase")
  expect_silent(mfsusie_plot(fit_s, effect_style = "errorbar"))
})

# --- lfsr bubble plot -------------------------------------------

test_that("mfsusie_plot_lfsr() runs on an M > 1 HMM-smoothed fit", {
  null_device()
  fit   <- build_multi_outcome_fit()
  fit_h <- mf_post_smooth(fit, method = "HMM")
  expect_silent(mfsusie_plot_lfsr(fit_h))
})

test_that("mfsusie_plot_lfsr() accepts truth on a multi-outcome fit", {
  null_device()
  fit   <- build_multi_outcome_fit()
  fit_h <- mf_post_smooth(fit, method = "HMM")
  # M = 2; truth must be length-M list per the documented API.
  truth_1 <- logical(32L); truth_1[14:18] <- TRUE
  truth_2 <- logical(64L); truth_2[20:28] <- TRUE
  expect_silent(mfsusie_plot_lfsr(fit_h, truth = list(truth_1, truth_2)))
  expect_error(mfsusie_plot_lfsr(fit_h, truth = truth_1),
               "length-M list when M > 1")
})

# --- show_grid_dots and lfsr_threshold sweeps -------------------

test_that("mfsusie_plot() honors show_grid_dots and various lfsr thresholds", {
  null_device()
  fit <- build_multi_outcome_fit()
  fit_h <- mf_post_smooth(fit, method = "HMM")
  expect_silent(mfsusie_plot(fit_h, show_grid_dots = TRUE,
                             lfsr_threshold = 0.001))
  expect_silent(mfsusie_plot(fit_h, lwd = 2.5, add_legend = FALSE))
})
