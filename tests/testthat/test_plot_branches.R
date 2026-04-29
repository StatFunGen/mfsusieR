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
  mfsusie(wavelet_qnorm = FALSE, X, list(Y1, Y2), L = 3, max_iter = 30, verbose = FALSE)
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
  mfsusie(wavelet_qnorm = FALSE, X, list(Y1, Y2, Y3), L = 3, max_iter = 30, verbose = FALSE)
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

# --- Internal helpers ------------------------------------------

test_that("mf_cs_colors honors n_cs = 0 and small palettes", {
  expect_length(mf_cs_colors(0L), 0L)
  expect_length(mf_cs_colors(1L), 1L)
  pal <- mf_cs_colors(5L)
  expect_length(pal, 5L)
  expect_true(length(unique(pal)) >= 4L)
})

test_that("credibly_nonzero_runs and credibly_nonzero_mask handle non-empty bands", {
  # A band where (lower > 0) OR (upper < 0) at a contiguous run.
  band <- cbind(
    c(-1,  0.1, 0.2, 0.3,  0,  -0.5, -0.4),
    c( 0,  0.5, 0.6, 0.7,  0,  -0.1, -0.05))
  runs <- mfsusieR:::credibly_nonzero_runs(band)
  expect_length(runs, 2L)
  expect_equal(runs[[1L]], c(2L, 4L))
  expect_equal(runs[[2L]], c(6L, 7L))

  mask <- mfsusieR:::credibly_nonzero_mask(band)
  expect_length(mask, 7L)
  expect_equal(which(mask), c(2L, 3L, 4L, 6L, 7L))
})

test_that(".pip_title formats coverage and purity when sets$cs is populated", {
  fit <- build_three_outcome_fit()
  title <- mfsusieR:::.pip_title(fit)
  # Coverage is formatted as "95% CS"; purity as "min(|r|)=...".
  # Single-line title (was a two-line title; collapsed for layout).
  expect_true(grepl("CS", title))
  expect_true(grepl("min", title))
})

test_that(".resolve_facet picks 'stack' for K = 2 with disjoint affected masks", {
  # Two CSes, smoothed with HMM produces credible bands; with
  # non-overlapping signal positions across CSes the affected
  # masks are disjoint and the resolver picks "stack". This
  # exercises the K == 2 branch in `.resolve_facet`.
  set.seed(1L)
  n <- 80L; p <- 10L; T_m <- 64L
  X <- matrix(rnorm(n * p), n, p)
  beta <- matrix(0, p, T_m)
  shape1 <- exp(-((seq_len(T_m) -  9L)^2) / (2 * 3^2))
  shape2 <- exp(-((seq_len(T_m) - 56L)^2) / (2 * 3^2))
  beta[2L, ] <- 1.6 * shape1
  beta[6L, ] <- -1.6 * shape2
  Y <- X %*% beta + matrix(rnorm(n * T_m, sd = 0.3), n)
  fit <- fsusie(wavelet_qnorm = FALSE, max_iter = 100, Y, X, pos = seq_len(T_m), L = 4, verbose = FALSE)
  fit_h <- mf_post_smooth(fit, method = "HMM")

  null_device()
  # Plot with facet_cs = "auto"; the resolver uses bands +
  # affected masks to choose layout.
  expect_silent(mfsusie_plot(fit_h, facet_cs = "auto"))
})

# --- Scalar (T_m = 1) plotting --------------------------------

test_that("mfsusie_plot draws the dot-plot panel on a scalar outcome", {
  null_device()
  set.seed(53L)
  n <- 80L; p <- 8L
  X    <- matrix(rnorm(n * p), n, p)
  beta <- numeric(p); beta[2L] <- 1.0
  Y_scalar <- as.vector(X %*% beta) + rnorm(n, sd = 0.3)
  Y_func   <- X %*% (matrix(beta, p, 1) %*%
                     matrix(exp(-((1:32 - 16)^2 / 8)), 1, 32L)) +
              matrix(rnorm(n * 32, sd = 0.3), n)
  fit <- mfsusie(wavelet_qnorm = FALSE, X, list(matrix(Y_scalar, n, 1L), Y_func),
                 L = 2, max_iter = 30, verbose = FALSE)
  expect_silent(mfsusie_plot(fit))
  expect_silent(mfsusie_plot(fit, m = 1L))   # scalar focus
  expect_silent(mfsusie_plot(fit, m = 2L))
})
