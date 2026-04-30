# Smoke tests for `mfsusie_plot()` and `mfsusie_plot_lfsr()`.
# We render to a NULL graphics device so the tests are silent and
# allocate no PNG files.

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

# Tiny single-outcome fit with two true effects. Reusable across
# tests in this file.
build_tiny_fit <- function(seed = 1L, T_m = 32L, p = 8L) {
  set.seed(seed)
  n <- 80L
  X <- matrix(rnorm(n * p), n, p)
  beta <- matrix(0, p, T_m)
  shape1 <- exp(-((seq_len(T_m) -  9L)^2) / (2 * 3^2))
  shape2 <- exp(-((seq_len(T_m) - 24L)^2) / (2 * 3^2))
  beta[2L, ] <- 1.6 * shape1
  beta[6L, ] <- -1.6 * shape2
  Y <- X %*% beta + matrix(rnorm(n * T_m, sd = 0.4), n)
  fsusie(wavelet_qnorm = FALSE, wavelet_standardize = FALSE, max_iter = 100, Y, X, pos = seq_len(T_m), L = 5, verbose = FALSE)
}

# ---- mfsusie_plot ---------------------------------------------------

test_that("mfsusie_plot() runs in band style on a raw fit", {
  null_device()
  fit <- build_tiny_fit()
  expect_silent(mfsusie_plot(fit))
})

test_that("mfsusie_plot() honors facet_cs = 'overlay' and 'stack'", {
  null_device()
  fit <- build_tiny_fit()
  expect_silent(mfsusie_plot(fit, facet_cs = "overlay"))
  expect_silent(mfsusie_plot(fit, facet_cs = "stack"))
})

test_that("mfsusie_plot() errors when errorbar is asked without bands", {
  null_device()
  fit <- build_tiny_fit()
  expect_error(mfsusie_plot(fit, effect_style = "errorbar"),
               "credible bands")
})

test_that("mfsusie_plot() runs in errorbar style after TI smoothing", {
  null_device()
  fit <- build_tiny_fit()
  fit_s <- mf_post_smooth(fit, method = "TI",
                          wavelet_filter = 1L,
                          wavelet_family = "DaubExPhase")
  expect_silent(mfsusie_plot(fit_s, effect_style = "errorbar"))
  expect_silent(mfsusie_plot(fit_s, effect_style = "errorbar",
                              facet_cs = "stack"))
})

test_that("mfsusie_plot() draws lfsr secondary axis after HMM smoothing", {
  null_device()
  fit <- build_tiny_fit()
  fit_h <- mf_post_smooth(fit, method = "HMM")
  expect_silent(mfsusie_plot(fit_h))
  expect_silent(mfsusie_plot(fit_h, lfsr_threshold = 0.05))
})

# ---- mfsusie_plot_lfsr ---------------------------------------------

test_that("mfsusie_plot_lfsr() runs on an HMM-smoothed fit, M = 1", {
  null_device()
  fit <- build_tiny_fit()
  fit_h <- mf_post_smooth(fit, method = "HMM")
  expect_silent(mfsusie_plot_lfsr(fit_h))
})

test_that("mfsusie_plot_lfsr() accepts truth (single bool vec, list of vecs)", {
  null_device()
  fit <- build_tiny_fit()
  fit_h <- mf_post_smooth(fit, method = "HMM")
  T_m <- length(fit$dwt_meta$pos[[1L]])
  m_a <- logical(T_m); m_a[6:12]  <- TRUE
  m_b <- logical(T_m); m_b[20:28] <- TRUE
  expect_silent(mfsusie_plot_lfsr(fit_h, truth = m_a))
  expect_silent(mfsusie_plot_lfsr(fit_h, truth = list(m_a, m_b)))
})

test_that("mfsusie_plot_lfsr() errors without HMM lfsr", {
  null_device()
  fit <- build_tiny_fit()
  expect_error(mfsusie_plot_lfsr(fit), "HMM-smoothed")
})

test_that("plot.mfsusie() forwards args", {
  null_device()
  fit <- build_tiny_fit()
  expect_silent(plot(fit))
  expect_silent(plot(fit, facet_cs = "overlay"))
})
