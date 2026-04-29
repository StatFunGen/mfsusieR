# Tests for mf_summarize_effects(): post-processing accessor that
# returns one row per (CS, outcome, contiguous run) triple where
# the credible band excludes zero.

make_smoothed_toy <- function(T_m = 32L, signal_idx = 3L,
                              seed = 1L) {
  set.seed(seed)
  n <- 80L; p <- 20L
  X <- matrix(rnorm(n * p), n)
  beta <- numeric(p); beta[signal_idx] <- 1.2
  shape <- exp(-((seq_len(T_m) - T_m / 2)^2) / (2 * 6^2))
  Y <- X %*% (matrix(beta, p, 1) %*% matrix(shape, 1, T_m)) +
         matrix(rnorm(n * T_m, sd = 0.3), n)
  fit <- fsusie(wavelet_qnorm = FALSE, Y, X, L = 1, max_iter = 30, verbose = FALSE)
  mf_post_smooth(fit, method = "TI")
}

test_that("mf_summarize_effects returns one row per (CS, outcome, run)", {
  fit_s <- make_smoothed_toy()
  out   <- mf_summarize_effects(fit_s)
  expect_s3_class(out, "data.frame")
  expect_named(out, c("cs_index", "outcome", "start", "end", "n_positions"))
  expect_true(all(vapply(out, is.integer, logical(1L))))
  expect_true(nrow(out) >= 1L)

  # Bounds: start <= end <= T_basis[1]; n_positions = end - start + 1.
  T_m <- fit_s$dwt_meta$T_basis[1L]
  expect_true(all(out$start >= 1L & out$start <= T_m))
  expect_true(all(out$end   >= 1L & out$end   <= T_m))
  expect_true(all(out$start <= out$end))
  expect_true(all(out$n_positions == out$end - out$start + 1L))
})

test_that("mf_summarize_effects errors when no smoothed bands are on the fit", {
  set.seed(1)
  n <- 60L; p <- 12L
  X <- matrix(rnorm(n * p), n)
  Y <- list(matrix(rnorm(n * 32L), n))
  fit <- mfsusie(wavelet_qnorm = FALSE, X, Y, L = 2, max_iter = 20, verbose = FALSE)
  expect_error(mf_summarize_effects(fit), "post-smoothed")
})

test_that("mf_summarize_effects respects the smoother priority and explicit smooth_method", {
  fit_s <- make_smoothed_toy()
  out_default  <- mf_summarize_effects(fit_s)
  out_explicit <- mf_summarize_effects(fit_s, smooth_method = "TI")
  expect_equal(out_default, out_explicit)
})

test_that("mf_summarize_effects on a multi-outcome fit returns rows per (m, l)", {
  set.seed(1)
  n <- 80L; p <- 12L
  X <- matrix(rnorm(n * p), n)
  T_per <- c(32L, 16L)
  beta <- numeric(p); beta[2] <- 1.5
  Y <- lapply(T_per, function(T_m) {
    shape <- exp(-((seq_len(T_m) - T_m / 2)^2) / (2 * 4^2))
    X %*% (matrix(beta, p, 1) %*% matrix(shape, 1, T_m)) +
      matrix(rnorm(n * T_m, sd = 0.3), n)
  })
  fit <- mfsusie(wavelet_qnorm = FALSE, X, Y, L = 1, max_iter = 30, verbose = FALSE)
  fit_s <- mf_post_smooth(fit, method = "TI")
  out <- mf_summarize_effects(fit_s)
  # Every row's outcome must be in 1..M.
  expect_true(all(out$outcome %in% seq_along(T_per)))
  # Every row's start/end fits within the corresponding T_basis[m].
  T_basis <- fit_s$dwt_meta$T_basis
  expect_true(all(out$end <= T_basis[out$outcome]))
})
