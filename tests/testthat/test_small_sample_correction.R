# Reference tests for `small_sample_correction = TRUE` on
# `mfsusie()` and `fsusie()`.
#
# 1. Default path is bit-identical (sanity: argument has no effect
#    when FALSE).
# 2. Johnson-t kernel runs to convergence on small n.
# 3. Structural: at small n with a known null variant, the Johnson-t
#    kernel produces PIPs <= the Wakefield kernel at the null
#    variant. The correction tightens variable selection.

build_small_n_sim <- function(seed = 23L, n = 40L, T_m = 32L,
                              p = 30L, signal_var = 5L) {
  set.seed(seed)
  X    <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p); beta[signal_var] <- 1
  Y    <- X %*% matrix(beta, ncol = 1L) %*% matrix(1, 1, T_m) +
          matrix(rnorm(n * T_m, sd = 0.5), n, T_m)
  list(X = X, Y = Y, p = p, T_m = T_m, n = n,
       signal_var = signal_var)
}

# --- 1. Default path unchanged -----------------------------------

test_that("small_sample_correction = FALSE leaves the fit bit-identical", {
  sim <- build_small_n_sim()
  fit_default <- fsusie(sim$Y, sim$X, L = 5,
                        verbose = FALSE)
  fit_explicit <- fsusie(sim$Y, sim$X, L = 5,
                         small_sample_correction = FALSE,
                         verbose = FALSE)
  expect_equal(fit_default$alpha, fit_explicit$alpha, tolerance = 0)
  expect_equal(fit_default$pip,   fit_explicit$pip,   tolerance = 0)
})

# --- 2. Johnson-t runs to convergence ----------------------------

test_that("small_sample_correction = TRUE runs to convergence on small n", {
  sim <- build_small_n_sim()
  fit_johnson <- fsusie(sim$Y, sim$X, L = 5,
                        small_sample_correction = TRUE,
                        verbose = FALSE)
  expect_true(fit_johnson$converged)
  expect_equal(length(fit_johnson$pip), sim$p)
})

# --- 3. Structural correction at null variants ------------------

test_that("Johnson-t reduces PIP at known-null variants vs Wakefield", {
  # n = 40 is the regime where the small-sample correction acts.
  # The signal variant should still be near 1 under both kernels;
  # null variants should have lower PIP under Johnson-t.
  sim <- build_small_n_sim(n = 40L, T_m = 32L, p = 30L)
  fit_wake <- fsusie(sim$Y, sim$X, L = 5,
                     small_sample_correction = FALSE,
                     verbose = FALSE)
  fit_john <- fsusie(sim$Y, sim$X, L = 5,
                     small_sample_correction = TRUE,
                     verbose = FALSE)
  # Signal variant: both kernels recover it (PIP > 0.9).
  expect_gt(fit_wake$pip[sim$signal_var], 0.9)
  expect_gt(fit_john$pip[sim$signal_var], 0.9)
  # Null variants: aggregate PIP at the non-signal variants is
  # not larger under Johnson-t. (Some individual variants may
  # tie or wiggle by EM noise; the aggregate sums out the noise.)
  null_idx  <- setdiff(seq_len(sim$p), sim$signal_var)
  expect_lte(sum(fit_john$pip[null_idx]),
             sum(fit_wake$pip[null_idx]) + 1e-6)
})

# --- 4. Fidelity vs upstream Johnson-t kernel -------------------

test_that("mixture_log_bf_per_scale_johnson matches fsusieR::log_BF with df", {
  skip_if_not_installed("fsusieR")
  set.seed(7L)
  n <- 50L; J <- 15L; T_m <- 64L
  data <- mfsusieR:::create_mf_individual(
    X = matrix(rnorm(n * J), n, J),
    Y = list(matrix(rnorm(n * T_m), n, T_m)),
    verbose = FALSE)

  sd_grid <- c(0, 0.5, 1.5)
  pi_grid <- c(0.6, 0.2, 0.2)
  g_norm  <- ashr::normalmix(pi = pi_grid,
                             mean = rep(0, length(pi_grid)),
                             sd   = sd_grid)
  rec <- list(fitted_g = g_norm); class(rec) <- "ash"
  G_prior <- rep(list(rec), length(data$scale_index[[1L]]))
  attr(G_prior, "class") <- "mixture_normal_per_scale"

  pw    <- data$xtx_diag
  D_1   <- data$D[[1L]]
  Bhat  <- crossprod(data$X, D_1) / pw
  Shat  <- sqrt(matrix(1 / pw, J, ncol(D_1)))

  df <- n - 1L

  ref <- fsusieR::log_BF(G_prior,
                         Bhat = Bhat, Shat = Shat,
                         lowc_wc = NULL,
                         indx_lst = data$scale_index[[1L]],
                         df = df)

  ours <- numeric(J)
  for (s in seq_along(data$scale_index[[1L]])) {
    idx <- data$scale_index[[1L]][[s]]
    ours <- ours + mfsusieR:::mixture_log_bf_per_scale_johnson(
      Bhat[, idx, drop = FALSE],
      Shat[, idx, drop = FALSE],
      sd_grid, pi_grid, V_scale = 1, df = df)
  }
  expect_equal(ours, ref, tolerance = 1e-12)
})

# --- 5. Argument validation -------------------------------------

test_that("small_sample_correction rejects non-logical input", {
  sim <- build_small_n_sim()
  expect_error(
    fsusie(sim$Y, sim$X, L = 5,
           small_sample_correction = "yes",
           verbose = FALSE),
    "TRUE` or `FALSE`"
  )
})
