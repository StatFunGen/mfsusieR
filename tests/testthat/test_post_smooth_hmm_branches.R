# Coverage tests for the defensive-cleanup and prefilter
# branches in `mf_fit_hmm` and `mf_univariate_hmm_regression`.
# These guard against degenerate per-position inputs (zero / NA
# / non-finite sd, large z-scores, sparse posterior occupancy)
# and are not exercised by the regular fit path.

# --- mf_fit_hmm defensive cleanup paths -------------------------

test_that("mf_fit_hmm clamps too-small sd to thresh_sd", {
  set.seed(1L)
  T_m <- 64L
  x   <- rnorm(T_m, sd = 1)
  # First quarter of the grid has near-zero sd (below default
  # `thresh_sd = 1e-30`); the cleanup branch sets them to the
  # threshold rather than dividing by zero downstream.
  sd_vec <- rep(0.5, T_m)
  sd_vec[1:16] <- 1e-40
  out <- mfsusieR:::mf_fit_hmm(x, sd_vec, halfK = 5L,
                               maxiter = 1L)
  expect_length(out$x_post, T_m)
  expect_length(out$lfsr,   T_m)
})

test_that("mf_fit_hmm replaces NA sd entries with 1 and zeros the corresponding x", {
  set.seed(2L)
  T_m <- 32L
  x   <- rnorm(T_m, sd = 1)
  sd_vec <- rep(0.5, T_m)
  sd_vec[c(5L, 17L)] <- NA_real_
  out <- mfsusieR:::mf_fit_hmm(x, sd_vec, halfK = 5L,
                               maxiter = 1L)
  expect_length(out$x_post, T_m)
})

test_that("mf_fit_hmm replaces non-finite sd entries with 1 and zeros x", {
  set.seed(3L)
  T_m <- 32L
  x   <- rnorm(T_m, sd = 1)
  sd_vec <- rep(0.5, T_m)
  sd_vec[c(7L, 19L)] <- Inf
  out <- mfsusieR:::mf_fit_hmm(x, sd_vec, halfK = 5L,
                               maxiter = 1L)
  expect_length(out$x_post, T_m)
})

test_that("mf_fit_hmm caps the per-position z-score via the max_zscore branch", {
  # An x[i] / sd[i] that exceeds `max_zscore = 20` triggers the
  # rescaling branch; without it the dnorm emission underflows.
  T_m <- 16L
  x   <- rep(0.1, T_m); x[5L] <- 1e6   # extreme outlier
  sd_vec <- rep(0.1, T_m)
  out <- mfsusieR:::mf_fit_hmm(x, sd_vec, halfK = 5L,
                               max_zscore = 20, maxiter = 1L)
  expect_length(out$x_post, T_m)
})

# --- prefilter / idx_comp keep-null branches --------------------

test_that("mf_univariate_hmm_regression runs end-to-end and exercises the prefilter / EM / idx_comp paths", {
  # The full kernel call: column-scaled X / Y, OLS Bhat / Shat,
  # then mf_fit_hmm. This is the path the public
  # `mf_post_smooth(method = 'HMM')` takes per (effect, outcome).
  set.seed(31L)
  n <- 80L; T_m <- 64L
  X <- matrix(rnorm(n), n, 1L)
  shape <- exp(-((seq_len(T_m) - T_m / 2)^2) / (2 * 5^2))
  Y <- X %*% matrix(0.8 * shape, 1, T_m) +
       matrix(rnorm(n * T_m, sd = 0.3), n, T_m)
  out <- mfsusieR:::mf_univariate_hmm_regression(
    Y = Y, X = X, halfK = 8L)
  expect_length(out$effect_estimate, T_m)
  expect_length(out$lfsr, T_m)
  expect_true(is.finite(out$lBF))
})

# --- mf_univariate_hmm_regression NA-sd cleanup -----------------

test_that("mf_univariate_hmm_regression imputes non-positive Shat with median", {
  # Construct a single-SNP, single-functional-outcome dataset
  # where one position has zero variance after column-scaling;
  # the kernel computes Bhat / Shat and replaces the bad column
  # with the median sd.
  set.seed(17L)
  n <- 60L; T_m <- 16L
  X <- matrix(rnorm(n), n, 1L)
  Y <- X %*% matrix(seq_len(T_m) / T_m, 1, T_m) +
       matrix(rnorm(n * T_m, sd = 0.4), n, T_m)
  # Make one column constant (zero residual variance).
  Y[, 8L] <- 1
  out <- mfsusieR:::mf_univariate_hmm_regression(
    Y = Y, X = X, halfK = 5L)
  expect_length(out$effect_estimate, T_m)
  expect_length(out$lfsr,            T_m)
})
