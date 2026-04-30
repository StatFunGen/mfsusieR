# `mf_post_smooth(method = "ash")` cycle-spinning + per-coefficient
# `ashr::ash` smoother. Counterpart to `method = "smash"` (smashr).

test_that("method = 'ash' works without smashr", {
  skip_if_not_installed("ashr")
  skip_if_not_installed("wavethresh")
  set.seed(1L)
  N <- 50L; J <- 30L; T_basis <- 16L
  X    <- matrix(rnorm(N * J), N, J)
  beta <- matrix(0, J, T_basis); beta[5L, ] <- 1.0
  Y1   <- X %*% beta + matrix(rnorm(N * T_basis, sd = 0.3), N, T_basis)
  fit  <- mfsusieR::mfsusie(X = X, Y = list(Y1),
                            L = 3L, max_iter = 10L, verbose = FALSE)
  fit  <- mfsusieR::mf_post_smooth(fit, method = "ash")
  payload <- fit$smoothed[["ash"]]
  expect_length(payload$effect_curves[[1L]], nrow(fit$alpha))
  for (l in seq_len(nrow(fit$alpha))) {
    expect_length(payload$effect_curves[[1L]][[l]], T_basis)
    expect_equal(dim(payload$credible_bands[[1L]][[l]]),
                 c(T_basis, 2L))
    expect_length(payload$lfsr_curves[[1L]][[l]], T_basis)
    expect_true(all(payload$lfsr_curves[[1L]][[l]] >= 0))
    expect_true(all(payload$lfsr_curves[[1L]][[l]] <= 0.5 + 1e-12))
  }
})

test_that("method = 'smash' and method = 'ash' produce the same payload shape", {
  skip_if_not_installed("smashr")
  set.seed(1L)
  N <- 50L; J <- 30L; T_basis <- 16L
  X    <- matrix(rnorm(N * J), N, J)
  beta <- matrix(0, J, T_basis); beta[5L, ] <- 1.0
  Y1   <- X %*% beta + matrix(rnorm(N * T_basis, sd = 0.3), N, T_basis)
  fit_a <- mfsusieR::mfsusie(X = X, Y = list(Y1),
                             L = 3L, max_iter = 10L, verbose = FALSE)
  fit_b <- fit_a
  fit_a <- mfsusieR::mf_post_smooth(fit_a, method = "ash")
  fit_b <- mfsusieR::mf_post_smooth(fit_b, method = "smash")

  pa <- fit_a$smoothed[["ash"]]
  pb <- fit_b$smoothed[["smash"]]
  expect_equal(length(pa$effect_curves[[1L]]),
               length(pb$effect_curves[[1L]]))
  for (l in seq_along(pa$effect_curves[[1L]])) {
    expect_equal(length(pa$effect_curves[[1L]][[l]]),
                 length(pb$effect_curves[[1L]][[l]]))
    expect_equal(dim(pa$credible_bands[[1L]][[l]]),
                 dim(pb$credible_bands[[1L]][[l]]))
  }
})

test_that("method = 'smash' errors when smashr is not installed", {
  if (requireNamespace("smashr", quietly = TRUE))
    skip("smashr is installed; cannot exercise the missing-dep error.")
  set.seed(1L)
  N <- 30L; J <- 20L; T_basis <- 8L
  X    <- matrix(rnorm(N * J), N, J)
  beta <- matrix(0, J, T_basis); beta[5L, ] <- 1.0
  Y1   <- X %*% beta + matrix(rnorm(N * T_basis, sd = 0.3), N, T_basis)
  fit  <- mfsusieR::mfsusie(X = X, Y = list(Y1),
                            L = 2L, max_iter = 5L, verbose = FALSE)
  expect_error(
    mfsusieR::mf_post_smooth(fit, method = "smash"),
    "smashr"
  )
})

test_that("mf_smash_ash returns the documented shape", {
  skip_if_not_installed("ashr")
  skip_if_not_installed("wavethresh")
  set.seed(2L)
  T_pos <- 60L
  truth <- c(rep(0, 20), rep(1, 20), rep(0, 20))
  noisy <- truth + rnorm(T_pos, sd = 0.5)
  out <- mfsusieR:::mf_smash_ash(noisy, noise_level = 0.5,
                                 n.shifts = 8L)
  expect_length(out$mu.est,     T_pos)
  expect_length(out$mu.est.var, T_pos)
  expect_true(all(out$mu.est.var >= 0))
})
