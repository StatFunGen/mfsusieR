# `mf_post_smooth(method = "TI", scaling = ...)` per-scale vs
# uniform shrinkage modes, plus `...` pass-through to ashr::ash.

test_that("TI uniform scaling differs from per_scale on the same fit", {
  skip_if_not_installed("ashr")
  skip_if_not_installed("wavethresh")
  set.seed(1L)
  N <- 60L; J <- 30L; T_basis <- 32L
  X    <- matrix(rnorm(N * J), N, J)
  beta <- matrix(0, J, T_basis)
  beta[5L, ] <- exp(-((seq_len(T_basis) - T_basis / 4)^2) / (2 * 3^2))
  Y1   <- X %*% beta + matrix(rnorm(N * T_basis, sd = 0.3), N, T_basis)
  base <- mfsusieR::mfsusie(X = X, Y = list(Y1),
                            L = 3L, max_iter = 10L, verbose = FALSE)

  fit_ps <- mfsusieR::mf_post_smooth(base, method = "TI",
                                     scaling = "per_scale")
  fit_un <- mfsusieR::mf_post_smooth(base, method = "TI",
                                     scaling = "uniform",
                                     overwrite_previous = TRUE)
  ec_ps <- fit_ps$smoothed[["TI"]]$effect_curves[[1L]]
  ec_un <- fit_un$smoothed[["TI"]]$effect_curves[[1L]]

  for (l in seq_along(ec_ps)) {
    expect_length(ec_un[[l]], length(ec_ps[[l]]))
  }
  # The two scaling modes apply different ash shrinkage and so
  # produce different position-space curves.
  any_diff <- any(vapply(seq_along(ec_ps), function(l)
    !isTRUE(all.equal(ec_ps[[l]], ec_un[[l]], tolerance = 1e-6)),
    logical(1L)))
  expect_true(any_diff)
})

test_that("TI `...` forwards nullweight through to ashr::ash", {
  skip_if_not_installed("ashr")
  skip_if_not_installed("wavethresh")
  set.seed(2L)
  N <- 60L; J <- 30L; T_basis <- 32L
  X    <- matrix(rnorm(N * J), N, J)
  beta <- matrix(0, J, T_basis); beta[5L, ] <- 1.0
  Y1   <- X %*% beta + matrix(rnorm(N * T_basis, sd = 0.3), N, T_basis)
  base <- mfsusieR::mfsusie(X = X, Y = list(Y1),
                            L = 3L, max_iter = 10L, verbose = FALSE)

  fit_default <- mfsusieR::mf_post_smooth(base, method = "TI")
  fit_strong  <- mfsusieR::mf_post_smooth(base, method = "TI",
                                          nullweight = 300,
                                          overwrite_previous = TRUE)
  ec_default <- fit_default$smoothed[["TI"]]$effect_curves[[1L]]
  ec_strong  <- fit_strong$smoothed[["TI"]]$effect_curves[[1L]]

  # `nullweight = 300` is much stronger than the per_scale default
  # (30) and shrinks the effect curve toward zero.
  any_diff <- any(vapply(seq_along(ec_default), function(l)
    !isTRUE(all.equal(ec_default[[l]], ec_strong[[l]], tolerance = 1e-12)),
    logical(1L)))
  expect_true(any_diff)
})

test_that("HMM `...` forwards nullweight through to ashr::ash", {
  skip_if_not_installed("ashr")
  set.seed(3L)
  N <- 60L; J <- 30L; T_basis <- 32L
  X    <- matrix(rnorm(N * J), N, J)
  beta <- matrix(0, J, T_basis); beta[5L, ] <- 1.0
  Y1   <- X %*% beta + matrix(rnorm(N * T_basis, sd = 0.3), N, T_basis)
  base <- mfsusieR::mfsusie(X = X, Y = list(Y1),
                            L = 3L, max_iter = 10L, verbose = FALSE)

  fit_default <- mfsusieR::mf_post_smooth(base, method = "HMM",
                                          halfK = 6L)
  fit_strong  <- suppressWarnings(
    mfsusieR::mf_post_smooth(base, method = "HMM",
                             halfK = 6L, nullweight = 300,
                             overwrite_previous = TRUE))
  ec_default <- fit_default$smoothed[["HMM"]]$effect_curves[[1L]]
  ec_strong  <- fit_strong$smoothed[["HMM"]]$effect_curves[[1L]]

  any_diff <- any(vapply(seq_along(ec_default), function(l)
    !isTRUE(all.equal(ec_default[[l]], ec_strong[[l]], tolerance = 1e-12)),
    logical(1L)))
  expect_true(any_diff)
})
