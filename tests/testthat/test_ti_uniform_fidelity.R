# C2 fidelity: `univariate_ti_regression(scaling = "uniform")` is
# a port of `fsusieR:::univariate_TI_regression_IS`. Both kernels
# are deterministic; bit-identity at machine-epsilon scale is
# required on every supported input shape.

test_that("univariate_ti_regression(uniform) matches fsusieR::univariate_TI_regression_IS at T = 64", {
  skip_if_no_fsusier()
  set.seed(101L)
  n   <- 60L
  T_m <- 64L
  X   <- matrix(rnorm(n), n, 1L)
  shape <- exp(-((seq_len(T_m) - T_m / 4)^2) / (2 * 4^2))
  Y   <- X %*% matrix(0.7 * shape, 1, T_m) +
         matrix(rnorm(n * T_m, sd = 0.5), n, T_m)

  alpha  <- 0.05
  z_crit <- stats::qnorm(1 - alpha / 2)
  ours <- mfsusieR:::univariate_ti_regression(
    Y_pos         = Y,
    x_eff         = as.numeric(X),
    filter_number = 1L,
    family        = "DaubExPhase",
    z_crit        = z_crit,
    scaling       = "uniform")
  ref  <- fsusieR:::univariate_TI_regression_IS(
    Y             = Y,
    X             = X,
    filter.number = 1L,
    family        = "DaubExPhase",
    alpha         = alpha)

  expect_equal(as.numeric(ours$effect_estimate),
               as.numeric(ref$effect_estimate),
               tolerance = 1e-14)
  expect_equal(as.numeric(ours$cred_band),
               as.numeric(ref$cred_band),
               tolerance = 1e-14)
})

test_that("univariate_ti_regression(uniform) matches fsusieR::univariate_TI_regression_IS at T = 128", {
  skip_if_no_fsusier()
  set.seed(102L)
  n   <- 80L
  T_m <- 128L
  X   <- matrix(rnorm(n), n, 1L)
  shape <- c(rep(0, 32), rep(1, 32), rep(0, 32), rep(-0.4, 32))
  Y   <- X %*% matrix(0.5 * shape, 1, T_m) +
         matrix(rnorm(n * T_m, sd = 0.4), n, T_m)

  alpha  <- 0.10
  z_crit <- stats::qnorm(1 - alpha / 2)
  ours <- mfsusieR:::univariate_ti_regression(
    Y_pos         = Y,
    x_eff         = as.numeric(X),
    filter_number = 1L,
    family        = "DaubExPhase",
    z_crit        = z_crit,
    scaling       = "uniform")
  ref  <- fsusieR:::univariate_TI_regression_IS(
    Y             = Y,
    X             = X,
    filter.number = 1L,
    family        = "DaubExPhase",
    alpha         = alpha)

  expect_equal(as.numeric(ours$effect_estimate),
               as.numeric(ref$effect_estimate),
               tolerance = 1e-14)
  expect_equal(as.numeric(ours$cred_band),
               as.numeric(ref$cred_band),
               tolerance = 1e-14)
})
