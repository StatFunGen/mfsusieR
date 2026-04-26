# Bit-identity tests for the HMM post-smoother. Tolerance = 0.

skip_if_no_fsusier_for_HMM <- function() {
  testthat::skip_if_not_installed("fsusieR")
  testthat::skip_if_not_installed("ashr")
}

test_that("mf_fit_hmm matches upstream fit_hmm bit-for-bit", {
  skip_if_no_fsusier_for_HMM()
  set.seed(7)
  T_pos <- 16L
  x   <- rnorm(T_pos, sd = 0.3)
  sd  <- runif(T_pos, 0.1, 0.4)
  ours <- mfsusieR:::mf_fit_hmm(x, sd, halfK = 5L)
  ref  <- fsusieR:::fit_hmm(x, sd, halfK = 5)
  expect_equal(ours$x_post, ref$x_post, tolerance = 0)
  expect_equal(ours$lfsr,   ref$lfsr,   tolerance = 0)
  expect_equal(ours$log_BF, ref$log_BF, tolerance = 0)
  expect_equal(ours$mu,     ref$mu,     tolerance = 0)
})

test_that("mf_fit_hmm matches across halfK / T_pos / random seeds", {
  skip_if_no_fsusier_for_HMM()
  cases <- expand.grid(halfK = c(5L, 10L, 20L),
                       T_pos = c(16L, 32L, 64L),
                       seed  = c(1L, 17L, 42L))
  for (i in seq_len(nrow(cases))) {
    set.seed(cases$seed[i])
    x  <- rnorm(cases$T_pos[i], sd = 0.4)
    sd <- runif(cases$T_pos[i], 0.1, 0.5)
    o  <- mfsusieR:::mf_fit_hmm(x, sd, halfK = cases$halfK[i])
    r  <- fsusieR:::fit_hmm(x, sd, halfK = cases$halfK[i])
    expect_equal(o$x_post, r$x_post, tolerance = 0,
                 info = sprintf("halfK=%d T=%d seed=%d",
                                cases$halfK[i], cases$T_pos[i], cases$seed[i]))
    expect_equal(o$lfsr, r$lfsr, tolerance = 0,
                 info = sprintf("halfK=%d T=%d seed=%d",
                                cases$halfK[i], cases$T_pos[i], cases$seed[i]))
    expect_equal(o$log_BF, r$log_BF, tolerance = 0)
  }
})

test_that("mf_univariate_hmm_regression matches upstream bit-for-bit", {
  skip_if_no_fsusier_for_HMM()
  set.seed(11)
  n     <- 60L
  T_pos <- 64L
  Y <- matrix(rnorm(n * T_pos), n)
  X <- matrix(rnorm(n), n, 1)
  ours <- mfsusieR:::mf_univariate_hmm_regression(Y, X, halfK = 20L)
  ref  <- fsusieR:::univariate_HMM_regression(Y, X, halfK = 20)
  expect_equal(ours$effect_estimate, as.numeric(ref$effect_estimate),
               tolerance = 0)
  expect_equal(ours$lfsr, as.numeric(ref$lfsr), tolerance = 0)
  expect_equal(ours$lBF,  ref$lBF,              tolerance = 0)
})

test_that("mf_univariate_hmm_regression matches across multiple simulated configs", {
  skip_if_no_fsusier_for_HMM()
  # Vary n, T, signal strength, and seed to cover several
  # operating points: pure null, sparse non-null, dense non-null.
  cases <- list(
    list(seed = 1L,  n = 50L,  T = 32L,  beta = 0,    sd = 0.3),
    list(seed = 23L, n = 100L, T = 64L,  beta = 1.2,  sd = 0.4),
    list(seed = 99L, n = 80L,  T = 128L, beta = -0.8, sd = 0.25),
    list(seed = 5L,  n = 60L,  T = 64L,  beta = 0.5,  sd = 0.5)
  )
  for (cs in cases) {
    set.seed(cs$seed)
    X     <- matrix(rnorm(cs$n), cs$n, 1)
    shape <- exp(-((seq_len(cs$T) - cs$T / 2)^2) / (2 * 4^2))
    Y     <- matrix(rnorm(cs$n * cs$T, sd = cs$sd), cs$n) +
             X %*% matrix(cs$beta * shape, 1, cs$T)
    ours <- mfsusieR:::mf_univariate_hmm_regression(Y, X, halfK = 15L)
    ref  <- fsusieR:::univariate_HMM_regression(Y, X, halfK = 15)
    expect_equal(ours$effect_estimate, as.numeric(ref$effect_estimate),
                 tolerance = 0,
                 info = sprintf("seed=%d n=%d T=%d beta=%g",
                                cs$seed, cs$n, cs$T, cs$beta))
    expect_equal(ours$lfsr, as.numeric(ref$lfsr), tolerance = 0,
                 info = sprintf("seed=%d", cs$seed))
    expect_equal(ours$lBF, ref$lBF, tolerance = 0,
                 info = sprintf("seed=%d", cs$seed))
  }
})
