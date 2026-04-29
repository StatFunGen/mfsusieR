# PR group 7e: `fsusie()` thin wrapper.
#
# Per design.md D8d, `fsusie(Y, X, pos, ...)` is a one-line wrapper
# around `mfsusie(X, list(Y), list(pos), ...)` matching
# `fsusieR::susiF`'s argument order. The wrapper adds NO numerics;
# bit-identical to mfsusie() on the equivalent inputs (tolerance 0
# is acceptable here -- they are literally the same call path).

# ---- smoke -----------------------------------------------------

test_that("fsusie() runs end-to-end on a single-modality fixture", {
  set.seed(mfsusier_test_seed())
  n <- 30; J <- 8; T <- 32
  X <- matrix(rnorm(n * J), nrow = n)
  beta <- numeric(J); beta[1] <- 1
  Y <- X %*% matrix(rep(beta, T), nrow = J) + matrix(rnorm(n * T, sd = 0.3), nrow = n)

  fit <- fsusie(wavelet_qnorm = FALSE, Y, X, L = 3, max_iter = 30, verbose = FALSE)

  expect_s3_class(fit, c("mfsusie", "susie"), exact = FALSE)
  expect_identical(ncol(fit$alpha), as.integer(J))
  expect_true(is.numeric(fit$pip))
  expect_identical(length(fit$pip), as.integer(J))
  expect_identical(length(fit$mu[[1]]), 1L)  # M = 1
})

# ---- drop-in compatibility with mfsusie() ---------------------

test_that("fsusie(Y, X, pos) produces a bit-identical fit to mfsusie(X, list(Y), list(pos))", {
  set.seed(mfsusier_test_seed())
  n <- 30; J <- 8; T <- 32
  X <- matrix(rnorm(n * J), nrow = n)
  beta <- numeric(J); beta[1] <- 1
  Y <- X %*% matrix(rep(beta, T), nrow = J) + matrix(rnorm(n * T, sd = 0.3), nrow = n)
  pos <- seq_len(T)

  fit_w <- fsusie(wavelet_qnorm = FALSE, Y, X, pos = pos, L = 3, max_iter = 30, verbose = FALSE)
  fit_m <- mfsusie(wavelet_qnorm = FALSE, X, list(Y), pos = list(pos), L = 3, max_iter = 30, verbose = FALSE)

  expect_equal(fit_w$alpha,        fit_m$alpha,        tolerance = 0)
  expect_equal(fit_w$mu,           fit_m$mu,           tolerance = 0)
  expect_equal(fit_w$mu2,          fit_m$mu2,          tolerance = 0)
  expect_equal(fit_w$pip,          fit_m$pip,          tolerance = 0)
  expect_equal(fit_w$elbo,         fit_m$elbo,         tolerance = 0)
  expect_equal(fit_w$sigma2,       fit_m$sigma2,       tolerance = 0)
  expect_identical(fit_w$niter,    fit_m$niter)
})

# ---- vector Y coerced to single-column matrix -----------------

test_that("fsusie() coerces a numeric vector Y to a single-column matrix", {
  set.seed(mfsusier_test_seed())
  n <- 30; J <- 6
  X <- matrix(rnorm(n * J), nrow = n)
  beta <- numeric(J); beta[1] <- 1
  y <- as.numeric(X %*% beta + rnorm(n, sd = 0.3))

  fit <- fsusie(wavelet_qnorm = FALSE, y, X, L = 3, max_iter = 30, verbose = FALSE)
  expect_s3_class(fit, "mfsusie")
  expect_identical(ncol(fit$alpha), as.integer(J))
})

# ---- forbidden multi-modality argument errors out -------------

test_that("fsusie() errors when the user passes a multi-modality argument", {
  set.seed(mfsusier_test_seed())
  X <- matrix(rnorm(30), nrow = 5)
  Y <- matrix(rnorm(40), nrow = 5)
  expect_error(
    fsusie(wavelet_qnorm = FALSE, max_iter = 100, Y, X, cross_outcome_prior = list()),
    "single-outcome wrapper"
  )
})

test_that("fsusie() errors on non-numeric Y", {
  X <- matrix(rnorm(30), nrow = 5)
  Y_bad <- letters[1:6]
  expect_error(fsusie(wavelet_qnorm = FALSE, max_iter = 100, Y_bad, X), "must be a numeric matrix")
})
