# Tests for the user-facing S3 methods on `mfsusie` fit objects:
# predict, coef, fitted, print, summary.

make_toy_fit <- function(n = 30, J = 8, T_per = c(64L, 32L),
                         L = 4, max_iter = 30,
                         seed = mfsusier_test_seed()) {
  set.seed(seed)
  X <- matrix(rnorm(n * J), nrow = n)
  beta <- numeric(J); beta[1] <- 1; beta[3] <- -0.5
  Y <- lapply(T_per, function(T_m) {
    eta <- X %*% beta
    matrix(rep(eta, T_m), nrow = n) + matrix(rnorm(n * T_m, sd = 0.3), nrow = n)
  })
  list(
    X = X, Y = Y,
    fit = mfsusie(X, Y, L = L, max_iter = max_iter, verbose = FALSE)
  )
}

# ---- fitted --------------------------------------------------

test_that("fitted.mfsusie returns a list of length M with n x T_m matrices", {
  toy <- make_toy_fit()
  out <- fitted(toy$fit)

  expect_true(is.list(out))
  expect_identical(length(out), 2L)
  for (m in seq_along(out)) {
    expect_true(is.matrix(out[[m]]))
    expect_identical(nrow(out[[m]]), nrow(toy$X))
    expect_identical(ncol(out[[m]]), ncol(toy$Y[[m]]))
    expect_false(any(is.na(out[[m]])))
  }
})

# ---- predict --------------------------------------------------

test_that("predict.mfsusie with newx returns predictions of shape n_new x T_m", {
  toy <- make_toy_fit()
  set.seed(2)
  newx <- matrix(rnorm(5 * ncol(toy$X)), nrow = 5)
  out <- predict(toy$fit, newx)

  expect_identical(length(out), 2L)
  for (m in seq_along(out)) {
    expect_identical(nrow(out[[m]]), 5L)
    expect_identical(ncol(out[[m]]), ncol(toy$Y[[m]]))
  }
})

test_that("predict.mfsusie with newx = NULL returns the training fitted curves", {
  toy <- make_toy_fit()
  out_pred    <- predict(toy$fit, newx = NULL)
  out_fitted  <- fitted(toy$fit)
  for (m in seq_along(out_pred)) {
    expect_equal(out_pred[[m]], out_fitted[[m]], tolerance = 1e-12)
  }
})

test_that("predict.mfsusie errors when newx is not a matrix or has wrong J", {
  toy <- make_toy_fit()
  expect_error(predict(toy$fit, newx = "bad"),    "must be a numeric matrix")
  expect_error(predict(toy$fit, newx = matrix(0, 5, 999)),
               "but the fit was trained")
})

# ---- coef --------------------------------------------------

test_that("coef.mfsusie returns L x T_m matrices per modality", {
  toy <- make_toy_fit()
  out <- coef(toy$fit)

  expect_identical(length(out), 2L)
  for (m in seq_along(out)) {
    expect_true(is.matrix(out[[m]]))
    expect_identical(nrow(out[[m]]), nrow(toy$fit$alpha))
    expect_identical(ncol(out[[m]]), ncol(toy$Y[[m]]))
    expect_false(any(is.na(out[[m]])))
  }
})

# ---- print --------------------------------------------------

test_that("print.mfsusie does not error and returns invisibly", {
  toy <- make_toy_fit()
  out <- expect_invisible(print(toy$fit))
  expect_identical(out, toy$fit)
})

test_that("print.mfsusie output contains key fit metadata", {
  toy <- make_toy_fit()
  txt <- capture.output(print(toy$fit))
  expect_true(any(grepl("p \\(predictors\\):", txt)))
  expect_true(any(grepl("L \\(effects\\):", txt)))
  expect_true(any(grepl("M \\(outcomes\\):", txt)))
  expect_true(any(grepl("iterations:", txt)))
  expect_true(any(grepl("top PIPs", txt)))
})

# ---- summary --------------------------------------------------

test_that("summary.mfsusie returns a list of class 'summary.mfsusie' with documented fields", {
  toy <- make_toy_fit()
  s <- summary(toy$fit)
  expect_s3_class(s, "summary.mfsusie")
  expect_named(s, c("n_effects", "n_variables", "n_outcomes", "T_basis",
                    "converged", "n_iter", "elbo_final", "pip", "cs"),
               ignore.order = TRUE)
  expect_identical(s$n_variables, ncol(toy$fit$alpha))
  expect_identical(s$n_effects, nrow(toy$fit$alpha))
  expect_identical(s$n_outcomes, 2L)
  expect_true(is.numeric(s$pip))
})

test_that("print(summary.mfsusie) does not error", {
  toy <- make_toy_fit()
  s <- summary(toy$fit)
  out <- expect_invisible(print(s))
  expect_identical(out, s)
})

# ---- predict + fitted on M = 1 ------------------------------

test_that("fitted.mfsusie + predict.mfsusie work for M = 1 single-modality fits", {
  toy <- make_toy_fit(T_per = 64L, L = 3)
  out_fit <- fitted(toy$fit)
  expect_identical(length(out_fit), 1L)
  expect_identical(dim(out_fit[[1]]), c(nrow(toy$X), 64L))

  set.seed(7)
  newx <- matrix(rnorm(3 * ncol(toy$X)), nrow = 3)
  out_pred <- predict(toy$fit, newx)
  expect_identical(dim(out_pred[[1]]), c(3L, 64L))
})
