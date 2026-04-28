# End-to-end smoke test for mfsusie(): does it run, return the
# documented shape, satisfy basic invariants?
#
# Per design.md D11d this is an "apple-to-orange" smoke test: no
# numerical comparison against any reference, just contract checks.
# Apple-to-apple fidelity tests against susieR (C1), fsusieR (C2),
# and mvf.susie.alpha (C3) live in PR groups 7c / 7d / 7e and depend
# on the public API landing first.

# ---- toy fixture --------------------------------------------------

make_toy_smoke <- function(n = 30, J = 8, T_per_outcome = c(64L, 32L),
                           seed = mfsusier_test_seed()) {
  set.seed(seed)
  X <- matrix(rnorm(n * J), nrow = n)
  beta <- numeric(J); beta[1] <- 1; beta[3] <- -0.5
  Y <- lapply(T_per_outcome, function(T_m) {
    eta <- X %*% beta
    matrix(rep(eta, T_m), nrow = n) + matrix(rnorm(n * T_m, sd = 0.3), nrow = n)
  })
  list(X = X, Y = Y, beta_true = beta)
}

# ---- Smoke test: fit returns the documented shape ----------------

test_that("mfsusie() runs end-to-end on a toy fixture and returns a fit of class c('mfsusie', 'susie')", {
  toy <- make_toy_smoke()

  fit <- mfsusie(toy$X, toy$Y, L = 5, max_iter = 30, verbose = FALSE)

  expect_s3_class(fit, c("mfsusie", "susie"), exact = FALSE)
  expect_true(is.list(fit))

  # Documented fit fields.
  expect_true(is.matrix(fit$alpha))
  expect_identical(ncol(fit$alpha), 8L)            # J = 8
  expect_true(nrow(fit$alpha) >= 1L)
  expect_true(nrow(fit$alpha) <= 5L)               # L max

  # mu / mu2: list[L] of list[M] of J x T_basis matrices.
  expect_true(is.list(fit$mu))
  expect_true(is.list(fit$mu2))
  for (l in seq_along(fit$mu)) {
    expect_true(is.list(fit$mu[[l]]))
    expect_identical(length(fit$mu[[l]]), 2L)      # M = 2
    for (m in seq_along(fit$mu[[l]])) {
      expect_identical(nrow(fit$mu[[l]][[m]]), 8L)
    }
  }

  # PIP, CS, ELBO, niter, converged.
  expect_true(is.numeric(fit$pip))
  expect_identical(length(fit$pip), 8L)
  expect_true(all(fit$pip >= 0 & fit$pip <= 1))

  expect_true(is.numeric(fit$elbo))
  expect_true(length(fit$elbo) >= 1L)

  expect_true(is.numeric(fit$niter) && length(fit$niter) == 1L)
  expect_true(is.logical(fit$converged) && length(fit$converged) == 1L)
})

# ---- ELBO monotonicity (per design.md D11d smoke checks) ---------

test_that("mfsusie() ELBO is non-decreasing across iterations (within numerical tolerance)", {
  toy <- make_toy_smoke()
  fit <- mfsusie(toy$X, toy$Y, L = 5, max_iter = 50, verbose = FALSE)

  # Skip the first ELBO diff: the initial sigma2 (var(Y) guess) gets
  # replaced by its first closed-form update, which is a one-shot
  # non-monotone step. Post-iter-1 the trajectory must be monotone
  # within numerical tolerance (mvsusieR pattern).
  if (length(fit$elbo) > 2L) {
    expect_true(all(diff(fit$elbo)[-1] >= -1e-6),
                info = paste("ELBO not monotone after iter 1:",
                             paste(round(diff(fit$elbo), 6), collapse = " ")))
  }
})

# ---- No NA / NaN in fit fields -----------------------------------

test_that("mfsusie() fit fields contain no NaN / NA in alpha, mu, mu2, pip, elbo", {
  toy <- make_toy_smoke()
  fit <- mfsusie(toy$X, toy$Y, L = 5, max_iter = 30, verbose = FALSE)

  expect_false(any(is.na(fit$alpha)))
  expect_false(any(is.nan(fit$alpha)))
  for (l in seq_along(fit$mu)) {
    for (m in seq_along(fit$mu[[l]])) {
      expect_false(any(is.na(fit$mu[[l]][[m]])))
      expect_false(any(is.na(fit$mu2[[l]][[m]])))
    }
  }
  expect_false(any(is.na(fit$pip)))
  expect_false(any(is.na(fit$elbo)))
})

# ---- Y_grid + X_eff are attached by default --------------------

test_that("mfsusie() attaches Y_grid + X_eff by default and can opt out", {
  toy <- make_toy_smoke()
  fit <- mfsusie(toy$X, toy$Y, L = 5, max_iter = 30, verbose = FALSE)
  expect_true(!is.null(fit$Y_grid))
  expect_identical(length(fit$Y_grid), length(toy$Y))
  for (m in seq_along(fit$Y_grid)) {
    expect_true(is.matrix(fit$Y_grid[[m]]))
    expect_identical(nrow(fit$Y_grid[[m]]), nrow(toy$X))
  }
  expect_identical(length(fit$X_eff), nrow(fit$alpha))
  for (l in seq_along(fit$X_eff)) {
    expect_identical(length(fit$X_eff[[l]]), nrow(toy$X))
  }

  # Opt-out drops both fields.
  fit_lite <- mfsusie(toy$X, toy$Y, L = 5, max_iter = 30,
                      attach_smoothing_inputs = FALSE,
                      verbose = FALSE)
  expect_null(fit_lite$Y_grid)
  expect_null(fit_lite$X_eff)
})

# ---- M = 1 single-modality functional path --------------------

test_that("mfsusie() works with M = 1 single-modality input", {
  toy <- make_toy_smoke(T_per_outcome = 64L)
  fit <- mfsusie(toy$X, toy$Y, L = 3, max_iter = 30, verbose = FALSE)
  expect_s3_class(fit, "mfsusie")
  expect_identical(length(fit$mu[[1]]), 1L)
  expect_true(is.numeric(fit$pip))
})

# ---- L_greedy passthrough --------------------

test_that("mfsusie(L_greedy = K) runs without error and respects the saturation bound", {
  toy <- make_toy_smoke()
  fit <- mfsusie(toy$X, toy$Y, L = 8, L_greedy = 2, greedy_lbf_cutoff = 0.1,
                 max_iter = 30, verbose = FALSE)
  # The greedy outer loop trims L; final number of effects must be in [L_greedy, L].
  expect_true(nrow(fit$alpha) >= 2L)
  expect_true(nrow(fit$alpha) <= 8L)
})
