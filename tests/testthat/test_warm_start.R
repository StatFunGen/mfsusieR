# Reference tests for the `model_init` warm-start argument on
# `mfsusie()` and `fsusie()`.
#
# 1. Convergence-shortcut: warm-starting from a converged fit
#    converges in a single iteration.
# 2. End-to-end: a warm-started fit on the same data matches a
#    cold fit at the contract tolerance.
# 3. Error path: model_init with mismatched L errors cleanly.

build_warm_sim <- function(seed = 17L, n = 100L, T_m = 64L,
                           p = 30L) {
  set.seed(seed)
  X    <- matrix(rnorm(n * p), n, p)
  # Two effects so cold-start takes >= 3 IBSS iters and the
  # warm-start contrast (warm < cold) is observable.
  beta <- rep(0, p); beta[5L] <- 1; beta[15L] <- 0.7
  Y    <- X %*% matrix(beta, ncol = 1L) %*% matrix(1, 1, T_m) +
          matrix(rnorm(n * T_m, sd = 0.5), n, T_m)
  list(X = X, Y = Y, p = p, T_m = T_m, n = n)
}

# --- 1. Convergence shortcut -------------------------------------

test_that("warm-started fit from a converged fit hits the convergence floor", {
  # The IBSS convergence check requires `abs(elbo[iter+1] -
  # elbo[iter]) < tol` and elbo[1] = -Inf is the sentinel, so
  # `niter = 2` is the absolute floor (one real ELBO comparison
  # against the seeded posterior). A warm-started fit from a
  # converged cold fit hits that floor; the cold fit takes
  # strictly more iterations.
  sim <- build_warm_sim()
  fit_cold <- fsusie(wavelet_qnorm = FALSE, wavelet_standardize = FALSE, max_iter = 100, sim$Y, sim$X, L = 5,
                     verbose = FALSE)
  expect_true(fit_cold$converged)
  # Cold-fit floor: at least 3 iters (PIP-difference convergence is
  # the default and on this fixture takes 3 iters; the prior ELBO-
  # default took 4+).
  expect_gte(fit_cold$niter, 3L)

  fit_warm <- fsusie(wavelet_qnorm = FALSE, wavelet_standardize = FALSE, max_iter = 100, sim$Y, sim$X, L = 5,
                     model_init = fit_cold,
                     verbose = FALSE)
  expect_true(fit_warm$converged)
  expect_equal(fit_warm$niter, 2L)
  expect_lt(fit_warm$niter, fit_cold$niter)
})

# --- 2. End-to-end agreement -------------------------------------

test_that("warm-started fit reaches the same posterior as a cold fit", {
  sim <- build_warm_sim()
  fit_cold <- fsusie(wavelet_qnorm = FALSE, wavelet_standardize = FALSE, max_iter = 100, sim$Y, sim$X, L = 5,
                     verbose = FALSE)
  fit_warm <- fsusie(wavelet_qnorm = FALSE, wavelet_standardize = FALSE, max_iter = 100, sim$Y, sim$X, L = 5,
                     model_init = fit_cold,
                     verbose = FALSE)
  # The warm fit re-runs at least one IBSS iteration and may
  # nudge the posterior by a few units of `tol`. Compare PIPs
  # within a tolerance commensurate with the IBSS tolerance.
  expect_equal(fit_warm$pip, fit_cold$pip, tolerance = 1e-3)
  expect_equal(length(fit_warm$sets$cs), length(fit_cold$sets$cs))
})

# --- 3. Mid-run resume ------------------------------------------

test_that("capping max_iter then resuming via model_init reaches convergence", {
  sim <- build_warm_sim()
  fit_partial <- fsusie(wavelet_qnorm = FALSE, wavelet_standardize = FALSE, sim$Y, sim$X, L = 5,
                        max_iter = 2,
                        verbose = FALSE)
  expect_false(isTRUE(fit_partial$converged))

  fit_resumed <- fsusie(wavelet_qnorm = FALSE, wavelet_standardize = FALSE, max_iter = 100, sim$Y, sim$X, L = 5,
                        model_init = fit_partial,
                        verbose = FALSE)
  expect_true(fit_resumed$converged)
})

# --- 4. L expansion ---------------------------------------------

test_that("model_init with smaller L is expanded to the requested L", {
  sim <- build_warm_sim()
  fit_small <- fsusie(wavelet_qnorm = FALSE, wavelet_standardize = FALSE, max_iter = 100, sim$Y, sim$X, L = 3,
                      verbose = FALSE)
  fit_grown <- fsusie(wavelet_qnorm = FALSE, wavelet_standardize = FALSE, max_iter = 100, sim$Y, sim$X, L = 5,
                      model_init = fit_small,
                      verbose = FALSE)
  expect_equal(nrow(fit_grown$alpha), 5L)
  expect_true(fit_grown$converged)
})

# --- 5. Error path: shrinking is not supported ------------------

test_that("model_init with larger L than requested errors cleanly", {
  sim <- build_warm_sim()
  fit_cold <- fsusie(wavelet_qnorm = FALSE, wavelet_standardize = FALSE, max_iter = 100, sim$Y, sim$X, L = 5,
                     verbose = FALSE)
  expect_error(
    fsusie(wavelet_qnorm = FALSE, wavelet_standardize = FALSE, max_iter = 100, sim$Y, sim$X, L = 3, model_init = fit_cold,
           verbose = FALSE),
    "shrinking is not supported"
  )
})
