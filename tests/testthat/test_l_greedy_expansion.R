# Regression test for the L_greedy warm-start expansion path
# (audit finding C-8). Compares a greedy fit (L = 10, starts at
# L_greedy = 5, expands by 5 to L = 10) against a cold-start fit
# at L = 10. The two paths should agree on the recovered PIPs and
# credible sets up to a tolerance commensurate with the warm-start
# numerics.

test_that("L_greedy expansion converges to a fit consistent with cold-start L = 10", {
  set.seed(11L)
  N <- 100L; J <- 60L; T_basis <- 16L
  X <- matrix(rnorm(N * J), N, J)
  beta <- matrix(0, J, T_basis)
  beta[5L, ]  <- 1.0
  beta[20L, ] <- 0.7
  beta[40L, ] <- 0.5
  Y <- X %*% beta + matrix(rnorm(N * T_basis, sd = 0.4), N, T_basis)

  fit_greedy <- mfsusieR::mfsusie(X = X, Y = list(Y),
                                  L = 10L, L_greedy = 5L,
                                  max_iter = 30L, verbose = FALSE)
  fit_cold   <- mfsusieR::mfsusie(X = X, Y = list(Y),
                                  L = 10L, L_greedy = NULL,
                                  max_iter = 30L, verbose = FALSE)

  # Top-PIP variants identified by both paths.
  pip_threshold <- 0.5
  signal_g <- which(fit_greedy$pip > pip_threshold)
  signal_c <- which(fit_cold$pip   > pip_threshold)
  expect_equal(signal_g, signal_c)

  # Per-variant PIPs: agreement on the dominant variants, with a
  # tolerance reflecting that the greedy expansion path adds
  # effects from a warm-start seed pi (audit finding C-8).
  expect_equal(fit_greedy$pip, fit_cold$pip, tolerance = 1e-3)
})
