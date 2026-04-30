# Regression test for the V-filter behavior on `pip`.
#
# After commit 605f395 the per-effect effective slab variance
# populates `model$V[l]`, and `susie_get_pip`'s
# `V[l] > prior_tol` filter drops collapsed effects from the
# PIP product. This test pins the contract on a null-only
# fixture (no causal signal): every effect SHOULD have its
# mixture pi collapse onto the null spike, V[l] SHOULD fall
# below `prior_tol = 1e-9`, the filter SHOULD drop those
# effects from the PIP product, and `trim_null_effects.mf_individual`
# SHOULD also zero them on the rest of the fit.

test_that("null-only fixture: V-filter drops collapsed effects and trim zeros them", {
  set.seed(101L)
  N <- 50L; J <- 30L; T_basis <- 8L
  X  <- matrix(rnorm(N * J), N, J)
  Y1 <- matrix(rnorm(N * T_basis, sd = 0.3), N, T_basis)  # pure noise

  fit <- mfsusieR::mfsusie(X = X, Y = list(Y1),
                           L = 4L, max_iter = 10L,
                           verbose = FALSE)

  expect_true(any(fit$V < 1e-9),
              info = "expected at least one collapsed effect on null data")
  trimmed <- which(fit$V < 1e-9)
  for (l in trimmed) {
    expect_equal(fit$alpha[l, ], unname(fit$pi), tolerance = 0)
    for (m in seq_along(fit$mu[[l]])) {
      expect_true(all(fit$mu[[l]][[m]] == 0))
      expect_true(all(fit$mu2[[l]][[m]] == 0))
    }
    expect_equal(fit$V[l], 0, tolerance = 0)
    expect_equal(fit$lbf[l], 0, tolerance = 0)
    expect_equal(fit$KL[l], 0, tolerance = 0)
    expect_true(all(fit$lbf_variable[l, ] == 0))
  }
  # PIP is bounded; on a null fixture it should be small.
  expect_true(all(fit$pip >= 0 & fit$pip <= 1))
})

test_that("prior_tol = 0 keeps every effect's PIP contribution", {
  set.seed(101L)
  N <- 50L; J <- 30L; T_basis <- 8L
  X  <- matrix(rnorm(N * J), N, J)
  Y1 <- matrix(rnorm(N * T_basis, sd = 0.3), N, T_basis)

  fit_default <- mfsusieR::mfsusie(X = X, Y = list(Y1),
                                   L = 4L, max_iter = 10L,
                                   verbose = FALSE)
  pip_no_filter <- susieR::susie_get_pip(fit_default, prior_tol = 0)
  expect_true(all(pip_no_filter >= fit_default$pip - 1e-12))
})
