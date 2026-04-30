# Per-variant credible-level lfsr stored in the smoother payload
# alongside `lfsr_curves`.

test_that("clfsr_curves is populated with the per-variant lfsr formula", {
  set.seed(1L)
  N <- 50L; J <- 30L; T_basis <- 8L
  X    <- matrix(rnorm(N * J), N, J)
  beta <- matrix(0, J, T_basis); beta[5L, ] <- 1.0; beta[15L, ] <- 0.7
  Y1   <- X %*% beta + matrix(rnorm(N * T_basis, sd = 0.3), N, T_basis)
  Y2   <- X %*% beta + matrix(rnorm(N * T_basis, sd = 0.3), N, T_basis)
  fit  <- mfsusieR::mfsusie(X = X, Y = list(Y1, Y2),
                            L = 4L, max_iter = 10L, verbose = FALSE)

  fit <- mfsusieR::mf_post_smooth(fit, method = "scalewise")
  smoothed <- fit$smoothed[["scalewise"]]
  L <- nrow(fit$alpha)
  M <- length(fit$mu[[1L]])

  expect_length(smoothed$clfsr_curves, M)
  for (m in seq_len(M)) {
    expect_length(smoothed$clfsr_curves[[m]], L)
    for (l in seq_len(L)) {
      mu_lm <- fit$mu[[l]][[m]]
      sd_lm <- sqrt(pmax(fit$mu2[[l]][[m]] - mu_lm^2, 0))
      expected <- mfsusieR:::lfsr_from_gaussian(mu_lm, sd_lm)
      expect_equal(smoothed$clfsr_curves[[m]][[l]], expected, tolerance = 0)
      expect_equal(dim(smoothed$clfsr_curves[[m]][[l]]), dim(mu_lm))
      expect_true(all(smoothed$clfsr_curves[[m]][[l]] >= 0))
      expect_true(all(smoothed$clfsr_curves[[m]][[l]] <= 0.5 + 1e-12))
    }
  }
})
