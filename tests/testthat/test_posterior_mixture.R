# Apple-to-apple comparison against:
#   mfsusieR pure-R oracles (R/reference_implementations.R) at <= 1e-12,
#   ashr::calc_logLR / ashr::postmean / ashr::postsd at <= 1e-12.
#
# Per design.md D14 (cpp11 kernel vs pure-R oracle apple-to-apple
# tolerance <= 1e-12).
#
# Manuscript references for the closed-form mixture-of-normals math:
#   methods/derivation.tex eq:post_f_mix
#   methods/derivation.tex eq:post_f2_mix

# ---- cpp11 vs pure-R oracle ------------------------------------------------

test_that("mixture_log_bf_per_scale: cpp11 matches pure-R oracle at <= 1e-12", {
  set.seed(mfsusier_test_seed())
  J <- 50; Ti <- 16; K <- 5
  B <- matrix(rnorm(J * Ti), nrow = J)
  S <- matrix(runif(J * Ti, 0.5, 1.5), nrow = J)
  sd_g <- c(0, sort(runif(K - 1, 0.1, 2)))
  pi_g <- c(0.6, rep(0.4 / (K - 1), K - 1))

  out_cpp <- mfsusieR:::mixture_log_bf_per_scale(B, S, sd_g, pi_g, V_scale = 1.5)
  out_R   <- mfsusieR:::mixture_log_bf_per_scale_R(B, S, sd_g, pi_g, V_scale = 1.5)
  expect_equal(out_cpp, out_R, tolerance = 1e-12)
})

test_that("mixture_posterior_per_scale: cpp11 matches pure-R oracle at <= 1e-12", {
  set.seed(mfsusier_test_seed())
  J <- 50; Ti <- 16; K <- 5
  B <- matrix(rnorm(J * Ti), nrow = J)
  S <- matrix(runif(J * Ti, 0.5, 1.5), nrow = J)
  sd_g <- c(0, sort(runif(K - 1, 0.1, 2)))
  pi_g <- c(0.6, rep(0.4 / (K - 1), K - 1))

  out_cpp <- mfsusieR:::mixture_posterior_per_scale(B, S, sd_g, pi_g, V_scale = 1.5)
  out_R   <- mfsusieR:::mixture_posterior_per_scale_R(B, S, sd_g, pi_g, V_scale = 1.5)
  expect_equal(out_cpp$pmean,  out_R$pmean,  tolerance = 1e-12)
  expect_equal(out_cpp$pmean2, out_R$pmean2, tolerance = 1e-12)
})

# ---- pure-R oracle vs ashr (sanity for the closed-form formula) -----------

test_that("mixture_log_bf_per_scale_R: closed form matches ashr::calc_logLR at <= 1e-12", {
  skip_if_not_installed("ashr")
  set.seed(mfsusier_test_seed())
  J <- 30; Ti <- 12; K <- 4
  B <- matrix(rnorm(J * Ti), nrow = J)
  S <- matrix(runif(J * Ti, 0.5, 1.5), nrow = J)
  sd_g <- c(0, 0.5, 1.0, 1.5)
  pi_g <- c(0.5, 0.2, 0.2, 0.1)
  V <- 1.7

  g <- ashr::normalmix(pi = pi_g, mean = rep(0, K), sd = sd_g * sqrt(V))
  ashr_lbf <- vapply(seq_len(J), function(j) {
    ashr::calc_logLR(g, ashr::set_data(B[j, ], S[j, ]))
  }, numeric(1))

  ours <- mfsusieR:::mixture_log_bf_per_scale_R(B, S, sd_g, pi_g, V_scale = V)
  expect_equal(ours, ashr_lbf, tolerance = 1e-12)
})

test_that("mixture_posterior_per_scale_R: closed form matches ashr::postmean / postsd at <= 1e-12", {
  skip_if_not_installed("ashr")
  set.seed(mfsusier_test_seed())
  J <- 30; Ti <- 12; K <- 4
  B <- matrix(rnorm(J * Ti), nrow = J)
  S <- matrix(runif(J * Ti, 0.5, 1.5), nrow = J)
  sd_g <- c(0, 0.5, 1.0, 1.5)
  pi_g <- c(0.5, 0.2, 0.2, 0.1)
  V <- 1.7

  g <- ashr::normalmix(pi = pi_g, mean = rep(0, K), sd = sd_g * sqrt(V))
  ashr_pm  <- matrix(0, J, Ti)
  ashr_psd <- matrix(0, J, Ti)
  for (j in seq_len(J)) {
    d <- ashr::set_data(B[j, ], S[j, ])
    ashr_pm[j, ]  <- ashr::postmean(g, d)
    ashr_psd[j, ] <- ashr::postsd(g, d)
  }

  ours <- mfsusieR:::mixture_posterior_per_scale_R(B, S, sd_g, pi_g, V_scale = V)
  expect_equal(ours$pmean,  ashr_pm, tolerance = 1e-12)
  expect_equal(ours$pmean2, ashr_psd^2 + ashr_pm^2, tolerance = 1e-12)
})

# ---- pure-null edge case (K = 1, sd = 0) --------------------------------

test_that("pure-null mixture (K = 1, sd_grid = 0) gives lbf = 0 and zero posterior", {
  # Exercises the `var_k[k] == 0.0` short-circuit branch in cpp11 and
  # the `if (sd_k_var == 0)` branch in the R oracle. The mixture is
  # entirely on the null component so log-BF and posterior moments are
  # exactly zero.
  set.seed(mfsusier_test_seed())
  J <- 12; Ti <- 8
  B <- matrix(rnorm(J * Ti), nrow = J)
  S <- matrix(runif(J * Ti, 0.5, 1.5), nrow = J)

  lbf <- mfsusieR:::mixture_log_bf_per_scale(B, S, sd_grid = 0, pi_grid = 1, V_scale = 1.5)
  expect_equal(lbf, rep(0, J), tolerance = 0)   # exact

  out <- mfsusieR:::mixture_posterior_per_scale(B, S, sd_grid = 0, pi_grid = 1, V_scale = 1.5)
  expect_equal(out$pmean,  matrix(0, J, Ti), tolerance = 0)
  expect_equal(out$pmean2, matrix(0, J, Ti), tolerance = 0)
})

# ---- susieR-degenerate single-component case (C1 sanity) -------------------

test_that("mixture_log_bf_per_scale degenerates to scalar log-BF when K = 1, no null", {
  # Single component, sd = 1, pi = 1; closed form is
  #   log_BF = -0.5 * sum_i [ log(1 + V/Shat^2) - V*Bhat^2 / (Shat^2 * (Shat^2 + V)) ]
  set.seed(mfsusier_test_seed())
  J <- 20; Ti <- 8
  B <- matrix(rnorm(J * Ti), nrow = J)
  S <- matrix(runif(J * Ti, 0.5, 1.5), nrow = J)
  V <- 2

  ours <- mfsusieR:::mixture_log_bf_per_scale(B, S, sd_grid = 1, pi_grid = 1, V_scale = V)
  Shat2 <- S^2
  per_pos <- -0.5 * (log(V + Shat2) - log(Shat2)
                     - B^2 * (1 / Shat2 - 1 / (V + Shat2)))
  expected <- rowSums(per_pos)
  expect_equal(ours, expected, tolerance = 1e-12)
})

test_that("mixture_posterior_per_scale degenerates to scalar Gaussian-conjugate when K = 1, no null", {
  # post_var = V * Shat^2 / (V + Shat^2)
  # post_mean = V / (V + Shat^2) * Bhat
  set.seed(mfsusier_test_seed())
  J <- 20; Ti <- 8
  B <- matrix(rnorm(J * Ti), nrow = J)
  S <- matrix(runif(J * Ti, 0.5, 1.5), nrow = J)
  V <- 2

  ours <- mfsusieR:::mixture_posterior_per_scale(B, S, sd_grid = 1, pi_grid = 1, V_scale = V)
  shrink <- V / (V + S^2)
  pm  <- shrink * B
  pv  <- shrink * S^2
  expect_equal(ours$pmean,  pm, tolerance = 1e-12)
  expect_equal(ours$pmean2, pv + pm^2, tolerance = 1e-12)
})
