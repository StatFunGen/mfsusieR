# Apple-to-apple cpp11 vs pure-R oracle bit-equivalence for the
# HMM smoother kernels and the full `mf_fit_hmm` driver.
# Tolerance <= 1e-12 per design.md D14.
#
# The kernel-level tests (forward / backward / xi) lock the cpp11
# implementations against the pure-R oracles in
# `R/reference_implementations.R`. The end-to-end test asserts
# that the production `mf_fit_hmm` (cpp11-backed forward / backward /
# xi loops, R-side ash refits) matches the pure-R reference
# `mf_fit_hmm_R` at <= 1e-12.

random_transition_matrix <- function(K, epsilon = 1e-2) {
  P <- matrix(0, K, K); diag(P) <- 0.5
  if (K > 1L) {
    P[1L, 2L:K] <- 0.5 / (K - 1L)
    P[2L:K, 1L] <- 0.5
  }
  P <- P + matrix(epsilon, K, K)
  if (K > 1L) {
    for (i in 2L:K) for (j in 2L:K) if (i != j) P[i, j] <- 0
  }
  P / rowSums(P)
}

random_emission_matrix <- function(T_pos, K) {
  matrix(runif(T_pos * K, 0.05, 1.0), T_pos, K)
}

# ---- kernel tests ----------------------------------------------------------

test_that("hmm_forward_cpp matches hmm_forward_R at <= 1e-12 (t1_normalize = TRUE)", {
  set.seed(mfsusier_test_seed())
  for (K in c(3L, 8L, 20L)) {
    for (T_pos in c(16L, 64L)) {
      emit <- random_emission_matrix(T_pos, K)
      P    <- random_transition_matrix(K)
      pi_v <- runif(K); pi_v <- pi_v / sum(pi_v)
      o_cpp <- mfsusieR:::hmm_forward_cpp(emit, P, pi_v, t1_normalize = TRUE)
      o_R   <- mfsusieR:::hmm_forward_R(  emit, P, pi_v, t1_normalize = TRUE)
      info_str <- sprintf("K=%d T=%d", K, T_pos)
      expect_equal(o_cpp$alpha_hat, o_R$alpha_hat, tolerance = 1e-12,
                   info = info_str)
      expect_equal(o_cpp$G_t,       o_R$G_t,       tolerance = 1e-12,
                   info = info_str)
    }
  }
})

test_that("hmm_forward_cpp matches hmm_forward_R at <= 1e-12 (t1_normalize = FALSE)", {
  # The pre-EM forward pass leaves alpha_hat[1, ] un-normalized
  # and G_t[1] = NA. Verify the cpp11 kernel matches the R oracle
  # in this regime, and that G_t[1] is NA in both.
  set.seed(mfsusier_test_seed())
  K     <- 6L
  T_pos <- 32L
  emit  <- random_emission_matrix(T_pos, K)
  P     <- random_transition_matrix(K)
  pi_v  <- runif(K); pi_v <- pi_v / sum(pi_v)
  o_cpp <- mfsusieR:::hmm_forward_cpp(emit, P, pi_v, t1_normalize = FALSE)
  o_R   <- mfsusieR:::hmm_forward_R(  emit, P, pi_v, t1_normalize = FALSE)
  expect_true(is.na(o_cpp$G_t[1L]))
  expect_true(is.na(o_R$G_t[1L]))
  expect_equal(o_cpp$alpha_hat,             o_R$alpha_hat,             tolerance = 1e-12)
  expect_equal(o_cpp$G_t[-1L],              o_R$G_t[-1L],              tolerance = 1e-12)
})

test_that("hmm_backward_cpp matches hmm_backward_R at <= 1e-12", {
  set.seed(mfsusier_test_seed())
  for (K in c(3L, 8L, 20L)) {
    for (T_pos in c(16L, 64L)) {
      emit <- random_emission_matrix(T_pos, K)
      P    <- random_transition_matrix(K)
      o_cpp <- mfsusieR:::hmm_backward_cpp(emit, P)
      o_R   <- mfsusieR:::hmm_backward_R(  emit, P)
      info_str <- sprintf("K=%d T=%d", K, T_pos)
      expect_equal(o_cpp$beta_hat,     o_R$beta_hat,     tolerance = 1e-12,
                   info = info_str)
      expect_equal(o_cpp$C_t[-T_pos],  o_R$C_t[-T_pos],  tolerance = 1e-12,
                   info = info_str)
      expect_true(is.na(o_cpp$C_t[T_pos]))
      expect_true(is.na(o_R$C_t[T_pos]))
    }
  }
})

test_that("hmm_xi_cpp matches hmm_xi_R at <= 1e-12", {
  set.seed(mfsusier_test_seed())
  for (K in c(3L, 8L, 20L)) {
    for (T_pos in c(16L, 64L)) {
      emit <- random_emission_matrix(T_pos, K)
      P    <- random_transition_matrix(K)
      pi_v <- runif(K); pi_v <- pi_v / sum(pi_v)
      fwd  <- mfsusieR:::hmm_forward_R( emit, P, pi_v, t1_normalize = TRUE)
      bwd  <- mfsusieR:::hmm_backward_R(emit, P)
      xi_cpp <- mfsusieR:::hmm_xi_cpp(fwd$alpha_hat, bwd$beta_hat, emit, P)
      xi_R   <- mfsusieR:::hmm_xi_R(  fwd$alpha_hat, bwd$beta_hat, emit, P)
      expect_equal(xi_cpp, xi_R, tolerance = 1e-12,
                   info = sprintf("K=%d T=%d", K, T_pos))
    }
  }
})

# ---- end-to-end driver test ------------------------------------------------

test_that("mf_fit_hmm (cpp11) matches mf_fit_hmm_R at <= 1e-12 across configs", {
  testthat::skip_if_not_installed("ashr")
  cases <- expand.grid(halfK = c(5L, 10L, 20L),
                       T_pos = c(32L, 64L, 128L),
                       seed  = c(1L, 17L, 42L))
  for (i in seq_len(nrow(cases))) {
    set.seed(cases$seed[i])
    x  <- rnorm(cases$T_pos[i], sd = 0.4)
    sd <- runif(cases$T_pos[i], 0.1, 0.5)
    o_cpp <- mfsusieR:::mf_fit_hmm(  x, sd, halfK = cases$halfK[i])
    o_R   <- mfsusieR:::mf_fit_hmm_R(x, sd, halfK = cases$halfK[i])
    info_str <- sprintf("halfK=%d T=%d seed=%d",
                        cases$halfK[i], cases$T_pos[i], cases$seed[i])
    expect_equal(o_cpp$x_post,  o_R$x_post,  tolerance = 1e-12, info = info_str)
    expect_equal(o_cpp$lfsr,    o_R$lfsr,    tolerance = 1e-12, info = info_str)
    expect_equal(o_cpp$log_BF,  o_R$log_BF,  tolerance = 1e-12, info = info_str)
    expect_equal(o_cpp$mu,      o_R$mu,      tolerance = 1e-12, info = info_str)
  }
})

test_that("mf_fit_hmm (cpp11) matches mf_fit_hmm_R on the bimodal manifestation case", {
  # Same construction as the Pattern A test in
  # `test_post_smooth_HMM.R`. Verifies that the cpp11 production
  # path exercises the `mu <- mu[idx_comp]` fix consistently with
  # the pure-R reference even when idx_comp is non-contiguous.
  testthat::skip_if_not_installed("ashr")
  set.seed(42)
  T_pos <- 200L
  x  <- c(rnorm(50, mean =  1.5, sd = 0.1),
          rnorm(50, mean = -1.5, sd = 0.1),
          rnorm(50, mean =  0,   sd = 0.1),
          rnorm(50, mean =  0.05,sd = 0.1))
  sd <- rep(0.1, T_pos)
  o_cpp <- mfsusieR:::mf_fit_hmm(  x, sd, halfK = 20L)
  o_R   <- mfsusieR:::mf_fit_hmm_R(x, sd, halfK = 20L)
  expect_equal(o_cpp$x_post, o_R$x_post, tolerance = 1e-12)
  expect_equal(o_cpp$lfsr,   o_R$lfsr,   tolerance = 1e-12)
  expect_equal(o_cpp$log_BF, o_R$log_BF, tolerance = 1e-12)
})
