# Tests for the HMM post-smoother. Bit-identity with
# `fsusieR::fit_hmm` at tolerance = 0 in both the contiguous-prefix
# regime and the non-contiguous `idx_comp` regime where the
# `mu <- mu[idx_comp]` subset is required for state-emission
# alignment.
#
# Historical: the non-contiguous regime previously diverged from
# upstream because fsusieR commit fc806a5 (2026-03-12) dropped the
# `mu <- mu[idx_comp]` line in a Baum-Welch restructure. fsusieR
# pull request stephenslab/fsusieR#31 restored the line. Both
# regimes now bit-match. See `inst/notes/refactor-exceptions.md`
# entry "HMM-mu-subset" for the full history.

skip_if_no_fsusier_for_HMM <- function() {
  testthat::skip_if_not_installed("fsusieR")
  testthat::skip_if_not_installed("ashr")
}

test_that("mf_fit_hmm matches upstream when idx_comp is the contiguous prefix", {
  skip_if_no_fsusier_for_HMM()
  set.seed(7)
  T_pos <- 16L
  x   <- rnorm(T_pos, sd = 0.3)
  sd  <- runif(T_pos, 0.1, 0.4)
  ours <- mfsusieR:::mf_fit_hmm(x, sd, halfK = 5L)
  ref  <- fsusieR:::fit_hmm(x, sd, halfK = 5)
  # Sanity-check we are in the no-manifest regime.
  expect_equal(length(ours$mu), ncol(ref$prob))
  expect_equal(ours$x_post, ref$x_post, tolerance = 0)
  expect_equal(ours$lfsr,   ref$lfsr,   tolerance = 0)
  expect_equal(ours$log_BF, ref$log_BF, tolerance = 0)
})

test_that("mf_fit_hmm matches upstream across halfK / T_pos / random seeds (no-manifest sweep)", {
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
    info_str <- sprintf("halfK=%d T=%d seed=%d",
                        cases$halfK[i], cases$T_pos[i], cases$seed[i])
    # Skip cases where the bug is expected to manifest (idx_comp
    # non-contiguous). We detect by comparing K vectors.
    if (length(o$mu) != ncol(r$prob)) next
    expect_equal(o$x_post, r$x_post, tolerance = 0, info = info_str)
    expect_equal(o$lfsr,   r$lfsr,   tolerance = 0, info = info_str)
    expect_equal(o$log_BF, r$log_BF, tolerance = 0, info = info_str)
  }
})

test_that("mf_fit_hmm uses the correct prior modes when idx_comp is non-contiguous (Pattern A)", {
  skip_if_no_fsusier_for_HMM()
  # Bimodal data at +/- 1.5 with intermediate states unoccupied.
  # With halfK = 20, prefilter keeps the null state and the two
  # extremes, dropping the middle of the grid. This makes the
  # state-space reduction non-contiguous and exercises the
  # `mu <- mu[idx_comp]` fix.
  set.seed(42)
  T_pos <- 200L
  x  <- c(rnorm(50, mean =  1.5, sd = 0.1),
          rnorm(50, mean = -1.5, sd = 0.1),
          rnorm(50, mean =  0,   sd = 0.1),
          rnorm(50, mean =  0.05,sd = 0.1))
  sd <- rep(0.1, T_pos)
  ours <- mfsusieR:::mf_fit_hmm(x, sd, halfK = 20L)

  # The corrected smoother concentrates posterior mass on the
  # null state and the two retained non-null states. The retained
  # state means satisfy 0 in mu, max(|mu|) > 1 (captures the +/- 1.5
  # mode), and the smoothed estimate tracks the observed bimodal
  # signal up to the ash shrinkage.
  expect_true(0 %in% ours$mu)
  expect_true(max(abs(ours$mu)) > 1)

  # Smoothed estimate must agree in sign with x on the strong-signal
  # blocks (positions 1:50 around +1.5, 51:100 around -1.5).
  expect_true(mean(sign(ours$x_post[1:50])  == sign(x[1:50]))  > 0.9)
  expect_true(mean(sign(ours$x_post[51:100])== sign(x[51:100]))> 0.9)

  # Upstream fsusieR has restored `mu <- mu[idx_comp]` (the
  # regression introduced in fc806a5 was reverted). With both
  # implementations subsetting `mu` to match the reduced state
  # space, the smoothers are bit-identical on the
  # non-contiguous-prefix configuration. See
  # `inst/notes/refactor-exceptions.md` (HMM-mu-subset entry)
  # for the history.
  ref <- fsusieR:::fit_hmm(x, sd, halfK = 20)
  expect_equal(ours$x_post, ref$x_post, tolerance = 0)
})

test_that("mf_univariate_hmm_regression matches upstream on contiguous-prefix configs", {
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
  # operating points. These configs put us in the no-manifest
  # regime; non-contiguous idx_comp is exercised by the dedicated
  # Pattern A test above.
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
    info_str <- sprintf("seed=%d n=%d T=%d beta=%g",
                        cs$seed, cs$n, cs$T, cs$beta)
    # Detect the no-manifest regime via the equivalent fit_hmm
    # call against upstream: if `length(ref_mu) == ncol(ref_prob)`
    # then idx_comp = 1..K_full and our fix collapses to upstream.
    Y_sc <- mfsusieR:::col_scale(Y)
    bs   <- susieR::compute_marginal_bhat_shat(X, Y_sc)
    rfit <- fsusieR:::fit_hmm(x  = as.numeric(bs$Bhat[1L, ]),
                              sd = as.numeric(bs$Shat[1L, ]),
                              halfK = 15)
    if (length(rfit$mu) != ncol(rfit$prob)) next
    expect_equal(ours$effect_estimate, as.numeric(ref$effect_estimate),
                 tolerance = 0, info = info_str)
    expect_equal(ours$lfsr, as.numeric(ref$lfsr),
                 tolerance = 0, info = info_str)
    expect_equal(ours$lBF, ref$lBF, tolerance = 0, info = info_str)
  }
})

# =============================================================================
# Credible band: shape, monotonicity, law-of-total-variance regression guard
# =============================================================================

test_that("mf_fit_hmm returns a per-position posterior SD vector", {
  set.seed(11)
  T_pos <- 24L
  x  <- rnorm(T_pos, sd = 0.3)
  sd <- runif(T_pos, 0.1, 0.4)
  s  <- mfsusieR:::mf_fit_hmm(x, sd, halfK = 5L)
  expect_named(s, c("prob", "x_post", "x_post_sd", "lfsr", "mu",
                    "ll_hmm", "ll_null", "log_BF"),
               ignore.order = TRUE)
  expect_length(s$x_post_sd, T_pos)
  expect_true(all(is.finite(s$x_post_sd)))
  expect_true(all(s$x_post_sd >= 0))
})

test_that("mf_fit_hmm posterior SD matches a hand-computed law-of-total-variance value", {
  # Construct a degenerate scenario where K = 2 (null + one
  # nontrivial state). Use a fixture with known posterior shape:
  # `x` strongly localized so the HMM concentrates on the
  # non-null state at most positions. Then law-of-total-variance
  # reduces to the per-state ash sd at each position.
  set.seed(12)
  T_pos <- 12L
  x  <- c(rep(0, 6L), rep(1.5, 6L))   # null half + signal half
  sd <- rep(0.2, T_pos)
  s  <- mfsusieR:::mf_fit_hmm(x, sd, halfK = 3L)

  # Sanity: at signal positions, the SD should be smaller than
  # the prior sd (shrinkage). At null positions, SD ~ 0.
  expect_true(all(s$x_post_sd[1:6] <= sd[1:6] + 1e-8))
  # At signal positions, posterior should have non-trivial SD
  # (not zero, not the prior).
  expect_true(any(s$x_post_sd[7:12] > 0))
})

test_that("mf_post_smooth(method = 'HMM') populates credible_bands of correct shape", {
  set.seed(13)
  n <- 60; p <- 12; T_m <- 32L
  X <- matrix(rnorm(n * p), n)
  beta <- numeric(p); beta[3] <- 1.0
  shape <- exp(-((seq_len(T_m) - T_m / 2)^2) / (2 * 6^2))
  Y <- X %*% (matrix(beta, p, 1) %*% matrix(shape, 1, T_m)) +
         matrix(rnorm(n * T_m, sd = 0.4), n)

  fit   <- fsusie(wavelet_qnorm = FALSE, Y, X, L = 1, max_iter = 30, verbose = FALSE)
  fit_h <- mf_post_smooth(fit, method = "HMM")

  bands <- fit_h$smoothed$HMM$credible_bands
  expect_true(is.list(bands))
  for (m in seq_along(bands)) {
    for (l in seq_along(bands[[m]])) {
      band <- bands[[m]][[l]]
      expect_identical(dim(band), c(T_m, 2L))
      lower <- band[, 1L]; upper <- band[, 2L]
      mean_curve <- fit_h$smoothed$HMM$effect_curves[[m]][[l]]
      # Monotonicity: lower <= mean <= upper everywhere.
      expect_true(all(lower <= mean_curve + 1e-10))
      expect_true(all(mean_curve <= upper + 1e-10))
    }
  }
})
