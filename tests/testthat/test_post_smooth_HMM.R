# Tests for the HMM post-smoother. Two regimes:
#
#  1. Non-manifestation cases (idx_comp = 1..K_full): assert bit-
#     identity with `fsusieR::fit_hmm` at tolerance = 0. Our port-
#     source fix (subsetting `mu` alongside `prob` and `P`) leaves
#     these cases numerically untouched because mu[k] == mu[idx_comp[k]]
#     when idx_comp is the contiguous prefix.
#
#  2. Manifestation cases (idx_comp non-contiguous): assert that
#     our smoother corresponds to the corrected HMM model
#     (state k <-> grid index idx_comp[k]). Upstream fsusieR has a
#     known regression here (commits 9f89333 -> fc806a5) and is
#     expected to differ. See `inst/notes/refactor-exceptions.md`
#     entry "HMM-mu-subset".
#
# Pattern A per refactor-discipline.md section 3.

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
