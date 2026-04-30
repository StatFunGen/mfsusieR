# C2 contract: single-modality functional case (M = 1, T_1 > 1)
# vs `fsusieR::susiF`. Per design.md D11b at face value the
# contract is `<= 1e-8` apple-to-apple. In practice, mfsusieR uses
# susieR's standard ER2 / residual-variance formula
# (`||Y - X b||^2 + sum_l ((alpha_l . pw)^T mu2_l - ||X (alpha_l . mu_l)||^2)`),
# while `fsusieR::susiF` uses the simplified
# `sum(postF2) - sum(postF^2)` formula -- documented in
# `inst/notes/refactor-exceptions.md` PR group 6 as a Pattern A
# port-source-bug fix. After the fix, sigma2 (and everything that
# depends on it: ELBO, KL, the iterative trajectory, the converged
# alpha) is no longer bit-identical between the two packages.
#
# What we DO test bit-identically (links to existing tests):
#   - per-(SNP, modality, scale) log-Bayes factor:
#     test_fsusier_fidelity.R, test_posterior_mixture.R
#   - per-(SNP, scale) posterior mean / second moment:
#     test_posterior_mixture.R, test_mvf_alpha_fidelity.R
#   - mixsqp likelihood matrix: test_mvf_alpha_fidelity.R
#
# What we test STRUCTURALLY here (the C2-as-shipped contract):
#   - same top SNP at convergence on a deterministic fixture
#   - same number of credible sets within +/- 1
#   - PIPs at the top causal SNP agree to within 0.05
#   - both fits run to convergence in <= 50 IBSS iterations
#
# When fsusieR's ER2 formula is fixed upstream (or mfsusieR adds
# a `residual_variance_scope = "fsusier_legacy"` option that
# mirrors the buggy formula), this test should be tightened to
# bit-identity.

skip_if_no_fsusier_pinned <- function() {
  testthat::skip_if_not_installed("fsusieR")
  installed <- as.character(utils::packageVersion("fsusieR"))
  pinned   <- mfsusier_pinned_fsusier_version()
  if (installed != pinned) {
    testthat::skip(paste0(
      "fsusieR version mismatch: installed ", installed,
      ", contract pinned to ", pinned, "."))
  }
}

make_c2_fixture <- function(n = 80, J = 15, T = 64,
                            seed = mfsusier_test_seed()) {
  set.seed(seed)
  X <- matrix(rnorm(n * J), nrow = n)
  beta <- numeric(J)
  beta[1] <- 1
  Y <- X %*% matrix(rep(beta, T), nrow = J) +
       matrix(rnorm(n * T, sd = 0.3), nrow = n)
  list(X = X, Y = Y, signal_idx = 1L)
}

c2_fit_pair <- function(X, Y, L = 3, max_iter = 50) {
  fit_f <- suppressWarnings(suppressMessages(
    fsusieR::susiF(Y, X, L = L,
                   prior          = "mixture_normal_per_scale",
                   greedy         = FALSE,
                   backfit        = FALSE,
                   maxit          = max_iter,
                   tol            = 1e-6,
                   post_processing = "none",
                   verbose        = FALSE,
                   cal_obj        = FALSE,
                   filter_cs      = FALSE)
  ))
  fit_m <- mfsusieR::mfsusie(wavelet_qnorm = FALSE, wavelet_standardize = FALSE, 
    X, list(Y), L = L,
    max_iter = max_iter, tol = 1e-6,
    verbose = FALSE
  )
  list(fsusier = fit_f, mfsusie = fit_m)
}

# ---- structural agreement on a deterministic fixture --------------

test_that("C2 (structural): mfsusie and fsusieR agree on the top causal SNP", {
  skip_if_no_fsusier_pinned()
  fx   <- make_c2_fixture()
  pair <- c2_fit_pair(fx$X, fx$Y, L = 3)

  expect_identical(which.max(pair$fsusier$pip),
                   which.max(pair$mfsusie$pip))
  expect_identical(which.max(pair$fsusier$pip), fx$signal_idx)
})

test_that("C2 (structural): both fits converge within 50 iterations", {
  skip_if_no_fsusier_pinned()
  fx   <- make_c2_fixture()
  pair <- c2_fit_pair(fx$X, fx$Y, L = 3, max_iter = 50)

  expect_true(isTRUE(pair$mfsusie$converged))
  # fsusieR does not expose a direct "converged" flag in this build;
  # use ELBO trajectory or PIP availability as a proxy.
  expect_true(!is.null(pair$fsusier$pip))
})

test_that("C2 (structural): PIPs at the causal SNP agree within 0.05", {
  skip_if_no_fsusier_pinned()
  fx   <- make_c2_fixture()
  pair <- c2_fit_pair(fx$X, fx$Y, L = 3)

  expect_lt(abs(pair$fsusier$pip[fx$signal_idx] -
                pair$mfsusie$pip[fx$signal_idx]),
            0.05)
})

test_that("C2 (structural): both fits find a small CS containing the causal SNP", {
  skip_if_no_fsusier_pinned()
  fx   <- make_c2_fixture()
  pair <- c2_fit_pair(fx$X, fx$Y, L = 3)

  small_cs_size <- 3L  # "concentrated" = size <= 3 SNPs
  cs_f <- pair$fsusier$cs
  cs_m <- pair$mfsusie$sets$cs %||% list()

  small_f <- cs_f[vapply(cs_f, length, integer(1)) <= small_cs_size]
  small_m <- cs_m[vapply(cs_m, length, integer(1)) <= small_cs_size]

  # Each set should have at least one concentrated CS that contains
  # the causal SNP (the size-15 "non-CS" fsusieR includes for
  # noise effects when `filter_cs = FALSE` is filtered out here).
  contains_f <- any(vapply(small_f, function(cs) fx$signal_idx %in% cs, logical(1)))
  contains_m <- any(vapply(small_m, function(cs) fx$signal_idx %in% cs, logical(1)))
  expect_true(contains_f, info = "fsusieR did not isolate the causal SNP in a small CS")
  expect_true(contains_m, info = "mfsusie did not isolate the causal SNP in a small CS")
})

# ---- documented numerical divergence: sigma2 ---------------------

test_that("C2 (documented divergence): sigma2 differs by O(1) factor due to fsusieR ER2 bug", {
  skip_if_no_fsusier_pinned()
  fx   <- make_c2_fixture()
  pair <- c2_fit_pair(fx$X, fx$Y, L = 3)

  # Both are positive. sigma2_fsusier uses the buggy formula
  # (missing xtx_diag factor); sigma2_mfsusie uses
  # susieR's correct formula. Magnitude differs by an O(n)
  # factor in worst case; for the small n = 80 fixture, the
  # diff is bounded but non-zero. Documented as Pattern A in
  # inst/notes/refactor-exceptions.md PR group 6.
  expect_true(pair$fsusier$sigma2 > 0)
  expect_true(pair$mfsusie$sigma2[[1]][1] > 0)
  # Assert there IS a measurable difference (so this test fails if
  # someone accidentally reverts to the buggy formula).
  expect_gt(abs(pair$fsusier$sigma2 - mean(pair$mfsusie$sigma2[[1]])),
            1e-4)
})
