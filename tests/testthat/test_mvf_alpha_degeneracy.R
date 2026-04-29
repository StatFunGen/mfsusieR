# C3 contract: end-to-end multi-modality functional case vs
# `mvf.susie.alpha::multfsusie`. Per design.md D11c the contract is
# `<= 1e-8` apple-to-apple, but bit-identity end-to-end is blocked
# by the ER2 / residual-variance formula divergence (Pattern A in
# inst/notes/refactor-exceptions.md PR group 6).
#
# This file gives the closest-to-end-to-end coverage we can offer
# while the ER2 bug stands upstream:
#
#   per-component bit-identity (already covered):
#     - test_mvf_alpha_fidelity.R: per-(SNP, modality) log-BF rows,
#       cross-modality combine, posterior moments, mixsqp L matrix
#       (all <= 1e-12)
#     - test_posterior_mixture.R: cpp11 vs R oracle (1e-15)
#
#   end-to-end structural (this file):
#     - same top SNP at convergence
#     - causal SNP isolated in a "small" CS on both sides
#     - PIP at the causal SNP within 0.05
#     - documented sigma2 divergence (Pattern A regression guard)

skip_if_no_mvf_alpha_pinned <- function() {
  testthat::skip_if_not_installed("mvf.susie.alpha")
  installed <- as.character(utils::packageVersion("mvf.susie.alpha"))
  pinned   <- mfsusier_pinned_mvf_version()
  if (installed != pinned) {
    testthat::skip(paste0(
      "mvf.susie.alpha version mismatch: installed ", installed,
      ", contract pinned to ", pinned, "."))
  }
}

make_c3_fixture <- function(n = 80, J = 12, T_per = c(32L, 64L),
                            seed = mfsusier_test_seed()) {
  set.seed(seed)
  X <- matrix(rnorm(n * J), nrow = n)
  beta <- numeric(J); beta[1] <- 1
  Y <- lapply(T_per, function(T_m) {
    X %*% matrix(rep(beta, T_m), nrow = J) +
      matrix(rnorm(n * T_m, sd = 0.3), nrow = n)
  })
  list(X = X, Y = Y, signal_idx = 1L)
}

c3_fit_pair <- function(X, Y, L = 3, max_iter = 30) {
  # mfsusie with `prior_variance_scope = "per_outcome"` to match
  # the (single-grid) mixture_normal prior that multfsusie supports.
  fit_m <- mfsusieR::mfsusie(wavelet_qnorm = FALSE, X, Y, L = L,
                             prior_variance_scope = "per_outcome",
                             max_iter = max_iter, tol = 1e-6,
                             verbose = FALSE)
  # multfsusie expects Y as a list with $Y_f (functional matrices)
  # and $Y_u (univariate; NULL here).
  Y_mvf <- list(Y_f = Y, Y_u = NULL)
  fit_mvf <- suppressMessages(suppressWarnings(
    mvf.susie.alpha::multfsusie(Y = Y_mvf, X = X, L = L,
                                prior          = "mixture_normal",
                                greedy         = FALSE,
                                backfit        = FALSE,
                                maxit          = max_iter,
                                tol            = 1e-6,
                                post_processing = "none",
                                verbose        = FALSE,
                                cal_obj        = FALSE,
                                filter_cs      = FALSE)
  ))
  list(mfsusie = fit_m, mvf = fit_mvf)
}

# ---- structural agreement on a deterministic fixture --------------

test_that("C3 (structural): mfsusie and multfsusie agree on the top causal SNP", {
  skip_if_no_mvf_alpha_pinned()
  fx   <- make_c3_fixture()
  pair <- c3_fit_pair(fx$X, fx$Y, L = 3)

  expect_identical(which.max(pair$mfsusie$pip),
                   which.max(pair$mvf$pip))
  expect_identical(which.max(pair$mfsusie$pip), fx$signal_idx)
})

test_that("C3 (structural): both fits isolate the causal SNP in a small CS", {
  skip_if_no_mvf_alpha_pinned()
  fx   <- make_c3_fixture()
  pair <- c3_fit_pair(fx$X, fx$Y, L = 3)

  small_cs_size <- 3L
  cs_m <- pair$mfsusie$sets$cs %||% list()
  cs_v <- pair$mvf$cs %||% list()

  small_m <- cs_m[vapply(cs_m, length, integer(1)) <= small_cs_size]
  small_v <- cs_v[vapply(cs_v, length, integer(1)) <= small_cs_size]

  contains_m <- any(vapply(small_m, function(cs) fx$signal_idx %in% cs, logical(1)))
  contains_v <- any(vapply(small_v, function(cs) fx$signal_idx %in% cs, logical(1)))
  expect_true(contains_m,
              info = "mfsusie did not isolate the causal SNP in a small CS")
  expect_true(contains_v,
              info = "multfsusie did not isolate the causal SNP in a small CS")
})

test_that("C3 (structural): PIP at causal SNP within 0.05 across mfsusie / multfsusie", {
  skip_if_no_mvf_alpha_pinned()
  fx   <- make_c3_fixture()
  pair <- c3_fit_pair(fx$X, fx$Y, L = 3)

  expect_lt(abs(pair$mfsusie$pip[fx$signal_idx] -
                pair$mvf$pip[fx$signal_idx]),
            0.05)
})

# ---- documented numerical divergence: sigma2 ---------------------

test_that("C3 (documented divergence): sigma2 differs because of mvf ER2 bug", {
  skip_if_no_mvf_alpha_pinned()
  fx   <- make_c3_fixture()
  pair <- c3_fit_pair(fx$X, fx$Y, L = 3)

  # multfsusie's sigma2 is a list with $sd_f (per modality scalar).
  s2_v <- pair$mvf$sigma2$sd_f
  s2_m <- vapply(pair$mfsusie$sigma2, function(s) mean(s), numeric(1))
  expect_true(all(s2_v > 0))
  expect_true(all(s2_m > 0))

  # Assert there IS a measurable difference (regression guard
  # against accidental revert to the buggy formula). Pattern A in
  # inst/notes/refactor-exceptions.md PR group 6.
  expect_gt(max(abs(s2_v - s2_m)), 1e-4)
})
