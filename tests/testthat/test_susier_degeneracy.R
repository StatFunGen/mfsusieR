# Apple-to-apple comparison against `susieR::susie` under the C1
# degeneracy contract from design.md D11a:
#
#   M = 1, T_1 = 1 (univariate trait, DWT short-circuits)
#   prior_variance_grid of length 1 (single mixture component)
#   null_prior_weight = 0 (no null component)
#   prior_variance_scope = "per_outcome"
#   residual_variance_scope = "per_outcome"
#   estimate_prior_method = "none" (V held fixed; equivalent to
#     susieR's `estimate_prior_variance = FALSE`)
#   L_greedy = NULL
#
# y is pre-scaled to unit variance on both sides so mfsusieR's
# internal sd-scaling becomes a no-op and both packages see
# numerically identical inputs. Design.md D11a sets the contract
# at `<= 1e-10`; we assert at `<= 1e-13` (3 orders of magnitude
# tighter than the contract floor) to catch any regression that
# introduces real drift. The remaining ~1e-13 gap is summation-
# order accumulation in `compute_kl` over n = 200 residuals,
# essentially the IEEE-754 double-precision floor for that
# operation. alpha, mu, sigma2, ELBO, niter all match
# bit-identically (~1e-16).
#
# Manuscript reference: methods/online_method.tex (mfsusieR
# reduces to susieR for the scalar SuSiE special case).

make_c1_fixture <- function(n = 200, J = 30, n_signals = 2,
                            seed = mfsusier_test_seed()) {
  set.seed(seed)
  X <- matrix(rnorm(n * J), nrow = n)
  beta <- numeric(J)
  signal_idx <- seq_len(n_signals)
  beta[signal_idx] <- c(1, -0.5, 0.3, 0.8)[seq_len(n_signals)]
  y_raw <- as.numeric(X %*% beta + rnorm(n, sd = 0.5))
  # Pre-scale y to unit variance so mfsusieR's univariate sd-scaling
  # is a no-op.
  y <- (y_raw - mean(y_raw)) / stats::sd(y_raw)
  list(X = X, y = y, signal_idx = signal_idx)
}

c1_fit_pair <- function(X, y, L = 5, V = 0.2) {
  fit_s <- susieR::susie(X, y, L = L,
                         scaled_prior_variance      = V,    # var(y) = 1, so V on both sides
                         estimate_prior_variance    = FALSE,
                         estimate_residual_variance = TRUE,
                         verbose                    = FALSE,
                         max_iter = 100, tol = 1e-8)
  # Match susieR's default convergence rule (ELBO) so iter counts and
  # final fits are bit-comparable. mfsusie's default is `"pip"` because
  # alpha is robust to mixsqp's GEM residual; for the susieR-degeneracy
  # contract we want apple-to-apple ELBO convergence.
  fit_m <- mfsusieR::mfsusie(wavelet_qnorm = FALSE, 
    X, list(matrix(y, ncol = 1)),
    L                        = L,
    prior_variance_grid      = V,
    prior_variance_scope     = "per_outcome",
    null_prior_weight        = 0,
    residual_variance_scope = "per_outcome",
    estimate_prior_variance = FALSE,
    convergence_method       = "elbo",
    L_greedy                 = NULL,
    max_iter = 100, tol = 1e-8,
    verbose = FALSE
  )
  list(susier = fit_s, mfsusie = fit_m)
}

# ---- Each L value: alpha, mu, mu2, lbf, KL, sigma2, elbo, pip ----

for (L_value in c(1, 5, 10)) {
  test_that(sprintf("C1: L = %d, alpha matches susieR at <= 1e-10", L_value), {
    fx <- make_c1_fixture()
    pair <- c1_fit_pair(fx$X, fx$y, L = L_value)
    expect_equal(pair$mfsusie$alpha, pair$susier$alpha, tolerance = 1e-13)
  })

  test_that(sprintf("C1: L = %d, sigma2 matches at <= 1e-10", L_value), {
    fx <- make_c1_fixture()
    pair <- c1_fit_pair(fx$X, fx$y, L = L_value)
    expect_equal(pair$mfsusie$sigma2[[1]], pair$susier$sigma2,
                 tolerance = 1e-13)
  })

  test_that(sprintf("C1: L = %d, pip matches at <= 1e-10", L_value), {
    fx <- make_c1_fixture()
    pair <- c1_fit_pair(fx$X, fx$y, L = L_value)
    expect_equal(unname(pair$mfsusie$pip), unname(pair$susier$pip),
                 tolerance = 1e-13)
  })

  test_that(sprintf("C1: L = %d, ELBO at convergence matches at <= 1e-10", L_value), {
    fx <- make_c1_fixture()
    pair <- c1_fit_pair(fx$X, fx$y, L = L_value)
    expect_equal(pair$mfsusie$elbo[length(pair$mfsusie$elbo)],
                 pair$susier$elbo[length(pair$susier$elbo)],
                 tolerance = 1e-13)
  })

  test_that(sprintf("C1: L = %d, niter matches", L_value), {
    fx <- make_c1_fixture()
    pair <- c1_fit_pair(fx$X, fx$y, L = L_value)
    expect_identical(pair$mfsusie$niter, pair$susier$niter)
  })

  test_that(sprintf("C1: L = %d, lbf per effect matches at <= 1e-10", L_value), {
    fx <- make_c1_fixture()
    pair <- c1_fit_pair(fx$X, fx$y, L = L_value)
    expect_equal(pair$mfsusie$lbf, pair$susier$lbf, tolerance = 1e-13)
  })

  test_that(sprintf("C1: L = %d, KL per effect matches at <= 1e-10", L_value), {
    fx <- make_c1_fixture()
    pair <- c1_fit_pair(fx$X, fx$y, L = L_value)
    expect_equal(pair$mfsusie$KL, pair$susier$KL, tolerance = 1e-13)
  })
}

# ---- Credible-set membership: exact ----

test_that("C1: credible-set membership exactly matches susieR (L = 5)", {
  fx <- make_c1_fixture()
  pair <- c1_fit_pair(fx$X, fx$y, L = 5)

  # Both fits should produce the same (unordered) set of CSes.
  cs_m <- if (!is.null(pair$mfsusie$sets$cs))
    lapply(pair$mfsusie$sets$cs, sort) else list()
  cs_s <- if (!is.null(pair$susier$sets$cs))
    lapply(pair$susier$sets$cs, sort) else list()

  expect_identical(length(cs_m), length(cs_s))
  for (i in seq_along(cs_m)) {
    expect_identical(cs_m[[i]], cs_s[[i]])
  }
})
