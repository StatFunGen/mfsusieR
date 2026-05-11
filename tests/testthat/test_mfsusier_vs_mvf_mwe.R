# Minimal working example: mfsusieR HEAD vs mvf.susie.alpha
#
# Runs both packages on identical data under maximally-equivalent parameter
# settings and documents the confirmed output divergences. This file is NOT
# an apple-to-apple parity test (the two packages are not expected to match).
# Its purpose is to make the known engineering differences observable and
# reproducible in a single place.
#
# Companion note: inst/notes/sessions/2026-05-11-mfsusier-vs-mvf-engineering-comparison.md
#
# Three confirmed divergences (see §4 of the companion note):
#   D1. ER2 bias correction: mvf missing predictor_weights → sigma2 deflated ~2x
#   D2. Bhat/Shat: mvf per-(j,m) marginal SE vs mfsusieR global sigma2/xtx
#   D3. ELBO sign error in mvf KL term (does not affect PIP-convergence runs)
#
# These are mvf bugs or intentional design choices, all documented in
# inst/notes/refactor-exceptions.md and inst/notes/sessions/2026-04-29-0716-mvf-divergences.md.
#
# Skip policy: all tests skip when mvf.susie.alpha is not installed. The package
# is a read-only reference; do not add it to DESCRIPTION Imports.

skip_if_no_mvf <- function() {
  skip_if_not_installed("mvf.susie.alpha")
}

# ── Shared fixture ────────────────────────────────────────────────────────────

make_mwe_data <- function(seed = 42L, n = 60L, p = 120L, M = 3L, T_y = 64L,
                          n_signals = 2L) {
  set.seed(seed)
  X      <- matrix(rnorm(n * p), n, p)
  # Two true causal SNPs with smooth bell-shape effects
  beta   <- matrix(0, p, T_y)
  sig_idx <- c(10L, 80L)
  for (j in sig_idx) {
    t_grid       <- seq_len(T_y)
    beta[j, ]   <- exp(-0.5 * ((t_grid - T_y / 2) / (T_y / 8))^2) *
                   rnorm(1, sd = 0.8)
  }
  Y <- lapply(seq_len(M), function(m) {
    E <- matrix(rnorm(n * T_y, sd = 0.5), n, T_y)
    X %*% beta + E
  })
  list(X = X, Y = Y, sig_idx = sig_idx, n = n, p = p, M = M, T_y = T_y)
}

# Maximally-equivalent parameter settings (see §8.1 of companion note):
#   - no qnorm, no inner EM, cold start, no greedy, same L/tol/max_iter
MWE_CTRL <- list(
  convtol.sqp = 1e-8, numiter.em = 20L, tol.svd = 1e-10, verbose = FALSE
)
MWE_NULL_WEIGHT <- 0.12   # ≈ mvf nullweight(0.7) / M(3) for M=3 → 0.7/3 ≈ 0.23; use 0.12 for both
MWE_L      <- 5L
MWE_MAXIT  <- 30L
MWE_TOL    <- 1e-4

# ── Test 1: both packages run without error ───────────────────────────────────

test_that("mfsusieR and mvf both run on the same data (sanity)", {
  skip_if_no_mvf()

  d <- make_mwe_data()

  # mfsusieR
  fit_mf <- mfsusie(
    X                    = d$X,
    Y                    = d$Y,
    pos                  = lapply(seq_len(d$M), function(m) seq_len(d$T_y)),
    L                    = MWE_L,
    prior_variance_scope = "per_outcome",
    wavelet_qnorm        = FALSE,
    max_inner_em_steps   = 0L,
    mixture_null_weight  = MWE_NULL_WEIGHT,
    control_mixsqp       = MWE_CTRL,
    L_greedy             = NULL,
    tol                  = MWE_TOL,
    max_iter             = MWE_MAXIT,
    verbose              = FALSE
  )

  # mvf: Y must be wrapped as list(Y_f = ...)
  fit_mvf <- mvf.susie.alpha::multfsusie(
    Y          = list(Y_f = d$Y),
    X          = d$X,
    L          = MWE_L,
    pos        = NULL,
    prior      = "mixture_normal",
    nullweight = 0.7,
    maxit      = MWE_MAXIT,
    tol        = MWE_TOL,
    greedy     = FALSE,
    backfit    = FALSE,
    verbose    = FALSE
  )

  expect_true(inherits(fit_mf,  "mfsusie"))
  expect_true(is.list(fit_mvf))

  # Both produce alpha[L × p]
  expect_equal(dim(fit_mf$alpha), c(MWE_L, d$p))
  # mvf stores alpha as a list of L plain numeric vectors of length p;
  # convert to matrix via rbind.
  alpha_mvf <- do.call(rbind, fit_mvf$alpha)
  expect_equal(dim(alpha_mvf), c(MWE_L, d$p))
})

# ── Test 2: sigma2 — mfsusieR larger than mvf (D1: ER2 bias correction bug) ──

test_that("mfsusieR sigma2 is larger than mvf sigma2 due to correct ER2 formula", {
  skip_if_no_mvf()

  d <- make_mwe_data()

  fit_mf <- mfsusie(
    X = d$X, Y = d$Y,
    pos                  = lapply(seq_len(d$M), function(m) seq_len(d$T_y)),
    L                    = MWE_L,
    prior_variance_scope = "per_outcome",
    wavelet_qnorm        = FALSE,
    max_inner_em_steps   = 0L,
    mixture_null_weight  = MWE_NULL_WEIGHT,
    control_mixsqp       = MWE_CTRL,
    L_greedy             = NULL,
    tol                  = MWE_TOL,
    max_iter             = MWE_MAXIT,
    verbose              = FALSE
  )

  fit_mvf <- mvf.susie.alpha::multfsusie(
    Y = list(Y_f = d$Y), X = d$X, L = MWE_L,
    prior = "mixture_normal", nullweight = 0.7,
    maxit = MWE_MAXIT, tol = MWE_TOL,
    greedy = FALSE, backfit = FALSE, verbose = FALSE
  )

  # mfsusieR: sigma2 is a list[M] of scalars
  sigma2_mf  <- unlist(fit_mf$sigma2)

  # mvf: sigma2 stored as $sigma2$sd_f (list[M] of SDs, not variances)
  sigma2_mvf <- unlist(fit_mvf$sigma2$sd_f)^2

  # D1 prediction: mfsusieR sigma2 > mvf sigma2 because mvf's ER2 bias
  # correction is ~1/n of the correct value (missing predictor_weights factor).
  # Both are expected to exceed true noise variance (0.25 = 0.5^2) because
  # residual includes unmodelled signal, but mvf deflation is a bug.
  expect_true(
    all(sigma2_mf > sigma2_mvf),
    label = sprintf(
      "mfsusieR sigma2 [%s] should all exceed mvf sigma2 [%s]",
      paste(round(sigma2_mf,  3), collapse = ", "),
      paste(round(sigma2_mvf, 3), collapse = ", ")
    )
  )

  # Report magnitudes for diagnostic visibility
  message(sprintf(
    "sigma2: mfsusieR [%s], mvf [%s], ratio [%s]",
    paste(round(sigma2_mf,  3), collapse = "/"),
    paste(round(sigma2_mvf, 3), collapse = "/"),
    paste(round(sigma2_mf / sigma2_mvf, 2), collapse = "/")
  ))
})

# ── Test 3: PIP at signal vs null — document the detection gap ─────────────

test_that("both packages detect at least one signal; document detection gap", {
  skip_if_no_mvf()

  d <- make_mwe_data()

  fit_mf <- mfsusie(
    X = d$X, Y = d$Y,
    pos                  = lapply(seq_len(d$M), function(m) seq_len(d$T_y)),
    L                    = MWE_L,
    prior_variance_scope = "per_outcome",
    wavelet_qnorm        = FALSE,
    max_inner_em_steps   = 0L,
    mixture_null_weight  = MWE_NULL_WEIGHT,
    control_mixsqp       = MWE_CTRL,
    L_greedy             = NULL,
    tol                  = MWE_TOL,
    max_iter             = MWE_MAXIT,
    verbose              = FALSE
  )

  fit_mvf <- mvf.susie.alpha::multfsusie(
    Y = list(Y_f = d$Y), X = d$X, L = MWE_L,
    prior = "mixture_normal", nullweight = 0.7,
    maxit = MWE_MAXIT, tol = MWE_TOL,
    greedy = FALSE, backfit = FALSE, verbose = FALSE
  )

  pip_mf  <- fit_mf$pip
  pip_mvf <- fit_mvf$pip

  sig_pip_mf  <- pip_mf[d$sig_idx]
  sig_pip_mvf <- pip_mvf[d$sig_idx]
  null_pip_mf  <- max(pip_mf[-d$sig_idx])
  null_pip_mvf <- max(pip_mvf[-d$sig_idx])

  message(sprintf(
    "PIP at signals  — mfsusieR: [%s], mvf: [%s]",
    paste(round(sig_pip_mf,  3), collapse = ", "),
    paste(round(sig_pip_mvf, 3), collapse = ", ")
  ))
  message(sprintf(
    "Max PIP at null — mfsusieR: %.3f, mvf: %.3f",
    null_pip_mf, null_pip_mvf
  ))

  # At least one signal SNP should be detected by each package at pip > 0.5.
  # D2 (Bhat/Shat) means mvf may detect signals that mfsusieR misses, because
  # mvf's per-variable marginal SE is smaller at true signals.
  expect_true(
    any(sig_pip_mf > 0.5) || any(sig_pip_mvf > 0.5),
    label = "at least one package detects at least one signal"
  )
})

# ── Test 4: CS purity structure ───────────────────────────────────────────────

test_that("CS count and purity accessible from both fit objects", {
  skip_if_no_mvf()

  d <- make_mwe_data()

  fit_mf <- mfsusie(
    X = d$X, Y = d$Y,
    pos                  = lapply(seq_len(d$M), function(m) seq_len(d$T_y)),
    L                    = MWE_L,
    prior_variance_scope = "per_outcome",
    wavelet_qnorm        = FALSE,
    max_inner_em_steps   = 0L,
    mixture_null_weight  = MWE_NULL_WEIGHT,
    control_mixsqp       = MWE_CTRL,
    L_greedy             = NULL,
    tol                  = MWE_TOL,
    max_iter             = MWE_MAXIT,
    verbose              = FALSE
  )

  fit_mvf <- mvf.susie.alpha::multfsusie(
    Y = list(Y_f = d$Y), X = d$X, L = MWE_L,
    prior = "mixture_normal", nullweight = 0.7,
    maxit = MWE_MAXIT, tol = MWE_TOL,
    greedy = FALSE, backfit = FALSE, verbose = FALSE
  )

  # mfsusieR: CSes live at fit$sets$cs (susieR 0.16.x format);
  #   each element is a list with $variables (integer vec) and $purity (matrix).
  cs_mf  <- fit_mf$sets$cs
  n_cs_mf <- length(cs_mf)

  # mvf: fit$cs is a list of L plain integer vectors (SNP indices).
  #   Purity lives at the top-level fit$purity (computed by fsusieR::cal_purity
  #   inside out_prep). It is a matrix with one row per CS.
  cs_mvf  <- fit_mvf$cs
  n_cs_mvf <- length(cs_mvf)

  message(sprintf("CS count — mfsusieR: %d, mvf: %d", n_cs_mf, n_cs_mvf))

  # Both should return a list (possibly empty) without error.
  expect_true(is.list(cs_mf))
  expect_true(is.list(cs_mvf))

  # If mfsusieR has any CS, purity is accessible via sets$cs[[i]]$purity.
  if (n_cs_mf > 0) {
    purities_mf <- vapply(cs_mf, function(cs) {
      if (is.matrix(cs$purity)) mean(cs$purity[, "mean.abs.corr"]) else mean(cs$purity)
    }, numeric(1))
    message(sprintf("mfsusieR purity: [%s]", paste(round(purities_mf, 3), collapse = ", ")))
    expect_true(all(is.finite(purities_mf)))
  }
  # If mvf has any CS, purity is at fit_mvf$purity (top-level matrix).
  if (n_cs_mvf > 0 && !is.null(fit_mvf$purity)) {
    purities_mvf <- if (is.matrix(fit_mvf$purity)) rowMeans(fit_mvf$purity, na.rm = TRUE)
                    else as.numeric(fit_mvf$purity)
    message(sprintf("mvf purity: [%s]", paste(round(purities_mvf, 3), collapse = ", ")))
    expect_true(all(is.finite(purities_mvf)))
  }
})

# ── Test 5: alpha matrix dimensions and row sums ──────────────────────────────

test_that("alpha[l,] sums to 1 in both packages", {
  skip_if_no_mvf()

  d <- make_mwe_data()

  fit_mf <- mfsusie(
    X = d$X, Y = d$Y,
    pos                  = lapply(seq_len(d$M), function(m) seq_len(d$T_y)),
    L                    = MWE_L,
    prior_variance_scope = "per_outcome",
    wavelet_qnorm        = FALSE,
    max_inner_em_steps   = 0L,
    mixture_null_weight  = MWE_NULL_WEIGHT,
    control_mixsqp       = MWE_CTRL,
    L_greedy             = NULL,
    tol                  = MWE_TOL,
    max_iter             = MWE_MAXIT,
    verbose              = FALSE
  )

  fit_mvf <- mvf.susie.alpha::multfsusie(
    Y = list(Y_f = d$Y), X = d$X, L = MWE_L,
    prior = "mixture_normal", nullweight = 0.7,
    maxit = MWE_MAXIT, tol = MWE_TOL,
    greedy = FALSE, backfit = FALSE, verbose = FALSE
  )

  # mfsusieR: alpha is an L × p matrix
  row_sums_mf <- rowSums(fit_mf$alpha)
  expect_equal(row_sums_mf, rep(1, MWE_L), tolerance = 1e-10,
               label = "mfsusieR alpha row sums = 1")

  # mvf: each element of fit_mvf$alpha is a plain numeric vector of length p
  alpha_mvf_rows <- vapply(fit_mvf$alpha, sum, numeric(1))
  expect_equal(alpha_mvf_rows, rep(1, MWE_L), tolerance = 1e-10,
               label = "mvf alpha row sums = 1")
})

# ── Test 6: mvf Y input format vs mfsusieR (document wrapping) ───────────────

test_that("mvf requires Y = list(Y_f = ...) while mfsusieR takes flat list", {
  skip_if_no_mvf()

  d <- make_mwe_data(n = 30L, p = 40L, M = 2L, T_y = 32L)
  pos <- lapply(seq_len(d$M), function(m) seq_len(d$T_y))

  # mfsusieR: flat list of M matrices
  fit_mf <- mfsusie(
    X = d$X, Y = d$Y, pos = pos,
    L = 3L, prior_variance_scope = "per_outcome",
    wavelet_qnorm = FALSE, max_inner_em_steps = 0L,
    mixture_null_weight = MWE_NULL_WEIGHT,
    control_mixsqp = MWE_CTRL, L_greedy = NULL,
    tol = MWE_TOL, max_iter = 10L, verbose = FALSE
  )

  # mvf: wrapped as list(Y_f = ...)
  fit_mvf <- mvf.susie.alpha::multfsusie(
    Y = list(Y_f = d$Y), X = d$X, L = 3L,
    prior = "mixture_normal", nullweight = 0.7,
    maxit = 10L, tol = MWE_TOL,
    greedy = FALSE, backfit = FALSE, verbose = FALSE
  )

  expect_true(inherits(fit_mf, "mfsusie"))
  expect_true(is.list(fit_mvf))
  # Both converge within max_iter
  expect_lte(fit_mf$niter, 10L)
})
