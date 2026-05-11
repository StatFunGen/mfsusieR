# Full-pipeline comparison: mfsusieR HEAD vs mvf.susie.alpha
#
# Covers: fit -> HMM post-smooth -> CS -> LFSR -> posthoc
# Not a parity test. Documents call signatures, output structure, and
# known divergences (D1-D3).
#
# D1. sigma2: mfsusieR > mvf because mvf ER2 bias correction is missing
#     predictor_weights (deflated ~2x).
# D2. Bhat/Shat: mvf uses per-variable marginal SE; mfsusieR uses global
#     sigma2/xtx. Different algorithm, not a bug.
# D3. ELBO KL sign error in mvf; does not affect PIP convergence.
#
# Companion: inst/notes/sessions/2026-05-11-mfsusier-vs-mvf-engineering-comparison.md

skip_if_no_mvf <- function() skip_if_not_installed("mvf.susie.alpha")

# ── Data and parameters ───────────────────────────────────────────────────────

make_mwe_data <- function(seed = 42L, n = 60L, p = 80L, M = 2L, T_y = 32L) {
  set.seed(seed)
  X    <- matrix(rnorm(n * p), n, p)
  beta <- matrix(0, p, T_y)
  for (j in c(10L, 60L)) {
    t    <- seq_len(T_y)
    beta[j, ] <- exp(-0.5 * ((t - T_y / 2) / (T_y / 8))^2) * rnorm(1, sd = 0.8)
  }
  Y <- lapply(seq_len(M), function(m)
    X %*% beta + matrix(rnorm(n * T_y, sd = 0.5), n, T_y))
  list(X = X, Y = Y, n = n, p = p, M = M, T_y = T_y, sig_idx = c(10L, 60L))
}

MWE_L    <- 5L
MWE_CTRL <- list(convtol.sqp = 1e-8, numiter.em = 20L, tol.svd = 1e-10, verbose = FALSE)

# ── Lazy fit cache ────────────────────────────────────────────────────────────
# Fits both packages once; shared across all tests below.

.cache <- new.env(parent = emptyenv())

.get_fits <- function() {
  if (isTRUE(.cache$ready)) return(.cache)
  skip_if_no_mvf()
  d <- make_mwe_data()
  pos <- lapply(seq_len(d$M), function(m) seq_len(d$T_y))

  # mfsusieR ─────────────────────────────────────────────────────────────────
  #   1. mfsusie(): returns fit with class "mfsusie"
  #   2. mf_post_smooth(method="HMM"): adds fit$smoothed$HMM with
  #        $effect_curves[[m]][[l]]  — length-T position-space curve
  #        $credible_bands[[m]][[l]] — T x 2 matrix [lower, upper]
  #        $lfsr_curves[[m]][[l]]    — length-T LFSR in [0,1]
  #      indexed by outcome m (1..M) and effect l (1..L, all L)
  #   3. susie_post_outcome_configuration(): reads fit$lbf_variable_outcome
  #      (L x p x M), enumerates 2^M configs per CS tuple

  fit_mf <- mfsusie(
    X = d$X, Y = d$Y, pos = pos,
    L                    = MWE_L,
    prior_variance_scope = "per_outcome",
    wavelet_qnorm        = FALSE,
    max_inner_em_steps   = 0L,
    mixture_null_weight  = 0.12,
    control_mixsqp       = MWE_CTRL,
    L_greedy             = NULL,
    tol = 1e-4, max_iter = 30L, verbose = FALSE
  )
  fit_mf <- mf_post_smooth(fit_mf, method = "HMM")
  pc_mf  <- susieR::susie_post_outcome_configuration(
    fit_mf, by = "outcome", method = "susiex"
  )

  # mvf ──────────────────────────────────────────────────────────────────────
  #   multfsusie() with post_processing="HMM" and posthoc=TRUE runs
  #   both inside out_prep(), trimming null effects before output.
  #   After trimming: length(fit$alpha) = length(fit$cs) = L_active (<= L).
  #   HMM output indexed by CS (l_cs = 1..L_active), modality (k = 1..M):
  #        fit$fitted_func[[l_cs]][[k]] — length-T effect curve
  #        fit$cred_band[[l_cs]][[k]]   — 2 x T matrix [lower; upper]  (rows, not cols)
  #        fit$lfsr[[l_cs]]$est_lfsr_functional[[k]] — length-T LFSR
  #   Posthoc at fit$posthoc[[l_cs]], called automatically.

  fit_mvf <- mvf.susie.alpha::multfsusie(
    Y               = list(Y_f = d$Y),
    X               = d$X,
    pos             = pos,
    L               = MWE_L,
    prior           = "mixture_normal",
    nullweight      = 0.7,
    maxit           = 30L,
    tol             = 1e-4,
    greedy          = FALSE,
    backfit         = FALSE,
    post_processing = "HMM",
    posthoc         = TRUE,
    verbose         = FALSE
  )

  .cache$d     <- d
  .cache$mf    <- fit_mf
  .cache$pc_mf <- pc_mf
  .cache$mvf   <- fit_mvf
  .cache$ready <- TRUE
  .cache
}

# ── Tests ─────────────────────────────────────────────────────────────────────

test_that("both pipelines complete without error", {
  c <- .get_fits()
  expect_true(inherits(c$mf, "mfsusie"))
  expect_true(is.list(c$mvf))
  expect_true(!is.null(c$mf$smoothed$HMM))
  expect_true(inherits(c$pc_mf, "susie_post_outcome_configuration"))
})

test_that("alpha: L x p matrix (mfsusieR) vs list of L_active vectors (mvf)", {
  c <- .get_fits(); d <- c$d
  # mfsusieR: always L x p, row sums = 1
  expect_equal(dim(c$mf$alpha), c(MWE_L, d$p))
  expect_equal(rowSums(c$mf$alpha), rep(1, MWE_L), tolerance = 1e-10)
  # mvf: out_prep trims null effects; length = L_active (<= L)
  L_active <- length(c$mvf$alpha)
  expect_true(L_active <= MWE_L)
  expect_equal(ncol(do.call(rbind, c$mvf$alpha)), d$p)
  expect_equal(vapply(c$mvf$alpha, sum, numeric(1)),
               rep(1, L_active), tolerance = 1e-10)
  message(sprintf("L_active: mvf=%d / L=%d", L_active, MWE_L))
})

test_that("PIP: both packages return a length-p vector", {
  c <- .get_fits(); d <- c$d
  expect_equal(length(c$mf$pip),  d$p)
  expect_equal(length(c$mvf$pip), d$p)
  sig <- d$sig_idx
  message(sprintf("PIP at signals: mfsusieR [%s], mvf [%s]",
    paste(round(c$mf$pip[sig],  3), collapse = "/"),
    paste(round(c$mvf$pip[sig], 3), collapse = "/")))
})

test_that("sigma2: mfsusieR > mvf (D1: correct vs deflated ER2)", {
  c <- .get_fits()
  # mfsusieR: list of M scalars
  s_mf  <- unlist(c$mf$sigma2)
  # mvf: $sigma2$sd_f stores SDs, not variances; square to compare
  s_mvf <- unlist(c$mvf$sigma2$sd_f)^2
  expect_true(all(s_mf > s_mvf))
  message(sprintf("sigma2 mfsusieR/mvf ratio: [%s]",
    paste(round(s_mf / s_mvf, 2), collapse = "/")))
})

test_that("CS: count and purity accessible from both fits", {
  c <- .get_fits()
  # mfsusieR: fit$sets$cs — named list of integer vectors (SNP indices)
  #           fit$sets$purity — data.frame, one row per CS
  cs_mf <- c$mf$sets$cs
  expect_true(is.list(cs_mf))
  if (length(cs_mf) > 0L) {
    expect_true(is.data.frame(c$mf$sets$purity) || is.matrix(c$mf$sets$purity))
    pur_mf <- c$mf$sets$purity[, "mean.abs.corr"]
    message(sprintf("mfsusieR: %d CS, mean purity [%s]",
      length(cs_mf), paste(round(pur_mf, 3), collapse = ", ")))
  }
  # mvf: fit$cs — list of L_active integer vectors (trimmed)
  #      fit$purity — matrix, one row per CS
  cs_mvf <- c$mvf$cs
  expect_true(is.list(cs_mvf))
  if (length(cs_mvf) > 0L && !is.null(c$mvf$purity)) {
    pur_mvf <- if (is.matrix(c$mvf$purity)) rowMeans(c$mvf$purity, na.rm = TRUE)
               else as.numeric(c$mvf$purity)
    message(sprintf("mvf: %d CS, mean purity [%s]",
      length(cs_mvf), paste(round(pur_mvf, 3), collapse = ", ")))
  }
  message(sprintf("CS count: mfsusieR=%d, mvf=%d",
    length(cs_mf), length(cs_mvf)))
})

test_that("HMM effect curves: shape and indexing", {
  c <- .get_fits(); d <- c$d
  hmm <- c$mf$smoothed$HMM

  # mfsusieR: [[outcome m]][[effect l]], all M x L slots, each length T_y
  expect_equal(length(hmm$effect_curves),        d$M)
  expect_equal(length(hmm$effect_curves[[1L]]), MWE_L)
  expect_equal(length(hmm$effect_curves[[1L]][[1L]]), d$T_y)

  # mfsusieR credible bands: T x 2 (lower, upper as columns)
  expect_equal(dim(hmm$credible_bands[[1L]][[1L]]), c(d$T_y, 2L))

  # mvf: [[CS index]][[modality k]], only L_active slots
  n_cs <- length(c$mvf$fitted_func)
  if (n_cs > 0L) {
    expect_equal(length(c$mvf$fitted_func[[1L]]), d$M)
    expect_equal(length(c$mvf$fitted_func[[1L]][[1L]]), d$T_y)
    # mvf credible bands: 2 x T (lower/upper as rows, opposite of mfsusieR)
    expect_equal(dim(c$mvf$cred_band[[1L]][[1L]]), c(2L, d$T_y))
  }
})

test_that("LFSR: values in [0,1]; indexing differs between packages", {
  c <- .get_fits(); d <- c$d

  # mfsusieR: $smoothed$HMM$lfsr_curves[[m]][[l]], all M x L effects
  lfsr_mf <- c$mf$smoothed$HMM$lfsr_curves
  expect_equal(length(lfsr_mf),        d$M)
  expect_equal(length(lfsr_mf[[1L]]), MWE_L)
  for (m in seq_len(d$M))
    for (l in seq_len(MWE_L))
      expect_true(all(lfsr_mf[[m]][[l]] >= 0 & lfsr_mf[[m]][[l]] <= 1))

  # mvf: $lfsr[[l_cs]]$est_lfsr_functional[[k]], only detected CS
  n_cs <- length(c$mvf$lfsr)
  if (n_cs > 0L) {
    for (l in seq_len(n_cs))
      for (k in seq_len(d$M)) {
        v <- c$mvf$lfsr[[l]]$est_lfsr_functional[[k]]
        expect_true(all(v >= 0 & v <= 1))
      }
  }
  message(sprintf(
    "LFSR slots: mfsusieR %dx%d=%d; mvf %d CS x %d outcomes=%d",
    d$M, MWE_L, d$M * MWE_L, n_cs, d$M, n_cs * d$M))
})

test_that("posthoc: explicit call (mfsusieR) vs automatic (mvf)", {
  c <- .get_fits(); d <- c$d

  # mfsusieR: lbf_variable_outcome (L x p x M) is stored during IBSS;
  #   user calls susie_post_outcome_configuration() after fit.
  expect_equal(dim(c$mf$lbf_variable_outcome), c(MWE_L, d$p, d$M))
  pc <- c$pc_mf
  expect_true(is.list(pc$susiex))
  if (length(pc$susiex) > 0L) {
    t1 <- pc$susiex[[1L]]
    expect_true(all(c("cs_indices", "logBF_trait", "configs",
                      "config_prob", "marginal_prob") %in% names(t1)))
    expect_equal(length(t1$logBF_trait), d$M)
    expect_equal(nrow(t1$configs), 2L^d$M)
  }

  # mvf: posthoc computed automatically inside out_prep(); stored at fit$posthoc.
  #   Length = L_active (same as fit$alpha, fit$cs).
  expect_true(is.list(c$mvf$posthoc))
  expect_true(length(c$mvf$posthoc) <= MWE_L)
  non_null <- Filter(Negate(is.null), c$mvf$posthoc)
  if (length(non_null) > 0L) {
    t1 <- non_null[[1L]]
    expect_true(all(c("logBF_trait", "posthoc", "active",
                      "configs", "config_prob") %in% names(t1)))
    expect_equal(nrow(t1$configs), 2L^length(t1$logBF_trait))
  }
  message(sprintf("posthoc susiex tuples: mfsusieR=%d, mvf=%d",
    length(pc$susiex), length(non_null)))
})
