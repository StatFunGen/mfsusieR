# Apple-to-apple comparison against mvf.susie.alpha at its pinned
# contract commit (see setup.R::mfsusier_pinned_mvf_sha). These
# tests verify per-component fidelity for the routines mvf reuses
# from fsusieR (log_BF, post_mat_mean, post_mat_sd, cal_zeta) plus
# the mvf-specific assemblies of those (the per-(SNP, modality)
# log-BF matrix that EM_pi_multsusie consumes).
#
# Per refactor-discipline section 3 / design.md D11d, tolerance is
# `<= 1e-12` for closed-form parity with the port-source path.
#
# Skip-on-bug policy: tests that would compare against a known
# port-source bug (e.g. mvf's get_ER2 formula -- see
# inst/notes/refactor-exceptions.md) are explicitly NOT here. See
# `test_variance_and_elbo.R` for the corrected ER2 formula tests.
#
# Manuscript references for the math being checked:
#   methods/derivation.tex eq:post_f_mix
#   methods/derivation.tex eq:post_f2_mix
#   methods/online_method.tex line 41 (cross-modality combine)

# ---- Per-modality log-BF matrix (the row-vector EM_pi_multsusie sums) ----

test_that("per-modality log-BF rows match what mvf.susie.alpha::log_BF builds via fsusieR::log_BF", {
  skip_if_pinned_mvf_mismatch()
  skip_if_no_fsusier()

  set.seed(mfsusier_test_seed())
  n <- 50; J <- 15
  data <- mfsusieR:::create_mf_individual(wavelet_qnorm = FALSE, wavelet_standardize = FALSE, 
    X = matrix(rnorm(n * J), n),
    Y = list(matrix(rnorm(n * 64), n),
             matrix(rnorm(n * 32), n)),
    verbose = FALSE)

  # Hand-built per-(modality) ash-shaped prior, identical at each scale
  # within a modality.
  sd_grid <- c(0, 0.5, 1.5)
  pi_grid <- c(0.6, 0.2, 0.2)
  build_per_outcome_g <- function(scale_idx_lst) {
    g_norm <- ashr::normalmix(pi = pi_grid,
                              mean = rep(0, length(pi_grid)),
                              sd   = sd_grid)
    rec <- list(fitted_g = g_norm); class(rec) <- "ash"
    out <- rep(list(rec), length(scale_idx_lst))
    attr(out, "class") <- "mixture_normal_per_scale"
    out
  }
  G_prior_outcomes <- lapply(data$scale_index, build_per_outcome_g)

  # Initial-state Bhat / Shat per modality (zero-init, sigma2 = 1).
  pw <- data$xtx_diag
  Bhat_per_m <- lapply(data$D, function(D_m) crossprod(data$X, D_m) / pw)
  Shat_per_m <- lapply(data$D, function(D_m) sqrt(matrix(1 / pw, J, ncol(D_m))))

  # Reference per-modality log-BF row k.
  ref_per_m <- vapply(seq_len(data$M), function(m) {
    fsusieR::log_BF(G_prior_outcomes[[m]],
                    Bhat = Bhat_per_m[[m]], Shat = Shat_per_m[[m]],
                    lowc_wc = NULL, indx_lst = data$scale_index[[m]])
  }, numeric(J))

  # mfsusieR per-modality log-BF (sum over scales of mixture_log_bf_per_scale).
  ours_per_m <- vapply(seq_len(data$M), function(m) {
    out <- numeric(J)
    for (s in seq_along(data$scale_index[[m]])) {
      idx <- data$scale_index[[m]][[s]]
      out <- out + mfsusieR:::mixture_log_bf_per_scale(
        Bhat_per_m[[m]][, idx, drop = FALSE],
        Shat_per_m[[m]][, idx, drop = FALSE],
        sd_grid, pi_grid)
    }
    out
  }, numeric(J))

  expect_equal(ours_per_m, ref_per_m, tolerance = 1e-12)
})

# ---- Cross-modality combine matches `apply(lBF_per_trait, 2, sum)` ----

test_that("combine_outcome_lbfs (independent default) matches mvf.susie.alpha's apply-sum across modalities", {
  skip_if_pinned_mvf_mismatch()

  set.seed(mfsusier_test_seed())
  J <- 20
  per_outcome_lbfs <- list(rnorm(J), rnorm(J), rnorm(J))
  ref <- apply(do.call(rbind, per_outcome_lbfs), 2, sum)
  ours <- mfsusieR:::combine_outcome_lbfs(
    mfsusieR:::cross_outcome_prior_independent(),
    per_outcome_lbfs,
    model_state = NULL)
  expect_equal(ours, ref, tolerance = 1e-15)   # FP summation order differs by ~ 1 ULP
})

# ---- mixture_posterior matches per-modality fsusieR::post_mat_* ----

test_that("mixture posterior matches fsusieR::post_mat_mean / post_mat_sd at <= 1e-12 for each (modality, scale)", {
  skip_if_pinned_mvf_mismatch()
  skip_if_no_fsusier()

  set.seed(mfsusier_test_seed())
  n <- 40; J <- 12
  data <- mfsusieR:::create_mf_individual(wavelet_qnorm = FALSE, wavelet_standardize = FALSE, 
    X = matrix(rnorm(n * J), n),
    Y = list(matrix(rnorm(n * 32), n)),
    verbose = FALSE)

  sd_grid <- c(0, 0.7, 2.0)
  pi_grid <- c(0.5, 0.3, 0.2)
  g_norm <- ashr::normalmix(pi = pi_grid,
                            mean = rep(0, length(pi_grid)),
                            sd   = sd_grid)
  rec <- list(fitted_g = g_norm); class(rec) <- "ash"
  G_prior <- rep(list(rec), length(data$scale_index[[1]]))
  attr(G_prior, "class") <- "mixture_normal_per_scale"

  pw <- data$xtx_diag
  Bhat <- crossprod(data$X, data$D[[1]]) / pw
  Shat <- sqrt(matrix(1 / pw, J, ncol(data$D[[1]])))

  for (s in seq_along(data$scale_index[[1]])) {
    idx <- data$scale_index[[1]][[s]]

    # fsusieR ref via element-wise ashr loop.
    ref_pm <- vapply(seq_len(J), function(j) {
      d <- ashr::set_data(Bhat[j, idx], Shat[j, idx])
      ashr::postmean(g_norm, d)
    }, numeric(length(idx)))
    ref_pm <- if (is.matrix(ref_pm)) t(ref_pm) else matrix(ref_pm, J, 1)

    ours <- mfsusieR:::mixture_posterior_per_scale(
      Bhat[, idx, drop = FALSE], Shat[, idx, drop = FALSE], sd_grid, pi_grid)
    expect_equal(ours$pmean, ref_pm, tolerance = 1e-12)
  }
})

# ---- mf_em_likelihood_per_scale matches fsusieR::cal_L_mixsq_s_per_scale ----

test_that("mf_em_likelihood_per_scale matches fsusieR::cal_L_mixsq_s_per_scale at <= 1e-12", {
  skip_if_pinned_mvf_mismatch()
  skip_if_no_fsusier()

  set.seed(mfsusier_test_seed())
  n <- 40; J <- 12
  data <- mfsusieR:::create_mf_individual(wavelet_qnorm = FALSE, wavelet_standardize = FALSE, 
    X = matrix(rnorm(n * J), n),
    Y = list(matrix(rnorm(n * 32), n)),
    verbose = FALSE)

  sd_grid <- c(0, 0.7, 2.0)
  g_norm  <- ashr::normalmix(pi = c(1, 0, 0),
                             mean = rep(0, length(sd_grid)),
                             sd   = sd_grid)
  rec <- list(fitted_g = g_norm); class(rec) <- "ash"
  G_prior <- rep(list(rec), length(data$scale_index[[1]]))
  attr(G_prior, "class") <- "mixture_normal_per_scale"

  pw   <- data$xtx_diag
  Bhat <- crossprod(data$X, data$D[[1]]) / pw
  Shat <- sqrt(matrix(1 / pw, J, ncol(data$D[[1]])))

  for (s in seq_along(data$scale_index[[1]])) {
    idx <- data$scale_index[[1]][[s]]
    ref_L <- fsusieR:::cal_L_mixsq_s_per_scale(G_prior, s, Bhat, Shat,
                                               indx_lst = data$scale_index[[1]])
    ours_L <- mfsusieR:::mf_em_likelihood_per_scale(
      Bhat[, idx, drop = FALSE], Shat[, idx, drop = FALSE], sd_grid)
    expect_equal(ours_L, ref_L, tolerance = 1e-12)
  }
})
