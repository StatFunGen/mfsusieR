# Per-effect S3 methods on `mf_individual` for the IBSS dispatch.
#
# Mirrors the susieR / mvsusieR pattern: compute_residuals computes
# the per-effect partial residual `R_l = Y - Xr_without_l` and
# caches `model$residuals = X^T R_l` for downstream SER stats.
# compute_ser_statistics turns the cached `XtR` into per-position
# `Bhat / Shat^2`. update_fitted_values folds the new posterior
# back into the running per-outcome fit `model$fitted[[m]]`.
#
# All three are per-outcome lapply over `seq_len(data$M)` plus
# matrix-level BLAS ops; no C++ needed at this layer.
#
# Manuscript references: methods/algorithms.tex eq:partial_resid;
# methods/derivation.tex eq:bhat_shat_per_outcome.

#' Per-effect partial residual on `mf_individual`
#'
#' For each outcome m, computes the contribution of effect `l`
#' to the running fit `Xb_lm = X %*% (alpha_l * mu_l[[m]])`,
#' removes it from `model$fitted[[m]]` to obtain the without-l
#' fit, and subtracts that from the wavelet matrix `data$D[[m]]`
#' to get the residual. Caches `XtR_m = X^T R_m` per outcome on
#' `model$residuals[[m]]`. Mirrors `compute_residuals.individual`
#' (susieR) and `compute_residuals.mv_individual` (mvsusieR).
#'
#' @keywords internal
#' @noRd
compute_residuals.mf_individual <- function(data, params, model, l, ...) {
  X <- data$X
  M <- data$M

  if (is.null(model$residuals))      model$residuals      <- vector("list", M)
  if (is.null(model$fitted_without_l)) model$fitted_without_l <- vector("list", M)
  if (is.null(model$raw_residuals))  model$raw_residuals  <- vector("list", M)

  alpha_l <- model$alpha[l, ]
  for (m in seq_len(M)) {
    # alpha_l . mu_l[[m]] is p x T_basis[m]; X is n x p. Xb is n x T_basis[m].
    b_lm  <- alpha_l * model$mu[[l]][[m]]
    Xb_lm <- X %*% b_lm

    Xr_without_lm <- model$fitted[[m]] - Xb_lm
    R_m  <- data$D[[m]] - Xr_without_lm
    XtR_m <- crossprod(X, R_m)

    model$residuals[[m]]        <- XtR_m
    model$fitted_without_l[[m]] <- Xr_without_lm
    model$raw_residuals[[m]]    <- R_m
  }
  model
}

#' Per-effect SER statistics on `mf_individual`
#'
#' Per outcome: `Bhat = X^T R / colSums(X^2)`,
#' `Shat^2 = sigma2 / colSums(X^2)` (single-effect-residual scalar
#' form). Mirrors `compute_ser_statistics.individual` (susieR);
#' the per-outcome lapply replaces the scalar T = 1 path.
#' `model$residuals` is set by the prior `compute_residuals` call
#' in the IBSS sweep; `data$xtx_diag` is cached at data
#' construction.
#'
#' Returned list shape (per outcome):
#'   `betahat[[m]]`  has shape `p x T_basis[m]` and
#'   `shat2[[m]]`    has shape `p x T_basis[m]`.
#' Plus the optim_init / optim_bounds / optim_scale fields the
#' susieR IBSS expects.
#'
#' @keywords internal
#' @noRd
compute_ser_statistics.mf_individual <- function(data, params, model, l, ...) {
  M <- data$M
  betahat <- vector("list", M)
  shat2   <- vector("list", M)
  for (m in seq_len(M)) {
    bs <- mf_per_outcome_bhat_shat(data, model, m)
    betahat[[m]] <- bs$bhat
    shat2[[m]]   <- bs$shat2
  }

  optim_init   <- log(max(c(unlist(betahat)^2 - unlist(shat2), 1),
                          na.rm = TRUE))
  optim_bounds <- c(-30, 15)
  optim_scale  <- "log"

  list(
    betahat      = betahat,
    shat2        = shat2,
    optim_init   = optim_init,
    optim_bounds = optim_bounds,
    optim_scale  = optim_scale
  )
}

# Per-outcome residual variance broadcast to per-position. Returns
# a length-T_basis[m] vector where each entry is the residual
# variance at that wavelet-coefficient position. Centralizes the
# scalar-vs-per-(scale, outcome) `model$sigma2` shape branch.
mf_sigma2_per_position <- function(data, model, m) {
  sigma2_m <- model$sigma2[[m]]
  T_pad    <- data$T_basis[m]
  if (length(sigma2_m) == 1L) {
    rep(sigma2_m, T_pad)
  } else {
    v <- numeric(T_pad)
    for (s in seq_along(sigma2_m)) {
      v[data$scale_index[[m]][[s]]] <- sigma2_m[s]
    }
    v
  }
}

# Per-outcome (Bhat, Shat) derivation from the cached residuals on
# `model$residuals[[m]]`. Shared by `compute_ser_statistics` and
# `calculate_posterior_moments` so the divide-by-pw and per-scale
# sigma2 broadcast lives in one place.
mf_per_outcome_bhat_shat <- function(data, model, m) {
  pw    <- data$xtx_diag
  XtR_m <- model$residuals[[m]]
  bhat_m <- XtR_m / pw
  sigma2_per_pos <- mf_sigma2_per_position(data, model, m)
  shat2_m <- outer(1 / pw, sigma2_per_pos)

  # Low-count mask: the IBSS treats flagged columns as
  # uninformative (Bhat = 0, Shat = 1) at every iteration.
  lowc <- data$lowc_idx[[m]]
  if (length(lowc) > 0L) {
    bhat_m [, lowc] <- 0
    shat2_m[, lowc] <- 1
  }

  list(bhat = bhat_m, shat2 = shat2_m)
}

#' Per-effect mixture-aware log-likelihood and weights on `mf_individual`
#'
#' For each outcome m and scale s, computes the per-(variable, scale)
#' log-Bayes factor via the cpp11 kernel
#' `mixture_log_bf_per_scale`. Sums across scales (within outcome)
#' to get a per-outcome log-BF vector of length p, combines across
#' outcomes via the S3 generic `combine_outcome_lbfs` (default
#' independence: `Reduce("+", ...)`), then computes alpha as the
#' softmax of `lbf_combined + log(model$pi)` and the SuSiE
#' `lbf_model = log(sum(pi * BF))` aggregate.
#'
#' Mirrors `loglik.individual` with K = 1, no NIG, and adds
#' the per-outcome / per-scale aggregation. Stores `model$alpha[l, ]`,
#' `model$lbf[l]`, `model$lbf_variable[l, ]` when `l` is non-NULL;
#' otherwise returns the aggregate `lbf_model` scalar (used by the
#' V-optim caller in `update_model_variance.mf_individual`).
#'
#' @references
#' Manuscript: methods/derivation.tex eq:post_alpha
#' Manuscript: methods/online_method.tex line 41 (cross-outcome combine)
#' @keywords internal
#' @noRd
loglik.mf_individual <- function(data, params, model, V, ser_stats, l = NULL, ...) {
  M <- data$M
  use_johnson <- isTRUE(params$small_sample_correction)
  df          <- params$small_sample_df
  outcome_lbfs <- vector("list", M)
  for (m in seq_len(M)) {
    Bhat_m <- ser_stats$betahat[[m]]
    Shat_m <- sqrt(ser_stats$shat2[[m]])
    G_m    <- model$G_prior[[m]]
    lbf_m  <- 0
    # Iterate G_prior groups directly; each entry carries the
    # column indices it covers (`$idx`). One uniform loop for
    # `per_outcome` (length 1) and `per_scale` (length S_m).
    for (s in seq_along(G_m)) {
      idx <- G_m[[s]]$idx
      g_s <- G_m[[s]]$fitted_g
      lbf_m <- lbf_m + if (use_johnson) {
        mixture_log_bf_per_scale_johnson(
          Bhat_m[, idx, drop = FALSE],
          Shat_m[, idx, drop = FALSE],
          g_s$sd, g_s$pi, V_scale = V, df = df)
      } else {
        mixture_log_bf_per_scale(
          Bhat_m[, idx, drop = FALSE],
          Shat_m[, idx, drop = FALSE],
          g_s$sd, g_s$pi, V_scale = V)
      }
    }
    outcome_lbfs[[m]] <- lbf_m
  }

  lbf_combined <- combine_outcome_lbfs(model$cross_outcome_combiner,
                                        outcome_lbfs, model)

  # Stable softmax with zero-pw handling. Mirrors susieR's
  # `lbf_stabilization` + `compute_posterior_weights` (7 lines vs 8;
  # inlined to keep the per-(outcome, position) shat2 logic explicit
  # rather than passing a marker vector to the susieR helper).
  zero_pw <- data$xtx_diag == 0
  if (any(zero_pw)) lbf_combined[zero_pw] <- 0
  lpo       <- lbf_combined + log(model$pi + sqrt(.Machine$double.eps))
  m_max     <- max(lpo)
  w         <- exp(lpo - m_max)
  alpha     <- w / sum(w)
  lbf_model <- m_max + log(sum(w))

  if (!is.null(l)) {
    model$alpha[l, ]        <- alpha
    model$lbf[l]            <- lbf_model
    model$lbf_variable[l, ] <- lbf_combined
    return(model)
  }
  lbf_model
}

#' Per-effect posterior moments on `mf_individual` (mixture-aware)
#'
#' For each outcome m and scale s, calls the cpp11 kernel
#' `mixture_posterior_per_scale` to compute per-(variable, position)
#' posterior mean and second moment under the per-scale mixture-of-
#' normals prior, then writes the result into the corresponding
#' columns of `model$mu[[l]][[m]]` and `model$mu2[[l]][[m]]`. Bhat
#' and Shat are re-derived from cached `model$residuals[[m]]` and
#' `data$xtx_diag` (mirrors `calculate_posterior_moments.individual`).
#'
#' Manuscript references:
#'   methods/derivation.tex eq:post_f_mix
#'   methods/derivation.tex eq:post_f2_mix
#'
#' @keywords internal
#' @noRd
calculate_posterior_moments.mf_individual <- function(data, params, model, V, l, ...) {
  p <- data$p
  for (m in seq_len(data$M)) {
    bs <- mf_per_outcome_bhat_shat(data, model, m)
    bhat_m <- bs$bhat
    shat_m <- sqrt(bs$shat2)

    G_m    <- model$G_prior[[m]]
    mu_lm  <- matrix(0, nrow = p, ncol = data$T_basis[m])
    mu2_lm <- matrix(0, nrow = p, ncol = data$T_basis[m])
    for (s in seq_along(G_m)) {
      idx <- G_m[[s]]$idx
      g_s <- G_m[[s]]$fitted_g
      out <- mixture_posterior_per_scale(
        bhat_m[, idx, drop = FALSE],
        shat_m[, idx, drop = FALSE],
        g_s$sd, g_s$pi, V_scale = V)
      mu_lm[, idx]  <- out$pmean
      mu2_lm[, idx] <- out$pmean2
    }
    model$mu[[l]][[m]]  <- mu_lm
    model$mu2[[l]][[m]] <- mu2_lm
  }
  model
}

#' Per-effect SER posterior expected log-likelihood
#'
#' Generalizes `SER_posterior_e_loglik.individual` to
#' per-outcome, per-(scale, outcome) residual variance. For each
#' outcome m,
#' `L_m = sum_t [ -0.5 * n * log(2 pi sigma2_t)
#'                - 0.5 / sigma2_t * (sum_n R[n,t]^2
#'                                     - 2 sum_n R[n,t] * (X * Eb)[n,t]
#'                                     + sum_j pw[j] * Eb2[j,t]) ]`,
#' where `sigma2_t` is broadcast from `model$sigma2[[m]]` via
#' `data$scale_index[[m]]`, `Eb_m = alpha_l * mu[[l]][[m]]`, and
#' `Eb2_m = alpha_l * mu2[[l]][[m]]`.
#'
#' Returns the scalar `sum_m L_m`.
#'
#' @references
#' Manuscript: methods/online_method.tex (per-effect ELBO contribution).
#' @keywords internal
#' @noRd
SER_posterior_e_loglik.mf_individual <- function(data, params, model, l, ...) {
  X       <- data$X
  alpha_l <- model$alpha[l, ]
  pw      <- data$xtx_diag
  out     <- 0
  for (m in seq_len(data$M)) {
    Eb_m  <- alpha_l * model$mu[[l]][[m]]
    Eb2_m <- alpha_l * model$mu2[[l]][[m]]
    R_m   <- model$raw_residuals[[m]]
    n     <- nrow(R_m)
    sigma2_per_pos <- mf_sigma2_per_position(data, model, m)

    XEb_m <- X %*% Eb_m
    sumR2_t    <- colSums(R_m * R_m)
    sumRXEb_t  <- colSums(R_m * XEb_m)
    sumpwEb2_t <- colSums(pw * Eb2_m)
    quad_t     <- sumR2_t - 2 * sumRXEb_t + sumpwEb2_t

    out <- out + sum(-0.5 * n * log(2 * pi * sigma2_per_pos) -
                     0.5 * quad_t / sigma2_per_pos)
  }
  out
}

#' Per-effect KL divergence on `mf_individual`
#'
#' Mirrors `compute_kl.individual`'s "Standard Gaussian KL"
#' branch generalized across outcomes and per-(scale, outcome)
#' residual variance:
#' `KL_l = -[ lbf_l + sum_m sum_n sum_t log dnorm(R_m[n, t]; 0, sqrt(sigma2_t)) ]
#'        + SER_posterior_e_loglik(l)`,
#' where `sigma2_t` comes from `mf_sigma2_per_position`.
#'
#' @references
#' Manuscript: methods/online_method.tex (per-effect ELBO contribution).
#' @keywords internal
#' @noRd
compute_kl.mf_individual <- function(data, params, model, l, ...) {
  L_null <- 0
  for (m in seq_len(data$M)) {
    R_m <- model$raw_residuals[[m]]
    n   <- nrow(R_m)
    sigma2_per_pos <- mf_sigma2_per_position(data, model, m)
    sumR2_t <- colSums(R_m * R_m)
    L_null  <- L_null + sum(-0.5 * n * log(2 * pi * sigma2_per_pos) -
                            0.5 * sumR2_t / sigma2_per_pos)
  }
  loglik_term <- model$lbf[l] + L_null
  # Generic dispatch (not the explicit `.mf_individual`) so any future
  # subclass of `mf_individual` can override.
  kl <- -loglik_term + SER_posterior_e_loglik(data, params, model, l)
  model$KL[l] <- kl
  model
}

#' Negative log-likelihood for the V-optim caller
#'
#' Mirrors `neg_loglik.individual`. Converts the optim
#' parameter from log-scale to V (per `compute_ser_statistics`'s
#' `optim_scale = "log"`) and returns `-loglik(...)`.
#'
#' @keywords internal
#' @noRd
neg_loglik.mf_individual <- function(data, params, model, V_param, ser_stats, ...) {
  V <- exp(V_param)
  -loglik.mf_individual(data, params, model, V, ser_stats)
}

#' Per-outcome expected squared residuals
#'
#' Returns the length-`T_basis[m]` vector of per-position
#' `E_q[||Y_m[, t] - sum_l X * b_lm[, t]||^2]` for the SER posterior:
#'
#' `ER2_m[t] = ||Y_m[, t] - sum_l X * (alpha_l * mu_lm)[, t]||^2
#'           + sum_l ((alpha_l * pw)^T * mu2_l_m[, t]
#'                    - ||X * (alpha_l * mu_lm[, t])||^2)`
#'
#' The first term reuses `model$fitted[[m]]`, the running per-outcome
#' fit maintained by `update_fitted_values`. Used by
#' `update_variance_components` and `Eloglik`.
#'
#' @references
#' Manuscript: methods/derivation.tex eq:ERSS.
#' @keywords internal
#' @noRd
mf_get_ER2_per_position <- function(data, model, m) {
  X     <- data$X
  pw    <- data$xtx_diag
  T_pad <- data$T_basis[m]

  # Cached running fit: model$fitted[[m]] = X %*% sum_l (alpha_l * mu_lm).
  res_m <- data$D[[m]] - model$fitted[[m]]
  rss_t <- colSums(res_m * res_m)

  # Bias correction summed over l.
  bias_t <- numeric(T_pad)
  for (l in seq_len(model$L)) {
    alpha_l <- model$alpha[l, ]
    mu_lm  <- model$mu[[l]][[m]]
    mu2_l_m <- model$mu2[[l]][[m]]
    XEb_lm <- X %*% (alpha_l * mu_lm)
    bias_t  <- bias_t + colSums((alpha_l * pw) * mu2_l_m) -
               colSums(XEb_lm * XEb_lm)
  }
  rss_t + bias_t
}

#' Update per-outcome (or per-(scale, outcome)) residual variance
#'
#' Computes `sigma2_m` per outcome using
#' `mf_get_ER2_per_position`, then aggregates to the shape selected
#' by `params$residual_variance_scope`:
#' - `"per_outcome"` (shared per outcome): `sigma2_m = sum_t ER2_m[t] / (n * T_pad[m])`
#' - `"per_scale"` (default): per scale `s`,
#'   `sigma2_{m,s} = sum_{t in idx_s} ER2_m[t] / (n * |idx_s|)`
#'
#' @references
#' Manuscript: methods/derivation.tex eq:residual_variance_per_scale.
#' @keywords internal
#' @noRd
update_variance_components.mf_individual <- function(data, params, model, ...) {
  method <- params$residual_variance_scope %||% "per_scale"
  if (!method %in% c("per_scale", "per_outcome")) {
    stop("`residual_variance_scope` must be one of 'per_scale' or 'per_outcome'.")
  }
  for (m in seq_len(data$M)) {
    er2_t <- mf_get_ER2_per_position(data, model, m)
    n     <- data$n
    if (method == "per_outcome") {
      model$sigma2[[m]] <- sum(er2_t) / (n * data$T_basis[m])
    } else {
      indx_m <- data$scale_index[[m]]
      sigma2_per_scale <- vapply(indx_m, function(idx) {
        sum(er2_t[idx]) / (n * length(idx))
      }, numeric(1))
      model$sigma2[[m]] <- sigma2_per_scale
    }
  }
  model
}

#' Per-effect prior update for `mf_individual` (mixsqp on `pi_V`)
#'
#' Overrides `optimize_prior_variance` for the
#' `mf_individual` data class. Per (outcome, scale), builds the
#' mixsqp likelihood matrix (`mf_em_likelihood_per_scale`) from the
#' per-outcome (Bhat, Shat) in `ser_stats`, runs the M step
#' (`mf_em_m_step_per_scale`) with `zeta = model$alpha[l, ]`, and
#' writes the new mixture weights into both
#' `model$G_prior[[m]][[s]]$fitted_g$pi` and `model$pi_V[[m]][s, ]`.
#'
#' The susieR-style scalar `V[l]` is held at 1 since the
#' mixture weights absorb the per-effect prior shape adaptation.
#' Returns `list(V = 1, model = updated_model)` per the susieR
#' generic contract.
#'
#' @references
#' Manuscript: methods/derivation.tex line 216
#' (factorized empirical-Bayes mixture-weight update).
#' @keywords internal
#' @noRd
optimize_prior_variance.mf_individual <- function(data, params, model, ser_stats,
                                                  l       = NULL,
                                                  alpha   = NULL,
                                                  moments = NULL,
                                                  V_init  = NULL) {
  if (is.null(l)) {
    stop("`optimize_prior_variance.mf_individual` requires the effect index `l`.")
  }
  # The joint per-effect SNP posterior `model$alpha[l, ]` is the
  # softmax of the joint log-Bayes-factor across all M outcomes
  # and S_m scales. Adding M outcomes increases the variance of
  # the per-SNP joint lbf by ~M, so the dominant alpha values
  # concentrate ~M-fold relative to the M = 1 case. The mixsqp
  # M-step's regularization-to-data balance is set by
  # `nullweight / max_alpha`; to hold this ratio fixed across M
  # we scale `mixsqp_null_penalty` by M. See
  # `inst/notes/cross-package-audit-derivations.md` section 1
  # for the full derivation.
  mixsqp_null_penalty <- (params$mixsqp_null_penalty %||% 0.1) *
                        max(1L, data$M)
  control    <- params$control_mixsqp %||% list()
  zeta_l     <- model$alpha[l, ]

  for (m in seq_len(data$M)) {
    bhat_m <- ser_stats$betahat[[m]]
    shat_m <- sqrt(ser_stats$shat2[[m]])
    G_m    <- model$G_prior[[m]]
    # Uniform group loop: each G_prior entry holds the column
    # indices it covers (`$idx`). For `per_outcome` there is one
    # group spanning every wavelet column (one mixsqp solve over
    # the full outcome); for `per_scale` there are S_m
    # groups, one per wavelet scale (S_m independent mixsqp solves).
    for (s in seq_along(G_m)) {
      idx     <- G_m[[s]]$idx
      sd_grid <- G_m[[s]]$fitted_g$sd
      L_mat <- mf_em_likelihood_per_scale(
        bhat_m[, idx, drop = FALSE],
        shat_m[, idx, drop = FALSE],
        sd_grid)
      new_pi <- mf_em_m_step_per_scale(
        L_mat, zeta_l, idx_size = length(idx),
        mixsqp_null_penalty = mixsqp_null_penalty,
        control_mixsqp = control)
      model$G_prior[[m]][[s]]$fitted_g$pi <- new_pi
      model$pi_V[[m]][s, ]                 <- new_pi
    }
  }
  list(V = 1, model = model)
}

#' Per-iteration residual variance + derived-quantity update for `mf_individual`
#'
#' Mirrors `update_model_variance.default`'s orchestration:
#' calls `update_variance_components` (refresh per-outcome or
#' per-(scale, outcome) sigma2), then `update_derived_quantities`
#' (refresh the running per-outcome fit). Bounds are not applied
#' per-element here because mfsusieR's sigma2 is a list of
#' length-`S_m` vectors per outcome; `params$residual_variance_lowerbound`

#' refinement.
#'
#' @keywords internal
#' @noRd
update_model_variance.mf_individual <- function(data, params, model, ...) {
  if (!isTRUE(params$estimate_residual_variance)) return(model)
  model <- update_variance_components(data, params, model)
  update_derived_quantities(data, params, model)
}

#' Recompute the running per-outcome fit `model$fitted[[m]]`
#'
#' Called by `update_model_variance` after a sigma2 update. The
#' fit is `X %*% sum_l alpha_l * mu_l[[m]]` per outcome. Mirrors
#' `update_derived_quantities.default`.
#'
#' @keywords internal
#' @noRd
update_derived_quantities.mf_individual <- function(data, params, model, ...) {
  for (m in seq_len(data$M)) {
    postF_m <- matrix(0, nrow = data$p, ncol = data$T_basis[m])
    for (l in seq_len(model$L)) {
      postF_m <- postF_m + model$alpha[l, ] * model$mu[[l]][[m]]
    }
    model$fitted[[m]] <- data$X %*% postF_m
  }
  model
}

#' Aggregate expected log-likelihood across outcomes
#'
#' Mirrors `Eloglik.individual`:
#' `Eloglik = sum_m sum_t -0.5 * n * log(2 pi sigma2_t) - 0.5 / sigma2_t * ER2_m[t]`
#'
#' @references
#' Manuscript: methods/online_method.tex (Eloglik aggregate).
#' @keywords internal
#' @noRd
Eloglik.mf_individual <- function(data, model, ...) {
  out <- 0
  for (m in seq_len(data$M)) {
    er2_t <- mf_get_ER2_per_position(data, model, m)
    s2_t  <- mf_sigma2_per_position(data, model, m)
    out   <- out + sum(-0.5 * data$n * log(2 * pi * s2_t) -
                       0.5 * er2_t / s2_t)
  }
  out
}

#' ELBO for an `mf_individual` fit
#'
#' Mirrors `get_objective.individual`:
#' `ELBO = Eloglik - sum_l KL_l`. The per-effect KL stored in
#' `model$KL[l]` (set by `compute_kl.mf_individual`) already
#' contains both the Gaussian component
#' (`KL(q(beta_l|gamma_l) || p(beta_l|gamma_l))`) and the
#' categorical mixture component
#' (`KL(q(gamma_l) || p(gamma_l))`), so no separate entropy
#' correction is needed here.
#'
#' @references
#' Manuscript: methods/online_method.tex (ELBO aggregate).
#' @keywords internal
#' @noRd
get_objective.mf_individual <- function(data, params, model, ...) {
  Eloglik.mf_individual(data, model) -
    if (all(is.na(model$KL))) 0 else sum(model$KL, na.rm = TRUE)
}

#' Fold the per-effect posterior into the running fit
#'
#' Per outcome, sets
#' `model$fitted[[m]] <- model$fitted_without_l[[m]] + X %*% (alpha_l * mu_lm)`,
#' where `mu_lm` is the per-effect, per-outcome posterior mean
#' (a `p x T_basis[m]` matrix). Called after
#' `calculate_posterior_moments` updates alpha_l, mu_l. Mirrors
#' the inline susie / mvsusie pattern of recomputing Xr from
#' Xr_without_l + the new effect contribution.
#'
#' @keywords internal
#' @noRd
update_fitted_values.mf_individual <- function(data, params, model, l, ...) {
  X <- data$X
  M <- data$M
  alpha_l <- model$alpha[l, ]
  for (m in seq_len(M)) {
    b_lm  <- alpha_l * model$mu[[l]][[m]]
    Xb_lm <- X %*% b_lm
    model$fitted[[m]] <- model$fitted_without_l[[m]] + Xb_lm
  }
  model
}
