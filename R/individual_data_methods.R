# Per-effect S3 methods on `mf_individual` for the IBSS dispatch.
#
# Mirrors the susieR / mvsusieR pattern: compute_residuals computes
# the per-effect partial residual `R_l = Y - Xr_without_l` and
# caches `model$residuals = X^T R_l` for downstream SER stats.
# compute_ser_statistics turns the cached `XtR` into per-position
# `Bhat / Shat^2`. update_fitted_values folds the new posterior
# back into the running per-modality fit `model$fitted[[m]]`.
#
# All three are per-modality lapply over `seq_len(data$M)` plus
# matrix-level BLAS ops; no C++ needed at this layer.
#
# Manuscript references: methods/algorithms.tex eq:partial_resid;
# methods/derivation.tex eq:bhat_shat_per_modality.

#' Per-effect partial residual on `mf_individual`
#'
#' For each modality m, computes the contribution of effect `l`
#' to the running fit `Xb_l_m = X %*% (alpha_l * mu_l[[m]])`,
#' removes it from `model$fitted[[m]]` to obtain the without-l
#' fit, and subtracts that from the wavelet matrix `data$D[[m]]`
#' to get the residual. Caches `XtR_m = X^T R_m` per modality on
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
    # alpha_l . mu_l[[m]] is J x T_padded[m]; X is n x J. Xb is n x T_padded[m].
    b_l_m  <- alpha_l * model$mu[[l]][[m]]
    Xb_l_m <- X %*% b_l_m

    Xr_without_l_m <- model$fitted[[m]] - Xb_l_m
    R_m  <- data$D[[m]] - Xr_without_l_m
    XtR_m <- crossprod(X, R_m)

    model$residuals[[m]]        <- XtR_m
    model$fitted_without_l[[m]] <- Xr_without_l_m
    model$raw_residuals[[m]]    <- R_m
  }
  model
}

#' Per-effect SER statistics on `mf_individual`
#'
#' Per modality: `Bhat = X^T R / colSums(X^2)`,
#' `Shat^2 = sigma2 / colSums(X^2)` (single-effect-residual scalar
#' form). Mirrors `compute_ser_statistics.individual` (susieR);
#' the per-modality lapply replaces the scalar T = 1 path.
#' `model$residuals` is set by the prior `compute_residuals` call
#' in the IBSS sweep; `data$predictor_weights` is cached at data
#' construction.
#'
#' Returned list shape (per modality):
#'   `betahat[[m]]`  has shape `J x T_padded[m]` and
#'   `shat2[[m]]`    has shape `J x T_padded[m]`.
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
    bs <- mf_per_modality_bhat_shat(data, model, m)
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

# Per-modality residual variance broadcast to per-position. Returns
# a length-T_padded[m] vector where each entry is the residual
# variance at that wavelet-coefficient position. Centralizes the
# scalar-vs-per-(scale, modality) `model$sigma2` shape branch.
mf_sigma2_per_position <- function(data, model, m) {
  sigma2_m <- model$sigma2[[m]]
  T_pad    <- data$T_padded[m]
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

# Per-modality (Bhat, Shat) derivation from the cached residuals on
# `model$residuals[[m]]`. Shared by `compute_ser_statistics` and
# `calculate_posterior_moments` so the divide-by-pw and per-scale
# sigma2 broadcast lives in one place.
mf_per_modality_bhat_shat <- function(data, model, m) {
  pw    <- data$predictor_weights
  XtR_m <- model$residuals[[m]]
  bhat_m <- XtR_m / pw
  sigma2_per_pos <- mf_sigma2_per_position(data, model, m)
  shat2_m <- outer(1 / pw, sigma2_per_pos)
  list(bhat = bhat_m, shat2 = shat2_m)
}

#' Per-effect mixture-aware log-likelihood and weights on `mf_individual`
#'
#' For each modality m and scale s, computes the per-(SNP, scale)
#' log-Bayes factor via the cpp11 kernel
#' `mixture_log_bf_per_scale`. Sums across scales (within modality)
#' to get a per-modality log-BF vector of length J, combines across
#' modalities via the S3 generic `combine_modality_lbfs` (default
#' independence: `Reduce("+", ...)`), then computes alpha as the
#' softmax of `lbf_combined + log(model$pi)` and the SuSiE
#' `lbf_model = log(sum(pi * BF))` aggregate.
#'
#' Mirrors `susieR::loglik.individual` with K = 1, no NIG, and adds
#' the per-modality / per-scale aggregation. Stores `model$alpha[l, ]`,
#' `model$lbf[l]`, `model$lbf_variable[l, ]` when `l` is non-NULL;
#' otherwise returns the aggregate `lbf_model` scalar (used by the
#' V-optim caller in `update_model_variance.mf_individual`, PR group 6).
#'
#' @references
#' Manuscript: methods/derivation.tex eq:post_alpha
#' Manuscript: methods/online_method.tex line 41 (cross-modality combine)
#' @keywords internal
#' @noRd
loglik.mf_individual <- function(data, params, model, V, ser_stats, l = NULL, ...) {
  M <- data$M
  modality_lbfs <- vector("list", M)
  for (m in seq_len(M)) {
    Bhat_m <- ser_stats$betahat[[m]]
    Shat_m <- sqrt(ser_stats$shat2[[m]])
    G_m    <- model$G_prior[[m]]
    indx_m <- data$scale_index[[m]]
    lbf_m  <- 0
    for (s in seq_along(indx_m)) {
      idx <- indx_m[[s]]
      g_s <- G_m[[s]]$fitted_g
      lbf_m <- lbf_m + mixture_log_bf_per_scale(
        Bhat_m[, idx, drop = FALSE],
        Shat_m[, idx, drop = FALSE],
        g_s$sd, g_s$pi, V_scale = V)
    }
    modality_lbfs[[m]] <- lbf_m
  }

  lbf_combined <- combine_modality_lbfs(model$cross_modality_combiner,
                                        modality_lbfs, model)

  # Stable softmax with zero-pw handling. Mirrors susieR's
  # `lbf_stabilization` + `compute_posterior_weights` (7 lines vs 8;
  # inlined to keep the per-(modality, position) shat2 logic explicit
  # rather than passing a marker vector to the susieR helper).
  zero_pw <- data$predictor_weights == 0
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
#' For each modality m and scale s, calls the cpp11 kernel
#' `mixture_posterior_per_scale` to compute per-(SNP, position)
#' posterior mean and second moment under the per-scale mixture-of-
#' normals prior, then writes the result into the corresponding
#' columns of `model$mu[[l]][[m]]` and `model$mu2[[l]][[m]]`. Bhat
#' and Shat are re-derived from cached `model$residuals[[m]]` and
#' `data$predictor_weights` (mirrors `susieR::calculate_posterior_moments.individual`).
#'
#' Manuscript references:
#'   methods/derivation.tex eq:post_f_mix
#'   methods/derivation.tex eq:post_f2_mix
#'
#' @keywords internal
#' @noRd
calculate_posterior_moments.mf_individual <- function(data, params, model, V, l, ...) {
  J <- data$J
  for (m in seq_len(data$M)) {
    bs <- mf_per_modality_bhat_shat(data, model, m)
    bhat_m <- bs$bhat
    shat_m <- sqrt(bs$shat2)

    G_m    <- model$G_prior[[m]]
    indx_m <- data$scale_index[[m]]
    mu_lm  <- matrix(0, nrow = J, ncol = data$T_padded[m])
    mu2_lm <- matrix(0, nrow = J, ncol = data$T_padded[m])
    for (s in seq_along(indx_m)) {
      idx <- indx_m[[s]]
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
#' Generalizes `susieR::SER_posterior_e_loglik.individual` to
#' per-modality, per-(scale, modality) residual variance. For each
#' modality m,
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
  pw      <- data$predictor_weights
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
#' Mirrors `susieR::compute_kl.individual`'s "Standard Gaussian KL"
#' branch generalized across modalities and per-(scale, modality)
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
#' Mirrors `susieR::neg_loglik.individual`. Converts the optim
#' parameter from log-scale to V (per `compute_ser_statistics`'s
#' `optim_scale = "log"`) and returns `-loglik(...)`.
#'
#' @keywords internal
#' @noRd
neg_loglik.mf_individual <- function(data, params, model, V_param, ser_stats, ...) {
  V <- exp(V_param)
  -loglik.mf_individual(data, params, model, V, ser_stats)
}

#' Per-modality expected squared residuals
#'
#' Returns the length-`T_padded[m]` vector of per-position
#' `E_q[||Y_m[, t] - sum_l X * b_l_m[, t]||^2]` for the SER posterior:
#'
#' `ER2_m[t] = ||Y_m[, t] - sum_l X * (alpha_l * mu_l_m)[, t]||^2
#'           + sum_l ((alpha_l * pw)^T * mu2_l_m[, t]
#'                    - ||X * (alpha_l * mu_l_m[, t])||^2)`
#'
#' The first term reuses `model$fitted[[m]]`, the running per-modality
#' fit maintained by `update_fitted_values`. Used by
#' `update_variance_components` and `Eloglik`.
#'
#' @references
#' Manuscript: methods/derivation.tex eq:ERSS.
#' @keywords internal
#' @noRd
mf_get_ER2_per_position <- function(data, model, m) {
  X     <- data$X
  pw    <- data$predictor_weights
  T_pad <- data$T_padded[m]

  # Cached running fit: model$fitted[[m]] = X %*% sum_l (alpha_l * mu_l_m).
  res_m <- data$D[[m]] - model$fitted[[m]]
  rss_t <- colSums(res_m * res_m)

  # Bias correction summed over l.
  bias_t <- numeric(T_pad)
  for (l in seq_len(model$L)) {
    alpha_l <- model$alpha[l, ]
    mu_l_m  <- model$mu[[l]][[m]]
    mu2_l_m <- model$mu2[[l]][[m]]
    XEb_l_m <- X %*% (alpha_l * mu_l_m)
    bias_t  <- bias_t + colSums((alpha_l * pw) * mu2_l_m) -
               colSums(XEb_l_m * XEb_l_m)
  }
  rss_t + bias_t
}

#' Update per-modality (or per-(scale, modality)) residual variance
#'
#' Computes `sigma2_m` per modality using
#' `mf_get_ER2_per_position`, then aggregates to the shape selected
#' by `params$residual_variance_method`:
#' - `"shared_per_modality"` (legacy): `sigma2_m = sum_t ER2_m[t] / (n * T_pad[m])`
#' - `"per_scale_modality"` (default): per scale `s`,
#'   `sigma2_{m,s} = sum_{t in idx_s} ER2_m[t] / (n * |idx_s|)`
#'
#' @references
#' Manuscript: methods/derivation.tex eq:residual_variance_per_scale.
#' @keywords internal
#' @noRd
update_variance_components.mf_individual <- function(data, params, model, ...) {
  method <- params$residual_variance_method %||% "per_scale_modality"
  if (!method %in% c("per_scale_modality", "shared_per_modality")) {
    stop("`residual_variance_method` must be one of 'per_scale_modality' or 'shared_per_modality'.")
  }
  for (m in seq_len(data$M)) {
    er2_t <- mf_get_ER2_per_position(data, model, m)
    n     <- data$n
    if (method == "shared_per_modality") {
      model$sigma2[[m]] <- sum(er2_t) / (n * data$T_padded[m])
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
#' Overrides `susieR::optimize_prior_variance` for the
#' `mf_individual` data class. Per (modality, scale), builds the
#' mixsqp likelihood matrix (`mf_em_likelihood_per_scale`) from the
#' per-modality (Bhat, Shat) in `ser_stats`, runs the M step
#' (`mf_em_m_step_per_scale`) with `zeta = model$alpha[l, ]`, and
#' writes the new mixture weights into both
#' `model$G_prior[[m]][[s]]$fitted_g$pi` and `model$pi_V[[m]][s, ]`.
#'
#' The susieR-style scalar `V[l]` is held at 1 here in v1 since the
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
  nullweight <- params$nullweight     %||% 0.7
  control    <- params$control_mixsqp %||% list()
  zeta_l     <- model$alpha[l, ]

  for (m in seq_len(data$M)) {
    bhat_m <- ser_stats$betahat[[m]]
    shat_m <- sqrt(ser_stats$shat2[[m]])
    indx_m <- data$scale_index[[m]]

    for (s in seq_along(indx_m)) {
      idx <- indx_m[[s]]
      sd_grid <- model$G_prior[[m]][[s]]$fitted_g$sd
      L_mat <- mf_em_likelihood_per_scale(
        bhat_m[, idx, drop = FALSE],
        shat_m[, idx, drop = FALSE],
        sd_grid)
      new_pi <- mf_em_m_step_per_scale(
        L_mat, zeta_l, idx_size = length(idx),
        nullweight = nullweight,
        control_mixsqp = control)
      model$G_prior[[m]][[s]]$fitted_g$pi <- new_pi
      model$pi_V[[m]][s, ]                 <- new_pi
    }
  }
  list(V = 1, model = model)
}

#' Per-iteration residual variance + derived-quantity update for `mf_individual`
#'
#' Mirrors `susieR::update_model_variance.default`'s orchestration:
#' calls `update_variance_components` (refresh per-modality or
#' per-(scale, modality) sigma2), then `update_derived_quantities`
#' (refresh the running per-modality fit). Bounds are not applied
#' per-element here because mfsusieR's sigma2 is a list of
#' length-`S_m` vectors per modality; `params$residual_variance_lowerbound`
#' / `_upperbound` semantics for that shape are deferred to a later
#' refinement.
#'
#' @keywords internal
#' @noRd
update_model_variance.mf_individual <- function(data, params, model, ...) {
  if (!isTRUE(params$estimate_residual_variance)) return(model)
  model <- update_variance_components(data, params, model)
  update_derived_quantities(data, params, model)
}

#' Recompute the running per-modality fit `model$fitted[[m]]`
#'
#' Called by `update_model_variance` after a sigma2 update. The
#' fit is `X %*% sum_l alpha_l * mu_l[[m]]` per modality. Mirrors
#' `susieR::update_derived_quantities.default`.
#'
#' @keywords internal
#' @noRd
update_derived_quantities.mf_individual <- function(data, params, model, ...) {
  for (m in seq_len(data$M)) {
    postF_m <- matrix(0, nrow = data$J, ncol = data$T_padded[m])
    for (l in seq_len(model$L)) {
      postF_m <- postF_m + model$alpha[l, ] * model$mu[[l]][[m]]
    }
    model$fitted[[m]] <- data$X %*% postF_m
  }
  model
}

#' Aggregate expected log-likelihood across modalities
#'
#' Mirrors `susieR::Eloglik.individual`:
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
#' Mirrors `susieR::get_objective.individual`:
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
#' Per modality, sets
#' `model$fitted[[m]] <- model$fitted_without_l[[m]] + X %*% (alpha_l * mu_l_m)`,
#' where `mu_l_m` is the per-effect, per-modality posterior mean
#' (a `J x T_padded[m]` matrix). Called after
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
    b_l_m  <- alpha_l * model$mu[[l]][[m]]
    Xb_l_m <- X %*% b_l_m
    model$fitted[[m]] <- model$fitted_without_l[[m]] + Xb_l_m
  }
  model
}
