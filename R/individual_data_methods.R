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

  # `model$residuals`, `fitted_without_l`, `raw_residuals` are
  # preallocated as length-M lists by
  # `initialize_susie_model.mf_individual` and overwritten in
  # place per effect.

  alpha_l <- model$alpha[l, ]
  for (m in seq_len(M)) {
    # alpha_l . mu_l[[m]] is p x T_basis[m]; X is n x p. Xb is n x T_basis[m].
    b_lm  <- alpha_l * model$mu[[l]][[m]]
    Xb_lm <- X %*% b_lm

    Xr_without_lm <- model$fitted[[m]] - Xb_lm
    R_m  <- data$D[[m]] - Xr_without_lm

    idx_m <- data$na_idx[[m]]
    XtR_m <- crossprod(X[idx_m, , drop = FALSE], R_m[idx_m, , drop = FALSE])

    model$residuals[[m]]        <- XtR_m
    model$fitted_without_l[[m]] <- Xr_without_lm
    model$raw_residuals[[m]]    <- R_m        # keep full n x T for reference
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

# Per-outcome (Bhat, Shat) from cached residuals. shat2_m only depends
# on xtx_diag (fixed) and sigma2[[m]] (changes at the per-iter sigma2
# update), so we read `model$iter_cache$shat2[[m]]` populated by
# `refresh_iter_cache` instead of redoing the outer product per (l, m).
mf_per_outcome_bhat_shat <- function(data, model, m) {
  pw    <- data$xtx_diag_list[[m]]
  XtR_m <- model$residuals[[m]]
  bhat_m <- XtR_m / pw

  cached_shat2 <- model$iter_cache$shat2[[m]]
  shat2_m <- if (!is.null(cached_shat2) &&
                 identical(dim(cached_shat2),
                           c(data$p, data$T_basis[m]))) {
    cached_shat2
  } else {
    # Fallback before iter_cache is populated (e.g., first SER call)
    # or in tests that bypass the IBSS orchestrator.
    sigma2_per_pos <- mf_sigma2_per_position(data, model, m)
    outer(1 / pw, sigma2_per_pos)
  }

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
loglik.mf_individual <- function(data, params, model, V, ser_stats, l = NULL,
                                  update_alpha = TRUE, ...) {
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
      b_sub <- Bhat_m[, idx, drop = FALSE]
      s_sub <- Shat_m[, idx, drop = FALSE]
      lbf_m <- lbf_m + if (inherits(g_s, "laplacemix")) {
        mixture_log_bf_laplace_per_scale(b_sub, s_sub,
                                          fitted_g = g_s, V_scale = V)
      } else if (use_johnson) {
        mixture_log_bf_per_scale_johnson(b_sub, s_sub,
                                          g_s$sd, g_s$pi,
                                          V_scale = V, df = df)
      } else {
        mixture_log_bf_per_scale(b_sub, s_sub,
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
    # `update_alpha = FALSE` is used by the post-iteration KL refresh
    # in `get_objective.mfsusie`: we want a fresh per-effect lbf[l]
    # against the iter-final pi_V, but the variational posterior's
    # alpha[l] is part of the q_t we are evaluating against and must
    # stay frozen during the refresh.
    if (update_alpha) model$alpha[l, ] <- alpha
    model$lbf[l]            <- lbf_model
    model$lbf_variable[l, ] <- lbf_combined

    # Stash ser_stats so the next call in susieR's SER orchestrator
    # -- calculate_posterior_moments, which does not receive
    # ser_stats -- can reuse it. Keyed by `l` to avoid stale reads.
    # refresh_lbf_kl calls with update_alpha = FALSE and never hits
    # CPM, so we skip the stash there.
    if (update_alpha) {
      model$ser_cache <- list(l = l,
                              betahat = ser_stats$betahat,
                              shat2   = ser_stats$shat2)
    }

    # Persist per-outcome BFs. The slot is pre-allocated in
    # `initialize_susie_model.mf_individual`. Each `outcome_lbfs[[m]]`
    # is a length-p vector of log BFs for outcome m, summed across
    # scales.
    for (m in seq_len(M)) {
      model$lbf_variable_outcome[l, , m] <- outcome_lbfs[[m]]
    }

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
  # Reuse the (bhat, shat2) stashed by the prior loglik.mf_individual
  # call (same SER orchestrator step) instead of recomputing.
  use_cache <- !is.null(model$ser_cache) &&
                identical(model$ser_cache$l, l)
  for (m in seq_len(data$M)) {
    if (use_cache) {
      bhat_m  <- model$ser_cache$betahat[[m]]
      shat2_m <- model$ser_cache$shat2[[m]]
    } else {
      bs      <- mf_per_outcome_bhat_shat(data, model, m)
      bhat_m  <- bs$bhat
      shat2_m <- bs$shat2
    }
    shat_m <- sqrt(shat2_m)

    G_m    <- model$G_prior[[m]]
    mu_lm  <- matrix(0, nrow = p, ncol = data$T_basis[m])
    mu2_lm <- matrix(0, nrow = p, ncol = data$T_basis[m])
    for (s in seq_along(G_m)) {
      idx   <- G_m[[s]]$idx
      g_s   <- G_m[[s]]$fitted_g
      b_sub <- bhat_m[, idx, drop = FALSE]
      s_sub <- shat_m[, idx, drop = FALSE]
      out <- if (inherits(g_s, "laplacemix")) {
        mixture_posterior_laplace_per_scale(b_sub, s_sub,
                                             fitted_g = g_s, V_scale = V)
      } else {
        mixture_posterior_per_scale(b_sub, s_sub,
                                     g_s$sd, g_s$pi, V_scale = V)
      }
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
  out     <- 0
  for (m in seq_len(data$M)) {
    idx_m <- data$na_idx[[m]]
    pw    <- data$xtx_diag_list[[m]]
    n     <- length(idx_m)

    Eb_m  <- alpha_l * model$mu[[l]][[m]]
    Eb2_m <- alpha_l * model$mu2[[l]][[m]]
    R_m   <- model$raw_residuals[[m]][idx_m, , drop = FALSE]
    sigma2_per_pos <- mf_sigma2_per_position(data, model, m)

    XEb_m <- (X %*% Eb_m)[idx_m, , drop = FALSE]
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
    idx_m <- data$na_idx[[m]]
    R_m   <- model$raw_residuals[[m]][idx_m, , drop = FALSE]
    n     <- length(idx_m)
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
#' Defensive override of susieR's default. The default
#' `neg_loglik.individual` calls `loglik.individual(...)` by
#' name (not via S3 dispatch), so if any future code path
#' invokes `neg_loglik(data, params, ...)` on `mf_individual`
#' data, the wrong (scalar-Gaussian) likelihood would fire.
#' This override forces `loglik.mf_individual`. Currently the
#' `optimize_prior_variance.mf_individual` path skips
#' `neg_loglik` entirely; this is a regression guard.
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
  pw    <- data$xtx_diag_list[[m]]
  idx_m <- data$na_idx[[m]]
  T_pad <- data$T_basis[m]

  res_m <- (data$D[[m]] - model$fitted[[m]])[idx_m, , drop = FALSE]
  rss_t <- colSums(res_m * res_m)

  # Bias correction: sum over observed rows only so the quadratic form
  # uses the same n_m observations as the RSS term above.
  bias_t <- numeric(T_pad)
  for (l in seq_len(model$L)) {
    alpha_l  <- model$alpha[l, ]
    mu_lm    <- model$mu[[l]][[m]]
    mu2_l_m  <- model$mu2[[l]][[m]]
    XEb_lm   <- (X %*% (alpha_l * mu_lm))[idx_m, , drop = FALSE]
    bias_t   <- bias_t + colSums((alpha_l * pw) * mu2_l_m) -
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
    n     <- length(data$na_idx[[m]])
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
  refresh_iter_cache.mf_individual(data, model)
}

# Recompute per-iter precomputes that depend on `sigma2` (and, for
# the mixsqp-backed prior classes, on the slab sd grid). Slots are
# invariant across the L-effect loop within one IBSS iteration;
# rebuilding them per effect would be L-fold redundant.
#
# Cache contents per outcome m:
#   $shat2[[m]]            -- p x T_basis[m] = outer(1/xtx_diag,
#                             sigma2_per_pos). Read by
#                             mf_per_outcome_bhat_shat() on every
#                             SER call (loglik + posterior moments)
#                             for every prior class. Always built.
#   $sdmat[[m]][[s]]       -- (p*|idx_s|) x K matrix used by
#                             mf_em_likelihood_per_scale, the
#                             mixsqp M-step kernel only. Skipped
#                             for the ebnm-backed prior classes
#                             (point_normal, point_laplace), which
#                             fit 2 scalars per (m, s) by direct
#                             optim and have no per-(j, t, k)
#                             aggregate to amortize.
#   $log_sdmat[[m]][[s]]   -- log of $sdmat[[m]][[s]] (mixsqp-only).
refresh_iter_cache.mf_individual <- function(data, model) {
  M  <- data$M
  is_mixsqp_prior <- inherits(
    model$G_prior[[1L]],
    c("mixture_normal", "mixture_normal_per_scale"))

  shat2_list     <- vector("list", M)
  sdmat_list     <- if (is_mixsqp_prior) vector("list", M) else NULL
  log_sdmat_list <- if (is_mixsqp_prior) vector("list", M) else NULL

  for (m in seq_len(M)) {
    pw <- data$xtx_diag_list[[m]]
    sigma2_per_pos <- mf_sigma2_per_position(data, model, m)
    shat2_m <- outer(1 / pw, sigma2_per_pos)
    shat2_list[[m]] <- shat2_m

    if (!is_mixsqp_prior) next

    G_m       <- model$G_prior[[m]]
    sdmat_m     <- vector("list", length(G_m))
    log_sdmat_m <- vector("list", length(G_m))
    for (s in seq_along(G_m)) {
      idx     <- G_m[[s]]$idx
      sd_grid <- G_m[[s]]$fitted_g$sd
      svec_s  <- as.numeric(sqrt(shat2_m[, idx, drop = FALSE]))
      sdmat   <- sqrt(outer(svec_s^2, sd_grid^2, "+"))
      sdmat_m[[s]]     <- sdmat
      log_sdmat_m[[s]] <- log(sdmat)
    }
    sdmat_list[[m]]     <- sdmat_m
    log_sdmat_list[[m]] <- log_sdmat_m
  }

  model$iter_cache <- list(
    shat2     = shat2_list,
    sdmat     = sdmat_list,
    log_sdmat = log_sdmat_list
  )
  model
}

#' Per-effect prior update for `mf_individual`
#'
#' Overrides `optimize_prior_variance` for the `mf_individual`
#' data class. Computes the alpha-thinned `(keep_idx, zeta_keep)`
#' shared by every solver and dispatches per (outcome, scale) on
#' the G_prior class:
#'
#' - `mixture_normal` / `mixture_normal_per_scale`: mixsqp via
#'   `.opv_mixsqp()`. Builds `mf_em_likelihood_per_scale()` then
#'   `mf_em_m_step_per_scale()` on the alpha-weighted
#'   `(|keep_idx| x |idx_s|)` rectangle of (Bhat, Shat); updates
#'   `fitted_g$pi`.
#' - `mixture_point_normal_per_scale`: ebnm via
#'   `.opv_ebnm_point_normal()`. Calls `ebnm::ebnm_point_normal()`
#'   on the unweighted `(|keep_idx| x |idx_s|)` rectangle of
#'   (Bhat, Shat); updates `fitted_g`.
#' - `mixture_point_laplace_per_scale`: ebnm via
#'   `.opv_ebnm_point_laplace()`. Calls `ebnm::ebnm_point_laplace()`
#'   on the same rectangle; updates `fitted_g`.
#'
#' Each helper mutates `model$G_prior[[m]][[s]]$fitted_g` and
#' `model$pi_V[[m]][s, ]`. The susieR-style scalar `V[l]` is held
#' at 1 since the mixture weights absorb the per-effect prior
#' shape adaptation. Returns `list(V = 1, model = updated_model)`
#' per the susieR generic contract.
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
  zeta_l <- model$alpha[l, ]

  # Adaptive variable subsetting: drop variables where the SuSiE
  # per-effect posterior alpha is below `alpha_thin_eps`. For the
  # mixsqp-backed paths this caps the L-matrix size; for the
  # ebnm-backed paths the kept variables form the multi-variable
  # observation set ebnm fits over (mirroring mixsqp's data
  # shape minus the per-row alpha weighting).
  alpha_eps <- params$alpha_thin_eps %||% 1e-6
  keep_idx  <- if (alpha_eps > 0) which(zeta_l > alpha_eps)
               else seq_along(zeta_l)
  if (length(keep_idx) < 1L) keep_idx <- which.max(zeta_l)
  zeta_keep <- zeta_l[keep_idx]

  G1      <- model$G_prior[[1L]]
  ebnm_fn <- switch(class(G1)[1L],
    mixture_point_normal_per_scale  = ebnm::ebnm_point_normal,
    mixture_point_laplace_per_scale = ebnm::ebnm_point_laplace,
    NULL)
  if (!is.null(ebnm_fn)) {
    model <- .opv_ebnm_point(data, params, model, ser_stats, keep_idx, ebnm_fn, l = l)
  } else if (inherits(G1, c("mixture_normal", "mixture_normal_per_scale"))) {
    model <- .opv_mixsqp(data, params, model, ser_stats, keep_idx, zeta_keep, l = l)
  } else {
    stop("Unknown prior class on G_prior[[1]]: ",
         paste(class(G1), collapse = ", "))
  }
  list(V = 1, model = model)
}

# =============================================================
# Per-effect EM step layout (mfsusieR vs standard susieR IBSS)
# =============================================================
#
# susieR scalar-V IBSS (per outer iter, per effect l):
#
#   1. partial residual r_l
#   2. compute_ser_statistics  (Bhat, Shat from r_l)
#   3. (if optim/uniroot/simple) optim V                 [pre hook]
#   4. loglik         -> lbf, alpha[l, ]
#   5. moments        -> mu[l], mu2[l]
#   6. KL             -> KL[l]
#   7. (if EM)        optim V using updated alpha+moments [post hook]
#
# Steps 4-7 run ONCE per effect per outer iter. For a scalar V
# prior this is enough: the V optim is fast and the (V, alpha)
# joint update converges in a couple of outer iters. susieR has
# no inner-EM concept.
#
# mfsusieR mixture-prior IBSS (this file, per outer iter, per
# effect l):
#
#   1. partial residual r_l
#   2. compute_ser_statistics  (per-outcome Bhat, Shat from r_l)
#   3. pre_loglik_prior_hook.mf_individual                [pre hook]
#        - swap fitted_g_per_effect[[l]] -> G_prior scratchpad
#          (so each effect sees its own persistent prior, not
#          another effect's last-written one)
#        - NO M-step here: at iter 1 alpha[l, ] is uniform 1/p
#          and mixsqp on uniform-weighted data collapses pi to
#          all-null, cascading lbf=0 through every effect.
#   4. loglik         -> lbf, alpha[l, ]
#   5. moments, KL
#   6. post_loglik_prior_hook.mf_individual               [post hook]
#        Inner EM loop, `params$inner_em_steps` M-step calls
#        (default 2; set 1 to disable). This control is local to
#        mfsusieR; susieR exposes no equivalent.
#        Per cycle:
#          a. M-step (mixsqp / ebnm) using current alpha[l, ]
#             -> writes new pi to G_prior + fitted_g_per_effect[[l]]
#          b. (skipped on the last cycle) re-run loglik / moments /
#             KL so alpha[l, ] reflects the new pi within this iter
#        With 2 cycles the per-effect (alpha, pi) end the iter
#        in lockstep.
#
# Big-picture intuition: scalar-V susieR is fine with one M-step
# per outer iter because V is one number and (V, alpha) settles
# quickly. mfsusieR's mixture pi is a K-dim simplex tightly
# coupled to alpha; the inner loop keeps the two in sync per
# effect per outer iter, so spurious effects collapse cleanly
# instead of slowly drifting on residual-noise correlations.
#
# Per-prior empirical note: per_scale_normal (ebnm) already
# converges in few outer iters even with a single inner step,
# because each scale's prior is fit independently and the
# pi-alpha drift is naturally decoupled. per_outcome (mixsqp on
# a single pooled prior per outcome) is the path that benefits
# most from the inner loop. The same code path applies to both;
# no special-casing needed.
#
# Cost: ~2x M-step + 1 extra loglik/moments/KL per effect per
# outer iter, recovered by ~2x fewer outer iters at the default
# inner_em_steps = 2. Higher values trade more inner work for
# fewer outer iters (diminishing returns past ~5).
# =============================================================

#' @keywords internal
#' @export
pre_loglik_prior_hook.mf_individual <- function(data, params, model, ser_stats,
                                                l, V_init) {
  if (!is.null(model$fitted_g_per_effect)) {
    fge_l <- model$fitted_g_per_effect[[l]]
    for (m in seq_len(data$M)) {
      G_m <- model$G_prior[[m]]
      for (s in seq_along(G_m)) {
        model$G_prior[[m]][[s]]$fitted_g <- fge_l[[m]][[s]]
      }
    }
  }
  list(V = 1, model = model)
}

#' @keywords internal
#' @export
post_loglik_prior_hook.mf_individual <- function(data, params, model, ser_stats,
                                                 l, V_init) {
  if (!isTRUE(params$estimate_prior_variance)) {
    return(list(V = 1, model = model))
  }
  # max_inner_em_steps caps EXTRA inner refinement cycles beyond
  # the always-on first M-step (0 = skip). Loop exits early when
  # the per-effect lbf change drops below `params$tol`, mirroring
  # the inner-loop convergence used by fsusieR's EM_pi.
  inner_cap <- max(0L, as.integer(params$max_inner_em_steps %||% 5L)) + 1L
  inner_tol <- params$tol %||% 1e-4
  loglik_fn  <- getFromNamespace("loglik",                       "susieR")
  cpm_fn     <- getFromNamespace("calculate_posterior_moments",  "susieR")
  ckl_fn     <- getFromNamespace("compute_kl",                   "susieR")
  get_alpha  <- getFromNamespace("get_alpha_l",                  "susieR")
  get_post   <- getFromNamespace("get_posterior_moments_l",      "susieR")

  for (k in seq_len(inner_cap)) {
    prev_lbf <- model$lbf[l]
    out <- optimize_prior_variance.mf_individual(
      data, params, model, ser_stats,
      l = l, alpha = get_alpha(model, l),
      moments = get_post(model, l), V_init = V_init)
    model <- out$model
    if (k == inner_cap) break
    # Re-run loglik / moments / KL against the freshly-M-stepped
    # G_prior so alpha[l, ] reflects the new pi before the next
    # M-step cycle.
    model <- loglik_fn(data, params, model, V = 1, ser_stats, l)
    model <- cpm_fn(data, params, model, V = 1, l)
    model <- ckl_fn(data, params, model, l)
    if (!is.na(prev_lbf) && abs(model$lbf[l] - prev_lbf) < inner_tol) break
  }
  list(V = 1, model = model)
}

#' mixsqp M-step on `pi_V` per (outcome, scale)
#'
#' Per (outcome, scale), builds the mixsqp likelihood matrix
#' (`mf_em_likelihood_per_scale`) from the per-outcome (Bhat, Shat)
#' in `ser_stats`, runs the M step (`mf_em_m_step_per_scale`)
#' with `zeta_keep = model$alpha[l, keep_idx]`, and writes the
#' new mixture weights into both
#' `model$G_prior[[m]][[s]]$fitted_g$pi` and `model$pi_V[[m]][s, ]`.
#'
#' @keywords internal
#' @noRd
.opv_mixsqp <- function(data, params, model, ser_stats,
                        keep_idx, zeta_keep, l = 1L) {
  # The joint per-effect variable posterior `model$alpha[l, ]` is
  # the softmax of the joint log-Bayes-factor across all M outcomes
  # and S_m scales. Adding M outcomes increases the variance of
  # the per-variable joint lbf by ~M, so the dominant alpha values
  # concentrate ~M-fold relative to the M = 1 case. The mixsqp
  # M-step's regularization-to-data balance is set by
  # `nullweight / max_alpha`; to hold this ratio fixed across M
  # we scale `mixture_null_weight` by M. See
  # `inst/notes/cross-package-audit-derivations.md` section 1
  # for the full derivation.
  mixture_null_weight <- (params$mixture_null_weight %||% 0.05) *
                        max(1L, data$M)
  control <- params$control_mixsqp %||% list()
  cache   <- model$iter_cache
  p_full  <- ncol(model$alpha)   # `model$alpha` is `L x p`

  # The cached sdmat / log_sdmat are laid out in column-major order
  # over `(j = 1..p_full, t in idx)`. Restricting to the rows that
  # correspond to `keep_idx` is the same index expression for both
  # caches, so we factor it once per `(m, s)`.
  cache_row_keep <- function(idx) {
    as.vector(outer(keep_idx, (seq_along(idx) - 1L) * p_full, "+"))
  }
  cache_subset <- function(cache_ms, idx) {
    if (is.null(cache_ms)) NULL
    else cache_ms[cache_row_keep(idx), , drop = FALSE]
  }
  for (m in seq_len(data$M)) {
    bhat_m <- ser_stats$betahat[[m]]
    shat_m <- sqrt(ser_stats$shat2[[m]])
    G_m    <- model$G_prior[[m]]
    cache_m_sdmat     <- if (is.null(cache)) NULL else cache$sdmat[[m]]
    cache_m_log_sdmat <- if (is.null(cache)) NULL else cache$log_sdmat[[m]]
    # Uniform group loop: each G_prior entry holds the column
    # indices it covers (`$idx`). For `per_outcome` there is one
    # group spanning every wavelet column (one mixsqp solve over
    # the full outcome); for `per_scale` there are S_m
    # groups, one per wavelet scale (S_m independent mixsqp solves).
    for (s in seq_along(G_m)) {
      idx     <- G_m[[s]]$idx
      sd_grid <- G_m[[s]]$fitted_g$sd
      bhat_sub <- bhat_m[keep_idx, idx, drop = FALSE]
      shat_sub <- shat_m[keep_idx, idx, drop = FALSE]
      sdmat_sub     <- cache_subset(cache_m_sdmat[[s]],     idx)
      log_sdmat_sub <- cache_subset(cache_m_log_sdmat[[s]], idx)
      L_mat <- mf_em_likelihood_per_scale(
        bhat_slice      = bhat_sub,
        shat_slice      = shat_sub,
        sd_grid         = sd_grid,
        sdmat_cache     = sdmat_sub,
        log_sdmat_cache = log_sdmat_sub)
      # Warm-start mixsqp from the previous (m, s) pi -- it changes
      # slowly between IBSS iters, cutting inner SQP iterations
      # from ~20 (cold) to ~1-3.
      pi_prev <- model$G_prior[[m]][[s]]$fitted_g$pi
      new_pi  <- mf_em_m_step_per_scale(
        L_mat, zeta_keep, idx_size = length(idx),
        mixture_null_weight = mixture_null_weight,
        control_mixsqp = control,
        pi_warm_start  = pi_prev)
      model$G_prior[[m]][[s]]$fitted_g$pi <- new_pi
      model$pi_V[[l]][[m]][s, ]            <- new_pi
      model$fitted_g_per_effect[[l]][[m]][[s]] <-
        model$G_prior[[m]][[s]]$fitted_g
    }
  }
  model
}

#' ebnm M-step on the per-(outcome, scale) point-* prior
#'
#' Per (outcome, scale), calls `ebnm_fn` on the alpha-thinned
#' `(|keep_idx| × |idx_s|)` (Bhat, Shat) rectangle. Writes
#' `fit$fitted_g` into `model$G_prior[[m]][[s]]$fitted_g` and
#' `fit$fitted_g$pi` into `model$pi_V[[m]][s, ]`. Mirrors
#' `.opv_mixsqp`'s data shape and per-(m, s) loop; the per-row
#' alpha weighting that mixsqp uses is dropped because ebnm
#' takes no observation weights.
#'
#' @keywords internal
#' @noRd
.opv_ebnm_point <- function(data, params, model, ser_stats,
                            keep_idx, ebnm_fn, l = 1L) {
  fix_g <- !isTRUE(params$estimate_prior_variance)
  for (m in seq_len(data$M)) {
    bhat_m <- ser_stats$betahat[[m]][keep_idx, , drop = FALSE]
    shat_m <- sqrt(ser_stats$shat2[[m]][keep_idx, , drop = FALSE])
    G_m    <- model$G_prior[[m]]
    for (s in seq_along(G_m)) {
      idx <- G_m[[s]]$idx
      # `mode = 0` locks the slab at zero (the SuSiE coefficient
      # is centered there). Skip when `fix_g = TRUE`; ebnm uses
      # `g_init$mean` and warns if both are passed.
      args <- list(x      = as.vector(bhat_m[, idx, drop = FALSE]),
                   s      = as.vector(shat_m[, idx, drop = FALSE]),
                   g_init = G_m[[s]]$fitted_g,
                   fix_g  = fix_g)
      if (!fix_g) args$mode <- 0
      fit <- do.call(ebnm_fn, args)
      model$G_prior[[m]][[s]]$fitted_g <- fit$fitted_g
      model$pi_V[[l]][[m]][s, ]         <- fit$fitted_g$pi
      model$fitted_g_per_effect[[l]][[m]][[s]] <- fit$fitted_g
    }
  }
  model
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

#' Refresh per-effect `lbf[l]` / `KL[l]` against the iter-final pi_V
#'
#' During an IBSS iteration, each effect's `model$lbf[l]` and
#' `model$KL[l]` are computed at the moment effect l was updated --
#' but mfsusieR's per-effect prior (`pi_V[m][s, ]`) is updated
#' between effects via mixsqp inside the SER step (see
#' `optimize_prior_variance.mf_individual`). When effect l+1 updates
#' pi_V, the previously-recorded `lbf[l]` / `KL[l]` for l < l+1 are
#' against a now-stale pi_V state. The reported per-iteration ELBO
#' `Eloglik(state_t) - sum(model$KL)` is therefore a hybrid quantity:
#' Eloglik against the iter-final state, `KL[l]` against
#' state-when-l-was-updated. The hybrid is empirically close to the
#' coherent variational free energy but is not strictly monotone.
#'
#' This refresh re-evaluates `lbf[l]` / `lbf_variable[l]` / `KL[l]`
#' for each effect l using the iter-final pi_V (and the iter-final
#' partial residual rebuilt from the frozen alpha / mu / mu2). The
#' variational posterior `alpha[l]` is preserved via the
#' `update_alpha = FALSE` flag on `loglik.mf_individual`. After the
#' refresh, `Eloglik(data, model) - sum(model$KL)` is the coherent
#' ELBO at iter end.
#'
#' Cost: ~5-10% per iter (an extra residual rebuild + SER stat
#' recomputation per effect). Worth it to make the ELBO trajectory
#' a defensible variational free energy rather than a hybrid.
#'
#' Reference: mvsusieR's `compute_multivariate_elbo_ss`
#' (`mvsusieR/R/sufficient_stats_methods.R`) follows the same
#' post-hoc-refresh pattern for the same reason.
#'
#' @keywords internal
#' @noRd
refresh_lbf_kl.mf_individual <- function(data, params, model) {
  L <- nrow(model$alpha)
  for (l in seq_len(L)) {
    # Rebuild R_l = Y - sum_{l' != l} X * alpha_{l'} * mu_{l'} from the
    # iter-final variational posterior.
    model     <- compute_residuals.mf_individual(data, params, model, l)
    ser_stats <- compute_ser_statistics.mf_individual(data, params, model, l)

    # Refresh lbf[l] / lbf_variable[l] (and lbf_variable_outcome[l, , m]
    # if attached) against the iter-final pi_V via the same loglik
    # pipeline IBSS uses, but freeze alpha[l] -- we are evaluating the
    # ELBO at the existing q_t, not advancing it.
    model <- loglik.mf_individual(data, params, model, V = 1L,
                                   ser_stats = ser_stats, l = l,
                                   update_alpha = FALSE)

    # Refresh KL[l] using the new lbf[l] and the new R_l.
    model <- compute_kl.mf_individual(data, params, model, l)
  }
  model
}

#' Coherent variational free energy for the `mfsusie` model class
#'
#' Overrides susieR's `get_objective.default`. On the ELBO
#' convergence path, refreshes per-effect `lbf[l]` / `KL[l]` against
#' the iter-final pi_V (see `refresh_lbf_kl.mf_individual`) so the
#' returned ELBO is a coherent variational free energy. On the PIP
#' convergence path the ELBO is unused by `check_convergence`, so
#' the refresh is wasted work and we short-circuit with `NA_real_`.
#'
#' @keywords internal
#' @noRd
get_objective.mfsusie <- function(data, params, model) {
  if (identical(params$convergence_method, "pip")) return(NA_real_)
  model <- refresh_lbf_kl.mf_individual(data, params, model)
  Eloglik(data, model) - sum(model$KL, na.rm = TRUE)
}

#' Class-specific extra fields stripped by cleanup_model.default
#'
#' susieR's `cleanup_model.default` strips a standard set of
#' temporary fields (`runtime`, `residuals`, `fitted_without_l`,
#' `residual_variance`, ...) and unions in whatever this generic
#' returns. mfsusieR adds `model$V` to the strip list because the
#' length-L scalar prior coefficient is held at 1 throughout IBSS:
#' the per-effect adaptation lives in the per-(scale, outcome)
#' mixture weights `pi_V[[m]]`, not in V. Keeping V on the fit
#' would mislead users into thinking it carries information.
#'
#' `model$fitted` (the wavelet-domain running fit, list of length
#' M) is intentionally NOT stripped: it is consumed by
#' `ibss_initialize.mf_individual` via `params$model_init$fitted`
#' to warm-start a follow-up call, saving one IBSS sweep.
#'
#' `raw_residuals` is stripped one rung up the dispatch chain by
#' `cleanup_model.individual`, so it does not appear here.
#'
#' @keywords internal
#' @noRd
cleanup_extra_fields.mf_individual <- function(data) {
  c("V")
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
    n_m <- length(data$na_idx[[m]])
    out <- out + sum(-0.5 * n_m * log(2 * pi * s2_t) -
                     0.5 * er2_t / s2_t)
  }
  out
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
