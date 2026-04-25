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
  pw <- data$predictor_weights

  # sigma2 may be:
  #   - a scalar per modality (legacy "shared_per_modality" mode); or
  #   - a length-S_m vector per modality (per-(scale, modality) default).
  # For per-position Shat in the IBSS, broadcast per scale_index.
  betahat <- vector("list", M)
  shat2   <- vector("list", M)
  for (m in seq_len(M)) {
    XtR_m <- model$residuals[[m]]
    betahat[[m]] <- XtR_m / pw

    sigma2_m <- model$sigma2[[m]]
    if (length(sigma2_m) == 1L) {
      shat2[[m]] <- matrix(sigma2_m / pw,
                           nrow = nrow(XtR_m), ncol = ncol(XtR_m))
    } else {
      # length-S_m vector: per scale. Broadcast via scale_index map.
      sigma2_per_pos <- numeric(ncol(XtR_m))
      for (s in seq_along(model$sigma2[[m]])) {
        sigma2_per_pos[data$scale_index[[m]][[s]]] <- sigma2_m[s]
      }
      shat2[[m]] <- outer(1 / pw, sigma2_per_pos)
    }
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
