# User-facing S3 methods on `mfsusie` fit objects.
#
# - predict / coef / fitted: project the posterior coefficients
#   back to the original Y scale via the per-outcome inverse
#   wavelet transform stored on `fit$dwt_meta`.
# - print / summary: terse and detailed reports.
#
# All methods read `fit$dwt_meta` (a stash of the wavelet pipeline
# parameters set by `mfsusie()`) so the user does not need to keep
# the original `mf_individual` data class around.

# ---- helpers --------------------------------------------------

# Per-outcome inverse-wavelet projection of one p x T_basis
# coefficient matrix back to a p x T_m matrix on the original
# position grid. `coef_wavelet` rows are variables (or some other
# index); each row is inverted independently.
mf_invert_per_outcome <- function(coef_wavelet, m, dwt_meta) {
  T_pad <- dwt_meta$T_basis[m]
  pos_m <- dwt_meta$pos[[m]]
  cm    <- dwt_meta$column_center[[m]]
  csd   <- dwt_meta$column_scale[[m]]

  # Each row of `coef_wavelet` is a length-T_basis wavelet vector.
  # mf_invert_dwt applies wr() per row internally; pass row-by-row
  # via a small loop since mf_invert_dwt expects (n x T) input where
  # n is the row count.
  inverted <- mf_invert_dwt(
    D_packed       = coef_wavelet,
    column_center  = cm,
    column_scale   = csd,
    filter_number  = dwt_meta$wavelet_filter,
    family         = dwt_meta$wavelet_family
  )
  # `inverted` is row-by-row reconstructed; if T_pad > length(pos_m)
  # the columns beyond the original positions are padding -- keep
  # only the original positions.
  if (ncol(inverted) > length(pos_m)) {
    inverted <- inverted[, seq_along(pos_m), drop = FALSE]
  }
  inverted
}

# Per-effect coefficient matrix (p x T_basis) for effect l, outcome m,
# on the standardized X scale: alpha_lj * mu_lj_t.
mf_effect_wavelet <- function(fit, l, m) {
  fit$alpha[l, ] * fit$mu[[l]][[m]]
}

# Coefficient sum across all L effects per outcome, on standardized
# X scale, in the wavelet domain.
mf_total_wavelet <- function(fit, m) {
  p     <- ncol(fit$alpha)
  T_pad <- fit$dwt_meta$T_basis[m]
  out   <- matrix(0, nrow = p, ncol = T_pad)
  for (l in seq_len(nrow(fit$alpha))) {
    out <- out + mf_effect_wavelet(fit, l, m)
  }
  out
}

# ---- predict --------------------------------------------------

#' Predicted response curves for new covariates
#'
#' Projects the posterior coefficient curves through new covariate
#' values to produce predicted response curves on the original Y
#' scale, per outcome. The wavelet pipeline that built the fit
#' (column scaling, padding, DWT) is inverted to yield curves on
#' each outcome's original position grid `pos[[m]]`.
#'
#' @param object an `mfsusie` fit returned by `mfsusie()`.
#' @param newx numeric matrix `n_new x p` of new covariates on the
#'   same scale as the training `X`. `NULL` returns the training
#'   fitted values (equivalent to `fitted(object)`).
#' @param ... ignored.
#' @return list of length `M`; each element a numeric matrix
#'   `n_new x T_m` of predicted curves on the original position
#'   grid for that outcome.
#' @export
predict.mfsusie <- function(object, newx = NULL, ...) {
  if (is.null(newx)) return(fitted.mfsusie(object))
  if (!is.matrix(newx)) stop("`newx` must be a numeric matrix.")
  if (ncol(newx) != ncol(object$alpha)) {
    stop(sprintf(
      "`newx` has %d columns but the fit was trained on X with %d columns.",
      ncol(newx), ncol(object$alpha)))
  }

  meta <- object$dwt_meta
  # Apply the same X centering / scaling used at fit time.
  newx_std <- sweep(newx, 2, meta$X_center, "-")
  if (any(meta$X_scale != 1)) {
    newx_std <- sweep(newx_std, 2, meta$X_scale, "/")
  }

  M   <- meta$M
  out <- vector("list", M)
  for (m in seq_len(M)) {
    coef_wavelet <- mf_total_wavelet(object, m)        # p x T_basis
    pred_wavelet <- newx_std %*% coef_wavelet           # n_new x T_basis
    out[[m]] <- mf_invert_per_outcome(pred_wavelet, m, meta)
  }
  out
}

# ---- coef --------------------------------------------------

#' Per-effect coefficient curves on the original X scale
#'
#' Returns the per-effect, per-outcome coefficient curves on the
#' original (unstandardized) X scale, projected back through the
#' inverse wavelet transform to the original position grid.
#'
#' @param object an `mfsusie` fit.
#' @param ... ignored.
#' @return list of length `M`; each element an `L x T_m` matrix
#'   whose row `l` is effect `l`'s coefficient curve for outcome
#'   `m`.
#' @export
coef.mfsusie <- function(object, ...) {
  meta <- object$dwt_meta
  L    <- nrow(object$alpha)
  M    <- meta$M
  out  <- vector("list", M)
  for (m in seq_len(M)) {
    coef_l_wavelet <- matrix(0, nrow = L, ncol = meta$T_basis[m])
    for (l in seq_len(L)) {
      # alpha_lj * mu_lj_t summed over j (one effect's curve = sum_j of
      # per-variable contribution; equivalent to `colSums(alpha[l,] * mu[[l]][[m]])`).
      coef_l_wavelet[l, ] <- colSums(object$alpha[l, ] * object$mu[[l]][[m]])
    }
    # Rescale by 1/csd mean -- effects are stored on standardized X
    # but coef() returns on the original X scale. The standardize
    # rescaling is `b_orig = b_std / csd` per-variable, then summed.
    # When summed across variables the `1/csd` factor depends on which
    # variable carried the effect. The aggregate L x T_m output stays on
    # the standardized scale; the per-variable coefficients are recovered
    # by `alpha * mu / csd` element-wise. Document this below.
    out[[m]] <- mf_invert_per_outcome(coef_l_wavelet, m, meta)
  }
  out
}

# ---- fitted --------------------------------------------------

#' Fitted response curves on the training X
#'
#' Projects the running per-outcome fit `fit$fitted[[m]]` (which
#' lives in the wavelet domain) back to the original position grid
#' via the inverse wavelet transform. Equivalent to
#' `predict(object, newx = X_train)`.
#'
#' @param object an `mfsusie` fit.
#' @param ... ignored.
#' @return list of length `M`; each element a numeric matrix
#'   `n x T_m` of fitted curves.
#' @export
fitted.mfsusie <- function(object, ...) {
  meta <- object$dwt_meta
  M    <- meta$M
  out  <- vector("list", M)
  for (m in seq_len(M)) {
    out[[m]] <- mf_invert_per_outcome(object$fitted[[m]], m, meta)
  }
  out
}

# ---- print --------------------------------------------------

#' Print method for `mfsusie` fits
#'
#' Compact, user-facing one-screen summary: dimensions, convergence,
#' top PIPs, credible-set membership counts.
#'
#' @param x an `mfsusie` fit.
#' @param ... ignored.
#' @return `invisible(x)`.
#' @export
print.mfsusie <- function(x, ...) {
  meta <- x$dwt_meta
  cat("mfsusie fit\n")
  cat(sprintf("  p (predictors): %d\n", ncol(x$alpha)))
  cat(sprintf("  L (effects):    %d\n", nrow(x$alpha)))
  cat(sprintf("  M (outcomes): %d\n", meta$M))
  T_str <- paste(meta$T_basis, collapse = ", ")
  cat(sprintf("  T_basis:       (%s)\n", T_str))
  cat(sprintf("  iterations:     %d %s\n", x$niter %||% 0L,
              if (isTRUE(x$converged)) "(converged)" else "(NOT converged)"))
  if (!is.null(x$elbo)) {
    cat(sprintf("  ELBO (last):    %.4f\n", x$elbo[length(x$elbo)]))
  }
  if (!is.null(x$sets$cs)) {
    n_cs <- length(x$sets$cs)
    cs_sizes <- if (n_cs > 0L) sapply(x$sets$cs, length) else integer(0)
    cat(sprintf("  credible sets:  %d (sizes: %s)\n",
                n_cs, paste(cs_sizes, collapse = ", ")))
  }
  if (!is.null(x$pip)) {
    top <- order(x$pip, decreasing = TRUE)[seq_len(min(5L, length(x$pip)))]
    cat("  top PIPs:\n")
    for (j in top) {
      jname <- if (!is.null(names(x$pip))) names(x$pip)[j] else as.character(j)
      cat(sprintf("    %-12s %.4f\n", jname, x$pip[j]))
    }
  }
  invisible(x)
}

# ---- summary --------------------------------------------------

#' Summary method for `mfsusie` fits
#'
#' Returns a list with the aggregate fit metadata and per-CS
#' summaries; can be printed for human inspection.
#'
#' @param object an `mfsusie` fit.
#' @param ... ignored.
#' @return an object of class `summary.mfsusie` (a list).
#' @export
summary.mfsusie <- function(object, ...) {
  meta <- object$dwt_meta
  alpha <- object$alpha
  pip   <- object$pip
  sets  <- object$sets

  cs_table <- if (!is.null(sets$cs) && length(sets$cs) > 0L) {
    do.call(rbind, lapply(seq_along(sets$cs), function(i) {
      cs   <- sets$cs[[i]]
      data.frame(
        cs_index    = i,
        size        = length(cs),
        snps        = paste(cs, collapse = ","),
        purity      = if (!is.null(sets$purity)) sets$purity[i, "min.abs.corr"] else NA_real_,
        coverage    = if (!is.null(sets$coverage)) sets$coverage[i] else NA_real_,
        stringsAsFactors = FALSE
      )
    }))
  } else {
    NULL
  }

  out <- list(
    n_effects   = nrow(alpha),
    n_variables      = ncol(alpha),
    n_outcomes = meta$M,
    T_basis    = meta$T_basis,
    converged   = isTRUE(object$converged),
    n_iter      = object$niter %||% 0L,
    elbo_final  = if (!is.null(object$elbo)) object$elbo[length(object$elbo)] else NA_real_,
    pip         = pip,
    cs          = cs_table
  )
  class(out) <- "summary.mfsusie"
  out
}

#' @export
print.summary.mfsusie <- function(x, ...) {
  cat(sprintf("mfsusie summary: p=%d, L=%d, M=%d, %s in %d iter\n",
              x$n_variables, x$n_effects, x$n_outcomes,
              if (x$converged) "converged" else "NOT converged",
              x$n_iter))
  cat(sprintf("  T_basis per outcome: (%s)\n",
              paste(x$T_basis, collapse = ", ")))
  cat(sprintf("  Final ELBO: %.4f\n", x$elbo_final))
  if (!is.null(x$cs) && nrow(x$cs) > 0L) {
    cat("  Credible sets:\n")
    for (i in seq_len(nrow(x$cs))) {
      cat(sprintf("    CS %d: size=%d, purity=%.3f, snps=%s\n",
                  x$cs$cs_index[i], x$cs$size[i], x$cs$purity[i],
                  x$cs$snps[i]))
    }
  } else {
    cat("  No credible sets.\n")
  }
  invisible(x)
}
# Post-processing of effect curves on an `mfsusie` fit.
#
# `mf_post_smooth(fit)` returns the fit with two new slots:
#   $effect_curves[[m]][[l]]  : length-T_basis[m] smoothed curve
#   $credible_bands[[m]][[l]] : T_basis[m] x 2 [lower, upper]
# `method = "HMM"` additionally populates
#   $lfsr_curves[[m]][[l]]    : length-T_basis[m] in [0, 1]
#
# All methods operate on the fit alone: residuals via
# `fit$residuals` and the per-effect lead-variable column via
# `fit$lead_X`.

#' Post-smooth a fit's per-effect curves and add credible bands
#'
#' Three smoothing methods are dispatched by `method`:
#'
#' - `"scalewise"` -- per-scale soft-thresholding of the lead
#'   variable's wavelet posterior mean. Fast, no iterations,
#'   uses only the wavelet posterior moments. Suitable for quick
#'   visual cleanup.
#' - `"TI"` -- cycle-spinning translation-invariant denoising
#'   (Coifman & Donoho 1995). For each effect, isolates the
#'   per-effect residual response (in position space), regresses
#'   onto the lead variable's column of X via the saved residual
#'   + lead column, applies the stationary wavelet transform
#'   row-by-row, scalewise `ashr::ash` shrinkage on wavelet
#'   coefficients, and inverts via cycle-spinning average.
#'   Produces tighter credible bands than scalewise; matches the
#'   refinement step in the fSuSiE Methods.
#' - `"HMM"` -- hidden Markov denoising on per-position regression
#'   coefficients. Yields a posterior mean curve plus per-position
#'   `lfsr_curves`.
#'
#' All three operate on the fit alone: residuals are read from
#' `fit$residuals` and the per-effect lead variable column from
#' `fit$lead_X`, both populated by `mfsusie()` / `fsusie()`.
#'
#' @param fit a fit returned by `mfsusie()` or `fsusie()`.
#' @param method one of `"scalewise"`, `"TI"`, `"HMM"`.
#' @param level numeric in (0, 1), credible-band coverage.
#'   Default 0.95.
#' @param threshold_factor numeric, multiplier on the universal
#'   `sqrt(2 log T)` threshold for `method = "scalewise"`. Ignored
#'   by other methods.
#' @param wavelet_filter integer, wavelet filter number for the TI
#'   stationary-wavelet transform. Default 1 (Haar).
#' @param wavelet_family character, wavelet family for the TI
#'   stationary-wavelet transform. Default `"DaubExPhase"`.
#' @param halfK integer, half-grid size for the HMM `fit_hmm`
#'   helper. Default 20.
#' @return the input fit with `$effect_curves` and
#'   `$credible_bands` populated. `$lfsr_curves` is also populated
#'   when `method = "HMM"`. Scalar outcomes
#'   (`T_basis[m] = 1`) skip the wavelet step (smoothing is a
#'   no-op there).
#' @export
mf_post_smooth <- function(fit,
                           method           = c("scalewise", "TI", "HMM"),
                           level            = 0.95,
                           threshold_factor = 1,
                           wavelet_filter   = 1L,
                           wavelet_family   = "DaubExPhase",
                           halfK            = 20L) {
  if (!inherits(fit, "mfsusie")) {
    stop("`fit` must be an `mfsusie` (or `fsusie`) fit object.")
  }
  if (is.null(fit$residuals) || is.null(fit$lead_X)) {
    stop("`fit$residuals` and `fit$lead_X` are required.")
  }
  method <- match.arg(method)
  if (level <= 0 || level >= 1) {
    stop("`level` must be in (0, 1).")
  }

  switch(method,
         "scalewise" = .post_smooth_scalewise(fit, level, threshold_factor),
         "TI"        = .post_smooth_ti(fit, level, wavelet_filter,
                                       wavelet_family),
         "HMM"       = .post_smooth_hmm(fit, level, halfK))
}

# ---- scalewise -----------------------------------------------------

.post_smooth_scalewise <- function(fit, level, threshold_factor) {
  z_crit <- stats::qnorm((1 + level) / 2)
  meta   <- fit$dwt_meta
  M      <- length(meta$T_basis)
  L      <- nrow(fit$alpha)

  effect_curves  <- vector("list", M)
  credible_bands <- vector("list", M)

  for (m in seq_len(M)) {
    T_m <- meta$T_basis[m]
    effect_curves[[m]]  <- vector("list", L)
    credible_bands[[m]] <- vector("list", L)

    for (l in seq_len(L)) {
      lead_l <- which.max(fit$alpha[l, ])
      mean_w <- fit$mu[[l]][[m]][lead_l, ]
      var_w  <- pmax(fit$mu2[[l]][[m]][lead_l, ] - mean_w^2, 0)

      shrunk_w <- if (T_m == 1L) mean_w else
        .scalewise_soft_threshold(mean_w, sd = sqrt(var_w),
          scale_index = meta$scale_index[[m]],
          T_padded    = T_m,
          factor      = threshold_factor)

      effect_curves[[m]][[l]] <-
        .invert_packed_curve(matrix(shrunk_w, nrow = 1L), meta, m)
      sd_pos <- abs(.invert_packed_curve(matrix(sqrt(var_w), nrow = 1L),
                                         meta, m))
      mean_pos <- effect_curves[[m]][[l]]
      credible_bands[[m]][[l]] <- cbind(mean_pos - z_crit * sd_pos,
                                        mean_pos + z_crit * sd_pos)
    }
  }
  fit$effect_curves  <- effect_curves
  fit$credible_bands <- credible_bands
  fit$lfsr_curves    <- NULL
  fit
}

# ---- TI: cycle-spinning translation-invariant wavelet denoising ----

.post_smooth_ti <- function(fit, level, wavelet_filter, wavelet_family) {
  z_crit <- stats::qnorm((1 + level) / 2)
  meta   <- fit$dwt_meta
  M      <- length(meta$T_basis)
  L      <- nrow(fit$alpha)

  effect_curves  <- vector("list", M)
  credible_bands <- vector("list", M)

  for (m in seq_len(M)) {
    T_m <- meta$T_basis[m]
    effect_curves[[m]]  <- vector("list", L)
    credible_bands[[m]] <- vector("list", L)
    if (T_m == 1L) {
      # Scalar outcome: TI is a no-op; fall back to scalewise.
      tmp <- .post_smooth_scalewise(fit, level, threshold_factor = 1)
      effect_curves[[m]]  <- tmp$effect_curves[[m]]
      credible_bands[[m]] <- tmp$credible_bands[[m]]
      next
    }

    Y_pos <- .iso_response_pos(fit, m)   # n x T_m: residual + sum_l lead_l*mu_l_lead
    for (l in seq_len(L)) {
      lead_l <- which.max(fit$alpha[l, ])
      mu_lead_w <- fit$mu[[l]][[m]][lead_l, ]
      # The "per-effect" position-space response: subtract every
      # other effect's lead-variable contribution from Y_pos.
      iso_pos <- Y_pos - .other_effects_pos(fit, m, exclude = l)
      x_lead  <- fit$lead_X[[l]]

      out <- .univariate_ti_regression(iso_pos, x_lead,
                                       wavelet_filter, wavelet_family,
                                       z_crit)
      effect_curves[[m]][[l]]  <- out$effect_estimate
      credible_bands[[m]][[l]] <- cbind(out$cred_band[2L, ],   # low
                                        out$cred_band[1L, ])   # up
      # Re-order to [lower, upper]:
      credible_bands[[m]][[l]] <- credible_bands[[m]][[l]][,
                                                          c(1L, 2L)]
    }
  }
  fit$effect_curves  <- effect_curves
  fit$credible_bands <- credible_bands
  fit$lfsr_curves    <- NULL
  fit
}

# Reconstruct the position-space response by inverting
# `D = residuals + fitted` (wavelet domain).
.iso_response_pos <- function(fit, m) {
  meta <- fit$dwt_meta
  # Wavelet-domain D = residuals + fitted.
  D_w <- fit$residuals[[m]] + fit$fitted[[m]]
  # Inverse-DWT row-by-row to get position-space Y.
  inv <- mf_invert_dwt(D_packed      = D_w,
                       column_center = meta$column_center[[m]],
                       column_scale  = meta$column_scale[[m]],
                       filter_number = meta$wavelet_filter,
                       family        = meta$wavelet_family)
  if (ncol(inv) > length(meta$pos[[m]])) {
    inv <- inv[, seq_along(meta$pos[[m]]), drop = FALSE]
  }
  inv
}

# Sum_{l != exclude} x_lead_l * inverse_DWT(mu_l[lead_l, ]).
.other_effects_pos <- function(fit, m, exclude) {
  meta <- fit$dwt_meta
  L    <- nrow(fit$alpha)
  T_pos <- length(meta$pos[[m]])
  out  <- matrix(0, nrow = length(fit$lead_X[[1]]), ncol = T_pos)
  for (l in seq_len(L)) {
    if (l == exclude) next
    lead_l <- which.max(fit$alpha[l, ])
    mu_w   <- fit$mu[[l]][[m]][lead_l, ]
    eff_pos <- .invert_packed_curve(matrix(mu_w, nrow = 1L), meta, m)
    if (length(eff_pos) > T_pos) eff_pos <- eff_pos[seq_len(T_pos)]
    out <- out + outer(fit$lead_X[[l]], eff_pos)
  }
  out
}

# Stationary-wavelet regression of one outcome's isolated
# response on the lead variable, with scalewise ash shrinkage of
# the wavelet coefficients and cycle-spinning average for the
# point estimate. Returns `effect_estimate` (length T_pos) and
# `cred_band` (2 x T_pos: row 1 "up", row 2 "low").
.univariate_ti_regression <- function(Y_pos, x_lead,
                                      filter_number, family,
                                      z_crit) {
  # Stationary wavelet transform of each row of Y_pos.
  T_pos <- ncol(Y_pos)
  dummy <- wavethresh::wd(Y_pos[1L, ], type = "station",
                          filter.number = filter_number,
                          family        = family)
  Y_f <- do.call(rbind, lapply(seq_len(nrow(Y_pos)), function(i)
    wavethresh::wd(Y_pos[i, ], type = "station",
                   filter.number = filter_number,
                   family        = family)$D))
  Y_c <- do.call(rbind, lapply(seq_len(nrow(Y_pos)), function(i)
    wavethresh::wd(Y_pos[i, ], type = "station",
                   filter.number = filter_number,
                   family        = family)$C))

  # Univariate regression of each wavelet column on x_lead.
  ss_x   <- sum(x_lead^2)
  bhat_d <- as.numeric(crossprod(x_lead, Y_f)) / ss_x
  bhat_c <- as.numeric(crossprod(x_lead, Y_c)) / ss_x
  resid_d <- Y_f - tcrossprod(x_lead, bhat_d)
  resid_c <- Y_c - tcrossprod(x_lead, bhat_c)
  shat_d <- sqrt(colMeans(resid_d^2) / ss_x)
  shat_c <- sqrt(colMeans(resid_c^2) / ss_x)

  # Scalewise ash shrinkage of D coefficients.
  fl  <- dummy$fl.dbase$first.last.d
  n   <- 2L^wavethresh::nlevelsWT(dummy)
  K   <- nrow(fl)
  wd_post <- numeric(length(bhat_d))
  wd_var  <- numeric(length(bhat_d))
  for (s in seq_len(K)) {
    first_s   <- fl[s, 1L]
    offset_s  <- fl[s, 3L]
    idx <- (offset_s + 1L - first_s):(offset_s + n - first_s)
    t_ash <- ashr::ash(bhat_d[idx], shat_d[idx],
                       nullweight = 30, mixcompdist = "normal")
    wd_post[idx] <- t_ash$result$PosteriorMean
    wd_var[idx]  <- t_ash$result$PosteriorSD^2
  }
  t_ash_c <- ashr::ash(bhat_c, shat_c,
                       nullweight = 3, mixcompdist = "normal")

  # Cycle-spinning average via av.basis.
  dummy$D <- wd_post
  dummy$C <- t_ash_c$result$PosteriorMean
  mywst <- wavethresh::convert(dummy)
  fitted_func <- wavethresh::av.basis(mywst,
                                      level = dummy$nlevels - 1L,
                                      ix1 = 0, ix2 = 1,
                                      filter = mywst$filter)

  # Exact pointwise variance via the squared-filter wd /
  # convert / AvBasis pipeline.
  var_wd       <- wd_variance(rep(0, T_pos),
                              filter.number = filter_number,
                              family        = family)
  var_wd$D     <- wd_var
  fitted_var   <- av_basis_variance(wst_variance(var_wd))
  fitted_sd    <- sqrt(pmax(fitted_var, 0))

  cred_band <- rbind(up  = fitted_func + z_crit * fitted_sd,
                     low = fitted_func - z_crit * fitted_sd)
  list(effect_estimate = fitted_func, cred_band = cred_band,
       fitted_sd = fitted_sd)
}

# ---- HMM denoising -------------------------------------------------

.post_smooth_hmm <- function(fit, level, halfK) {
  z_crit <- stats::qnorm((1 + level) / 2)
  meta   <- fit$dwt_meta
  M      <- length(meta$T_basis)
  L      <- nrow(fit$alpha)

  effect_curves  <- vector("list", M)
  credible_bands <- vector("list", M)
  lfsr_curves    <- vector("list", M)

  for (m in seq_len(M)) {
    T_m <- meta$T_basis[m]
    effect_curves[[m]]  <- vector("list", L)
    credible_bands[[m]] <- vector("list", L)
    lfsr_curves[[m]]    <- vector("list", L)

    if (T_m == 1L) {
      tmp <- .post_smooth_scalewise(fit, level, threshold_factor = 1)
      effect_curves[[m]]  <- tmp$effect_curves[[m]]
      credible_bands[[m]] <- tmp$credible_bands[[m]]
      lfsr_curves[[m]]    <- replicate(L, NULL, simplify = FALSE)
      next
    }

    Y_pos <- .iso_response_pos(fit, m)
    for (l in seq_len(L)) {
      iso_pos <- Y_pos - .other_effects_pos(fit, m, exclude = l)
      x_lead  <- fit$lead_X[[l]]
      out     <- .univariate_hmm_regression(iso_pos, x_lead, halfK,
                                            z_crit)
      effect_curves[[m]][[l]]  <- out$effect_estimate
      credible_bands[[m]][[l]] <- out$cred_band
      lfsr_curves[[m]][[l]]    <- out$lfsr
    }
  }
  fit$effect_curves  <- effect_curves
  fit$credible_bands <- credible_bands
  fit$lfsr_curves    <- lfsr_curves
  fit
}

.univariate_hmm_regression <- function(Y_pos, x_lead, halfK, z_crit) {
  ss_x  <- sum(x_lead^2)
  bhat  <- as.numeric(crossprod(x_lead, Y_pos)) / ss_x
  resid <- Y_pos - tcrossprod(x_lead, bhat)
  shat  <- sqrt(colMeans(resid^2) / ss_x)

  bad <- !is.finite(shat) | shat <= 0
  if (any(bad)) {
    bhat[bad] <- 0
    shat[bad] <- median(shat[!bad], na.rm = TRUE)
    if (!is.finite(median(shat[!bad], na.rm = TRUE))) shat[bad] <- 1
  }

  s <- .fit_hmm(bhat, shat, halfK = halfK)
  cred_band <- cbind(s$x_post - z_crit * s$x_sd,
                     s$x_post + z_crit * s$x_sd)
  list(effect_estimate = s$x_post,
       cred_band       = cred_band,
       lfsr            = s$lfsr)
}

# Forward-backward over a discrete grid of effect-size values
# (-halfK*max_sd, +halfK*max_sd) with Gaussian emission and a
# uniform prior across grid points. Returns the smoothed
# posterior mean, posterior SD, and per-position lfsr.
.fit_hmm <- function(x, sd, halfK = 20L) {
  T_pos <- length(x)
  max_sd <- max(sd, na.rm = TRUE)
  if (!is.finite(max_sd) || max_sd <= 0) max_sd <- 1
  grid <- seq(-halfK, halfK, length.out = 2L * halfK + 1L) * max_sd

  # Emission log-prob: T x K
  emit <- vapply(seq_along(grid),
                 function(k) stats::dnorm(x, mean = grid[k], sd = sd,
                                          log = TRUE),
                 numeric(T_pos))

  # Forward-backward with uniform transition (i.e. iid posterior
  # smoothed by the Gaussian mixture prior over the grid).
  emit_max <- apply(emit, 1L, max)
  w <- exp(emit - emit_max)
  w <- w / rowSums(w)

  x_post <- as.numeric(w %*% grid)
  x2_post <- as.numeric(w %*% grid^2)
  x_var <- pmax(x2_post - x_post^2, 0)
  x_sd  <- sqrt(x_var)

  # lfsr: P(sign != sign(x_post)) under the smoothed posterior.
  pos_w <- w[, grid >  0, drop = FALSE]
  neg_w <- w[, grid <  0, drop = FALSE]
  zero_w <- w[, grid == 0, drop = FALSE]
  p_pos  <- rowSums(pos_w)
  p_neg  <- rowSums(neg_w)
  p_zero <- rowSums(zero_w)
  lfsr <- ifelse(x_post > 0,
                 p_neg + p_zero,
                 ifelse(x_post < 0,
                        p_pos + p_zero,
                        1))
  list(x_post = x_post, x_sd = x_sd, lfsr = lfsr)
}

# --- helpers --------------------------------------------------------

# Soft-threshold the packed wavelet vector per scale at
# `factor * sd * sqrt(2 log T)`.
.scalewise_soft_threshold <- function(coef_vec, sd, scale_index,
                                      T_padded, factor) {
  out <- coef_vec
  T_eff <- T_padded
  for (idx in scale_index) {
    if (length(idx) == 0L) next
    sigma <- mean(sd[idx], na.rm = TRUE)
    if (!is.finite(sigma) || sigma == 0) next
    thr <- factor * sigma * sqrt(2 * log(T_eff))
    x <- out[idx]
    out[idx] <- sign(x) * pmax(abs(x) - thr, 0)
  }
  out
}

# Use the existing mf_invert_dwt machinery to bring a length-T_padded
# packed wavelet vector back to position space.
.invert_packed_curve <- function(D_packed_row, meta, m) {
  inverted <- mf_invert_dwt(
    D_packed      = D_packed_row,
    column_center = rep(0, ncol(D_packed_row)),
    column_scale  = rep(1, ncol(D_packed_row)),
    filter_number = meta$wavelet_filter %||% 10L,
    family        = meta$wavelet_family %||% "DaubLeAsymm")
  as.numeric(inverted)
}
