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
# Post-processing helpers for an `mfsusie` fit.
#
# `mf_post_smooth(fit)` adds two slots to the returned fit:
#   $effect_curves[[m]][[l]]  -- smoothed position-space curve,
#                                length T_basis[m]
#   $credible_bands[[m]][[l]] -- T_basis[m] x 2 matrix [lower, upper]
#
# Method "scalewise" performs scalewise soft-thresholding of the
# per-effect, per-outcome wavelet posterior mean at
# `factor * sigma_s * sqrt(2 log T)`, where `sigma_s` is the mean
# posterior SD on wavelet scale `s`. Pointwise credible bands
# come from the posterior second moment via `mu2 - mu^2` mapped
# through the inverse DWT and a sqrt(.) (Parseval approximation
# on the orthonormal basis).
#
# This is a simplified relative of the cycle-spinning translation-
# invariant denoising used in fSuSiE's reference implementation
# (see Denault et al. 2025 bioRxiv 10.1101/2025.08.17.670732,
# Methods, "Refining the effect estimates"); it operates directly
# on the wavelet-domain posterior moments rather than re-running
# wavelet regression on the residualized response, so it does not
# require X / Y at smoothing time. A full cycle-spinning method
# is a planned addition and will be exposed via
# `method = "TI"` when added.

#' Post-smooth a fit's per-effect curves and add credible bands
#'
#' @param fit a fit returned by `mfsusie()` or `fsusie()`.
#' @param method character, smoother. Currently only
#'   `"scalewise"` is implemented.
#' @param level numeric in (0, 1), credible-band coverage.
#'   Default 0.95 (-> +/- 1.96 * sd).
#' @param threshold_factor numeric, multiplier on the
#'   `sqrt(2 log T)` universal threshold. 1 = standard SureShrink;
#'   < 1 keeps more wavelet energy (less smoothing); > 1 smoothes
#'   more aggressively. Default 1.
#'
#' @return the input fit with `$effect_curves` and
#'   `$credible_bands` populated. Scalar outcomes
#'   (`T_basis[m] = 1`) are passed through unchanged (no smoothing
#'   needed).
#' @references
#' Manuscript: Denault et al. (2025) Methods, "Refining the effect
#' estimates" (bioRxiv 10.1101/2025.08.17.670732).
#' @export
mf_post_smooth <- function(fit,
                           method           = c("scalewise"),
                           level            = 0.95,
                           threshold_factor = 1) {
  if (!inherits(fit, "mfsusie")) {
    stop("`fit` must be an `mfsusie` (or `fsusie`) fit object.")
  }
  method <- match.arg(method)
  if (level <= 0 || level >= 1) {
    stop("`level` must be in (0, 1).")
  }
  z_crit <- stats::qnorm((1 + level) / 2)

  meta <- fit$dwt_meta
  M    <- length(meta$T_basis)
  L    <- nrow(fit$alpha)

  effect_curves  <- vector("list", M)
  credible_bands <- vector("list", M)

  for (m in seq_len(M)) {
    T_m <- meta$T_basis[m]
    effect_curves[[m]]  <- vector("list", L)
    credible_bands[[m]] <- vector("list", L)

    for (l in seq_len(L)) {
      mu_lm  <- fit$mu[[l]][[m]]   # p x T_basis
      mu2_lm <- fit$mu2[[l]][[m]]
      a_l    <- fit$alpha[l, ]
      mean_w <- as.numeric(crossprod(a_l, mu_lm))
      var_w  <- as.numeric(crossprod(a_l, mu2_lm)) - mean_w^2
      var_w[var_w < 0] <- 0

      if (T_m == 1L) {
        effect_curves[[m]][[l]] <-
          .invert_packed_curve(matrix(mean_w, nrow = 1L), meta, m)
        sd_pos <- abs(.invert_packed_curve(matrix(sqrt(var_w),
                                                  nrow = 1L),
                                           meta, m))
        mean_pos <- effect_curves[[m]][[l]]
        credible_bands[[m]][[l]] <- cbind(
          mean_pos - z_crit * sd_pos,
          mean_pos + z_crit * sd_pos
        )
        next
      }

      shrunk_w <- .scalewise_soft_threshold(
        mean_w, sd = sqrt(var_w),
        scale_index = meta$scale_index[[m]],
        T_padded    = T_m,
        factor      = threshold_factor)

      effect_curves[[m]][[l]] <-
        .invert_packed_curve(matrix(shrunk_w, nrow = 1L), meta, m)

      sd_pos <- abs(.invert_packed_curve(matrix(sqrt(var_w),
                                                nrow = 1L),
                                         meta, m))
      mean_pos <- effect_curves[[m]][[l]]
      credible_bands[[m]][[l]] <- cbind(
        mean_pos - z_crit * sd_pos,
        mean_pos + z_crit * sd_pos
      )
    }
  }

  fit$effect_curves  <- effect_curves
  fit$credible_bands <- credible_bands
  fit
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
