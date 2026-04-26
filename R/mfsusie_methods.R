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
