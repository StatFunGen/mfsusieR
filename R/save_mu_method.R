# `save_mu_method` storage policy.
#
# Three modes shrink fit$mu and fit$mu2 in different ways at finalize:
#   - "complete":        keep p x T_basis[m] per-(effect, outcome). Default.
#                        The only mode that supports model_init warm-starts,
#                        predict.mfsusie(newx), and the per-variant lfsr toggle.
#   - "alpha_collapsed": replace each p x T matrix by the alpha-weighted
#                        1 x T summary `alpha[l, ] %*% mu_full`. Storage is
#                        factor-p smaller; coef.mfsusie / mf_post_smooth /
#                        summary / print / alpha-aggregated plots remain
#                        numerically equivalent to "complete". A separate
#                        coef_wavelet[[l]][[m]] (1 x T) is precomputed for the
#                        raw-X coef path, since per-j csd_X scaling cannot be
#                        recovered after alpha-collapse.
#   - "lead":            keep only the lead variable per effect, where
#                        j* = which.max(alpha[l, ]). mu[[l]][[m]] becomes
#                        mu_full[l, j*, ] (1 x T) and fit$top_index[l] = j*.
#                        Cheap single-variable summary, biased toward the lead.
#
# The trim runs in mfsusie() after susie_workhorse + dwt_meta attachment.
# An `attr(fit, "save_mu_method")` records the choice so downstream consumers
# (and mf_thin) can dispatch.

# ---- internal: consumer helpers ----------------------------------

# Returns the mode string for a fit. NULL fits or fits without the attribute
# are treated as "complete" for backward compatibility with serialized
# pre-feature objects.
mf_save_mu_method <- function(fit) {
  attr(fit, "save_mu_method") %||% "complete"
}

# Returns the alpha-aggregated wavelet-domain (mean_w, var_w) for effect l,
# outcome m. Length-T_basis[m] vectors.
#
# - complete: mean_w = alpha %*% mu, var_w = alpha %*% mu2 - mean_w^2
# - alpha_collapsed: mean_w = mu (1 x T), var_w = mu2 - mean_w^2
#   (math identity: alpha %*% mu_full = stored mu;
#    alpha %*% mu2_full - (alpha %*% mu_full)^2 = stored mu2 - stored mu^2)
# - lead: mean_w = mu (lead row), var_w = mu2 - mean_w^2 = posterior variance
#   of lead variant only. Biased; documented.
get_effect_wavelet_moments <- function(fit, l, m) {
  mu_lm  <- fit$mu[[l]][[m]]
  mu2_lm <- fit$mu2[[l]][[m]]
  if (NROW(mu_lm) == 1L) {
    mean_w <- as.numeric(mu_lm)
    var_w  <- as.numeric(mu2_lm) - mean_w^2
  } else {
    alpha_l <- fit$alpha[l, ]
    mean_w  <- as.numeric(alpha_l %*% mu_lm)
    var_w   <- as.numeric(alpha_l %*% mu2_lm) - mean_w^2
  }
  list(mean_w = mean_w, var_w = pmax(var_w, 0))
}

# Returns the per-effect raw-X wavelet-domain coefficient curve (length
# T_basis[m]) for effect l, outcome m. Used by coef.mfsusie.
#
# - complete: sum_j alpha[l, j] * mu[l, j, ] / csd_X[j], standard formula.
# - alpha_collapsed: read fit$coef_wavelet[[l]][[m]] directly (precomputed
#   at finalize because per-j csd_X scaling cannot be recovered post-collapse).
# - lead: mu[l, j*, ] / csd_X[j*], cheap and biased.
get_coef_wavelet_curve <- function(fit, l, m) {
  mode    <- mf_save_mu_method(fit)
  X_scale <- fit$dwt_meta$X_scale
  if (mode == "complete") {
    mu_raw_X <- sweep(fit$mu[[l]][[m]], 1L, X_scale, "/")
    return(as.numeric(fit$alpha[l, ] %*% mu_raw_X))
  }
  if (mode == "alpha_collapsed") {
    return(as.numeric(fit$coef_wavelet[[l]][[m]]))
  }
  if (mode == "lead") {
    j_star <- fit$top_index[l]
    return(as.numeric(fit$mu[[l]][[m]]) / X_scale[j_star])
  }
  stop(sprintf("Unknown save_mu_method on fit: %s", mode))
}

# Logical predicate. TRUE iff the fit carries per-variant mu/mu2 (so
# per-variant lfsr, predict(newx), and model_init can run).
mf_has_per_variant_mu <- function(fit) {
  mf_save_mu_method(fit) == "complete"
}

# Internal stop helper for the three call sites that require complete mu/mu2.
# Keeps wording consistent.
stop_save_mu_method_combo <- function(operation, mode) {
  stop(sprintf(
    "%s requires save_mu_method = \"complete\"; got \"%s\". %s",
    operation, mode,
    "Refit with save_mu_method = \"complete\", or call mf_thin() only after running predict / per-variant lfsr / model_init."), call. = FALSE)
}

# ---- internal: finalize-time trim --------------------------------

# Apply the storage policy to a fresh fit returned by susie_workhorse +
# dwt_meta attachment. Sets attr(fit, "save_mu_method") in all cases. Returns
# the modified fit.
mf_apply_save_mu_method <- function(fit, mode) {
  if (mode == "complete") {
    attr(fit, "save_mu_method") <- "complete"
    return(fit)
  }
  if (!(mode %in% c("alpha_collapsed", "lead"))) {
    stop(sprintf("Unknown save_mu_method: %s", mode))
  }

  L <- nrow(fit$alpha)
  M <- fit$dwt_meta$M
  X_scale <- fit$dwt_meta$X_scale

  if (mode == "alpha_collapsed") {
    coef_wavelet <- vector("list", L)
    for (l in seq_len(L)) {
      coef_wavelet[[l]] <- vector("list", M)
      alpha_l <- fit$alpha[l, ]
      for (m in seq_len(M)) {
        mu_lm  <- fit$mu[[l]][[m]]
        mu2_lm <- fit$mu2[[l]][[m]]
        # Standardised-X alpha-collapsed mu / mu2 (1 x T). Used by post_smooth.
        mu_collapsed  <- matrix(alpha_l %*% mu_lm,  nrow = 1L)
        mu2_collapsed <- matrix(alpha_l %*% mu2_lm, nrow = 1L)
        # Raw-X alpha-collapsed coefficient (1 x T). Used by coef.mfsusie.
        # Per-j csd_X division must happen before the alpha-collapse, so we
        # cannot recover this from mu_collapsed alone.
        mu_raw_X <- sweep(mu_lm, 1L, X_scale, "/")
        coef_lm  <- matrix(alpha_l %*% mu_raw_X, nrow = 1L)
        fit$mu[[l]][[m]]      <- mu_collapsed
        fit$mu2[[l]][[m]]     <- mu2_collapsed
        coef_wavelet[[l]][[m]] <- coef_lm
      }
    }
    fit$coef_wavelet <- coef_wavelet
    attr(fit, "save_mu_method") <- "alpha_collapsed"
    return(fit)
  }

  # mode == "lead"
  top_index <- integer(L)
  for (l in seq_len(L)) {
    j_star <- which.max(fit$alpha[l, ])
    top_index[l] <- j_star
    for (m in seq_len(M)) {
      fit$mu[[l]][[m]]  <- fit$mu[[l]][[m]][j_star, , drop = FALSE]
      fit$mu2[[l]][[m]] <- fit$mu2[[l]][[m]][j_star, , drop = FALSE]
    }
  }
  fit$top_index <- top_index
  attr(fit, "save_mu_method") <- "lead"
  fit
}

# ---- exported: post-fit thinning ---------------------------------

#' Thin an `mfsusie` fit by replacing per-variant posterior moments
#'
#' Applies the same storage trim as `mfsusie(save_mu_method = method)` but
#' as a post-fit operation, returning a new fit object. Lets callers keep a
#' complete fit for warm-starts and a thinned copy for distribution / saving.
#'
#' @param fit an `mfsusie` fit object. Must currently be in `"complete"`
#'   storage; thinning a fit that is already 1D is rejected.
#' @param method one of `"alpha_collapsed"` or `"lead"`. See `mfsusie()`'s
#'   `save_mu_method` documentation for the per-mode semantics.
#' @return a new `mfsusie` fit with the chosen storage shape. The original
#'   `fit` is not modified.
#' @details Both 1D modes drop information that the SER-step IBSS update
#'   needs, so the returned fit cannot be passed back as `model_init`. Use
#'   the thinned copy for `coef`, `mf_post_smooth`, `summary`, `print`, and
#'   alpha-aggregated plots; keep the complete fit for `predict(newx)`,
#'   per-variant lfsr, and warm-start checkpoints.
#' @examples
#' \dontrun{
#'   fit_full <- mfsusie(X, Y)                              # complete
#'   saveRDS(fit_full, "checkpoint.rds")                    # warm-start ready
#'   fit_thin <- mf_thin(fit_full, method = "alpha_collapsed")
#'   saveRDS(fit_thin, "for_distribution.rds")              # ~factor p smaller
#' }
#' @export
mf_thin <- function(fit, method = c("alpha_collapsed", "lead")) {
  if (!inherits(fit, "mfsusie")) {
    stop("`fit` must be an `mfsusie` object.")
  }
  method  <- match.arg(method)
  current <- mf_save_mu_method(fit)
  if (current != "complete") {
    stop(sprintf(
      "mf_thin() requires a fit with save_mu_method = \"complete\"; got \"%s\". A fit can only be thinned once.",
      current), call. = FALSE)
  }
  mf_apply_save_mu_method(fit, method)
}
