# Miscellaneous internal utilities (internal).
#
# Covers the `save_mu_method` storage policy: helpers for inspecting and
# applying the three mu/mu2 storage modes ("complete", "alpha_collapsed",
# "lead"), plus the post-fit `mf_thin()` trimmer.

# ---- save_mu_method helpers ------------------------------------------

mf_save_mu_method <- function(fit) {
  attr(fit, "save_mu_method") %||% "complete"
}

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

mf_has_per_variant_mu <- function(fit) {
  mf_save_mu_method(fit) == "complete"
}

stop_save_mu_method_combo <- function(operation, mode) {
  stop(sprintf(
    "%s requires save_mu_method = \"complete\"; got \"%s\". %s",
    operation, mode,
    "Refit with save_mu_method = \"complete\", or call mf_thin() only after running predict / per-variant lfsr / model_init."), call. = FALSE)
}

# ---- finalize-time trim -----------------------------------------------

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
        mu_collapsed  <- matrix(alpha_l %*% mu_lm,  nrow = 1L)
        mu2_collapsed <- matrix(alpha_l %*% mu2_lm, nrow = 1L)
        mu_raw_X <- sweep(mu_lm, 1L, X_scale, "/")
        coef_lm  <- matrix(alpha_l %*% mu_raw_X, nrow = 1L)
        fit$mu[[l]][[m]]       <- mu_collapsed
        fit$mu2[[l]][[m]]      <- mu2_collapsed
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

# ---- post-fit thinning ------------------------------------------------

#' @keywords internal
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
