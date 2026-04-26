# `mfsusie` model class.
#
# The model object returned by `mfsusie()` and threaded through
# `susie_workhorse` via the per-iteration S3 dispatch. The
# initialize method `initialize_susie_model.mf_individual` is the
# constructor; `ibss_fit` calls it once at the start of the
# IBSS loop with the data class (`mf_individual`) and the params
# list.
#
# Fields are sized at construction; the IBSS loop mutates them
# in-place per iteration. Finalize (`ibss_finalize.mf_individual`,
# PR group 7) attaches the post-fit fields (`pip`, `cs`, `elbo`,
# `niter`, `converged`, `residuals`).
#
# @references
# Manuscript: methods/online_method.tex
# (per-modality wavelet decomposition + IBSS dispatch).

#' Initialize the `mfsusie` model from an `mf_individual` data class
#'
#' @param data an `mf_individual` object from `create_mf_individual`.
#' @param params named list of fit-time parameters. Required fields:
#'   `L` (effect-count upper bound), `prior_weights` (length-p
#'   variable-selection prior, defaults to uniform `1/p` if NULL),
#'   `prior` (an `mf_prior_scale_mixture` object containing
#'   `V_grid`, `pi`, `null_prior_weight`), and
#'   `residual_variance` (initial sigma2; either a length-M list
#'   of scalars for legacy mode or a length-M list of length-`S_m`
#'   vectors for per-(scale, modality) mode).
#' @param var_y per-modality variance of the observed response
#'   curves, used as a fallback when `params$residual_variance` is
#'   NULL.
#' @param ... ignored.
#' @return list of class `c("mfsusie", "susie")` with the IBSS-
#'   iterated fields pre-allocated.
#' @keywords internal
#' @noRd
initialize_susie_model.mf_individual <- function(data, params, var_y, ...) {
  L <- as.integer(params$L)
  p <- data$p
  M <- data$M
  T_padded <- data$T_padded

  # Per-effect posterior moments. mu[[l]][[m]] and mu2[[l]][[m]] are
  # p x T_padded[m] matrices; the nested-list shape lets ragged T_m
  # coexist without padding to a common T.
  mu  <- vector("list", L)
  mu2 <- vector("list", L)
  for (l in seq_len(L)) {
    mu[[l]]  <- vector("list", M)
    mu2[[l]] <- vector("list", M)
    for (m in seq_len(M)) {
      mu[[l]][[m]]  <- matrix(0, nrow = p, ncol = T_padded[m])
      mu2[[l]][[m]] <- matrix(0, nrow = p, ncol = T_padded[m])
    }
  }

  # Prior structure. V_grid: list[M] of length-K vectors (or list[M]
  # of S_m x K matrices for `per_scale_modality`). pi_V: list[M] of
  # S_m x K mixture-weight matrices. null_prior_weight: scalar.
  prior <- params$prior
  V_grid             <- if (is.null(prior)) NULL else prior$V_grid
  pi_V               <- if (is.null(prior)) NULL else prior$pi
  null_prior_weight  <- if (is.null(prior)) 0   else prior$null_prior_weight
  G_prior            <- if (is.null(prior)) NULL else prior$G_prior

  # Cross-modality combiner. Defaults to the trivial independence
  # combiner.
  cross_modality_combiner <- if (is.null(params$cross_modality_prior)) {
    cross_modality_prior_independent()
  } else {
    params$cross_modality_prior
  }

  # Residual variance. `sigma2` is a list[M]: scalar entries under
  # `shared_per_modality`, length-S_m vectors under
  # `per_scale_modality` (the default). Falls back to `var_y` when
  # `params$residual_variance` is NULL.
  sigma2 <- params$residual_variance
  if (is.null(sigma2)) {
    sigma2 <- var_y
  }

  # V[l] scales the entire mixture grid for effect l. The IBSS loop
  # updates it via `update_model_variance.mf_individual`.
  V <- rep(1, L)

  # Variable-selection prior over the p predictors.
  prior_weights <- params$prior_weights
  if (is.null(prior_weights)) {
    prior_weights <- rep(1 / p, p)
  }

  # Per-modality fitted values (running sum of effect curves).
  fitted_values <- vector("list", M)
  for (m in seq_len(M)) {
    fitted_values[[m]] <- matrix(0, nrow = data$n, ncol = T_padded[m])
  }

  model <- list(
    alpha             = matrix(1 / p, nrow = L, ncol = p),
    mu                = mu,
    mu2               = mu2,
    V                 = V,
    V_grid            = V_grid,
    pi_V              = pi_V,
    null_prior_weight = null_prior_weight,
    G_prior           = G_prior,
    cross_modality_combiner = cross_modality_combiner,
    KL                = rep(NA_real_, L),
    lbf               = rep(NA_real_, L),
    lbf_variable      = matrix(NA_real_, nrow = L, ncol = p),
    sigma2            = sigma2,
    pi                = prior_weights,
    fitted            = fitted_values,
    intercept         = rep(0, M),
    L                 = L,
    p                 = p,
    M                 = M
  )
  class(model) <- c("mfsusie", "susie")
  model
}


# =============================================================================
# Per-effect getter / setter accessors on the `mfsusie` model
#
# Used by the IBSS loop and by user-facing inspection helpers.
# Registered on the susieR namespace by `.onLoad`.
# =============================================================================

#' @keywords internal
#' @noRd
get_alpha_l.mfsusie <- function(model, l) {
  model$alpha[l, ]
}

#' @keywords internal
#' @noRd
get_posterior_moments_l.mfsusie <- function(model, l) {
  list(
    mu  = model$mu[[l]],
    mu2 = model$mu2[[l]]
  )
}

#' Per-effect posterior mean curves (alpha-weighted mu)
#'
#' Returns a list of length M where element m is a `p x T_padded[m]`
#' matrix; row j of element m is `alpha[l, j] * mu[[l]][[m]][j, ]`,
#' the contribution of SNP j to the posterior mean of effect l on
#' modality m.
#' @keywords internal
#' @noRd
get_posterior_mean_l.mfsusie <- function(model, l) {
  alpha_l <- model$alpha[l, ]
  lapply(model$mu[[l]], function(mu_lm) sweep(mu_lm, 1, alpha_l, "*"))
}

#' Posterior mean curves summed over all L effects
#'
#' Returns a list of length M where element m is a `p x T_padded[m]`
#' matrix giving `sum_l alpha[l, j] * mu[[l]][[m]][j, ]`. Used by
#' `get_posterior_mean_sum` callers in the IBSS dispatch.
#' @keywords internal
#' @noRd
get_posterior_mean_sum.mfsusie <- function(model) {
  L <- model$L
  M <- model$M
  out <- vector("list", M)
  for (m in seq_len(M)) {
    out[[m]] <- matrix(0, nrow = nrow(model$mu[[1]][[m]]),
                          ncol = ncol(model$mu[[1]][[m]]))
  }
  for (l in seq_len(L)) {
    pm_l <- get_posterior_mean_l.mfsusie(model, l)
    for (m in seq_len(M)) {
      out[[m]] <- out[[m]] + pm_l[[m]]
    }
  }
  out
}

#' @keywords internal
#' @noRd
get_prior_variance_l.mfsusie <- function(model, l) {
  model$V[l]
}

#' @keywords internal
#' @noRd
set_prior_variance_l.mfsusie <- function(model, l, V) {
  model$V[l] <- V
  model
}
