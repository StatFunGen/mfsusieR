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
# attaches the post-fit fields (`pip`, `cs`, `elbo`,
# `niter`, `converged`).
#
# @references
# Manuscript: methods/online_method.tex
# (per-outcome wavelet decomposition + IBSS dispatch).

#' Initialize the `mfsusie` model from an `mf_individual` data class
#'
#' @param data an `mf_individual` object from `create_mf_individual`.
#' @param params named list of fit-time parameters. Required fields:
#'   `L` (effect-count upper bound), `prior_weights` (length-p
#'   variable-selection prior, defaults to uniform `1/p` if NULL),
#'   `prior` (an `mf_prior_scale_mixture` object containing
#'   `V_grid`, `pi`, `null_prior_weight`), and
#'   `residual_variance` (initial sigma2; either a length-M list
#'   of scalars for the per-outcome mode or a length-M list of
#'   length-`S_m` vectors for the per-(scale, outcome) mode).
#' @param var_y per-outcome variance of the observed response
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
  T_basis <- data$T_basis

  # Per-effect posterior moments. mu[[l]][[m]] and mu2[[l]][[m]] are
  # p x T_basis[m] matrices; the nested-list shape lets ragged T_m
  # coexist without padding to a common T.
  mu  <- vector("list", L)
  mu2 <- vector("list", L)
  for (l in seq_len(L)) {
    mu[[l]]  <- vector("list", M)
    mu2[[l]] <- vector("list", M)
    for (m in seq_len(M)) {
      mu[[l]][[m]]  <- matrix(0, nrow = p, ncol = T_basis[m])
      mu2[[l]][[m]] <- matrix(0, nrow = p, ncol = T_basis[m])
    }
  }

  # Prior structure. V_grid: list[M] of length-K vectors (or list[M]
  # of S_m x K matrices for `per_scale`). pi_V: list[M] of
  # S_m x K mixture-weight matrices. null_prior_weight: scalar.
  prior <- params$prior
  V_grid             <- if (is.null(prior)) NULL else prior$V_grid
  pi_V               <- if (is.null(prior)) NULL else prior$pi
  null_prior_weight  <- if (is.null(prior)) 0   else prior$null_prior_weight
  G_prior            <- if (is.null(prior)) NULL else prior$G_prior

  # Cross-outcome combiner. Defaults to the trivial independence
  # combiner.
  cross_outcome_combiner <- if (is.null(params$cross_outcome_prior)) {
    cross_outcome_prior_independent()
  } else {
    params$cross_outcome_prior
  }

  # Residual variance. `sigma2` is a list[M]: scalar entries under
  # `per_outcome`, length-S_m vectors under
  # `per_scale` (the default). Falls back to `var_y` when
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

  # Per-outcome fitted values (running sum of effect curves).
  fitted_values <- vector("list", M)
  for (m in seq_len(M)) {
    fitted_values[[m]] <- matrix(0, nrow = data$n, ncol = T_basis[m])
  }

  # `lbf_variable_outcome` is an L x p x M array of per-(effect, variant,
  # outcome) log Bayes factors, populated by `loglik.mf_individual`
  # during the IBSS sweep. Always allocated; consumed by
  # `susie_post_outcome_configuration(fit, by = "outcome")`.
  lbf_variable_outcome <- array(NA_real_, dim = c(L, p, M))

  model <- list(
    alpha             = matrix(1 / p, nrow = L, ncol = p),
    mu                = mu,
    mu2               = mu2,
    V                 = V,
    V_grid            = V_grid,
    pi_V              = pi_V,
    null_prior_weight = null_prior_weight,
    G_prior           = G_prior,
    cross_outcome_combiner = cross_outcome_combiner,
    KL                = rep(NA_real_, L),
    lbf               = rep(NA_real_, L),
    lbf_variable      = matrix(NA_real_, nrow = L, ncol = p),
    lbf_variable_outcome       = lbf_variable_outcome,
    sigma2            = sigma2,
    pi                = prior_weights,
    fitted            = fitted_values,
    # Per-outcome SER scratch slots, populated in-place by
    # `compute_residuals.mf_individual` per effect.
    residuals         = vector("list", M),
    fitted_without_l  = vector("list", M),
    raw_residuals     = vector("list", M),
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
get_posterior_moments_l.mfsusie <- function(model, l) {
  list(
    mu  = model$mu[[l]],
    mu2 = model$mu2[[l]]
  )
}

#' Per-iteration mixture-null-mass summary for the verbose row
#'
#' Appends `pi_null=[min, max]` to the verbose tabular row during
#' IBSS, where the values are the smallest and largest null-component
#' weight across all (outcome, scale) cells of `model$pi_V`. The same
#' quantity is summarized by `summary.mfsusie()` after convergence;
#' surfacing it per-iteration lets the user watch the mixture prior
#' concentrate on the null component as the IBSS sweep progresses.
#'
#' @keywords internal
#' @noRd
format_extra_diag.mfsusie <- function(model) {
  if (is.null(model$pi_V)) return("")
  null_mass <- unlist(
    lapply(model$pi_V, function(piVm) piVm[, 1L]),
    use.names = FALSE)
  if (length(null_mass) == 0L) return("")
  sprintf("pi_null=[%.2f, %.2f]", min(null_mass), max(null_mass))
}

#' Compact verbose-row formatting for the list-of-vectors `sigma2`
#'
#' susieR's default `format_sigma2_summary` formats a single scalar
#' via `sprintf("%.4f", model$sigma2)`; mfsusieR's `sigma2` is a
#' length-`M` list of either scalars (per-outcome scope) or
#' length-`S_m` vectors (per-(scale, outcome) scope), so the default
#' would crash. We collapse to a min/median/max summary for readable
#' fixed-width verbose output.
#'
#' @keywords internal
#' @noRd
format_sigma2_summary.mfsusie <- function(model) {
  vals <- unlist(model$sigma2, use.names = FALSE)
  if (length(vals) == 0L) return("(empty)")
  if (length(vals) <= 5L) {
    paste0("[", paste(sprintf("%.3f", vals), collapse = ", "), "]")
  } else {
    sprintf("[%.3f, %.3f, %.3f]",
            min(vals), stats::median(vals), max(vals))
  }
}

#' Per-effect posterior mean curves (alpha-weighted mu)
#'
#' Returns a list of length M where element m is a `p x T_basis[m]`
#' matrix; row j of element m is `alpha[l, j] * mu[[l]][[m]][j, ]`,
#' the contribution of variable j to the posterior mean of effect l on
#' outcome m.
#' @keywords internal
#' @noRd
get_posterior_mean_l.mfsusie <- function(model, l) {
  alpha_l <- model$alpha[l, ]
  lapply(model$mu[[l]], function(mu_lm) sweep(mu_lm, 1, alpha_l, "*"))
}

#' Posterior mean curves summed over all L effects
#'
#' Returns a list of length M where element m is a `p x T_basis[m]`
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

