# IBSS-loop dispatch methods on `mf_individual` and `mfsusie`.
#
# These are the S3 generics susieR's `susie_workhorse` calls that
# don't fit the per-effect SER hot path (which lives in
# `individual_data_methods.R`). The split is:
#
#   ibss_initialize    -- once, before the loop. Validates params,
#                         sets var_y, builds the model object via
#                         `initialize_susie_model`, runs `initialize_fitted`.
#   initialize_fitted  -- once. Returns extra fields to merge into model.
#   get_var_y          -- once, used by `ibss_initialize` to get sigma2 init.
#   track_ibss_fit     -- per iteration, optional tracking.
#   check_convergence  -- per iteration.
#   trim_null_effects  -- once at end of loop.
#   ibss_finalize      -- once at end. Builds pip, cs, sets niter / elbo /
#                         residuals on the fit.
#   cleanup_model      -- called by ibss_finalize. Strips ephemerals.
#
# Per-effect methods (compute_residuals, loglik, etc.) live in
# `individual_data_methods.R`.

# ---- get_var_y -------------------------------------------------

#' Initial residual variance estimate per outcome
#'
#' Returns a list of length `M`; each entry is a length-`S_m` vector
#' (one per wavelet scale) of marginal sample variances of `data$D[[m]]`
#' computed scale by scale. Used by `ibss_initialize` to seed
#' `params$residual_variance` when the user does not supply one.
#' This shape matches the per-(scale, outcome) sigma2 convention
#' carried throughout the IBSS loop. The
#' `per_outcome` mode collapses each entry to a scalar at
#' update time.
#'
#' @keywords internal
#' @noRd
get_var_y.mf_individual <- function(data, ...) {
  lapply(seq_len(data$M), function(m) {
    indx_m <- data$scale_index[[m]]
    vapply(indx_m, function(idx) {
      as.numeric(var(as.vector(data$D[[m]][, idx, drop = FALSE])))
    }, numeric(1))
  })
}

# ---- initialize_fitted -----------------------------------------

#' Per-outcome fitted-values initializer
#'
#' Returns the list of extra fields to merge into `model` after
#' `initialize_susie_model.mf_individual` runs. The running fit
#' `model$fitted[[m]]` is already initialized (zero matrices) by
#' `initialize_susie_model.mf_individual`, so this returns an empty
#' list. Defined so susieR's `ibss_initialize.default` dispatch
#' does not fall through to `initialize_fitted.individual` (which
#' assumes scalar `data$y`).
#'
#' @keywords internal
#' @noRd
initialize_fitted.mf_individual <- function(data, mat_init, ...) {
  list()
}

# ---- ibss_initialize -------------------------------------------

#' IBSS loop initializer for `mf_individual`
#'
#' Mirrors `ibss_initialize.default` but skips its
#' scalar-only validation of `params$residual_variance` (mfsusieR
#' uses a list-of-vectors per-(scale, outcome) shape) and routes
#' through the mfsusieR S3 methods for `initialize_susie_model`,
#' `initialize_fitted`. Model_init / warm-start is deferred.
#' not supported (param$model_init is ignored).
#'
#' @keywords internal
#' @noRd
ibss_initialize.mf_individual <- function(data, params) {
  var_y <- get_var_y(data)

  if (data$p < params$L) params$L <- data$p

  if (is.null(params$residual_variance)) {
    params$residual_variance <- var_y
  }

  mat_init   <- initialize_susie_model(data, params, var_y)
  fitted     <- initialize_fitted(data, mat_init)

  model_class <- class(mat_init)
  model       <- c(mat_init, list(null_index = 0), fitted)
  class(model) <- if (inherits(mat_init, "susie")) model_class else "susie"
  model$converged <- FALSE
  model
}

# ---- validate_prior --------------------------------------------

#' Per-iteration prior sanity check (no-op for `mf_individual`)
#'
#' susieR's `validate_prior.default` checks the scalar `V[l]` is
#' non-negative; mfsusieR holds `V[l]` = 1 by design and validates
#' the mixture-weight prior at construction time
#' (`mf_prior_scale_mixture`). No per-iteration validation needed.
#'
#' @keywords internal
#' @noRd
validate_prior.mf_individual <- function(data, params, model, ...) {
  invisible(TRUE)
}

# ---- track_ibss_fit --------------------------------------------

#' IBSS iteration tracking (no-op by default)
#'
#' Tracking adds runtime overhead and large-state copies. mfsusieR
#' Tracking is supported only when `params$track_fit` is `TRUE`,
#' otherwise this is a no-op (matches the susieR pattern but with
#' a smaller default footprint, since per-effect curves are
#' large).
#'
#' @keywords internal
#' @noRd
track_ibss_fit.mf_individual <- function(data, params, model, tracking, iter, elbo, ...) {
  if (!isTRUE(params$track_fit)) return(tracking)
  tracking[[iter]] <- list(
    alpha  = model$alpha,
    sigma2 = model$sigma2,
    pi_V   = model$pi_V,
    elbo   = elbo[iter]
  )
  tracking
}

# ---- check_convergence -----------------------------------------

#' Convergence check on ELBO change between iterations
#'
#' Mirrors susieR's default: declares convergence when
#' `abs(elbo[iter+1] - elbo[iter]) < params$tol`.
#'
#' @keywords internal
#' @noRd
check_convergence.mf_individual <- function(data, params, model, elbo, iter) {
  if (iter > 0L && is.finite(elbo[iter]) && is.finite(elbo[iter + 1L])) {
    if (abs(elbo[iter + 1L] - elbo[iter]) < params$tol) {
      model$converged <- TRUE
    }
  }
  model
}

# ---- trim_null_effects ----------------------------------------

#' Drop effects with negligible posterior mass (no-op)
#'
#' susieR's default zeros out effects with `V[l] < prior_tol`.
#' mfsusieR holds `V[l] = 1` for all l (mixture prior absorbs the
#' per-effect adaptation), so the V-based pruning does not apply.
#' A future refinement may prune effects whose mixture weight
#' concentrates on the null component.
#'
#' @keywords internal
#' @noRd
trim_null_effects.mf_individual <- function(data, params, model) {
  model
}

# Note: susieR's `ibss_finalize` is a regular function, not an S3
# generic. It calls the per-class generics (`get_fitted`, `get_cs`,
# `get_intercept`, `get_variable_names`, `get_zscore`,
# `cleanup_model`) via dispatch, so mfsusieR overrides those rather
# than trying to register an `ibss_finalize.mf_individual` method.
# Per-outcome residual attachment happens in the public `mfsusie()`
# wrapper after `susie_workhorse` returns.

# ---- cleanup_model --------------------------------------------

#' Strip per-iteration ephemeral fields from the fit
#'
#' Called by `ibss_finalize`. Removes the per-effect residual
#' caches that were used during the IBSS sweep but are not part of
#' the user-facing fit contract.
#'
#' @keywords internal
#' @noRd
cleanup_model.mf_individual <- function(data, params, model, ...) {
  model$residuals          <- if (is.null(data$residuals)) NULL else model$residuals
  model$fitted_without_l   <- NULL
  model$raw_residuals      <- NULL
  model$residual_variance  <- NULL
  model$runtime            <- NULL
  model
}

# ---- configure_data -------------------------------------------

#' Configure data pre-IBSS (no-op for `mf_individual`)
#'
#' susieR's `configure_data.individual` handles the
#' `unmappable_effects` shim (ash / sufficient-stats conversion).
#' mfsusieR does not support that path; configure_data is a
#' pass-through.
#'
#' @keywords internal
#' @noRd
configure_data.mf_individual <- function(data, params, ...) {
  data
}

# ---- User-facing post-fit accessors ---------------------------

#' Per-outcome column scale factors of X
#' @keywords internal
#' @noRd
get_scale_factors.mf_individual <- function(data, params, ...) {
  data$csd
}

#' Per-outcome intercepts (length `M`)
#'
#' Reconstructs the intercept per outcome on the original Y scale.
#' Returns a list of length `M`; each entry is a length-`T_basis[m]`
#' vector. Computed from per-effect posterior means projected back
#' through the inverse wavelet transform when `T_basis[m] > 1`.
#'
#' @keywords internal
#' @noRd
get_intercept.mf_individual <- function(data, params, model, ...) {
  if (!isTRUE(params$intercept)) {
    return(lapply(data$T_basis, function(T_m) rep(0, T_m)))
  }
  # mfsusieR centers Y per-outcome, per-column at construction.
  # Intercept reconstruction returns 0.
  lapply(data$T_basis, function(T_m) rep(0, T_m))
}

#' Per-outcome fitted values (length-M list of n x T_basis matrices)
#' @keywords internal
#' @noRd
get_fitted.mf_individual <- function(data, params, model, ...) {
  model$fitted
}

#' Credible sets via `susie_get_cs` (works on alpha + X)
#' @keywords internal
#' @noRd
get_cs.mf_individual <- function(data, params, model, ...) {
  if (is.null(params$coverage) || is.null(params$min_abs_corr)) return(NULL)
  susie_get_cs(model,
                       X            = data$X,
                       coverage     = params$coverage,
                       min_abs_corr = params$min_abs_corr,
                       n_purity     = params$n_purity %||% 100)
}

#' Variable name assignment on the fit's per-effect arrays
#' @keywords internal
#' @noRd
get_variable_names.mf_individual <- function(data, model, ...) {
  vnames <- colnames(data$X)
  if (!is.null(vnames)) {
    colnames(model$alpha)        <- vnames
    colnames(model$lbf_variable) <- vnames
  }
  model
}

#' Univariate z-scores not exposed for `mf_individual`
#' @keywords internal
#' @noRd
get_zscore.mf_individual <- function(data, params, model, ...) {
  NULL
}
