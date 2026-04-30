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
      v <- as.numeric(var(as.vector(data$D[[m]][, idx, drop = FALSE]),
                          na.rm = TRUE))
      if (!is.finite(v) || v == 0) 1 else v
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

# ---- expand_model_init_to_L -----------------------------------

#' Expand a warm-start model to a new L
#'
#' When the L_greedy path in `susie_workhorse` grows L between
#' rounds, the supplied `model_init` from round `k` has `L_prev`
#' effects but round `k + 1` requests `L_new >= L_prev`. This
#' helper appends zero-state effects to `alpha`, `mu`, `mu2`,
#' `V`, `KL`, and `lbf` so the warm-loaded fields have the
#' requested L. Mirrors `prune_single_effects` from susieR but
#' for the mfsusieR list-of-list `mu` / `mu2` shape.
#'
#' Errors when `L_prev > L_new` (no shrinking; the L_greedy
#' loop only grows L).
#'
#' @keywords internal
#' @noRd
expand_model_init_to_L <- function(mi, L_new, p, M, T_basis) {
  if (is.null(mi$alpha)) {
    stop("model_init must carry an `alpha` matrix.")
  }
  L_prev <- nrow(mi$alpha)
  if (L_prev > L_new) {
    stop(sprintf(
      "model_init has L = %d but the call requested L = %d; shrinking is not supported.",
      L_prev, L_new))
  }
  if (L_prev == L_new) return(mi)

  L_diff <- L_new - L_prev

  mi$alpha <- rbind(mi$alpha,
                    matrix(1 / p, nrow = L_diff, ncol = p))

  if (!is.null(mi$mu)) {
    extra_mu <- vector("list", L_diff)
    for (l in seq_len(L_diff)) {
      extra_mu[[l]] <- vector("list", M)
      for (m in seq_len(M)) {
        extra_mu[[l]][[m]] <- matrix(0, nrow = p, ncol = T_basis[m])
      }
    }
    mi$mu <- c(mi$mu, extra_mu)
  }
  if (!is.null(mi$mu2)) {
    extra_mu2 <- vector("list", L_diff)
    for (l in seq_len(L_diff)) {
      extra_mu2[[l]] <- vector("list", M)
      for (m in seq_len(M)) {
        extra_mu2[[l]][[m]] <- matrix(0, nrow = p, ncol = T_basis[m])
      }
    }
    mi$mu2 <- c(mi$mu2, extra_mu2)
  }
  if (!is.null(mi$V)) mi$V <- c(mi$V, rep(1, L_diff))
  if (!is.null(mi$KL))  mi$KL  <- c(mi$KL,  rep(NA_real_, L_diff))
  if (!is.null(mi$lbf)) mi$lbf <- c(mi$lbf, rep(NA_real_, L_diff))

  # Extend per-effect prior storage. When growing L, append L_diff
  # cold-start copies of effect 1's prior state. Detects (and
  # upgrades) the legacy non-per-effect pi_V shape from older
  # warm-starts.
  if (!is.null(mi$pi_V) && !is.null(mi$pi_V[[1L]]) &&
      is.matrix(mi$pi_V[[1L]])) {
    # Legacy shape: list[M] of matrix. Wrap into list[L_prev] of that.
    mi$pi_V <- lapply(seq_len(L_prev), function(.) mi$pi_V)
  }
  if (!is.null(mi$pi_V)) {
    seed_pi <- mi$pi_V[[1L]]
    mi$pi_V <- c(mi$pi_V, lapply(seq_len(L_diff), function(.) seed_pi))
  }
  if (!is.null(mi$fitted_g_per_effect)) {
    seed_fge <- mi$fitted_g_per_effect[[1L]]
    mi$fitted_g_per_effect <- c(mi$fitted_g_per_effect,
      lapply(seq_len(L_diff), function(.) seed_fge))
  }
  mi
}

# ---- ibss_initialize -------------------------------------------

#' IBSS loop initializer for `mf_individual`
#'
#' Mirrors `ibss_initialize.default` but skips its
#' scalar-only validation of `params$residual_variance` (mfsusieR
#' uses a list-of-vectors per-(scale, outcome) shape) and routes
#' through the mfsusieR S3 methods for `initialize_susie_model`,
#' `initialize_fitted`. When `params$model_init` is non-NULL, the
#' working model is seeded from the supplied fit's `alpha`, `mu`,
#' `mu2`, `KL`, `lbf`, `V`, `pi_V`, `sigma2`, `fitted`, and
#' `intercept` so the IBSS loop resumes from the prior posterior
#' state rather than the zero state.
#'
#' @keywords internal
#' @noRd
ibss_initialize.mf_individual <- function(data, params) {
  var_y <- get_var_y(data)

  if (data$p < params$L) params$L <- data$p

  if (is.null(params$residual_variance)) {
    params$residual_variance <- var_y
  }

  mat_init <- initialize_susie_model(data, params, var_y)

  if (!is.null(params$model_init)) {
    mi <- expand_model_init_to_L(params$model_init, params$L,
                                 data$p, data$M, data$T_basis)
    warm_fields <- c("alpha", "mu", "mu2", "V", "V_grid",
                     "pi_V", "fitted_g_per_effect",
                     "G_prior", "sigma2", "fitted",
                     "intercept")
    for (f in warm_fields) {
      if (!is.null(mi[[f]])) mat_init[[f]] <- mi[[f]]
    }
    mat_init$KL  <- rep(NA_real_, params$L)
    mat_init$lbf <- rep(NA_real_, params$L)
  }

  fitted     <- initialize_fitted(data, mat_init)

  model_class <- class(mat_init)
  model       <- c(mat_init, list(null_index = 0), fitted)
  class(model) <- if (inherits(mat_init, "susie")) model_class else "susie"
  model$converged <- FALSE
  # Populate the per-iter precompute cache so the first SER
  # iteration has cached `shat2` (and, on mixsqp paths, `sdmat`
  # / `log_sdmat`). Re-populated by
  # `update_variance_components.mf_individual` after each sigma2
  # update.
  if (inherits(data, "mf_individual")) {
    model <- refresh_iter_cache.mf_individual(data, model)
  }
  model
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

# ---- trim_null_effects ----------------------------------------

#' Drop effects with negligible effective slab variance
#'
#' Zeros out `alpha[l, ]`, `mu[[l]][[m]]`, `mu2[[l]][[m]]`,
#' `lbf[l]`, `KL[l]`, `lbf_variable[l, ]`, and
#' `lbf_variable_outcome[l, , ]` for every effect `l` with
#' `model$V[l] < params$prior_tol`. Mirrors
#' `susieR::trim_null_effects.default` over mfsusieR's
#' list-of-list `mu` / `mu2` shape so that effects dropped by
#' `susie_get_pip`'s `V > prior_tol` filter are also zeroed
#' on the rest of the fit (no double-counting).
#'
#' @keywords internal
#' @noRd
trim_null_effects.mf_individual <- function(data, params, model) {
  null_idx <- which(model$V < (params$prior_tol %||% 1e-9))
  if (length(null_idx) == 0L) return(model)

  model$V[null_idx]              <- 0
  model$alpha[null_idx, ]        <- rep(model$pi, each = length(null_idx))
  model$lbf[null_idx]            <- 0
  model$KL[null_idx]             <- 0
  model$lbf_variable[null_idx, ] <- 0
  if (!is.null(model$lbf_variable_outcome)) {
    model$lbf_variable_outcome[null_idx, , ] <- 0
  }
  for (l in null_idx) {
    for (m in seq_len(data$M)) {
      model$mu[[l]][[m]][]  <- 0
      model$mu2[[l]][[m]][] <- 0
    }
  }
  model
}

# Note: susieR's `ibss_finalize` is a regular function, not an S3
# generic. It calls the per-class generics (`get_fitted`, `get_cs`,
# `get_intercept`, `get_variable_names`, `get_zscore`,
# `cleanup_model`) via dispatch, so mfsusieR overrides those rather
# than trying to register an `ibss_finalize.mf_individual` method.
# Per-outcome residual attachment happens in the public `mfsusie()`
# wrapper after `susie_workhorse` returns.

# ---- User-facing post-fit accessors ---------------------------

#' Per-outcome column scale factors of X
#' @keywords internal
#' @noRd
get_scale_factors.mf_individual <- function(data, params, ...) {
  data$csd
}

#' Per-outcome intercepts (length `M`)
#'
#' Returns a list of length `M`; each entry is a length-`T_basis[m]`
#' zero vector. mfsusieR centers Y per-outcome, per-column at
#' construction, so the wavelet-domain intercept is identically
#' zero. The on-original-Y-scale reconstruction (per
#' `predict.mfsusie`) folds the per-column center back in via
#' `dwt_meta$column_center[[m]]` rather than through this method.
#'
#' @keywords internal
#' @noRd
get_intercept.mf_individual <- function(data, params, model, ...) {
  lapply(data$T_basis, function(T_m) rep(0, T_m))
}

#' Per-outcome fitted values (length-M list of n x T_basis matrices)
#' @keywords internal
#' @noRd
get_fitted.mf_individual <- function(data, params, model, ...) {
  model$fitted
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
#'
#' Override the `.individual` default (which would crash on the
#' list-shaped `data$D`). mfsusieR does not surface per-variable
#' z-scores; use `fit$lbf_variable` for per-variable evidence.
#'
#' @keywords internal
#' @noRd
get_zscore.mf_individual <- function(data, params, model, ...) {
  NULL
}
