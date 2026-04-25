# `mfsusie` model class.
#
# The model object returned by `mfsusie()` and threaded through
# `susieR::susie_workhorse` via the per-iteration S3 dispatch. The
# initialize method `initialize_susie_model.mf_individual` is the
# constructor; `susieR::ibss_fit` calls it once at the start of the
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
#'   `L` (effect-count upper bound), `prior_weights` (length-J
#'   variable-selection prior, defaults to uniform `1/J` if NULL),
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
  J <- data$J
  M <- data$M
  T_padded <- data$T_padded

  # ---- Per-effect posterior moments (nested list per design.md D2) ----
  # mu[[l]][[m]] and mu2[[l]][[m]] are J x T_padded[m] matrices.
  # The nested-list shape lets ragged T_m coexist without padding to
  # a common T across modalities.
  mu  <- vector("list", L)
  mu2 <- vector("list", L)
  for (l in seq_len(L)) {
    mu[[l]]  <- vector("list", M)
    mu2[[l]] <- vector("list", M)
    for (m in seq_len(M)) {
      mu[[l]][[m]]  <- matrix(0, nrow = J, ncol = T_padded[m])
      mu2[[l]][[m]] <- matrix(0, nrow = J, ncol = T_padded[m])
    }
  }

  # ---- Prior structure (per design.md D6) ----
  # V_grid: list[M] of length-K vectors (or list[M] of S_m x K
  #   matrices when prior_variance_scope = "per_scale_modality").
  # pi_V:   list[M] of S_m x K mixture-weight matrices.
  # null_prior_weight: scalar (default 2 per design.md D5/D7).
  prior <- params$prior
  V_grid             <- if (is.null(prior)) NULL else prior$V_grid
  pi_V               <- if (is.null(prior)) NULL else prior$pi
  null_prior_weight  <- if (is.null(prior)) 0   else prior$null_prior_weight

  # ---- Residual variance (per design.md D2) ----
  # `sigma2` is either a list[M] of scalars (legacy
  # `shared_per_modality` mode) or a list[M] of length-S_m vectors
  # (per-(scale, modality) default). Falls back to `var_y` if
  # `params$residual_variance` is NULL.
  sigma2 <- params$residual_variance
  if (is.null(sigma2)) {
    sigma2 <- var_y
  }

  # ---- Per-effect scalar scaling of the prior ----
  # V[l] scales the entire mixture grid for effect l; mvsusieR uses
  # the same pattern. Initial value 1 (no scaling); the IBSS loop
  # updates V[l] via `update_model_variance.mf_individual`.
  V <- rep(1, L)

  # ---- Variable-selection prior over the J SNPs ----
  prior_weights <- params$prior_weights
  if (is.null(prior_weights)) {
    prior_weights <- rep(1 / J, J)
  }

  # ---- Per-modality fitted values (running sum of effect curves) ----
  fitted_values <- vector("list", M)
  for (m in seq_len(M)) {
    fitted_values[[m]] <- matrix(0, nrow = data$n, ncol = T_padded[m])
  }

  model <- list(
    alpha             = matrix(1 / J, nrow = L, ncol = J),
    mu                = mu,
    mu2               = mu2,
    V                 = V,
    V_grid            = V_grid,
    pi_V              = pi_V,
    null_prior_weight = null_prior_weight,
    KL                = rep(NA_real_, L),
    lbf               = rep(NA_real_, L),
    lbf_variable      = matrix(NA_real_, nrow = L, ncol = J),
    sigma2            = sigma2,
    pi                = prior_weights,
    fitted            = fitted_values,
    intercept         = rep(0, M),
    L                 = L,
    J                 = J,
    M                 = M
  )
  class(model) <- c("mfsusie", "susie")
  model
}
