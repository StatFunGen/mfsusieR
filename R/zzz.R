# Package load hook + namespace plumbing.
#
# Registers S3 methods on susieR's internal generics so that
# `susieR::ibss_fit` dispatches to the mfsusieR per-iteration methods
# when the data class inherits from `mf_individual`. Mirrors the
# pattern in `mvsusieR/R/zzz.R` (registering on `mv_individual`).
#
# Caches susieR internal helpers (`lbf_stabilization`,
# `compute_posterior_weights`, `warning_message`) as package-level
# bindings so the methods can call them without `susieR:::`.
#
# Numerical routines do not live here.

# Null-coalescing operator fallback for R < 4.4.0 (mvsusieR pattern).
if (!exists("%||%", baseenv())) {
  `%||%` <- function(x, y) if (is.null(x)) y else x
}

# Cached susieR internals - populated by .onLoad()
warning_message            <- NULL
# Cached susieR S3 generics so the per-effect / per-iteration methods
# can call the generic (and pick up future subclass overrides) without
# going through `susieR:::`.
SER_posterior_e_loglik     <- NULL
update_variance_components <- NULL
update_derived_quantities  <- NULL
get_var_y                  <- NULL
initialize_susie_model     <- NULL
initialize_fitted          <- NULL
get_cs                     <- NULL

#' @useDynLib mfsusieR, .registration = TRUE
#' @keywords internal
#' @noRd
.onLoad <- function(libname, pkgname) {
  susie_ns <- asNamespace("susieR")
  pkg_ns   <- asNamespace(pkgname)

  # Cache susieR internals used by `loglik.mf_individual` and
  # forthcoming PR-7 finalize methods.
  for (fn in c("warning_message",
               "SER_posterior_e_loglik",
               "update_variance_components",
               "update_derived_quantities",
               "get_var_y",
               "initialize_susie_model",
               "initialize_fitted",
               "get_cs")) {
    assign(fn, get(fn, envir = susie_ns), envir = pkg_ns)
  }

  # Register S3 methods for the susieR generics dispatched on
  # `mf_individual`. New per-iteration methods land here as PR groups
  # 6, 6b, 7 add them; the list grows.
  mf_generics <- c(
    # per-effect SER step
    "compute_residuals",
    "compute_ser_statistics",
    "optimize_prior_variance",
    "loglik",
    "neg_loglik",
    "calculate_posterior_moments",
    "compute_kl",
    "SER_posterior_e_loglik",
    "update_fitted_values",
    # per-iteration
    "update_variance_components",
    "update_derived_quantities",
    "update_model_variance",
    "Eloglik",
    "get_objective",
    "track_ibss_fit",
    "check_convergence",
    "validate_prior",
    "trim_null_effects",
    # one-shot init / finalize
    "configure_data",
    "get_var_y",
    "initialize_susie_model",
    "initialize_fitted",
    "ibss_initialize",
    "cleanup_model",
    # post-fit accessors
    "get_scale_factors",
    "get_intercept",
    "get_fitted",
    "get_cs",
    "get_variable_names",
    "get_zscore"
  )
  for (g in mf_generics) {
    method_fn <- get(paste0(g, ".mf_individual"), envir = pkg_ns)
    registerS3method(g, "mf_individual", method_fn, envir = susie_ns)
  }

  # Register S3 methods for the `mfsusie` model class (per-effect
  # accessors mirroring mvsusieR's mvsusie methods).
  mfsusie_generics <- c(
    "get_alpha_l",
    "get_posterior_mean_l",
    "get_posterior_mean_sum",
    "get_posterior_moments_l",
    "get_prior_variance_l",
    "set_prior_variance_l"
  )
  for (g in mfsusie_generics) {
    method_fn <- get(paste0(g, ".mfsusie"), envir = pkg_ns)
    registerS3method(g, "mfsusie", method_fn, envir = susie_ns)
  }
}
