# Package load hook + namespace plumbing. Registers S3 methods on
# susieR's per-iteration generics so that `ibss_fit` dispatches to
# mfsusieR methods when the data class inherits from `mf_individual`.
# Caches susieR internals as package-level bindings so methods call
# them bare. No numerics here.

# Null-coalescing operator fallback for R < 4.4.0.
if (!exists("%||%", baseenv())) {
  `%||%` <- function(x, y) if (is.null(x)) y else x
}

# Cached susieR internals (populated by .onLoad).
warning_message            <- NULL
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

  # Cache susieR internals used by per-iteration methods.
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
  # `mf_individual`.
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
    "trim_null_effects",
    # one-shot init / finalize
    "get_var_y",
    "initialize_susie_model",
    "initialize_fitted",
    "ibss_initialize",
    "cleanup_model",
    # post-fit accessors
    "get_scale_factors",
    "get_intercept",
    "get_fitted",
    "get_variable_names",
    "get_zscore"
  )
  for (g in mf_generics) {
    method_fn <- get(paste0(g, ".mf_individual"), envir = pkg_ns)
    registerS3method(g, "mf_individual", method_fn, envir = susie_ns)
  }

  # Register S3 accessors for the `mfsusie` model class.
  mfsusie_generics <- c(
    "get_posterior_mean_l",
    "get_posterior_mean_sum",
    "get_posterior_moments_l"
  )
  for (g in mfsusie_generics) {
    method_fn <- get(paste0(g, ".mfsusie"), envir = pkg_ns)
    registerS3method(g, "mfsusie", method_fn, envir = susie_ns)
  }
}
