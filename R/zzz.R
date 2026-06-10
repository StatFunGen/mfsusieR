# Package load hook + namespace plumbing. Registers S3 methods on
# susieR's per-iteration generics so that `ibss_fit` dispatches to
# mfsusieR methods when the data class inherits from `mf_individual`.
# Caches susieR internals as package-level bindings so methods call
# them bare. No numerics here.

# Null-coalescing operator fallback for R < 4.4.0.
if (!exists("%||%", baseenv())) {
  `%||%` <- function(x, y) if (is.null(x)) y else x
}

# Cached susieR + ebnm internals (populated by .onLoad).
warning_message            <- NULL
SER_posterior_e_loglik     <- NULL
update_variance_components <- NULL
update_derived_quantities  <- NULL
get_var_y                  <- NULL
initialize_susie_model     <- NULL
initialize_fitted          <- NULL
get_cs                     <- NULL
make_track_snapshot        <- NULL
vloglik_point_laplace      <- NULL

#' @keywords internal
#' @noRd
.onLoad <- function(libname, pkgname) {
  susie_ns <- asNamespace("susieR")
  ebnm_ns  <- asNamespace("ebnm")
  pkg_ns   <- asNamespace(pkgname)

  # Cache susieR internals used by per-iteration methods.
  for (fn in c("warning_message",
               "SER_posterior_e_loglik",
               "update_variance_components",
               "update_derived_quantities",
               "get_var_y",
               "initialize_susie_model",
               "initialize_fitted",
               "get_cs",
               "make_track_snapshot")) {
    assign(fn, get(fn, envir = susie_ns), envir = pkg_ns)
  }
  assign("vloglik_point_laplace",
         get("vloglik_point_laplace", envir = ebnm_ns),
         envir = pkg_ns)

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
    "track_ibss_fit",
    "trim_null_effects",
    # one-shot init / finalize
    "get_var_y",
    "initialize_susie_model",
    "initialize_fitted",
    "ibss_initialize",
    "cleanup_extra_fields",
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
    "get_posterior_moments_l",
    # Override susieR's `get_objective.default` to first refresh per-
    # effect lbf[l] / KL[l] against the iter-final pi_V before the
    # ELBO read. Without the refresh, the per-iteration ELBO is a
    # hybrid quantity (Eloglik at iter-final state, KL[l] at state-
    # when-l-was-updated). See `refresh_lbf_kl.mf_individual`.
    "get_objective",
    # Verbose-row formatting hook used by check_convergence.default
    # so the per-iteration line shows the list-of-vectors sigma2 in a
    # compact, fixed-width form rather than the default scalar
    # `sprintf("%.4f", model$sigma2)` (which would crash here).
    "format_sigma2_summary",
    # Extra diagnostic appended to the verbose row: per-iteration
    # mixture-null-mass summary across (outcome, scale).
    "format_extra_diag"
  )
  for (g in mfsusie_generics) {
    method_fn <- get(paste0(g, ".mfsusie"), envir = pkg_ns)
    registerS3method(g, "mfsusie", method_fn, envir = susie_ns)
  }
}
