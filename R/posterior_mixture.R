# Per-(SNP, scale) mixture-of-normals helpers, R-wrapper layer.
#
# The kernels themselves live in `src/posterior_mixture.cpp` (cpp11);
# this file is the thin wrapper that input-validates and forwards.
# Pure-R oracles for the same kernels live in
# `R/reference_implementations.R` and are tested at <= 1e-12 against
# the cpp11 output.
#
# Per design.md D14 (Phase 3 cpp11 stance).
#
# Manuscript references:
#   methods/derivation.tex eq:post_f_mix
#   methods/derivation.tex eq:post_f2_mix

#' Per-(SNP, scale) log-Bayes factor under a mixture-of-normals prior
#'
#' Forwards to the cpp11 kernel
#' `mixture_log_bf_per_scale_cpp` in `src/posterior_mixture.cpp`. The
#' R reference oracle `mixture_log_bf_per_scale_R` in
#' `R/reference_implementations.R` matches this at `<= 1e-12`.
#'
#' @param bhat_slice p x |idx_s| matrix of marginal effect estimates
#'   for one modality, one scale.
#' @param shat_slice p x |idx_s| matrix of standard errors.
#' @param sd_grid length-K vector of mixture-component standard
#'   deviations (`sd_k = 0` for the null component).
#' @param pi_grid length-K vector of mixture-component weights.
#' @param V_scale numeric scalar, per-effect prior-variance multiplier
#'   (so component variance is `V_scale * sd_k^2`).
#' @return numeric vector of length p.
#' @keywords internal
#' @noRd
mixture_log_bf_per_scale <- function(bhat_slice,
                                     shat_slice,
                                     sd_grid,
                                     pi_grid,
                                     V_scale = 1) {
  if (!is.matrix(bhat_slice) || !is.matrix(shat_slice)) {
    stop("`bhat_slice` and `shat_slice` must be matrices.")
  }
  mixture_log_bf_per_scale_cpp(bhat_slice, shat_slice,
                               as.numeric(sd_grid),
                               as.numeric(pi_grid),
                               as.numeric(V_scale))
}

#' Per-(SNP, position) posterior moments under a mixture-of-normals prior
#'
#' Forwards to the cpp11 kernel
#' `mixture_posterior_per_scale_cpp`. Returns the p x |idx_s|
#' posterior mean and second moment per (SNP, position).
#'
#' @inheritParams mixture_log_bf_per_scale
#' @return list with elements `pmean` and `pmean2`, each a p x |idx_s|
#'   matrix.
#' @keywords internal
#' @noRd
mixture_posterior_per_scale <- function(bhat_slice,
                                        shat_slice,
                                        sd_grid,
                                        pi_grid,
                                        V_scale = 1) {
  if (!is.matrix(bhat_slice) || !is.matrix(shat_slice)) {
    stop("`bhat_slice` and `shat_slice` must be matrices.")
  }
  mixture_posterior_per_scale_cpp(bhat_slice, shat_slice,
                                  as.numeric(sd_grid),
                                  as.numeric(pi_grid),
                                  as.numeric(V_scale))
}
