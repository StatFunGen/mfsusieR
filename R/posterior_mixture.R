# Per-(variable, scale) mixture-of-normals helpers, R-wrapper layer.
#
# The kernels themselves live in `src/posterior_mixture.cpp` (cpp11);
# this file is the thin wrapper that input-validates and forwards.
# Pure-R oracles for the same kernels live in
# `R/reference_implementations.R` and are tested at <= 1e-12 against
# the cpp11 output.
#
# Compiled cpp11 kernels for the per-effect SER hot path.
#
# Manuscript references:
#   methods/derivation.tex eq:post_f_mix
#   methods/derivation.tex eq:post_f2_mix

#' Per-(variable, scale) log-Bayes factor under a mixture-of-normals prior
#'
#' Forwards to the cpp11 kernel
#' `mixture_log_bf_per_scale_cpp` in `src/posterior_mixture.cpp`. The
#' R reference oracle `mixture_log_bf_per_scale_R` in
#' `R/reference_implementations.R` matches this at `<= 1e-12`.
#'
#' @param bhat_slice p x |idx_s| matrix of marginal effect estimates
#'   for one outcome, one scale.
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

#' Per-(variable, scale) log-Bayes factor under a mixture-of-Student-t prior
#'
#' R-only kernel that mirrors `mixture_log_bf_per_scale` but
#' replaces each Normal density with a scaled Student-t (location
#' 0, scale `sqrt(V_scale * sd_k^2 + Shat^2)`, degrees of freedom
#' `df`). Selected by `small_sample_correction = TRUE` on
#' `mfsusie()` / `fsusie()` when the per-outcome sample size is
#' small enough that the Wakefield Normal approximation
#' under-propagates residual-variance uncertainty into the
#' per-variable Bayes factor. Posterior moments
#' (`mixture_posterior_per_scale`) are unchanged: the correction
#' acts on variable selection probabilities, not on the
#' conditional posterior given inclusion.
#'
#' @inheritParams mixture_log_bf_per_scale
#' @param df positive integer, degrees of freedom for the
#'   Student-t marginal. Conventionally `n - 1` per outcome.
#' @return numeric vector of length p.
#' @keywords internal
#' @noRd
mixture_log_bf_per_scale_johnson <- function(bhat_slice,
                                             shat_slice,
                                             sd_grid,
                                             pi_grid,
                                             V_scale = 1,
                                             df) {
  if (!is.matrix(bhat_slice) || !is.matrix(shat_slice)) {
    stop("`bhat_slice` and `shat_slice` must be matrices.")
  }
  if (missing(df) || !is.numeric(df) || length(df) != 1L ||
      !is.finite(df) || df <= 0) {
    stop("`df` must be a positive scalar.")
  }
  shat_safe <- shat_slice
  shat_safe[shat_safe <= 0] <- 1e-32
  dims <- dim(bhat_slice)

  K       <- length(sd_grid)
  pi_grid <- as.numeric(pi_grid)
  sd_grid <- as.numeric(sd_grid)
  V_scale <- as.numeric(V_scale)

  # `LaplacesDemon::dstp` does not preserve matrix shape; flatten,
  # compute, reshape.
  bhat_vec  <- as.vector(bhat_slice)
  shat2_vec <- as.vector(shat_safe^2)

  # Per-component density mixture under the alternative.
  tt <- numeric(length(bhat_vec))
  for (k in seq_len(K)) {
    sd_k_var <- V_scale * sd_grid[k]^2
    tt       <- tt + pi_grid[k] *
      LaplacesDemon::dstp(bhat_vec,
                          tau = 1 / (sd_k_var + shat2_vec),
                          nu  = df)
  }
  log_dens_null <- LaplacesDemon::dstp(bhat_vec,
                                       tau = 1 / shat2_vec,
                                       nu  = df, log = TRUE)
  per_cell <- matrix(log(tt) - log_dens_null,
                     nrow = dims[1L], ncol = dims[2L])

  # Match the upstream kernel: per-(variable) summed across the
  # |idx_s| columns, returning a length-p vector.
  rowSums(per_cell)
}

#' Per-(variable, position) posterior moments under a mixture-of-normals prior
#'
#' Forwards to the cpp11 kernel
#' `mixture_posterior_per_scale_cpp`. Returns the p x |idx_s|
#' posterior mean and second moment per (variable, position).
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
