# Mixture-weight EM helpers for the per-(outcome, scale) prior update.
#
# Two pure-R helpers that build the (p * |idx_s| + 1) x K likelihood
# matrix consumed by mixsqp and run mixsqp at one (outcome, scale)
# pair. The orchestration (looping over outcomes and scales,
# reading per-effect zeta from `model$alpha[l, ]`, etc.) happens in
# `update_model_variance.mf_individual` (R/individual_data_methods.R).
#
# Manuscript reference: methods/derivation.tex line 216 (factorized
# empirical-Bayes mixture-weight update).

#' Build the mixsqp likelihood matrix at one (outcome, scale)
#'
#' Returns a `(p * |idx_s| + 1)` x `K` matrix when
#' `is_ebmvfr = FALSE`, or a `(p * |idx_s|)` x `K` matrix when
#' `is_ebmvfr = TRUE`. The penultimate `(p * |idx_s|)` rows are
#' `dnorm(Bhat_jt; 0, sqrt(sd_k^2 + Shat_jt^2)) / sd`,
#' one per `(j, t)` pair flattened in column-major order. When
#' `is_ebmvfr = FALSE` the leading row is a null-component
#' penalty pseudo-observation (`c(100, 0, ..., 0)`) that nudges
#' the mixsqp solution toward the null component when data is
#' uninformative; this matches the SuSiE / mfsusie path.
#' `is_ebmvfr = TRUE` omits the penalty row; the EB
#' wavelet-domain regression in `mf_adjust_for_covariates()`
#' uses this branch.
#'
#' @param bhat_slice p x |idx_s| numeric matrix.
#' @param shat_slice p x |idx_s| numeric matrix.
#' @param sd_grid length-K numeric vector of mixture-component standard
#'   deviations (the null component has `sd = 0`).
#' @param is_ebmvfr logical, drop the penalty row when TRUE.
#'   Default FALSE preserves the SuSiE / mfsusie path.
#' @return numeric matrix; `(p * |idx_s| + 1) x K` when
#'   `is_ebmvfr = FALSE`, `(p * |idx_s|) x K` otherwise.
#' @references
#' Manuscript: methods/derivation.tex line 216
#' (factorized empirical-Bayes mixture-weight update).
#' @importFrom stats dnorm median
#' @keywords internal
#' @noRd
mf_em_likelihood_per_scale <- function(bhat_slice, shat_slice, sd_grid,
                                       is_ebmvfr = FALSE,
                                       sdmat_cache = NULL,
                                       log_sdmat_cache = NULL) {
  bvec <- as.numeric(bhat_slice)
  K    <- length(sd_grid)

  if (!is.null(sdmat_cache) && !is.null(log_sdmat_cache)) {
    # Loop-invariant fast path: caller pre-built sdmat once per
    # (m, s) per IBSS iter and cached it on `model$em_cache`.
    sdmat     <- sdmat_cache
    log_sdmat <- log_sdmat_cache
  } else {
    svec <- as.numeric(shat_slice)
    sdmat     <- sqrt(outer(svec^2, sd_grid^2, "+"))
    log_sdmat <- log(sdmat)
  }
  log_L <- mf_em_log_likelihood_per_scale_cpp(bvec, sdmat, log_sdmat)

  # NA imputation for ill-conditioned Shat: replace each NA
  # column with the column median so mixsqp never sees an NA
  # row weight.
  # Common case has no NAs; skip the imputation entirely. When NAs
  # are present, replace each per-column with the column median in
  # a vectorised pass instead of iterating R-level over columns.
  nas <- is.na(log_L)
  if (any(nas)) {
    col_med <- apply(log_L, 2L, median, na.rm = TRUE)
    log_L[nas] <- col_med[col(log_L)[nas]]
  }
  L <- exp(log_L)

  if (is_ebmvfr) L else rbind(c(100, rep(0, K - 1)), L)
}

#' Run mixsqp at one (outcome, scale)
#'
#' For `is_ebmvfr = FALSE` (default; mfsusie path):
#' builds the weight vector
#' `w = (mixsqp_null_penalty * idx_size, zeta_repeated)`,
#' prepended with the null-component penalty pseudo-observation.
#'
#' For `is_ebmvfr = TRUE` (covariate-adjustment path): omits the
#' penalty prepend; the weight vector is just
#' `w = rep(1, p * idx_size)` (each (j, t) observation
#' counts equally) and the L matrix has no penalty row.
#'
#' Both branches collapse to exact null when more than
#' `1 - tol_null_prior` of the mass lands on the null
#' component.
#'
#' @param L mixsqp likelihood matrix from
#'   `mf_em_likelihood_per_scale` (with the matching
#'   `is_ebmvfr` flag).
#' @param zeta length-p variable-level posterior weights (the per-effect
#'   alpha for the SER step that triggered this update).
#'   Ignored when `is_ebmvfr = TRUE`.
#' @param idx_size integer, |idx_s|, positions in this scale.
#' @param mixsqp_null_penalty numeric, null-component penalty weight.
#' @param init_pi0_w numeric, starting null-component mass.
#' @param tol_null_prior numeric, threshold below which non-null
#'   mass collapses to zero.
#' @param control_mixsqp optional named list of mixsqp control args.
#' @param is_ebmvfr logical, run the covariate-adjustment-flavored
#'   M-step (no penalty prepend, uniform weights). Default FALSE.
#' @return length-K numeric vector of mixture proportions summing to 1.
#' @references
#' Manuscript: methods/derivation.tex line 216
#' (factorized empirical-Bayes mixture-weight update).
#' @keywords internal
#' @noRd
mf_em_m_step_per_scale <- function(L, zeta, idx_size,
                                   mixsqp_null_penalty = 0.1,
                                   init_pi0_w     = 0.5,
                                   tol_null_prior = 0.001,
                                   control_mixsqp = NULL,
                                   is_ebmvfr      = FALSE) {
  K <- ncol(L)
  w <- if (is_ebmvfr) {
    # Covariate-adjustment branch: zeta is uniform across the
    # K covariate-effect cells, repeated `idx_size` times. No
    # null-component penalty row.
    rep(zeta, idx_size)
  } else {
    c(mixsqp_null_penalty * idx_size, rep(zeta, idx_size))
  }
  ctrl <- if (is.null(control_mixsqp)) list(verbose = FALSE) else {
    if (is.null(control_mixsqp$verbose)) control_mixsqp$verbose <- FALSE
    control_mixsqp
  }
  out <- mixsqp(L, w,
                x0      = c(init_pi0_w, rep(1e-6, K - 1)),
                log     = FALSE,
                control = ctrl)$x
  if (out[1] > 1 - tol_null_prior) {
    out    <- numeric(K)
    out[1] <- 1
  }
  out
}
