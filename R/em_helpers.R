# Mixture-weight EM helpers for the per-(modality, scale) prior update.
#
# Two pure-R helpers that build the (p * |idx_s| + 1) x K likelihood
# matrix consumed by mixsqp and run mixsqp at one (modality, scale)
# pair. The orchestration (looping over modalities and scales,
# reading per-effect zeta from `model$alpha[l, ]`, etc.) happens in
# `update_model_variance.mf_individual` (R/individual_data_methods.R).
#
# Manuscript reference: methods/derivation.tex line 216 (factorized
# empirical-Bayes mixture-weight update).

#' Build the mixsqp likelihood matrix at one (modality, scale)
#'
#' Returns a `(p * |idx_s| + 1)` x `K` matrix where:
#' row 1 is the null-component penalty pseudo-observation
#' (`c(100, 0, ..., 0)`); rows 2..(p*|idx_s|+1) are
#' `dnorm(Bhat_jt; 0, sqrt(sd_k^2 + Shat_jt^2)) / sd`,
#' one per `(j, t)` pair flattened in column-major order.
#' The first row nudges the mixsqp solution toward the null
#' component when data is uninformative.
#'
#' @param bhat_slice p x |idx_s| numeric matrix.
#' @param shat_slice p x |idx_s| numeric matrix.
#' @param sd_grid length-K numeric vector of mixture-component standard
#'   deviations (the null component has `sd = 0`).
#' @return `(p * |idx_s| + 1)` x `K` numeric matrix.
#' @references
#' Manuscript: methods/derivation.tex line 216
#' (factorized empirical-Bayes mixture-weight update).
#' @importFrom stats dnorm median
#' @keywords internal
#' @noRd
mf_em_likelihood_per_scale <- function(bhat_slice, shat_slice, sd_grid) {
  bvec <- as.numeric(bhat_slice)
  svec <- as.numeric(shat_slice)
  K    <- length(sd_grid)

  sdmat <- sqrt(outer(svec^2, sd_grid^2, "+"))
  log_L <- dnorm(outer(bvec, rep(0, K), "-") / sdmat, log = TRUE) - log(sdmat)

  # NA imputation for ill-conditioned Shat (port-fidelity).
  log_L <- apply(log_L, 2, function(col) {
    na <- is.na(col)
    if (any(na)) col[na] <- median(col, na.rm = TRUE)
    col
  })
  L <- exp(log_L)

  rbind(c(100, rep(0, K - 1)), L)
}

#' Run mixsqp at one (modality, scale)
#'
#' Builds the weight vector
#' `w = (nullweight * idx_size, zeta_repeated)`, calls `mixsqp`,
#' and collapses to exact null when more than `1 - tol_null_prior` of
#' the mass lands on the null component.
#'
#' @param L mixsqp likelihood matrix from
#'   `mf_em_likelihood_per_scale`.
#' @param zeta length-p SNP-level posterior weights (the per-effect
#'   alpha for the SER step that triggered this update).
#' @param idx_size integer, |idx_s|, positions in this scale.
#' @param nullweight numeric, null-component penalty weight.
#' @param init_pi0_w numeric, starting null-component mass.
#' @param tol_null_prior numeric, threshold below which non-null
#'   mass collapses to zero.
#' @param control_mixsqp optional named list of mixsqp control args.
#' @return length-K numeric vector of mixture proportions summing to 1.
#' @references
#' Manuscript: methods/derivation.tex line 216
#' (factorized empirical-Bayes mixture-weight update).
#' @keywords internal
#' @noRd
mf_em_m_step_per_scale <- function(L, zeta, idx_size,
                                   nullweight     = 0.7,
                                   init_pi0_w     = 0.5,
                                   tol_null_prior = 0.001,
                                   control_mixsqp = NULL) {
  w <- c(nullweight * idx_size, rep(zeta, idx_size))
  K <- ncol(L)
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
