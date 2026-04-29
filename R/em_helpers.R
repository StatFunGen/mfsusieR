# Mixture-weight EM helpers for the per-(outcome, scale) prior update.
#
# Two pure-R helpers, called by the orchestrator in
# `optimize_prior_variance.mf_individual` (R/individual_data_methods.R):
#   * `mf_em_likelihood_per_scale`  -- build the mixsqp L matrix
#   * `mf_em_m_step_per_scale`      -- run mixsqp (warm-started)
#
# Manuscript reference: methods/derivation.tex line 216 (factorized
# empirical-Bayes mixture-weight update).

# All-null pi vector (1, 0, ..., 0); used by mixsqp's `tol_null_prior`
# collapse path.
mf_em_pi_null <- function(K) {
  out <- numeric(K)
  out[1L] <- 1
  out
}

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
    # (m, s) per IBSS iter and cached it on `model$iter_cache`.
    sdmat     <- sdmat_cache
    log_sdmat <- log_sdmat_cache
  } else {
    svec <- as.numeric(shat_slice)
    sdmat     <- sqrt(outer(svec^2, sd_grid^2, "+"))
    log_sdmat <- log(sdmat)
  }
  # log_L[r, k] = log dnorm(bvec[r]; 0, sdmat[r, k])
  #             = -0.5 log(2 pi) - log_sdmat[r, k] - 0.5 (bvec[r] / sdmat[r, k])^2
  # Vectorized via broadcast: bvec/sdmat is (N x K), log_sdmat is (N x K).
  log_L <- -0.5 * log(2 * pi) - log_sdmat - 0.5 * (bvec / sdmat)^2

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
#' `w = (mixture_null_weight * idx_size, zeta_repeated)`,
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
#' @param mixture_null_weight numeric, null-component penalty weight.
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
                                   mixture_null_weight = 0.05,
                                   init_pi0_w     = 0.5,
                                   tol_null_prior = 0.001,
                                   control_mixsqp = NULL,
                                   is_ebmvfr      = FALSE,
                                   pi_warm_start  = NULL) {
  K <- ncol(L)
  w <- if (is_ebmvfr) {
    # Covariate-adjustment branch: zeta is uniform across the
    # K covariate-effect cells, repeated `idx_size` times. No
    # null-component penalty row.
    rep(zeta, idx_size)
  } else {
    c(mixture_null_weight * idx_size, rep(zeta, idx_size))
  }

  # mixsqp defaults tuned for warm-started inner-loop use. tol.svd = 0
  # skips an irlba SVD that only helps low-rank L (ours is full rank
  # by construction). Caller-supplied `control_mixsqp` (from
  # `params$control_mixsqp`) overrides any of these.
  default_ctrl <- list(verbose      = FALSE,
                       tol.svd      = 0,
                       numiter.em   = 10L,
                       convtol.sqp  = 1e-6)
  ctrl <- if (is.null(control_mixsqp)) default_ctrl else {
    out_ctrl <- default_ctrl
    out_ctrl[names(control_mixsqp)] <- control_mixsqp
    out_ctrl
  }

  cold_x0 <- function(K_eff) c(init_pi0_w, rep(1e-6, K_eff - 1L))
  warm_x0 <- function(pi_prev, valid_idx, K_eff) {
    if (is.null(pi_prev) || length(pi_prev) != K) return(cold_x0(K_eff))
    # Floor near-zero entries so mixsqp's interior-point start is well-defined.
    x <- pmax(pi_prev[valid_idx], 1e-6)
    x / sum(x)
  }

  # Filter all-zero columns of L before handing to mixsqp. Each such
  # column corresponds to a mixture-component standard deviation that
  # has zero data support at this (outcome, scale, l) -- e.g., a
  # very fine-scale sd_k against an alpha-thinned variant subset
  # whose density at that scale is below numerical precision. mixsqp
  # is correct on these columns (it sets their weight to 0 and emits
  # a warning), but the warning leaks to user output once per
  # `optimize_prior_variance` call -- which fires L * S_m times per
  # IBSS iteration. We drop the all-zero columns up front, run
  # mixsqp on the survivors, and pad the solution back with zeros
  # in the dropped slots. Numerically identical to mixsqp's own
  # all-zeros handling, but quiet.
  zero_col <- colSums(abs(L)) == 0
  if (any(zero_col)) {
    keep_cols <- which(!zero_col)
    L_keep    <- L[, keep_cols, drop = FALSE]
    K_keep    <- length(keep_cols)
    if (K_keep == 0L) return(mf_em_pi_null(K))
    x0_keep  <- warm_x0(pi_warm_start, keep_cols, K_keep)
    res      <- mixsqp(L_keep, w, x0 = x0_keep, log = FALSE,
                       control = ctrl)
    out <- numeric(K)
    out[keep_cols] <- res$x
  } else {
    res <- mixsqp(L, w,
                  x0      = warm_x0(pi_warm_start, seq_len(K), K),
                  log     = FALSE,
                  control = ctrl)
    out <- res$x
  }
  status_ok <- is.null(res$status) ||
                identical(res$status, "converged to optimal solution") ||
                identical(res$status, "converged")
  if (!status_ok) {
    warning_message(sprintf(
      "mixsqp prior-update did not converge (status: \"%s\"). If this repeats, pass `control_mixsqp = list(numiter.em = 20, convtol.sqp = 1e-8, tol.svd = 1e-10)` to mfsusie() / fsusie() to use mixsqp's stock control.",
      as.character(res$status)),
      style = "hint")
  }
  if (out[1] > 1 - tol_null_prior) out <- mf_em_pi_null(K)
  out
}
