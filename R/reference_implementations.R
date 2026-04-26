# Pure-R reference implementations of functions with cpp11 counterparts.
# Used as the oracle for cpp11 unit tests (per design.md D14).
# Mirrors the mvsusieR pattern (`mvsusieR/R/reference_implementations.R`).
#
# These functions SHALL stay literal-equivalent to the cpp11 kernels:
# the unit tests in `tests/testthat/test_posterior_mixture.R` assert
# agreement at `<= 1e-12` on randomized inputs.
#
# The naming convention is `<kernel>_R` for the R oracle and
# `<kernel>` (no suffix) for the cpp11 kernel, matching mvsusieR.

# =============================================================================
# Mixture-of-normals per-(SNP, scale) helpers (oracle for src/posterior_mixture.cpp)
# =============================================================================

#' Per-(SNP, scale) log-Bayes factor: pure-R oracle
#'
#' Closed-form log-BF for the per-scale mixture-of-normals prior.
#' For SNP j and the scale's |idx_s| positions, computes
#' `sum_i logSumExp_k(log(pi_k) + log dnorm(Bhat[j,i]; 0, V*sd_k^2 + Shat[j,i]^2)
#'                                    - log dnorm(Bhat[j,i]; 0, Shat[j,i]^2))`.
#' The null component (`sd_k = 0`) contributes `log(pi_k)` per position.
#'
#' Used in tests as the `<= 1e-12` oracle for the cpp11 implementation
#' `mixture_log_bf_per_scale`.
#'
#' @param bhat_slice p x |idx_s| matrix.
#' @param shat_slice p x |idx_s| matrix.
#' @param sd_grid length-K vector of mixture-component standard deviations.
#' @param pi_grid length-K vector of mixture-component weights.
#' @param V_scale numeric, per-effect prior-variance multiplier.
#' @return numeric vector of length p.
#' @keywords internal
#' @noRd
mixture_log_bf_per_scale_R <- function(bhat_slice,
                                       shat_slice,
                                       sd_grid,
                                       pi_grid,
                                       V_scale = 1) {
  K <- length(sd_grid)
  log_pi <- log(pi_grid)
  Shat2 <- shat_slice^2
  log_dens_null <- -0.5 * log(2 * pi * Shat2) - 0.5 * bhat_slice^2 / Shat2

  per_k <- vector("list", K)
  for (k in seq_len(K)) {
    sd_k_var <- sd_grid[k]^2 * V_scale
    per_k[[k]] <- if (sd_k_var == 0) {
      matrix(log_pi[k], nrow = nrow(bhat_slice), ncol = ncol(bhat_slice))
    } else {
      var_alt <- sd_k_var + Shat2
      log_dens_alt <- -0.5 * log(2 * pi * var_alt) - 0.5 * bhat_slice^2 / var_alt
      log_pi[k] + log_dens_alt - log_dens_null
    }
  }

  m <- per_k[[1]]
  if (K > 1) for (k in 2:K) m <- pmax(m, per_k[[k]])
  s <- exp(per_k[[1]] - m)
  if (K > 1) for (k in 2:K) s <- s + exp(per_k[[k]] - m)
  rowSums(m + log(s))
}

#' Per-(SNP, position) posterior moments: pure-R oracle
#'
#' Closed-form posterior mean and second moment under the per-scale
#' mixture-of-normals prior, treating each position as conditionally
#' independent given the per-component variance.
#'
#' Per-component closed form (per `(j, i)`):
#' `var_alt = V*sd_k^2 + Shat^2`,
#' `shrink = V*sd_k^2 / var_alt`,
#' `post_mean_k = shrink * Bhat`,
#' `post_var_k = shrink * Shat^2`,
#' `weight_k = softmax_k( log(pi_k) + log dnorm_alt - log dnorm_null )`.
#' Mixture-collapsed: `E\[beta\] = sum_k weight_k * post_mean_k`,
#' `E\[beta^2\] = sum_k weight_k * (post_var_k + post_mean_k^2)`.
#'
#' Used in tests as the `<= 1e-12` oracle for the cpp11 implementation
#' `mixture_posterior_per_scale`.
#'
#' @inheritParams mixture_log_bf_per_scale_R
#' @return list with `pmean` and `pmean2`, both p x |idx_s| matrices.
#' @keywords internal
#' @noRd
mixture_posterior_per_scale_R <- function(bhat_slice,
                                          shat_slice,
                                          sd_grid,
                                          pi_grid,
                                          V_scale = 1) {
  K <- length(sd_grid)
  p <- nrow(bhat_slice)
  T_idx <- ncol(bhat_slice)
  log_pi <- log(pi_grid)
  Shat2 <- shat_slice^2
  log_dens_null <- -0.5 * log(2 * pi * Shat2) - 0.5 * bhat_slice^2 / Shat2

  per_k_log_w <- vector("list", K)
  per_k_pmean <- vector("list", K)
  per_k_pvar  <- vector("list", K)
  for (k in seq_len(K)) {
    sd_k_var <- sd_grid[k]^2 * V_scale
    if (sd_k_var == 0) {
      per_k_log_w[[k]] <- matrix(log_pi[k], nrow = p, ncol = T_idx)
      per_k_pmean[[k]] <- matrix(0, nrow = p, ncol = T_idx)
      per_k_pvar[[k]]  <- matrix(0, nrow = p, ncol = T_idx)
    } else {
      var_alt <- sd_k_var + Shat2
      log_dens_alt <- -0.5 * log(2 * pi * var_alt) - 0.5 * bhat_slice^2 / var_alt
      per_k_log_w[[k]] <- log_pi[k] + log_dens_alt - log_dens_null
      shrink <- sd_k_var / var_alt
      per_k_pmean[[k]] <- shrink * bhat_slice
      per_k_pvar[[k]]  <- shrink * Shat2
    }
  }

  m_max <- per_k_log_w[[1]]
  if (K > 1) for (k in 2:K) m_max <- pmax(m_max, per_k_log_w[[k]])
  e_sum <- exp(per_k_log_w[[1]] - m_max)
  if (K > 1) for (k in 2:K) e_sum <- e_sum + exp(per_k_log_w[[k]] - m_max)

  pmean  <- matrix(0, nrow = p, ncol = T_idx)
  pmean2 <- matrix(0, nrow = p, ncol = T_idx)
  for (k in seq_len(K)) {
    w_k <- exp(per_k_log_w[[k]] - m_max) / e_sum
    pmean  <- pmean  + w_k * per_k_pmean[[k]]
    pmean2 <- pmean2 + w_k * (per_k_pvar[[k]] + per_k_pmean[[k]]^2)
  }
  list(pmean = pmean, pmean2 = pmean2)
}
