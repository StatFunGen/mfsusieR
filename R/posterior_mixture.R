# Per-(variable, scale) mixture-of-normals helpers (vectorized R).
#
# Closed-form mixture log-BF and posterior moments for the per-effect
# SER step. Pure vectorized R + BLAS; no cpp dependency. Numerical
# parity with the prior cpp implementation guaranteed at <= 1e-12 by
# `tests/testthat/test_posterior_mixture.R`.
#
# Manuscript references:
#   methods/derivation.tex eq:post_f_mix
#   methods/derivation.tex eq:post_f2_mix

#' Per-(variable, scale) log-Bayes factor under a mixture-of-normals prior
#'
#' Closed-form log Bayes factor for the per-scale
#' mixture-of-normals prior. For variable j and the scale's
#' `|idx_s|` positions, computes
#' `sum_i logSumExp_k( log(pi_k) +
#'                     log dnorm(Bhat[j,i]; 0, V*sd_k^2 + Shat[j,i]^2)
#'                     - log dnorm(Bhat[j,i]; 0, Shat[j,i]^2) )`.
#' The null component (`sd_k = 0`) contributes `log(pi_k)` per
#' position. Vectorized via per-component matrices and `pmax` /
#' `exp` / `rowSums` aggregation.
#'
#' @param bhat_slice p x |idx_s| matrix of marginal effect estimates
#'   for one outcome, one scale.
#' @param shat_slice p x |idx_s| matrix of standard errors.
#' @param sd_grid length-K vector of mixture-component standard
#'   deviations (`sd_k = 0` for the null component).
#' @param pi_grid length-K vector of mixture-component weights.
#' @param V_scale numeric scalar, per-effect prior-variance multiplier
#'   (component variance is `V_scale * sd_k^2`).
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
  K <- length(sd_grid)
  log_pi <- log(pi_grid)
  Shat2  <- shat_slice^2
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

#' Per-(variable, scale) log-Bayes factor under a mixture-of-Student-t prior
#'
#' Replaces each Normal density in `mixture_log_bf_per_scale` with a
#' scaled Student-t (location 0, scale `sqrt(V_scale * sd_k^2 + Shat^2)`,
#' degrees of freedom `df`). Selected by `small_sample_correction = TRUE`
#' on `mfsusie()` / `fsusie()` when the per-outcome sample size is
#' small enough that the Wakefield Normal approximation
#' under-propagates residual-variance uncertainty into the
#' per-variable Bayes factor. Posterior moments
#' (`mixture_posterior_per_scale`) are unchanged: the correction
#' acts on variable selection probabilities, not on the conditional
#' posterior given inclusion.
#'
#' @inheritParams mixture_log_bf_per_scale
#' @param df positive integer, degrees of freedom for the Student-t
#'   marginal. Conventionally `n - 1` per outcome.
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
  rowSums(per_cell)
}

#' Per-(variable, position) posterior moments under a mixture-of-normals prior
#'
#' Closed-form posterior mean and second moment under the per-scale
#' mixture-of-normals prior, treating each position as conditionally
#' independent given the per-component variance. Per-component closed
#' form (per `(j, i)`):
#' `var_alt = V*sd_k^2 + Shat^2`,
#' `shrink  = V*sd_k^2 / var_alt`,
#' `post_mean_k = shrink * Bhat`,
#' `post_var_k  = shrink * Shat^2`,
#' `weight_k = softmax_k( log(pi_k) + log dnorm_alt - log dnorm_null )`.
#' Mixture-collapsed: `E[beta] = sum_k weight_k * post_mean_k`,
#' `E[beta^2] = sum_k weight_k * (post_var_k + post_mean_k^2)`.
#'
#' @inheritParams mixture_log_bf_per_scale
#' @return list with elements `pmean` and `pmean2`, each a
#'   p x |idx_s| matrix.
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
  K <- length(sd_grid)
  p <- nrow(bhat_slice)
  T_idx <- ncol(bhat_slice)
  log_pi <- log(pi_grid)
  Shat2  <- shat_slice^2
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

#' Per-(variable, scale) log-Bayes factor under a per-scale point-Laplace prior
#'
#' Per (j, t) computes the marginal log-likelihood under the
#' point-Laplace prior via `ebnm:::vloglik_point_laplace()`,
#' subtracts the null marginal log-likelihood, and sums across
#' the scale's `|idx_s|` positions. The Normal-Laplace
#' convolution + spike mixture is ebnm's authoritative formula;
#' calling its internal vector helper avoids a math fork.
#'
#' @inheritParams mixture_log_bf_per_scale
#' @param fitted_g a `laplacemix` record returned by
#'   `ebnm::ebnm_point_laplace()`: a list with `pi`, `mean`, and
#'   `scale` of length 2, where `scale[1] == 0` (spike) and
#'   `scale[2]` is the Laplace slab scale `lambda`.
#' @return numeric vector of length p.
#' @keywords internal
#' @noRd
mixture_log_bf_laplace_per_scale <- function(bhat_slice,
                                              shat_slice,
                                              fitted_g,
                                              V_scale = 1) {
  if (!is.matrix(bhat_slice) || !is.matrix(shat_slice)) {
    stop("`bhat_slice` and `shat_slice` must be matrices.")
  }
  # V_scale acts on the slab variance. Laplace variance is
  # 2 * lambda^2; rescale lambda by sqrt(V_scale) so the
  # post-rescale variance equals V_scale * 2 * lambda^2.
  lambda <- fitted_g$scale[2L] * sqrt(V_scale)
  pi_0   <- fitted_g$pi[1L]
  pi_1   <- 1 - pi_0
  mu     <- fitted_g$mean[2L]

  if (pi_1 == 0 || lambda == 0) {
    # No slab mass; sum_t log(pi_0) per j.
    return(rep(ncol(bhat_slice) * log(pi_0), nrow(bhat_slice)))
  }

  bhat_vec <- as.vector(bhat_slice)
  shat_vec <- as.vector(shat_slice)
  log_dens_mix <- ebnm:::vloglik_point_laplace(
    x  = bhat_vec,
    s  = shat_vec,
    w  = pi_1,
    a  = 1 / lambda,
    mu = mu
  )
  log_dens_null <- dnorm(bhat_vec, mean = 0, sd = shat_vec, log = TRUE)
  rowSums(matrix(log_dens_mix - log_dens_null,
                 nrow = nrow(bhat_slice), ncol = ncol(bhat_slice)))
}

#' Per-(variable, position) posterior moments under a per-scale point-Laplace prior
#'
#' Closed-form posterior mean and second moment under the per-
#' scale point-Laplace prior. Delegates to
#' `ebnm::ebnm_point_laplace(g_init = fitted_g, fix_g = TRUE)` to
#' get the per-observation `(posterior$mean, posterior$second_moment)`
#' over the flattened `(j, t)` rectangle, then reshapes back to
#' `p x |idx_s|` matrices. ebnm holds `g` fixed via `fix_g` and
#' returns the analytic posterior under the Laplace + Normal
#' convolution.
#'
#' @inheritParams mixture_log_bf_laplace_per_scale
#' @return list with `pmean` and `pmean2`, each a `p x |idx_s|` matrix.
#' @keywords internal
#' @noRd
mixture_posterior_laplace_per_scale <- function(bhat_slice,
                                                 shat_slice,
                                                 fitted_g,
                                                 V_scale = 1) {
  if (!is.matrix(bhat_slice) || !is.matrix(shat_slice)) {
    stop("`bhat_slice` and `shat_slice` must be matrices.")
  }
  p     <- nrow(bhat_slice)
  T_idx <- ncol(bhat_slice)
  # V-rescale the slab as in the lbf kernel.
  g_eff <- fitted_g
  g_eff$scale[2L] <- fitted_g$scale[2L] * sqrt(V_scale)

  fit <- ebnm::ebnm_point_laplace(
    x      = as.vector(bhat_slice),
    s      = as.vector(shat_slice),
    g_init = g_eff,
    fix_g  = TRUE,
    output = c("posterior_mean", "posterior_second_moment")
  )
  list(pmean  = matrix(fit$posterior$mean,          nrow = p, ncol = T_idx),
       pmean2 = matrix(fit$posterior$second_moment, nrow = p, ncol = T_idx))
}
