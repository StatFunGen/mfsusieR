# Pure-R reference implementations of functions with cpp11 counterparts.
# Pure-R oracle for cpp11 unit tests.
# Mirrors the mvsusieR pattern (`mvsusieR/R/reference_implementations.R`).
#
# These functions SHALL stay literal-equivalent to the cpp11 kernels:
# the unit tests in `tests/testthat/test_posterior_mixture.R` assert
# agreement at `<= 1e-12` on randomized inputs.
#
# The naming convention is `<kernel>_R` for the R oracle and
# `<kernel>` (no suffix) for the cpp11 kernel, matching mvsusieR.

# =============================================================================
# Mixture-of-normals per-(variable, scale) helpers (oracle for src/posterior_mixture.cpp)
# =============================================================================

#' Per-(variable, scale) log-Bayes factor: pure-R oracle
#'
#' Closed-form log-BF for the per-scale mixture-of-normals prior.
#' For variable j and the scale's |idx_s| positions, computes
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

#' Per-(variable, position) posterior moments: pure-R oracle
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

# =============================================================================
# HMM forward / backward / xi kernels (oracle for src/post_smooth_hmm.cpp)
# =============================================================================

#' HMM forward pass on a precomputed emission matrix: pure-R oracle
#'
#' Runs the scaled forward recursion
#'   alpha_tilde[t+1, ] = (alpha_hat[t, ] %*% P) * emit[t+1, ]
#'   G_t[t+1]           = sum(alpha_tilde[t+1, ])
#'   alpha_hat[t+1, ]   = alpha_tilde[t+1, ] / G_t[t+1]
#' with t1 initialization controlled by `t1_normalize`. When `FALSE`,
#' alpha_hat[1, ] is `pi * emit[1, ]` un-normalized and `G_t[1]` is
#' NA, matching the upstream `fsusieR::fit_hmm` convention for the
#' pre-EM warm-up forward pass. When `TRUE`, alpha_hat[1, ] is
#' normalized and `G_t[1] = sum(pi * emit[1, ])`, matching the
#' upstream EM-iteration forward pass.
#'
#' Used in tests as the `<= 1e-12` oracle for `hmm_forward_cpp`.
#'
#' @param emit T_pos x K matrix of emission densities.
#' @param P K x K transition matrix.
#' @param pi length-K initial-state distribution.
#' @param t1_normalize logical; see Details above.
#' @return list with `alpha_hat` (T_pos x K) and `G_t` (length T_pos).
#' @keywords internal
#' @noRd
hmm_forward_R <- function(emit, P, pi, t1_normalize) {
  T_pos <- nrow(emit)
  K     <- ncol(emit)
  alpha_hat <- matrix(NA_real_, T_pos, K)
  G_t       <- rep(NA_real_, T_pos)

  alpha_tilde_1 <- pi * emit[1L, ]
  if (t1_normalize) {
    G_t[1L]          <- sum(alpha_tilde_1)
    alpha_hat[1L, ]  <- alpha_tilde_1 / G_t[1L]
  } else {
    alpha_hat[1L, ]  <- alpha_tilde_1
  }

  for (t in seq_len(T_pos - 1L)) {
    m              <- alpha_hat[t, ] %*% P
    alpha_tilde    <- as.numeric(m) * emit[t + 1L, ]
    G_t[t + 1L]    <- sum(alpha_tilde)
    alpha_hat[t + 1L, ] <- alpha_tilde / G_t[t + 1L]
  }
  list(alpha_hat = alpha_hat, G_t = G_t)
}

#' HMM backward pass on a precomputed emission matrix: pure-R oracle
#'
#' Runs the scaled backward recursion
#'   beta_tilde[t, k] = sum_j P[k, j] * beta_hat[t+1, j] * emit[t+1, j]
#'   C_t[t]           = max(beta_tilde[t, ])
#'   beta_hat[t, ]    = beta_tilde[t, ] / C_t[t]
#' with `beta_hat[T_pos, ] = 1` and `C_t[T_pos] = NA`.
#'
#' Used in tests as the `<= 1e-12` oracle for `hmm_backward_cpp`.
#'
#' @inheritParams hmm_forward_R
#' @return list with `beta_hat` (T_pos x K) and `C_t` (length T_pos).
#' @keywords internal
#' @noRd
hmm_backward_R <- function(emit, P) {
  T_pos <- nrow(emit)
  K     <- ncol(emit)
  beta_hat <- matrix(NA_real_, T_pos, K)
  C_t      <- rep(NA_real_, T_pos)
  beta_hat[T_pos, ] <- 1
  for (t in (T_pos - 1L):1L) {
    weight       <- beta_hat[t + 1L, ] * emit[t + 1L, ]
    beta_tilde_t <- rowSums(sweep(P, 2L, weight, "*"))
    C_t[t]       <- max(beta_tilde_t)
    beta_hat[t, ] <- beta_tilde_t / C_t[t]
  }
  list(beta_hat = beta_hat, C_t = C_t)
}

#' HMM xi accumulator on precomputed alpha / beta / emission: pure-R oracle
#'
#' Computes
#'   xi = sum_t (outer(alpha_hat[t, ], beta_hat[t+1, ] * emit[t+1, ]) * P) / sum_xi_t
#' for t = 1..T_pos - 1. Skips t-slabs whose pre-normalization sum
#' is zero, matching the cpp11 implementation.
#'
#' Used in tests as the `<= 1e-12` oracle for `hmm_xi_cpp`.
#'
#' @param alpha_hat T_pos x K matrix from `hmm_forward_R`.
#' @param beta_hat  T_pos x K matrix from `hmm_backward_R`.
#' @param emit      T_pos x K emission matrix.
#' @param P         K x K transition matrix.
#' @return K x K xi matrix.
#' @keywords internal
#' @noRd
hmm_xi_R <- function(alpha_hat, beta_hat, emit, P) {
  T_pos <- nrow(alpha_hat)
  K     <- ncol(alpha_hat)
  xi <- matrix(0, K, K)
  for (t in seq_len(T_pos - 1L)) {
    xi_t <- outer(alpha_hat[t, ],
                  beta_hat[t + 1L, ] * emit[t + 1L, ]) * P
    s <- sum(xi_t)
    if (s > 0) xi <- xi + xi_t / s
  }
  xi
}

# =============================================================================
# Pure-R reference smoother for the HMM post-processor
# =============================================================================

#' Forward-backward EM smoother on a hub-and-spoke effect-size grid: pure-R oracle
#'
#' This is the pure-R reference implementation of `mf_fit_hmm`, used
#' in tests as the `<= 1e-12` oracle for the cpp11-backed production
#' version. It runs the same algorithm as `mf_fit_hmm` but does the
#' forward / backward / xi loops in pure R rather than via the
#' kernels in `src/post_smooth_hmm.cpp`. The pre-EM forward pass uses
#' the upstream `fsusieR::fit_hmm` t1 convention (alpha_hat[1, ]
#' un-normalized, G_t[1] = NA); the EM-iteration forward pass uses
#' the normalized t1 convention. See `hmm_forward_R`.
#'
#' @inheritParams mf_fit_hmm
#' @return Same shape as `mf_fit_hmm`.
#' @keywords internal
#' @noRd
mf_fit_hmm_R <- function(x, sd,
                         halfK            = 100L,
                         mult             = 3,
                         thresh           = 1e-5,
                         prefilter        = TRUE,
                         thresh_prefilter = 1e-30,
                         maxiter          = 3L,
                         max_zscore       = 20,
                         thresh_sd        = 1e-30,
                         epsilon          = 1e-2) {
  too_small <- which(sd < thresh_sd)
  if (length(too_small) > 0L) sd[too_small] <- thresh_sd
  bad <- which(is.na(sd))
  if (length(bad) > 0L) { x[bad] <- 0; sd[bad] <- 1 }
  bad <- which(!is.finite(sd))
  if (length(bad) > 0L) { x[bad] <- 0; sd[bad] <- 1 }
  big_z <- which(abs(x / sd) > max_zscore)
  if (length(big_z) > 0L) sd[big_z] <- abs(x[big_z]) / max_zscore

  K <- 2L * halfK - 1L
  X <- x
  T_pos <- length(X)

  pos <- seq(0, 1, length.out = halfK)
  mu  <- (pos^(1 / mult)) * 1.5 * max(abs(X))
  mu  <- c(mu, -mu[-1L])

  if (prefilter) {
    avg_post <- colMeans(do.call(rbind, lapply(seq_len(T_pos), function(i) {
      e <- dnorm(X[i], mean = mu, sd = sd[i])
      e / sum(e)
    })), na.rm = TRUE)
    keep <- which(avg_post > thresh_prefilter)
    if (!(1L %in% keep)) keep <- c(1L, keep)
    mu <- mu[keep]
    K  <- length(mu)
  }

  P <- matrix(0, K, K); diag(P) <- 0.5
  if (K > 1L) {
    P[1L, 2L:K] <- 0.5 / (K - 1L)
    P[2L:K, 1L] <- 0.5
  }
  P <- P + matrix(epsilon, K, K)
  if (K > 1L) {
    for (i in 2L:K) for (j in 2L:K) if (i != j) P[i, j] <- 0
  }
  P <- P / rowSums(P)

  pi_init <- rep(epsilon, K)
  pi_init[1L] <- if (K > 1L) 1 - sum(pi_init[-1L]) else 1

  emit_dn <- vapply(seq_len(T_pos),
                    function(t) dnorm(X[t], mean = mu, sd = sd[t]),
                    numeric(K))
  emit_dn <- if (K == 1L) matrix(emit_dn, ncol = 1L) else t(emit_dn)

  fwd <- hmm_forward_R(emit_dn, P, pi_init, t1_normalize = FALSE)
  bwd <- hmm_backward_R(emit_dn, P)
  alpha_hat <- fwd$alpha_hat
  beta_hat  <- bwd$beta_hat

  ab   <- alpha_hat * beta_hat
  prob <- ab / rowSums(ab)

  xi <- hmm_xi_R(alpha_hat, beta_hat, emit_dn, P)
  if (K > 1L) {
    for (i in 2L:K) for (j in 2L:K) if (i != j) xi[i, j] <- 0
  }
  rs <- rowSums(xi); rs[rs == 0] <- 1
  P  <- xi / rs
  P[P < epsilon & P > 0] <- epsilon
  P  <- P / rowSums(P)

  idx_comp <- which(colMeans(prob) > thresh)
  if (!(1L %in% idx_comp)) idx_comp <- c(1L, idx_comp)

  ash_obj <- vector("list", length(idx_comp))
  x_post  <- numeric(T_pos)
  for (i in 2L:length(idx_comp)) {
    mu_ash <- mu[idx_comp[i]]
    weight <- prob[, idx_comp[i]]
    ash_obj[[i]] <- ash(x, sd, weight = weight, mode = mu_ash,
                        mixcompdist = "normal")
    x_post <- x_post + weight * ash_obj[[i]]$result$PosteriorMean
  }
  prob <- prob[, idx_comp, drop = FALSE]
  K    <- length(idx_comp)
  P    <- P[idx_comp, idx_comp, drop = FALSE]
  P    <- P / rowSums(P)
  mu   <- mu[idx_comp]   # see refactor-exceptions.md HMM-mu-subset

  build_emit_iter <- function() {
    out <- matrix(NA_real_, T_pos, K)
    for (t in seq_len(T_pos)) {
      data_t <- set_data(X[t], sd[t])
      out[t, 1L] <- dnorm(X[t], mean = 0, sd = sd[t])
      if (K > 1L) for (k in 2L:K)
        out[t, k] <- exp(calc_loglik(ash_obj[[k]], data_t))
    }
    out
  }

  iter <- 1L
  G_t  <- NULL
  while (iter < maxiter) {
    pi_iter <- prob[1L, ]
    pi_iter[pi_iter < epsilon] <- epsilon
    pi_iter <- pi_iter / sum(pi_iter)

    emit_iter <- build_emit_iter()
    fwd       <- hmm_forward_R(emit_iter, P, pi_iter, t1_normalize = TRUE)
    bwd       <- hmm_backward_R(emit_iter, P)
    alpha_hat <- fwd$alpha_hat
    beta_hat  <- bwd$beta_hat
    G_t       <- fwd$G_t

    ab   <- alpha_hat * beta_hat
    prob <- ab / rowSums(ab)

    ash_obj <- vector("list", K)
    x_post  <- numeric(T_pos)
    for (k in 2L:K) {
      ash_obj[[k]] <- ash(x, sd, weight = prob[, k], mode = mu[k],
                          mixcompdist = "normal")
      x_post <- x_post + prob[, k] * ash_obj[[k]]$result$PosteriorMean
    }

    ab   <- alpha_hat * beta_hat
    prob <- ab / rowSums(ab)
    xi   <- hmm_xi_R(alpha_hat, beta_hat, emit_iter, P)
    if (K > 1L) {
      for (i in 2L:K) for (j in 2L:K) if (i != j) xi[i, j] <- 0
    }
    rs <- rowSums(xi); rs[rs == 0] <- 1
    P  <- xi / rs
    P[P < epsilon & P > 0] <- epsilon
    P  <- P / rowSums(P)

    iter <- iter + 1L
  }

  lfsr_est <- prob[, 1L]
  if (K > 1L) for (k in 2L:K)
    lfsr_est <- lfsr_est + prob[, k] * ash_obj[[k]]$result$lfsr

  ll_hmm  <- sum(log(G_t))
  ll_null <- sum(dnorm(X, mean = 0, sd = sd, log = TRUE))
  log_BF  <- ll_hmm - ll_null

  list(prob = prob, x_post = x_post, lfsr = lfsr_est, mu = mu,
       ll_hmm = ll_hmm, ll_null = ll_null, log_BF = log_BF)
}
