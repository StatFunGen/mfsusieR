# HMM-based wavelet-free post-smoother for `mfsusie` / `fsusie`.
#
# `mf_fit_hmm()` is the forward-backward EM smoother on a discrete
# grid of effect-size states with hub-and-spoke transitions and
# ashr refinement at each iteration. The univariate caller uses
# `susieR::compute_marginal_bhat_shat` for marginal Bhat / Shat.

#' Forward-backward EM smoother on a hub-and-spoke effect-size grid
#'
#' @param x numeric vector of position-wise effect estimates.
#' @param sd numeric vector of position-wise standard errors,
#'   same length as `x`.
#' @param halfK integer, half of the grid size. The full grid is
#'   `2 * halfK - 1` after symmetric mirroring around zero.
#' @param mult numeric, exponent on the position-wise grid
#'   spacing. Default 3.
#' @param thresh numeric, posterior-mass cutoff for keeping a
#'   non-null state in the second EM iteration.
#' @param prefilter logical; if `TRUE`, drop grid states with
#'   negligible average posterior weight before EM.
#' @param thresh_prefilter numeric, prefilter weight threshold.
#' @param maxiter integer, EM iteration cap.
#' @param max_zscore numeric, z-score floor on `abs(x/sd)` to
#'   avoid emission underflow.
#' @param thresh_sd numeric, lower clamp on `sd`.
#' @param epsilon numeric, transition-matrix smoothing constant.
#' @return a list with elements
#'   `prob` (T x K posterior state probabilities),
#'   `x_post` (length-T smoothed effect estimate),
#'   `lfsr` (length-T local false sign rate),
#'   `mu` (length-K state means), `ll_hmm`, `ll_null`, `log_BF`.
#' @keywords internal
#' @noRd
mf_fit_hmm <- function(x, sd,
                       halfK            = 100L,
                       mult             = 3,
                       thresh           = 1e-5,
                       prefilter        = TRUE,
                       thresh_prefilter = 1e-30,
                       maxiter          = 3L,
                       max_zscore       = 20,
                       thresh_sd        = 1e-30,
                       epsilon          = 1e-2,
                       ...) {
  # `...` collects extra args forwarded to every internal
  # `ashr::ash` call (e.g., `nullweight`, `mixcompdist`).
  ash_extras <- utils::modifyList(list(mixcompdist = "normal"),
                                  list(...))
  # Defensive sd / x cleanup.
  too_small <- which(sd < thresh_sd)
  if (length(too_small) > 0L) sd[too_small] <- thresh_sd
  bad <- which(is.na(sd))
  if (length(bad) > 0L) {
    x [bad] <- 0
    sd[bad] <- 1
  }
  bad <- which(!is.finite(sd))
  if (length(bad) > 0L) {
    x [bad] <- 0
    sd[bad] <- 1
  }
  big_z <- which(abs(x / sd) > max_zscore)
  if (length(big_z) > 0L) {
    sd[big_z] <- abs(x[big_z]) / max_zscore
  }

  K <- 2L * halfK - 1L
  X <- x
  T_pos <- length(X)

  # Effect-size grid: symmetric around zero with non-linear spacing.
  pos <- seq(0, 1, length.out = halfK)
  mu  <- (pos^(1 / mult)) * 1.5 * max(abs(X))
  mu  <- c(mu, -mu[-1L])

  # Prefilter rarely-occupied grid points.
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

  # Hub-and-spoke transition matrix.
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

  # Pre-EM forward / backward / occupancy on the dnorm emission.
  emit_dn <- function(t) dnorm(X[t], mean = mu, sd = sd[t])

  alpha_hat   <- matrix(NA_real_, T_pos, K)
  alpha_tilde <- matrix(NA_real_, T_pos, K)
  G_t         <- rep(NA_real_, T_pos)
  alpha_hat[1L, ]   <- pi_init * emit_dn(1L)
  alpha_tilde[1L, ] <- alpha_hat[1L, ]
  for (t in seq_len(T_pos - 1L)) {
    m <- alpha_hat[t, ] %*% P
    alpha_tilde[t + 1L, ] <- m * emit_dn(t + 1L)
    G_t[t + 1L] <- sum(alpha_tilde[t + 1L, ])
    alpha_hat[t + 1L, ] <- alpha_tilde[t + 1L, ] / G_t[t + 1L]
  }

  beta_hat   <- matrix(NA_real_, T_pos, K)
  beta_tilde <- matrix(NA_real_, T_pos, K)
  C_t        <- rep(NA_real_, T_pos)
  beta_hat[T_pos, ]   <- 1
  beta_tilde[T_pos, ] <- 1
  for (t in (T_pos - 1L):1L) {
    e <- emit_dn(t + 1L)
    beta_tilde[t, ] <- rowSums(sweep(P, 2L, beta_hat[t + 1L, ] * e, "*"))
    C_t[t] <- max(beta_tilde[t, ])
    beta_hat[t, ] <- beta_tilde[t, ] / C_t[t]
  }

  ab   <- alpha_hat * beta_hat
  prob <- ab / rowSums(ab)

  # Pre-EM transition update + structural-zero enforcement.
  xi <- matrix(0, K, K)
  for (t in seq_len(T_pos - 1L)) {
    xi_t <- outer(alpha_hat[t, ],
                  beta_hat[t + 1L, ] * emit_dn(t + 1L)) * P
    xi <- xi + xi_t / sum(xi_t)
  }
  if (K > 1L) {
    for (i in 2L:K) for (j in 2L:K) if (i != j) xi[i, j] <- 0
  }
  rs <- rowSums(xi); rs[rs == 0] <- 1
  P <- xi / rs
  P[P < epsilon & P > 0] <- epsilon
  P <- P / rowSums(P)

  # Pick states whose average posterior occupancy clears `thresh`,
  # always keeping the null state.
  idx_comp <- which(colMeans(prob) > thresh)
  if (!(1L %in% idx_comp)) idx_comp <- c(1L, idx_comp)

  # First ash refinement seeds x_post and x_post_sd. The
  # second moment is accumulated by law of total variance over
  # the state mixture; the null state contributes zero.
  ash_obj   <- vector("list", length(idx_comp))
  x_post    <- numeric(T_pos)
  x_post_m2 <- numeric(T_pos)
  for (i in 2L:length(idx_comp)) {
    mu_ash <- mu[idx_comp[i]]
    weight <- prob[, idx_comp[i]]
    ash_obj[[i]] <- do.call(ash,
                            c(list(x, sd, weight = weight, mode = mu_ash),
                              ash_extras))
    pm_i  <- ash_obj[[i]]$result$PosteriorMean
    psd_i <- ash_obj[[i]]$result$PosteriorSD
    x_post    <- x_post    + weight * pm_i
    x_post_m2 <- x_post_m2 + weight * (psd_i^2 + pm_i^2)
  }
  x_post_sd <- sqrt(pmax(x_post_m2 - x_post^2, 0))
  prob <- prob[, idx_comp, drop = FALSE]
  K    <- length(idx_comp)
  P    <- P[idx_comp, idx_comp, drop = FALSE]
  P    <- P / rowSums(P)
  # State k (k = 1..K) corresponds to original grid index
  # `idx_comp[k]`; subset `mu` so the EM emission `mu[k]` aligns
  # with state k. See `inst/notes/refactor-exceptions.md`
  # (HMM-mu-subset).
  mu   <- mu[idx_comp]

  # T_pos x K emission matrix: column 1 is the null dnorm at
  # zero, columns 2..K are the per-state ash marginal densities
  # via `calc_vloglik` (one S4 dispatch per state, vector input).
  build_emit_iter <- function() {
    out <- matrix(NA_real_, T_pos, K)
    out[, 1L] <- dnorm(X, mean = 0, sd = sd)
    if (K > 1L) {
      data_all <- set_data(X, sd)
      for (k in 2L:K)
        out[, k] <- exp(calc_vloglik(ash_obj[[k]]$fitted_g, data_all))
    }
    out
  }

  iter <- 1L
  while (iter < maxiter) {
    emit_mat <- build_emit_iter()

    alpha_hat   <- matrix(NA_real_, T_pos, K)
    alpha_tilde <- matrix(NA_real_, T_pos, K)
    G_t         <- rep(NA_real_, T_pos)

    pi_iter <- prob[1L, ]
    pi_iter[pi_iter < epsilon] <- epsilon
    pi_iter <- pi_iter / sum(pi_iter)
    alpha_tilde[1L, ] <- pi_iter * emit_mat[1L, ]
    G_t[1L] <- sum(alpha_tilde[1L, ])
    alpha_hat[1L, ] <- alpha_tilde[1L, ] / G_t[1L]

    for (t in seq_len(T_pos - 1L)) {
      m <- alpha_hat[t, ] %*% P
      alpha_tilde[t + 1L, ] <- m * emit_mat[t + 1L, ]
      G_t[t + 1L] <- sum(alpha_tilde[t + 1L, ])
      alpha_hat[t + 1L, ] <- alpha_tilde[t + 1L, ] / G_t[t + 1L]
    }

    beta_hat   <- matrix(NA_real_, T_pos, K)
    beta_tilde <- matrix(NA_real_, T_pos, K)
    C_t        <- rep(NA_real_, T_pos)
    beta_hat[T_pos, ]   <- 1
    beta_tilde[T_pos, ] <- 1
    for (t in (T_pos - 1L):1L) {
      e <- emit_mat[t + 1L, ]
      beta_tilde[t, ] <- rowSums(sweep(P, 2L, beta_hat[t + 1L, ] * e, "*"))
      C_t[t] <- max(beta_tilde[t, ])
      beta_hat[t, ] <- beta_tilde[t, ] / C_t[t]
    }

    ab   <- alpha_hat * beta_hat
    prob <- ab / rowSums(ab)

    ash_obj <- vector("list", K)
    x_post  <- numeric(T_pos)
    # Per-position posterior second moment for the credible band
    # via law of total variance. State 1 (null) has point mass at 0,
    # contributes 0 to both mean and second moment. For k >= 2 the
    # state's ash refinement gives PosteriorMean and PosteriorSD;
    # second moment per state is sd^2 + mean^2.
    x_post_m2 <- numeric(T_pos)
    for (k in 2L:K) {
      ash_obj[[k]] <- do.call(ash,
                              c(list(x, sd, weight = prob[, k], mode = mu[k]),
                                ash_extras))
      pm_k  <- ash_obj[[k]]$result$PosteriorMean
      psd_k <- ash_obj[[k]]$result$PosteriorSD
      x_post    <- x_post    + prob[, k] * pm_k
      x_post_m2 <- x_post_m2 + prob[, k] * (psd_k^2 + pm_k^2)
    }
    x_post_var <- pmax(x_post_m2 - x_post^2, 0)
    x_post_sd  <- sqrt(x_post_var)

    # Baum-Welch transition update. Re-emit against the freshly
    # refit `ash_obj`: the xi accumulator depends on post-refit
    # ash means / weights, not the pre-refit emissions used in
    # forward / backward.
    ab   <- alpha_hat * beta_hat
    prob <- ab / rowSums(ab)
    emit_mat <- build_emit_iter()
    xi   <- matrix(0, K, K)
    for (t in seq_len(T_pos - 1L)) {
      xi_t <- outer(alpha_hat[t, ],
                    beta_hat[t + 1L, ] * emit_mat[t + 1L, ]) * P
      xi <- xi + xi_t / sum(xi_t)
    }
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

  list(prob     = prob,
       x_post   = x_post,
       x_post_sd = x_post_sd,
       lfsr     = lfsr_est,
       mu       = mu,
       ll_hmm   = ll_hmm,
       ll_null  = ll_null,
       log_BF   = log_BF)
}

#' Univariate HMM regression of a position-space response on one
#' covariate column
#'
#' Used by `mf_post_smooth(method = "HMM")`. Per-position
#' marginal regression via `compute_marginal_bhat_shat()` followed
#' by HMM-mixture shrinkage; returns the position-space effect,
#' the local false sign rate, and the underlying log-BF.
#' @keywords internal
#' @noRd
mf_univariate_hmm_regression <- function(Y, X, halfK = 20L, ...) {
  Y <- as.matrix(Y)
  X <- as.matrix(X)
  stopifnot(ncol(X) == 1L)

  Y     <- col_scale(Y)
  csd_Y <- attr(Y, "scaled:scale")
  X     <- col_scale(X)
  csd_X <- attr(X, "scaled:scale")

  res <- compute_marginal_bhat_shat(matrix(X[, 1L], ncol = 1L), Y)
  est <- as.numeric(res$Bhat[1L, ])
  sds <- as.numeric(res$Shat[1L, ])

  bad <- is.na(sds) | sds <= 0
  if (any(bad)) {
    est[bad] <- 0
    sds[bad] <- median(sds[!bad], na.rm = TRUE)
  }

  s <- mf_fit_hmm(x = est, sd = sds, halfK = halfK, ...)
  list(effect_estimate = s$x_post    * csd_Y / csd_X,
       effect_sd       = s$x_post_sd * csd_Y / csd_X,
       lfsr            = s$lfsr,
       lBF             = s$log_BF)
}
