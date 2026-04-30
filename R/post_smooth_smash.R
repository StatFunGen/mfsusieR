# Post-smoother kernel + per-position smash routine. Two flavors:
#   "ash"    -- cycle-spinning + per-coefficient `ashr::ash` on the
#               wavelet decomposition of the per-position OLS
#               estimate, with the per-position OLS Shat used as
#               the noise level. No `smashr` dependency.
#   "smashr" -- `smashr::smash.gaus` on the per-position OLS
#               estimate. Requires the `smashr` Suggests package;
#               the `requireNamespace` gate lives in `mf_post_smooth`.

.post_smooth_smash <- function(fit, level, flavor = "ash") {
  z_crit <- stats::qnorm((1 + level) / 2)
  alpha  <- 1 - level
  Y_pos_cache <- vector("list", length(fit$dwt_meta$T_basis))
  kernel <- function(fit, l, m, T_m, level) {
    if (is.null(Y_pos_cache[[m]]))
      Y_pos_cache[[m]] <<- .iso_response_pos(fit, m)
    iso <- Y_pos_cache[[m]] - .other_effects_pos(fit, m, exclude = l)
    out <- univariate_smash_regression(iso, fit$X_eff[[l]], alpha,
                                       flavor = flavor)
    # Recover sd from the credible-band half-width:
    # (upper - lower) / 2 = z_crit * sd.
    sd_pos <- (out$cred_band[1L, ] - out$cred_band[2L, ]) /
              (2 * z_crit)
    list(effect_estimate = out$effect_estimate,
         credible_band   = cbind(out$cred_band[2L, ], out$cred_band[1L, ]),
         lfsr            = lfsr_from_gaussian(out$effect_estimate,
                                              sd_pos))
  }
  .smoother_loop(fit, level, kernel,
                 method_name = "smash", kind = "wavelet")
}

#' Per-position smash regression of a multi-position response
#'
#' For a single regressor `X` (length `n`) and a multi-position
#' response `Y` (`n x T`), compute the per-position OLS estimate
#' of `Y` on `X` and pass it to one of two empirical-Bayes wavelet
#' shrinkers selected by `flavor`. Returns the smoothed
#' position-space effect estimate and a `(1 - alpha)` credible
#' band derived from the posterior variance.
#'
#' Used internally by `mf_post_smooth(method = "smash")`; also
#' useful as a standalone post-processing utility on a single
#' variable and a single functional outcome.
#'
#' @param Y numeric matrix (`n x T`) or numeric vector of length
#'   `n` (treated as `n x 1`). Per-position response.
#' @param X numeric matrix (`n x 1`) or numeric vector of length
#'   `n`. Single regressor.
#' @param alpha numeric in `(0, 1)`. Credible-band level is
#'   `1 - alpha`. Default `0.05`.
#' @param flavor `"ash"` (default) runs cycle-spinning +
#'   per-coefficient `ashr::ash`; `"smashr"` calls
#'   `smashr::smash.gaus` and requires the `smashr` Suggests
#'   package.
#' @return A list with components `effect_estimate` (length `T`
#'   numeric) and `cred_band` (`2 x T` matrix with rows
#'   `c("up", "low")`).
#' @examples
#' \donttest{
#' set.seed(1L)
#' X <- matrix(rnorm(60), 60, 1L)
#' Y <- X %*% matrix(c(rep(0, 24), rep(1, 16), rep(0, 24)),
#'                   1L, 64L) + matrix(rnorm(60 * 64, sd = 0.5), 60)
#' out <- univariate_smash_regression(Y, X)
#' plot(out$effect_estimate, type = "l")
#' }
#' @export
univariate_smash_regression <- function(Y, X, alpha = 0.05,
                                        flavor = c("ash", "smashr")) {
  Y <- as.matrix(Y)
  X <- as.matrix(X)
  stopifnot(ncol(X) == 1L)
  flavor <- match.arg(flavor)

  coeff <- stats::qnorm(1 - alpha / 2)

  # Column-scale Y and X via the package helper, which mean-
  # centers and sd-rescales each column and stashes the per-
  # column scale on the result.
  Y <- col_scale(Y)
  csd_Y <- attr(Y, "scaled:scale")
  X <- col_scale(X)
  csd_X <- attr(X, "scaled:scale")

  # Per-position OLS estimate of Y on X via the susieR helper.
  # Bhat / Shat are J x T matrices; with J = 1 the SER kernel
  # collapses to a single-regressor path.
  res <- compute_marginal_bhat_shat(X, Y)
  est <- as.numeric(res$Bhat[1L, ])
  sds <- as.numeric(res$Shat[1L, ])

  # Defensive fix: drop NA / non-positive sds.
  bad <- is.na(sds) | sds <= 0
  if (any(bad)) {
    est[bad] <- 0
    sds[bad] <- stats::median(sds[!bad], na.rm = TRUE)
  }

  s <- if (flavor == "smashr") {
    smashr::smash.gaus(x = est, sigma = sds,
                       ashparam = list(optmethod = "mixVBEM"),
                       post.var = TRUE)
  } else {
    mf_smash_ash(noisy_signal = est, noise_level = sds)
  }

  # Undo X scaling.
  fitted_func <- s$mu.est * csd_Y / csd_X
  fitted_var  <- s$mu.est.var * (csd_Y / csd_X)^2

  up  <- fitted_func + coeff * sqrt(fitted_var)
  low <- fitted_func - coeff * sqrt(fitted_var)
  cred_band <- rbind(up = up, low = low)

  list(effect_estimate = fitted_func, cred_band = cred_band)
}

# Cycle-spinning + per-coefficient `ashr::ash` shrinkage on the
# wavelet decomposition of `noisy_signal`, given a per-position
# noise sd `noise_level`. Pads to a power-of-2 by reflection,
# double-reflects for boundary handling, averages the inverse
# transforms across `n.shifts` cyclic shifts. Pure dependencies:
# `wavethresh` + `ashr`.
mf_smash_ash <- function(noisy_signal, noise_level = 1,
                         n.shifts = 50L,
                         filter.number = 1L,
                         family = "DaubExPhase") {
  x <- as.numeric(noisy_signal)
  n <- length(x)
  sds <- if (length(noise_level) == 1L) rep(noise_level, n) else noise_level
  if (length(sds) != n)
    stop("`noise_level` must be length 1 or length(noisy_signal).")

  # Pad to next power-of-2 by reflection.
  if ((log2(n) %% 1) != 0) {
    next_pow2 <- 2L^ceiling(log2(n))
    extra     <- next_pow2 - n
    x_padded   <- c(x, rev(x[seq_len(extra)]))
    sds_padded <- c(sds, rev(sds[seq_len(extra)]))
  } else {
    x_padded   <- x
    sds_padded <- sds
  }

  # Double-reflect for boundary handling.
  x_reflect <- c(x_padded, rev(x_padded))
  s_reflect <- c(sds_padded, rev(sds_padded))
  n_padded  <- length(x_padded)
  pos_orig        <- seq_len(n)
  pos_orig_padded <- (n_padded + 1L):(n_padded + n)
  n_r <- length(x_reflect)
  k   <- floor(n / n.shifts)

  W <- wavethresh::GenW(n = n_r,
                        filter.number = filter.number,
                        family = family)
  wavelet_var <- as.numeric(W^2 %*% (s_reflect^2))

  est     <- vector("list", n.shifts)
  est_var <- vector("list", n.shifts)
  idx_wave <- gen_wavelet_indx(log2(n_r))

  for (i in seq_len(n.shifts)) {
    shifted_x <- c(x_reflect[(i * k + 1L):n_r], x_reflect[seq_len(i * k)])
    wd_shifted <- wavethresh::wd(shifted_x, filter.number = filter.number,
                                 family = family)
    d  <- rep(0, length(wd_shifted$D))
    d2 <- rep(0, length(wd_shifted$D) + 1L)  # extra slot for the C coef
    for (s_idx in seq_len(length(idx_wave) - 1L)) {
      ix <- idx_wave[[s_idx]]
      t_ash <- ashr::ash(wd_shifted$D[ix],
                         sqrt(wavelet_var[ix]),
                         nullweight = 300)
      d[ix]  <- t_ash$result$PosteriorMean
      d2[ix] <- t_ash$result$PosteriorSD^2
    }
    wd_shifted$D <- d
    est[[i]]     <- wavethresh::wr(wd_shifted)
    est_var[[i]] <- as.numeric(W^2 %*% d2)
  }

  recover <- function(shifted_x, i, k) {
    n_s <- length(shifted_x)
    shift_back <- (n_s - i * k) %% n_s
    c(shifted_x[(shift_back + 1L):n_s], shifted_x[seq_len(shift_back)])
  }
  est_f   <- do.call(rbind, lapply(seq_along(est),
                                   function(i) recover(est[[i]], i, k)))
  est_v_f <- do.call(rbind, lapply(seq_along(est_var),
                                   function(i) recover(est_var[[i]], i, k)))

  est_orig   <- matrix(0, n.shifts, n)
  est_v_orig <- matrix(0, n.shifts, n)
  for (i in seq_len(n.shifts)) {
    est_orig[i, ]   <- 0.5 * (est_f[i, pos_orig] +
                              rev(est_f[i, pos_orig_padded]))
    est_v_orig[i, ] <- 0.5 * (est_v_f[i, pos_orig] +
                              rev(est_v_f[i, pos_orig_padded]))
  }

  list(mu.est     = colMeans(est_orig),
       mu.est.var = colMeans(est_v_orig))
}
