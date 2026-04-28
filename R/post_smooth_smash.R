# Post-smoother that delegates to smashr::smash.gaus on a per-
# effect, per-outcome basis.
#
# Per-effect, per-outcome:
#   1. Per-effect, per-outcome, isolate the position-space
#      response by subtracting the alpha-weighted contribution
#      of every other effect.
#   2. Compute the per-position OLS estimate `(Bhat, Shat)` of
#      that response on the per-effect alpha-weighted X
#      aggregate, after column-scaling both sides.
#   3. Pass `(Bhat, Shat)` to `smashr::smash.gaus(post.var = TRUE)`.
#   4. Undo the column scaling and assemble the position-space
#      effect estimate and credible band.
#
# `smashr` is a Suggests dependency; the dispatch in
# `mf_post_smooth` checks `requireNamespace` before reaching this
# kernel.

.post_smooth_smash <- function(fit, level) {
  z_crit <- stats::qnorm((1 + level) / 2)
  meta   <- fit$dwt_meta
  M      <- length(meta$T_basis)
  L      <- nrow(fit$alpha)
  alpha  <- 1 - level

  effect_curves  <- vector("list", M)
  credible_bands <- vector("list", M)
  lfsr_curves    <- vector("list", M)

  for (m in seq_len(M)) {
    T_m <- meta$T_basis[m]
    effect_curves[[m]]  <- vector("list", L)
    credible_bands[[m]] <- vector("list", L)
    lfsr_curves[[m]]    <- vector("list", L)
    if (T_m == 1L) {
      # Scalar outcome: smash is a no-op; fall back to scalewise.
      tmp <- .post_smooth_scalewise(fit, level, threshold_factor = 1)
      effect_curves[[m]]  <- tmp$effect_curves[[m]]
      credible_bands[[m]] <- tmp$credible_bands[[m]]
      lfsr_curves[[m]]    <- tmp$lfsr_curves[[m]]
      next
    }

    Y_pos <- .iso_response_pos(fit, m)
    for (l in seq_len(L)) {
      iso_pos <- Y_pos - .other_effects_pos(fit, m, exclude = l)
      x_eff   <- fit$X_eff[[l]]

      out <- univariate_smash_regression(iso_pos, x_eff, alpha)
      effect_curves[[m]][[l]]  <- out$effect_estimate
      credible_bands[[m]][[l]] <- cbind(out$cred_band[2L, ],
                                        out$cred_band[1L, ])
      # Recover sd from the credible-band half-width: (upper -
      # lower) / 2 = z_crit * sd, so sd = (upper - lower) /
      # (2 * z_crit). Then lfsr = pnorm(-|mean| / sd).
      sd_pos <- (out$cred_band[1L, ] - out$cred_band[2L, ]) /
                (2 * z_crit)
      lfsr_curves[[m]][[l]] <- gaussian_lfsr(out$effect_estimate,
                                              sd_pos)
    }
  }
  list(effect_curves  = effect_curves,
       credible_bands = credible_bands,
       lfsr_curves    = lfsr_curves)
}

#' Per-position smash regression of a multi-position response
#'
#' For a single regressor `X` (length `n`) and a multi-position
#' response `Y` (`n x T`), compute the per-position OLS estimate
#' of `Y` on `X` and pass it to `smashr::smash.gaus` for
#' empirical-Bayes wavelet shrinkage. Returns the smoothed
#' position-space effect estimate and a `(1 - alpha)` credible
#' band derived from the smash posterior variance.
#'
#' Used internally by `mf_post_smooth(method = "smash")`; also
#' useful as a standalone post-processing utility on a single
#' SNP and a single functional outcome.
#'
#' @param Y numeric matrix (`n x T`) or numeric vector of length
#'   `n` (treated as `n x 1`). Per-position response.
#' @param X numeric matrix (`n x 1`) or numeric vector of length
#'   `n`. Single regressor.
#' @param alpha numeric in `(0, 1)`. Credible-band level is
#'   `1 - alpha`. Default `0.05`.
#' @return A list with components `effect_estimate` (length `T`
#'   numeric) and `cred_band` (`2 x T` matrix with rows
#'   `c("up", "low")`).
#' @examples
#' \donttest{
#' if (requireNamespace("smashr", quietly = TRUE)) {
#'   set.seed(1L)
#'   X <- matrix(rnorm(60), 60, 1L)
#'   Y <- X %*% matrix(c(rep(0, 24), rep(1, 16), rep(0, 24)),
#'                     1L, 64L) + matrix(rnorm(60 * 64, sd = 0.5), 60)
#'   out <- univariate_smash_regression(Y, X)
#'   plot(out$effect_estimate, type = "l")
#' }
#' }
#' @export
univariate_smash_regression <- function(Y, X, alpha = 0.05) {
  Y <- as.matrix(Y)
  X <- as.matrix(X)
  stopifnot(ncol(X) == 1L)

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

  s <- smashr::smash.gaus(
    x        = est,
    sigma    = sds,
    ashparam = list(optmethod = "mixVBEM"),
    post.var = TRUE
  )

  # Undo X scaling.
  fitted_func <- s$mu.est * csd_Y / csd_X
  fitted_var  <- s$mu.est.var * (csd_Y / csd_X)^2

  up  <- fitted_func + coeff * sqrt(fitted_var)
  low <- fitted_func - coeff * sqrt(fitted_var)
  cred_band <- rbind(up = up, low = low)

  list(effect_estimate = fitted_func, cred_band = cred_band)
}

