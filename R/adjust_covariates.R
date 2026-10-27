# Public covariate-adjustment utility for functional fine-mapping.
#
# `mf_adjust_for_covariates(Y, Z, X = NULL, method = ...)` removes
# the effect of covariates Z from a functional response Y prior to
# fine-mapping. Two methods:
#
#   - "wavelet_eb" (default): wavelet-domain empirical-Bayes
#     regression with a per-scale mixture-of-normals prior.
#   - "ols": closed-form Frisch-Waugh-Lovell residualization,
#     `(I - Z (Z'Z)^{-1} Z') Y`. The no-shrinkage limit of the
#     wavelet-EB path.
#
# When `X` is supplied the function additionally returns
# `X_adjusted = (I - Z (Z'Z)^{-1} Z') X`, which is the FWL
# correction for genotype-covariate correlation.
#
# Manuscript references:
#   methods/derivation.tex line 216 (factorized empirical-Bayes
#     mixture-weight update; reused for the per-scale prior here).
#   methods/online_method.tex eq:wavelet_residualization (covariate
#     adjustment in the wavelet domain).

#' Adjust a functional response for covariate effects
#'
#' Returns `Y_adjusted = Y - Z %*% fitted_func` where `fitted_func`
#' is a `K x T` matrix of per-position covariate effects fit by one
#' of two methods:
#'
#' - `"wavelet_eb"` (default): wavelet-domain empirical-Bayes
#'   regression with a per-scale mixture-of-normals prior. Suitable
#'   when the covariate effects are smooth functions of position.
#' - `"ols"`: closed-form Frisch-Waugh-Lovell residualization. The
#'   no-shrinkage limit of the wavelet-EB path; suitable when the
#'   covariate effects can vary arbitrarily across positions.
#'
#' When `X` is supplied, the function also returns
#' `X_adjusted = (I - Z (Z'Z)^{-1} Z') X` to handle covariate-
#' genotype correlation. Upstream covariate-adjustment workflows
#' typically do not residualize `X`; the option is provided for
#' downstream fine-mapping users who need it.
#'
#' @param Y `n x T` numeric matrix; the functional response. `T`
#'   must be a power of two for `method = "wavelet_eb"`. With
#'   `method = "ols"`, rows with any missing value in `Y`, `Z` (or
#'   `X`) are left out of the fit and returned as `NA`; the
#'   `"wavelet_eb"` method requires complete data.
#' @param Z `n x K` numeric matrix; covariates to adjust away.
#' @param X optional `n x p` numeric matrix; the genotype matrix
#'   to FWL-residualize. Default `NULL`.
#' @param method `"wavelet_eb"` (default) or `"ols"`.
#' @param wavelet_basis_order integer, selects the wavelet basis
#'   member within `wavelet_family` (number of vanishing moments
#'   for Daubechies families). Forwarded to `wavethresh::wd`'s
#'   `filter.number`. Default `10L`.
#' @param wavelet_family character, the `family` argument to
#'   `wd`. Default `"DaubExPhase"`.
#' @param max_iter integer, maximum outer iterations for the
#'   wavelet-EB path. Default `50`.
#' @param tol numeric, convergence tolerance on the per-scale
#'   prior weights `||pi - pi_prev|| / log(K) < tol`.
#'   Default `1e-3`.
#' @param init_pi0_w,control_mixsqp Legacy optimizer arguments
#'   retained for call compatibility. The wavelet-EB path now uses
#'   an exact responsibility M-step (Kim et al.) instead of mixsqp,
#'   so these arguments have no effect.
#' @param grid_mult numeric, multiplier for the ash mixture grid.
#'   Default `sqrt(2)`.
#' @param wavelet_magnitude_cutoff non-negative numeric. Wavelet
#'   columns whose `median(|column|)` is at or below this cutoff
#'   are zeroed and treated as uninformative (`Bhat = 0`,
#'   `Shat = 1`) at every outer iteration. Default `0`.
#' @return A named list with:
#'   \describe{
#'     \item{`Y_adjusted`}{`n x T` matrix, the residualized response.}
#'     \item{`X_adjusted`}{`n x p` matrix when `X` was supplied,
#'       else `NULL`.}
#'     \item{`fitted_func`}{`K x T` matrix of fitted covariate
#'       effects.}
#'     \item{`sigma2`}{numeric, residual variance estimate.}
#'     \item{`niter`}{integer, number of outer iterations
#'       (wavelet-EB only).}
#'     \item{`converged`}{logical, convergence status (wavelet-EB
#'       only).}
#'     \item{`method`}{character, the method used.}
#'   }
#' @references
#' Manuscript: methods/online_method.tex eq:wavelet_residualization
#' (covariate adjustment in the wavelet domain).
#' @export
mf_adjust_for_covariates <- function(Y, Z, X = NULL,
                                     method = c("wavelet_eb", "ols"),
                                     wavelet_basis_order = 10L,
                                     wavelet_family        = "DaubLeAsymm",
                                     max_iter              = 50L,
                                     tol                   = 1e-3,
                                     init_pi0_w            = 1,
                                     control_mixsqp        = list(verbose = FALSE,
                                                                  eps = 1e-6,
                                                                  numiter.em = 4L),
                                     grid_mult                = sqrt(2),
                                     wavelet_magnitude_cutoff = 0) {
  method <- match.arg(method)

  if (!is.matrix(Y) || !is.numeric(Y))
    stop("`Y` must be a numeric matrix.")
  if (!is.matrix(Z) || !is.numeric(Z))
    stop("`Z` must be a numeric matrix.")
  if (nrow(Y) != nrow(Z))
    stop("`Y` and `Z` must have the same number of rows.")
  if (!is.null(X)) {
    if (!is.matrix(X) || !is.numeric(X))
      stop("`X` must be a numeric matrix.")
    if (nrow(X) != nrow(Y))
      stop("`X` and `Y` must have the same number of rows.")
  }
  if (!is.numeric(wavelet_magnitude_cutoff) || length(wavelet_magnitude_cutoff) != 1L ||
      wavelet_magnitude_cutoff < 0)
    stop("`wavelet_magnitude_cutoff` must be a non-negative scalar.")

  if (method == "ols") {
    return(mf_residualize_ols(Y, Z, X))
  }
  if (anyNA(Y) || anyNA(Z) || (!is.null(X) && anyNA(X)))
    stop("`method = \"wavelet_eb\"` does not support missing values; ",
         "use `method = \"ols\"` or drop incomplete rows first.")

  mf_residualize_wavelet_eb(Y, Z, X = X,
                            wavelet_basis_order   = wavelet_basis_order,
                            wavelet_family        = wavelet_family,
                            max_iter              = max_iter,
                            tol                   = tol,
                            init_pi0_w            = init_pi0_w,
                            control_mixsqp        = control_mixsqp,
                            grid_mult             = grid_mult,
                            wavelet_magnitude_cutoff = wavelet_magnitude_cutoff)
}

#' OLS Frisch-Waugh-Lovell residualization
#'
#' Closed-form `Y_adjusted = (I - Z (Z'Z)^{-1} Z') Y`. Fast and
#' shrinkage-free; the no-shrinkage limit of
#' `mf_residualize_wavelet_eb()`. Internal helper; users should
#' call `mf_adjust_for_covariates(method = "ols")`.
#'
#' @inheritParams mf_adjust_for_covariates
#' @return Named list with `Y_adjusted`, `X_adjusted`,
#'   `fitted_func`, `method`. `niter` and `converged` are absent
#'   (the operation is closed-form).
#' @keywords internal
#' @noRd
mf_residualize_ols <- function(Y, Z, X = NULL) {
  # Fit on rows complete in Y and Z (and X when supplied); other rows
  # come back as NA, matching how `mfsusie()` drops rows with missing
  # Y via `complete.cases()`. Without this, one NA row turns every
  # entry of `Y_adjusted` into NA.
  obs <- complete.cases(Y) & complete.cases(Z)
  if (!is.null(X)) obs <- obs & complete.cases(X)
  if (!any(obs)) stop("No rows are complete in `Y` and `Z`.")
  if (all(obs)) {
    ZtZ_inv_Zt <- solve(crossprod(Z), t(Z))   # K x n
    fitted_func <- ZtZ_inv_Zt %*% Y           # K x T
    Y_adjusted  <- Y - Z %*% fitted_func
    X_adjusted  <- if (is.null(X)) NULL else X - Z %*% (ZtZ_inv_Zt %*% X)
  } else {
    Z_obs       <- Z[obs, , drop = FALSE]
    ZtZ_inv_Zt  <- solve(crossprod(Z_obs), t(Z_obs))
    fitted_func <- ZtZ_inv_Zt %*% Y[obs, , drop = FALSE]
    Y_adjusted  <- Y
    Y_adjusted[]     <- NA_real_
    Y_adjusted[obs, ] <- Y[obs, , drop = FALSE] - Z_obs %*% fitted_func
    X_adjusted  <- if (is.null(X)) NULL else {
      X_adj <- X
      X_adj[] <- NA_real_
      X_adj[obs, ] <- X[obs, , drop = FALSE] -
        Z_obs %*% (ZtZ_inv_Zt %*% X[obs, , drop = FALSE])
      X_adj
    }
  }
  list(
    Y_adjusted  = Y_adjusted,
    X_adjusted  = X_adjusted,
    fitted_func = fitted_func,
    method      = "ols"
  )
}

#' Wavelet-domain empirical-Bayes covariate residualization
#'
#' Inner workhorse for `method = "wavelet_eb"`. Implements the
#' per-coefficient empirical-Bayes regression of a functional
#' response `Y` on a covariate matrix `Z` (optionally controlling
#' for `X` via the Frisch-Waugh-Lovell construction). Internal
#' helper; users should call
#' `mf_adjust_for_covariates(method = "wavelet_eb")`.
#'
#' @inheritParams mf_adjust_for_covariates
#' @return Same named-list shape as `mf_adjust_for_covariates()`.
#' @keywords internal
#' @noRd
mf_residualize_wavelet_eb <- function(Y, Z, X = NULL,
                                      wavelet_basis_order = 10L,
                                      wavelet_family        = "DaubLeAsymm",
                                      max_iter              = 50L,
                                      tol                   = 1e-3,
                                      init_pi0_w            = 1,
                                      control_mixsqp        = list(verbose = FALSE,
                                                                   eps = 1e-6,
                                                                   numiter.em = 4L),
                                      grid_mult                = sqrt(2),
                                      wavelet_magnitude_cutoff = 0) {
  n     <- nrow(Y)
  T_pos <- ncol(Y)

  if (!is_wholenumber(log2(T_pos)))
    stop("`mf_residualize_wavelet_eb` requires `ncol(Y)` to be a power of two.",
         " For unevenly-spaced grids, remap with `remap_data()` first.")

  Y_org <- Y
  X_org <- X

  # Step 1: center+scale Z, center Y.
  Z_scaled <- col_scale(Z, scale = TRUE)
  csd_Z    <- attr(Z_scaled, "scaled:scale")
  Y_cent   <- col_scale(Y, scale = FALSE)

  # Step 2: row-wise DWT.
  W       <- dwt_matrix(Y_cent,
                        filter_number = wavelet_basis_order,
                        family        = wavelet_family,
                        max_scale     = log2(T_pos))
  Y_wd    <- cbind(W$D, W$C)        # n x T_basis
  T_pad   <- ncol(Y_wd)             # equals T_pos
  lev_res <- log2(T_pad)
  indx_lst <- gen_wavelet_indx(lev_res)

  # Low-count columns excluded from fitting and prior learning (matching
  # fsusieR EBmvFR behaviour: prior_groups excludes lowc_wc).
  lowc_idx <- mf_low_count_indices(Y_wd, threshold = wavelet_magnitude_cutoff)
  active_per_scale <- lapply(indx_lst, function(idx_s) {
    if (length(lowc_idx) > 0L) setdiff(idx_s, lowc_idx) else idx_s
  })

  # Step 3: prior init via the existing helper. Pass Y_wd without
  # zeroing lowc columns — the prior grid is set from the full data
  # distribution, which is unchanged.
  prior_obj <- init_scale_mixture_prior_default(
    Y_m            = Y_wd,
    X              = Z_scaled,
    prior_class    = "mixture_normal_per_scale",
    groups         = indx_lst,
    grid_multiplier = grid_mult)
  G_prior <- prior_obj$G_prior

  # Step 4: state init.
  K      <- ncol(Z_scaled)
  pw     <- colSums(Z_scaled^2)   # n-normalised predictor norms, reused each iter
  fitted_wc  <- matrix(0, K, T_pad)
  fitted_wc2 <- matrix(0, K, T_pad)
  sigma2     <- mean(apply(Y_cent, 2, var))   # initial guess (informational only)

  pi_hist <- list(.flatten_pi_per_scale(G_prior, indx_lst))

  # Step 5: outer loop. Gauss-Seidel coordinate descent over
  # covariates; exact responsibility M-step (Kim et al.) after the
  # sweep. The E-step and M-step are fused: `.wavelet_eb_normal_means`
  # returns posterior moments AND per-component responsibilities in
  # one pass, removing the need for a separate mixsqp sub-problem.
  check     <- 3 * tol
  iter      <- 1L
  converged <- FALSE
  while (check > tol && iter < max_iter) {
    # 5a: Gauss-Seidel sweep; accumulate responsibilities across j.
    mixture_counts <- lapply(seq_along(indx_lst), function(s) {
      matrix(0, K, length(G_prior[[s]]$fitted_g$pi))
    })
    for (j in seq_len(K)) {
      fitted_wc[j, ] <- 0        # exclude covariate j from residual
      partial_resid <- Y_wd - Z_scaled %*% fitted_wc
      bs <- mf_compute_marginal_bhat_shat(
        X = matrix(Z_scaled[, j], ncol = 1L),
        Y = partial_resid)
      bhat_j <- as.numeric(bs$Bhat)
      shat_j <- pmax(as.numeric(bs$Shat), 1e-32)
      for (s in seq_along(indx_lst)) {
        act_s <- active_per_scale[[s]]
        if (!length(act_s)) next
        post <- .wavelet_eb_normal_means(bhat_j[act_s], shat_j[act_s],
                                         G_prior[[s]]$fitted_g)
        fitted_wc [j, act_s] <- post$mean
        fitted_wc2[j, act_s] <- post$var
        mixture_counts[[s]][j, ] <- post$counts
      }
    }

    # 5b: residual variance on active columns only (matches fsusieR
    # EBmvFR: estimate_residual_variance uses active_wc).
    active_all <- unlist(active_per_scale)
    if (length(active_all)) {
      resid_act <- Y_wd[, active_all, drop = FALSE] -
                   (Z_scaled %*% fitted_wc)[, active_all, drop = FALSE]
      er2    <- sum(resid_act^2) +
                sum(pw * rowSums(fitted_wc2[, active_all, drop = FALSE]))
      sigma2 <- er2 / (n * length(active_all))
    }

    # 5c: exact EM M-step — average posterior responsibilities across
    # covariates and positions (Kim et al.; no mixsqp sub-problem).
    for (s in seq_along(indx_lst)) {
      totals <- colSums(mixture_counts[[s]])
      if (sum(totals) == 0) next
      G_prior[[s]]$fitted_g$pi <- totals / sum(totals)
    }

    # 5d: convergence on flattened pi (sum |delta| / log(K_total)).
    flat_now  <- .flatten_pi_per_scale(G_prior, indx_lst)
    flat_prev <- pi_hist[[length(pi_hist)]]
    pi_hist[[length(pi_hist) + 1L]] <- flat_now
    check <- sum(abs(flat_now - flat_prev)) / log(length(flat_now))

    iter <- iter + 1L
  }
  converged <- (check <= tol)

  # Step 6: per-covariate inverse DWT to get fitted_func in
  # position space. Each row of fitted_wc is a wavelet-coefficient
  # vector; mf_invert_dwt expects a (n x T) matrix where each row
  # is the packed (D, C) representation. We invert each covariate
  # independently and stack.
  fitted_func <- matrix(0, K, T_pos)
  for (j in seq_len(K)) {
    fc_row <- mf_invert_dwt(
      D_packed      = matrix(fitted_wc[j, ], nrow = 1L),
      column_center = rep(0, T_pos),
      column_scale  = rep(1, T_pos),
      filter_number = wavelet_basis_order,
      family        = wavelet_family)
    fitted_func[j, ] <- as.numeric(fc_row) / csd_Z[j]
  }

  # Step 7: output adjustment.
  Y_adjusted <- Y_org - Z %*% fitted_func
  X_adjusted <- if (is.null(X_org)) NULL else {
    ZtZ_inv_Zt <- solve(crossprod(Z), t(Z))
    X_org - Z %*% (ZtZ_inv_Zt %*% X_org)
  }

  list(
    Y_adjusted  = Y_adjusted,
    X_adjusted  = X_adjusted,
    fitted_func = fitted_func,
    sigma2      = sigma2,
    niter       = iter,
    converged   = converged,
    method      = "wavelet_eb"
  )
}

# Internal: compute posterior moments and mixture-component responsibilities
# for one covariate at one wavelet scale, under a zero-centred normal mixture
# prior. Implements the same computation as fsusieR's `.ebmvfr_normal_means`.
#
# Returns a list with:
#   mean   -- posterior mean, length |bhat|
#   var    -- posterior variance, length |bhat|
#   counts -- colSums of responsibilities, length K_mix (for exact M-step)
.wavelet_eb_normal_means <- function(bhat, shat, g) {
  keep <- which(g$pi > 0)
  if (!length(keep)) {
    return(list(mean   = numeric(length(bhat)),
                var    = numeric(length(bhat)),
                counts = numeric(length(g$pi))))
  }
  prior_var  <- g$sd[keep]^2
  se2        <- shat^2
  # total_var[t, k] = se2[t] + prior_var[k]
  total_var  <- outer(se2, prior_var, "+")
  log_weight <- sweep(-0.5 * (log(2 * pi * total_var) + bhat^2 / total_var),
                      2L, log(g$pi[keep]), "+")
  # Subtract row max before exp to prevent underflow (same result up to
  # normalisation).
  offset     <- matrixStats::rowMaxs(log_weight)
  weight     <- exp(log_weight - offset)
  weight     <- weight / rowSums(weight)    # responsibilities: T x length(keep)
  # Posterior mean and variance under the mixture.
  shrink     <- sweep(1 / total_var, 2L, prior_var, "*")   # T x K_keep
  comp_mean  <- bhat * shrink
  post_mean  <- rowSums(weight * comp_mean)
  post_var   <- rowSums(weight * (se2 * shrink + (comp_mean - post_mean)^2))
  # Counts = sum of responsibilities over positions, per mixture component.
  counts        <- numeric(length(g$pi))
  counts[keep]  <- colSums(weight)
  list(mean = post_mean, var = post_var, counts = counts)
}

# Internal: flatten the per-scale pi vectors of `G_prior` into a
# single concatenated numeric vector for convergence tracking.
.flatten_pi_per_scale <- function(G_prior, indx_lst) {
  unlist(lapply(seq_along(indx_lst), function(s) {
    G_prior[[s]]$fitted_g$pi
  }))
}
