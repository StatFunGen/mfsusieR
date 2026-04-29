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
#'   must be a power of two for `method = "wavelet_eb"`.
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
#' @param null_prior_weight numeric, penalty on the null component
#'   (passed through to the EB prior init). Default `10`.
#' @param init_pi0_w numeric, starting null-component mass for
#'   mixsqp. Default `1`.
#' @param control_mixsqp optional named list of `mixsqp` control
#'   args. Default `list(verbose = FALSE, eps = 1e-6,
#'   numiter.em = 4L)`.
#' @param grid_mult numeric, multiplier for the ash mixture grid.
#'   Default `sqrt(2)`.
#' @param wavelet_magnitude_cutoff non-negative numeric. Wavelet
#'   columns whose `median(|column|)` is at or below this cutoff
#'   are zeroed and treated as uninformative (`Bhat = 0`,
#'   `Shat = 1`) at every outer iteration. Default `0`.
#' @param wavelet_qnorm logical. When `TRUE` (default), applies a
#'   column-wise rank-based normal quantile transform to the
#'   wavelet-domain response before the EB regression loop. The
#'   returned `Y_adjusted` is on the original position scale so
#'   downstream `fsusie(Y_adjusted, X)` calls operate on the same
#'   units as the input `Y`.
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
#' @importFrom ashr postmean postsd set_data get_fitted_g
#' @export
mf_adjust_for_covariates <- function(Y, Z, X = NULL,
                                     method = c("wavelet_eb", "ols"),
                                     wavelet_basis_order = 10L,
                                     wavelet_family        = "DaubLeAsymm",
                                     max_iter              = 50L,
                                     tol                   = 1e-3,
                                     null_prior_weight     = 10,
                                     init_pi0_w            = 1,
                                     control_mixsqp        = list(verbose = FALSE,
                                                                  eps = 1e-6,
                                                                  numiter.em = 4L),
                                     grid_mult                = sqrt(2),
                                     wavelet_magnitude_cutoff = 0,
                                     wavelet_qnorm            = TRUE) {
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
  if (!is.logical(wavelet_qnorm) || length(wavelet_qnorm) != 1L)
    stop("`wavelet_qnorm` must be `TRUE` or `FALSE`.")

  if (method == "ols") {
    return(mf_residualize_ols(Y, Z, X))
  }

  mf_residualize_wavelet_eb(Y, Z, X = X,
                            wavelet_basis_order = wavelet_basis_order,
                            wavelet_family        = wavelet_family,
                            max_iter              = max_iter,
                            tol                   = tol,
                            null_prior_weight     = null_prior_weight,
                            init_pi0_w            = init_pi0_w,
                            control_mixsqp        = control_mixsqp,
                            grid_mult             = grid_mult,
                            wavelet_magnitude_cutoff      = wavelet_magnitude_cutoff,
                            wavelet_qnorm         = wavelet_qnorm)
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
  ZtZ_inv_Zt <- solve(crossprod(Z), t(Z))   # K x n
  fitted_func <- ZtZ_inv_Zt %*% Y           # K x T
  Y_adjusted  <- Y - Z %*% fitted_func
  X_adjusted  <- if (is.null(X)) NULL else X - Z %*% (ZtZ_inv_Zt %*% X)
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
                                      null_prior_weight     = 10,
                                      init_pi0_w            = 1,
                                      control_mixsqp        = list(verbose = FALSE,
                                                                   eps = 1e-6,
                                                                   numiter.em = 4L),
                                      grid_mult                = sqrt(2),
                                      wavelet_magnitude_cutoff = 0,
                                      wavelet_qnorm            = TRUE) {
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

  # Optional preprocessing on the wavelet-domain matrix.
  # `wavelet_qnorm` is applied first; the median-based low-count
  # mask is computed on the (possibly transformed) coefficients
  # so flagged columns are treated as uninformative for the rest
  # of the run.
  if (wavelet_qnorm) {
    Y_wd <- mf_quantile_normalize(Y_wd)
  }
  lowc_idx <- mf_low_count_indices(Y_wd, threshold = wavelet_magnitude_cutoff)
  if (length(lowc_idx) > 0L) {
    Y_wd[, lowc_idx] <- 0
  }

  # Step 3: prior init via the existing helper.
  prior_obj <- init_scale_mixture_prior_default(
    Y_m            = Y_wd,
    X              = Z_scaled,
    prior_class    = "mixture_normal_per_scale",
    groups         = indx_lst,
    grid_multiplier = grid_mult)
  G_prior <- prior_obj$G_prior

  # Step 4: state init.
  K <- ncol(Z_scaled)
  fitted_wc  <- matrix(0, K, T_pad)
  fitted_wc2 <- matrix(0, K, T_pad)
  MLE_wc     <- matrix(0, K, T_pad)
  MLE_wc2    <- matrix(0, K, T_pad)
  sigma2     <- mean(apply(Y_cent, 2, var))   # initial guess

  pi_hist <- list(.flatten_pi_per_scale(G_prior, indx_lst))

  # Step 5: outer loop. Per-covariate Gauss-Seidel coordinate
  # descent inside; per-scale prior-weight EM after the j sweep.
  check     <- 3 * tol
  iter      <- 1L
  converged <- FALSE
  while (check > tol && iter < max_iter) {
    # 5a: Gauss-Seidel sweep over covariates.
    for (j in seq_len(K)) {
      fitted_wc[j, ] <- 0
      partial_resid <- Y_wd - Z_scaled %*% fitted_wc
      bs <- compute_marginal_bhat_shat(
        X = matrix(Z_scaled[, j], ncol = 1L),
        Y = partial_resid)
      bhat_j <- as.numeric(bs$Bhat)
      shat_j <- pmax(as.numeric(bs$Shat), 1e-32)
      # Low-count mask: flagged columns are treated as
      # uninformative (Bhat = 0, Shat = 1).
      if (length(lowc_idx) > 0L) {
        bhat_j[lowc_idx] <- 0
        shat_j[lowc_idx] <- 1
      }
      MLE_wc [j, ] <- bhat_j
      # MLE_wc2 stores the per-coefficient second moment Shat^2
      # so `sqrt(MLE_wc2) = Shat` recovers the standard error
      # for the L_mixsq call below.
      MLE_wc2[j, ] <- shat_j^2

      # ash posterior per scale, in wavethresh ordering.
      pm <- numeric(T_pad)
      ps <- numeric(T_pad)
      for (s in seq_along(indx_lst)) {
        idx_s <- indx_lst[[s]]
        m     <- G_prior[[s]]
        data_s <- ashr::set_data(bhat_j[idx_s], shat_j[idx_s])
        g <- ashr::get_fitted_g(m)
        pm[idx_s] <- ashr::postmean(g, data_s)
        ps[idx_s] <- ashr::postsd  (g, data_s)
      }
      fitted_wc [j, ] <- pm
      fitted_wc2[j, ] <- ps^2
    }

    # 5b: residual variance update. The formula is intentionally
    # an approximation here (sigma2 is downstream-inert because
    # `compute_marginal_bhat_shat` ignores it); see
    # `inst/notes/refactor-exceptions.md` for the full discussion.
    # The stored sigma2 field is informational only and does not
    # affect `Y_adjusted`.
    # Posterior-expected squared-error decomposition under the
    # factorized variational posterior:
    #   E[||Y - X*beta||^2]
    #     = ||Y - X*mu||^2 + sum_j pw[j] * sum_t Var(beta[j,t])
    # where mu = E[beta] (= fitted_wc), pw[j] = colSums(X^2),
    # and Var(beta[j,t]) = fitted_wc2[j,t] (set to postsd^2
    # by the ash posterior block above).
    pw     <- colSums(Z_scaled^2)
    rss    <- sum((Y_wd - Z_scaled %*% fitted_wc)^2)
    er2    <- rss + sum(pw * rowSums(fitted_wc2))
    sigma2 <- er2 / (n * T_pad)

    # 5c: prior-weight update via mixsqp (is_ebmvfr = TRUE).
    new_G_prior <- G_prior
    for (s in seq_along(indx_lst)) {
      idx_s  <- indx_lst[[s]]
      sd_grid <- G_prior[[s]]$fitted_g$sd
      L <- mf_em_likelihood_per_scale(
        bhat_slice = MLE_wc [, idx_s, drop = FALSE],
        shat_slice = sqrt(MLE_wc2[, idx_s, drop = FALSE]),
        sd_grid    = sd_grid,
        is_ebmvfr  = TRUE)
      new_pi_s <- mf_em_m_step_per_scale(
        L              = L,
        zeta           = rep(1 / K, K),
        idx_size       = length(idx_s),
        mixsqp_null_penalty = null_prior_weight,
        init_pi0_w     = init_pi0_w,
        tol_null_prior = 0.001,
        control_mixsqp = control_mixsqp,
        is_ebmvfr      = TRUE)
      new_G_prior[[s]]$fitted_g$pi <- new_pi_s
    }
    G_prior <- new_G_prior

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

# Internal: flatten the per-scale pi vectors of `G_prior` into a
# single concatenated numeric vector for convergence tracking.
.flatten_pi_per_scale <- function(G_prior, indx_lst) {
  unlist(lapply(seq_along(indx_lst), function(s) {
    G_prior[[s]]$fitted_g$pi
  }))
}
