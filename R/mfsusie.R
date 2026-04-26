# Public entry point for mfsusieR.
#
# `mfsusie()` is a thin orchestrator over susieR's IBSS workhorse:
#   1. construct the `mf_individual` data class (wavelet decomposition,
#      X centering / scaling, prior init);
#   2. assemble the `params` list expected by `susie_workhorse`;
#   3. call `susie_workhorse(data, params)`, which dispatches
#      every per-effect / per-iteration step to the `.mf_individual` /
#      `.mfsusie` S3 methods registered by `.onLoad`;
#   4. return the resulting `mfsusie` fit.
#
# The wavelet-specific work happens at step 1 (data prep) and step 4
# (post-processors, in a separate PR group). The IBSS loop itself
# is susieR's, end-to-end.

#' Multi-functional, multi-modality SuSiE
#'
#' Fit the multi-functional, multi-modality Sum of Single Effects
#' (mfSuSiE) regression model. Each modality `m` carries a per-sample
#' functional response `Y_m` of length `T_m`; `T_m` may differ across
#' modalities, and either `T_m = 1` (univariate response) or `T_m > 1`
#' (functional response on a regular grid) is supported. mfsusieR runs
#' a per-modality wavelet decomposition, then dispatches into susieR's
#' IBSS workhorse with mfsusieR-specific S3 methods that handle the
#' per-(scale, modality) mixture-of-normals prior and residual
#' variance.
#'
#' @param X numeric matrix `n x p` of covariates (e.g., genotype dosages).
#' @param Y list of length `M`; each element a numeric matrix `n x T_m`
#'   (or a length-n vector when `T_m = 1`).
#' @param pos optional list of length `M`; each element a numeric vector
#'   of length `T_m` recording the sampling positions for that
#'   modality. Defaults to `seq_len(T_m)` per modality.
#' @param L integer, maximum number of single effects to fit.
#' @param prior_variance_grid optional length-K vector of mixture
#'   variances for the scale-mixture-of-normals prior. When `NULL`,
#'   the data-driven path
#'   (`compute_marginal_bhat_shat` + `ash`) selects the
#'   grid per modality.
#' @param prior_variance_scope `"per_scale_modality"` (default) or
#'   `"per_modality"`. Controls whether the mixture-weight matrix
#'   `pi_V[[m]]` has a per-scale row dimension or collapses across
#'   scales.
#' @param null_prior_weight numeric, weight on the null prior
#'   component. Default 2 per the manuscript and design.md D5.
#' @param cross_modality_prior optional cross-modality combiner
#'   object. Defaults to the trivial independence combiner
#'   (`cross_modality_prior_independent()`).
#' @param prior_weights optional length-p numeric vector, the
#'   variable-selection prior. Defaults to uniform `1/p`.
#' @param residual_variance optional list of length `M`, initial
#'   residual variance per modality. Defaults to the per-modality
#'   sample variance.
#' @param residual_variance_method `"per_scale_modality"` (default)
#'   or `"shared_per_modality"`. Controls the sigma2 update shape.
#' @param standardize logical, scale `X` columns to unit variance.
#' @param intercept logical, center `X` columns to mean zero.
#' @param save_residuals logical, store per-modality residuals on
#'   the fit (default `TRUE`).
#' @param max_iter integer, maximum IBSS iterations.
#' @param tol numeric, ELBO change tolerance for convergence.
#' @param coverage numeric in (0, 1), credible-set coverage.
#' @param min_abs_corr numeric, minimum SNP-to-SNP correlation
#'   inside a credible set (CS purity threshold).
#' @param L_greedy integer or `NULL`. When non-`NULL`, run susieR's
#'   greedy outer loop, growing the effect count from `L_greedy` up
#'   to `L` in linear steps until the fit saturates
#'   (`min(lbf) < lbf_min`).
#' @param lbf_min numeric saturation threshold for the greedy outer
#'   loop. Default 0.1.
#' @param mixture_weight_method one of `"mixsqp"` (default) or
#'   `"none"`. `"mixsqp"` runs the per-effect empirical-Bayes
#'   update of the mixture weights `pi_V[[m]]` per
#'   (modality, scale) using the `mixsqp` solver. `"none"` holds
#'   the prior fixed at the initial `prior_variance_grid` /
#'   `null_prior_weight`; required by the C1 (susieR) degeneracy
#'   contract.
#' @param verbose logical.
#' @param track_fit logical, retain a per-iteration tracking list on
#'   the fit. Default `FALSE`.
#' @param max_padded_log2 integer, log2 cap on the post-remap grid
#'   length per modality. Default 10.
#' @param wavelet_filter_number integer; see
#'   `filter.select`.
#' @param wavelet_family character; see `wd`.
#' @param control_mixsqp optional named list of `mixsqp` control
#'   arguments forwarded to the per-(modality, scale) M-step.
#' @param nullweight numeric, mixsqp null-component penalty weight.
#'   Default 0.7.
#'
#' @return A list of class `c("mfsusie", "susie")` carrying:
#' \describe{
#'   \item{`alpha`}{`L x p` posterior inclusion probabilities.}
#'   \item{`mu`, `mu2`}{`list[L]` of `list[M]` of `p x T_padded[m]`
#'     matrices: per-effect, per-modality posterior mean and second
#'     moment in the wavelet domain.}
#'   \item{`pi_V`}{`list[M]` of `S_m x K` mixture-weight matrices.}
#'   \item{`G_prior`}{`list[M]` of `list[S_m]` ash-shaped prior
#'     records (mutated by the IBSS loop).}
#'   \item{`sigma2`}{`list[M]` of either scalar (legacy mode) or
#'     length-`S_m` per-(scale, modality) residual variances.}
#'   \item{`elbo`}{numeric vector of ELBO values per iteration.}
#'   \item{`niter`, `converged`}{IBSS termination metadata.}
#'   \item{`pip`}{length-p posterior inclusion probabilities.}
#'   \item{`sets`}{credible sets via `susie_get_cs`.}
#'   \item{`fitted`}{`list[M]` of running per-modality fits in the
#'     wavelet domain.}
#'   \item{`residuals`}{`list[M]` per-modality residuals, when
#'     `save_residuals = TRUE`.}
#' }
#'
#' @references
#' Manuscript: methods/online_method.tex (mfsusieR algorithm).
#' @export
mfsusie <- function(X, Y,
                    pos                       = NULL,
                    L                         = 10,
                    prior_variance_grid       = NULL,
                    prior_variance_scope      = c("per_scale_modality",
                                                  "per_modality"),
                    null_prior_weight         = 2,
                    cross_modality_prior      = NULL,
                    prior_weights             = NULL,
                    residual_variance         = NULL,
                    residual_variance_method  = c("per_scale_modality",
                                                  "shared_per_modality"),
                    standardize               = TRUE,
                    intercept                 = TRUE,
                    save_residuals            = TRUE,
                    max_iter                  = 100,
                    tol                       = 1e-3,
                    coverage                  = 0.95,
                    min_abs_corr              = 0.5,
                    L_greedy                  = NULL,
                    lbf_min                   = 0.1,
                    mixture_weight_method     = c("mixsqp", "none"),
                    verbose                   = FALSE,
                    track_fit                 = FALSE,
                    max_padded_log2           = 10,
                    wavelet_filter_number     = 10,
                    wavelet_family            = "DaubLeAsymm",
                    control_mixsqp            = NULL,
                    nullweight                = 0.7) {
  prior_variance_scope     <- match.arg(prior_variance_scope)
  residual_variance_method <- match.arg(residual_variance_method)
  mixture_weight_method    <- match.arg(mixture_weight_method)
  # Translate the public choice to susieR's internal vocabulary.
  # susieR's `single_effect_regression.default` skips the
  # per-effect prior update entirely when `estimate_prior_method
  # == "none"`; "mixsqp" routes through our mfsusieR override
  # `optimize_prior_variance.mf_individual`, which runs the
  # mixsqp M step on `pi_V` per (modality, scale).
  estimate_prior_method <- if (mixture_weight_method == "mixsqp") "optim" else "none"

  # The wavelet basis needs at least 4 sampled positions per curve.
  # 1 position routes to the scalar (degenerate) path; 2 or 3
  # positions sit in a no-man's-land that fits per-position
  # scalars without smoothing benefit.
  ncol_Y <- vapply(Y,
    function(Ym) if (is.matrix(Ym)) ncol(Ym) else length(Ym),
    integer(1))
  short_Y <- which(ncol_Y > 1L & ncol_Y <= 3L)
  if (length(short_Y) > 0L) {
    detail <- paste(sprintf("Y[[%d]] has %d columns",
                            short_Y, ncol_Y[short_Y]),
                    collapse = "; ")
    warning_message(sprintf(
      "%s. Each Y[[m]] is a response matrix with one column per sampled position; the wavelet basis needs at least 4 columns to add power. Either collapse the response to a single column (treat it as a scalar phenotype) or measure more positions before calling mfsusie().",
      detail), style = "hint")
  }

  # 1. Construct the data class.
  data <- create_mf_individual(
    X                     = X,
    Y                     = Y,
    pos                   = pos,
    max_padded_log2       = max_padded_log2,
    wavelet_filter_number = wavelet_filter_number,
    wavelet_family        = wavelet_family,
    standardize           = standardize,
    intercept             = intercept,
    save_residuals        = save_residuals,
    verbose               = verbose
  )

  # 2. Build the prior (scale-mixture-of-normals + cross-modality
  #    combiner). The G_prior carried on the fit is mutated per
  #    effect by `optimize_prior_variance.mf_individual`.
  prior <- mf_prior_scale_mixture(
    data,
    prior_variance_grid  = prior_variance_grid,
    prior_variance_scope = prior_variance_scope,
    null_prior_weight    = null_prior_weight
  )

  # 3. Assemble params for `susie_workhorse`.
  params <- list(
    L                          = L,
    prior_weights              = prior_weights,
    prior                      = prior,
    cross_modality_prior       = cross_modality_prior,
    residual_variance          = residual_variance,
    residual_variance_method   = residual_variance_method,
    estimate_residual_variance = TRUE,
    estimate_prior_variance    = (estimate_prior_method != "none"),
    estimate_prior_method      = estimate_prior_method,   # forwarded to single_effect_regression scaffolding
    check_null_threshold       = 0,
    prior_tol                  = 1e-9,
    max_iter                   = max_iter,
    tol                        = tol,
    coverage                   = coverage,
    min_abs_corr               = min_abs_corr,
    n_purity                   = 100,
    intercept                  = intercept,
    standardize                = standardize,
    track_fit                  = track_fit,
    verbose                    = verbose,
    nullweight                 = nullweight,
    control_mixsqp             = control_mixsqp,
    L_greedy                   = L_greedy,
    lbf_min                    = lbf_min,
    refine                     = FALSE,
    unmappable_effects         = "none",
    residual_variance_lowerbound = 0,
    residual_variance_upperbound = Inf
  )

  # 4. Run the susieR workhorse. All per-effect and per-iteration
  #    work dispatches to the .mf_individual / .mfsusie S3 methods
  #    registered by `.onLoad`.
  fit <- susie_workhorse(data, params)

  # 5. Attach per-modality wavelet-domain residuals (D - fitted)
  #    when save_residuals = TRUE. susieR's `ibss_finalize` is not
  #    a generic, so we do this in the wrapper rather than via S3.
  if (isTRUE(save_residuals)) {
    fit$residuals <- vector("list", data$M)
    for (m in seq_len(data$M)) {
      fit$residuals[[m]] <- data$D[[m]] - fit$fitted[[m]]
    }
  }

  # 6. Stash the wavelet pipeline metadata + X scaling on the fit so
  #    `predict.mfsusie`, `coef.mfsusie`, `fitted.mfsusie` can
  #    project posterior coefficients back to the original Y scale
  #    without the user passing `data` again.
  fit$mf_meta <- list(
    M               = data$M,
    pos             = data$pos,
    T_padded        = data$T_padded,
    scale_index     = data$scale_index,
    csd_X           = data$csd_X,
    X_center        = data$wavelet_meta$X_center,
    X_scale         = data$wavelet_meta$X_scale,
    column_center   = data$wavelet_meta$column_center,
    column_scale    = data$wavelet_meta$column_scale,
    wavelet_filter  = data$wavelet_meta$filter_number,
    wavelet_family  = data$wavelet_meta$family
  )

  fit
}
