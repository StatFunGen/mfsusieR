# `mfsusie()` is multi-outcome functional regression using susieR's IBSS workhorse:
#   1. construct the `mf_individual` data class (wavelet decomposition,
#      X centering / scaling, prior init);
#   2. assemble the `params` list expected by `susie_workhorse`;
#   3. call `susie_workhorse(data, params)`, which dispatches
#      every per-effect / per-iteration step to the `.mf_individual` /
#      `.mfsusie` S3 methods registered by `.onLoad`;
#   4. return the resulting `mfsusie` fit.
#
# The wavelet-specific work happens at step 1 (data prep) and step 4
# (post-processors, separate API). 

#' Multi-functional, multi-outcome SuSiE
#'
#' Fit the multi-functional, multi-outcome Sum of Single Effects
#' (mfSuSiE) regression model. Each outcome `m` carries a per-sample
#' functional response `Y_m` of length `T_m`; `T_m` may differ across
#' outcomes, and either `T_m = 1` (univariate response) or `T_m > 1`
#' (functional response on a regular grid) is supported; for the latter
#' mfsusieR runs a per-outcome wavelet decomposition, then dispatches into susieR's
#' IBSS workhorse with mfsusieR-specific S3 methods that handle the
#' per-(scale, outcome) mixture-of-normals prior and residual
#' variance.
#'
#' @param X numeric matrix `n x p` of covariates (e.g., genotype dosages).
#' @param Y list of length `M`; each element a numeric matrix `n x T_m`
#'   (or a length-n vector when `T_m = 1`).
#' @param pos optional list of length `M`; each element a numeric vector
#'   of length `T_m` recording the sampling positions for that
#'   outcome. Defaults to `seq_len(T_m)` per outcome.
#' @param L integer, maximum number of single effects to fit.
#'   Default `20`. With the default `L_greedy = 5`, the IBSS
#'   workhorse grows the effect count from 5 toward this cap in
#'   steps of 5 until the fit saturates (`min(lbf) <
#'   greedy_lbf_cutoff`); set `L_greedy = NULL` to fit the full
#'   `L` effects in one round.
#' @param prior_variance_grid optional length-K vector of mixture
#'   variances for the scale-mixture-of-normals prior. When `NULL`,
#'   the data-driven path
#'   (`compute_marginal_bhat_shat` + `ash`) selects the
#'   grid per outcome.
#' @param prior_variance_scope `"per_outcome"` (default) or
#'   `"per_scale"`. Controls whether the mixture-weight matrix
#'   `pi_V[[m]]` collapses across scales (one weight vector per
#'   outcome) or keeps a separate row per wavelet scale. The
#'   per-outcome scope is the default because it is roughly
#'   `S_m`-fold cheaper per IBSS iter (M mixsqp calls per iter
#'   instead of M*S_m); switch to `"per_scale"` only when you
#'   need scale-specific mixture weights for power on shape-
#'   varying signals.
#' @param null_prior_init numeric in `[0, 1]`, initial `pi[null]`
#'   for the scale-mixture prior on `b_l` at IBSS iter 1. The
#'   EM M-step overwrites it within a few iterations. Default
#'   `0`.
#' @param cross_outcome_prior optional object that controls how the
#'   per-outcome log-Bayes factors are combined into the joint
#'   log-Bayes factor used by the SER step. Each IBSS effect first
#'   computes a length-`p` log-BF per outcome `m` (one entry per
#'   variable, summed over wavelet positions); these are then
#'   reduced to a single length-`p` joint log-BF before the
#'   posterior `alpha` is taken as a softmax over variables.
#'   Defaults to outcome independence
#'   (`cross_outcome_prior_independent()`), under which the joint
#'   log-BF is the elementwise sum of the per-outcome log-BFs;
#'   equivalently, the joint BF is the product of per-outcome BFs.
#'   The argument is an extension point for non-independence
#'   combiners (e.g. modality-covariance priors) registered as S3
#'   methods on `combine_outcome_lbfs`; only the independence
#'   combiner is shipped, so most users leave this `NULL`.
#' @param prior_weights optional length-p numeric vector, the
#'   variable-selection prior. Defaults to uniform `1/p`.
#' @param residual_variance optional list of length `M`, initial
#'   residual variance per outcome. Defaults to the per-outcome
#'   sample variance.
#' @param residual_variance_scope `"per_outcome"` (default) or
#'   `"per_scale"`. Controls the sigma2 update shape: a single
#'   scalar per outcome (default) or one scalar per wavelet scale
#'   per outcome.
#' @param standardize logical, scale `X` columns to unit variance.
#' @param intercept logical, center `X` columns to mean zero.
#' @param max_iter integer, maximum IBSS iterations. Default `50`.
#' @param tol numeric, ELBO change tolerance for convergence.
#' @param coverage numeric in (0, 1), credible-set coverage.
#' @param min_abs_corr numeric, minimum variable-to-variable correlation
#'   inside a credible set (CS purity threshold).
#' @param L_greedy integer or `NULL`. When non-`NULL` (default `5`),
#'   run susieR's greedy outer loop, growing the effect count from
#'   `L_greedy` up to `L` in linear steps until the fit saturates
#'   (`min(lbf) < greedy_lbf_cutoff`). Set `NULL` to fit the full
#'   `L` effects in one round.
#' @param greedy_lbf_cutoff numeric saturation threshold for the
#'   greedy outer loop. Default 0.1.
#' @param estimate_prior_variance logical. When `TRUE` (default),
#'   run the per-effect empirical-Bayes update of the mixture
#'   weights `pi_V[[m]]` per (outcome, scale) using the `mixsqp`
#'   solver. When `FALSE`, hold the prior fixed at the initial
#'   `prior_variance_grid` / `null_prior_init`; useful when
#'   collapsing `mfsusie()` to `susieR::susie()` for sanity
#'   checks.
#' @param convergence_method one of `"pip"` (default) or
#'   `"elbo"`. `"pip"` stops when
#'   `max(abs(prev_alpha - alpha)) < tol`, with stall detection
#'   via `pip_stall_window`. `"elbo"` stops when the per-iter
#'   ELBO change is below `tol`.
#' @param pip_stall_window integer, number of consecutive
#'   iterations without PIP improvement after which `"pip"`
#'   convergence declares done. Default `5`.
#' @param estimate_residual_variance logical. When `TRUE`
#'   (default), update `sigma2` per IBSS iteration via
#'   `update_variance_components`. When `FALSE`, hold `sigma2`
#'   at the initial value across iterations; useful for
#'   sensitivity analysis.
#' @param verbose logical, default `TRUE`.
#' @param track_fit logical, retain a per-iteration tracking list on
#'   the fit. Default `FALSE`.
#' @param max_padded_log2 integer, log2 cap on the post-remap grid
#'   length per outcome. Default 10.
#' @param wavelet_basis_order integer; selects the wavelet basis
#'   member within `wavelet_family`. For Daubechies families this
#'   equals the number of vanishing moments (higher = smoother,
#'   longer-support filter). Forwarded to `wavethresh::wd`'s
#'   `filter.number`; see also `filter.select`. Default `10`.
#' @param wavelet_family character; selects the wavelet family.
#'   Forwarded to `wavethresh::wd`'s `family`. Default
#'   `"DaubLeAsymm"` (Daubechies least-asymmetric, a.k.a. Symmlet).
#' @param wavelet_magnitude_cutoff non-negative numeric. After
#'   the wavelet decomposition, each wavelet column `t` of the
#'   response has its median absolute value `median(|Y_wd[, t]|)`
#'   compared against this cutoff; columns at or below are zeroed
#'   in the data and treated as `Bhat = 0`, `Shat = 1` at every
#'   IBSS iteration. Useful for sparse-coverage assays where
#'   high-frequency wavelet coefficients are dominated by zeros.
#'   Default `0`; only columns with a strictly-zero absolute-value
#'   median are masked.
#' @param wavelet_qnorm logical. When `TRUE` (default), applies a
#'   column-wise rank-based normal quantile transform to the
#'   wavelet-domain response before the IBSS loop. Useful for
#'   non-Gaussian wavelet coefficients arising from heavy-tailed
#'   assays.
#' @param control_mixsqp optional named list of `mixsqp` control
#'   arguments forwarded to the per-(outcome, scale) M-step.
#' @param mixture_null_weight numeric in `[0, 1]`, ratio of
#'   null-pseudo-weight to data weight in the mixsqp M-step.
#'   Regularizes the EB mixture prior toward null. `0` is
#'   unregularized MLE. Internally scaled by `M` so the
#'   null:data balance is invariant to outcome count; single-
#'   outcome fits are unchanged. Default `0.05`.
#' @param mixsqp_alpha_eps numeric, threshold below which a
#'   variable's per-effect posterior `alpha[l, j]` is dropped
#'   from the mixsqp M-step input. Set to `0` to use every
#'   variable. Default `5e-5`.
#' @param small_sample_correction logical. When `TRUE`, replaces
#'   the per-variable Wakefield Normal marginal Bayes factor in
#'   the SER step with a Johnson 2005 scaled Student-t marginal
#'   with `df = n - 1` degrees of freedom. Useful for small
#'   sample sizes where the Wakefield approximation
#'   under-propagates residual-variance uncertainty into the
#'   per-variable BF, inflating PIPs at null variants. The
#'   correction acts on variable selection probabilities only;
#'   posterior moments given inclusion are unchanged. Default
#'   `FALSE`.
#' @param model_init optional `mfsusie` fit object from a prior
#'   call. When supplied, the IBSS loop is seeded from the
#'   supplied `alpha`, `mu`, `mu2`, `KL`, `lbf`, `V`, `pi_V`,
#'   `sigma2`, `fitted`, and `intercept` rather than the cold-
#'   start zero state. The supplied fit must have the same `L`
#'   as the requested `L`; otherwise the call errors. Useful
#'   for resuming a long-running fit after a per-call iteration
#'   budget. Default `NULL` (cold start).
#' @param attach_smoothing_inputs logical. When `TRUE` (default),
#'   the fit carries `Y_grid` (post-remap position-space Y) and
#'   `X_eff` (per-effect alpha-weighted aggregate of X) so that
#'   `mf_post_smooth(fit)` runs without re-passing the data.
#'   Set `FALSE` to drop these and call
#'   `mf_post_smooth(fit, X = X, Y = Y, ...)` instead; useful
#'   when sharing fits where the per-individual data should
#'   not travel with the fit.
#'
#' @return A list of class `c("mfsusie", "susie")` carrying:
#' \describe{
#'   \item{`alpha`}{`L x p` posterior inclusion probabilities.}
#'   \item{`mu`, `mu2`}{`list[L]` of `list[M]` of `p x T_basis[m]`
#'     matrices: per-effect, per-outcome posterior mean and second
#'     moment in the wavelet domain.}
#'   \item{`pi_V`}{`list[M]` of `S_m x K` mixture-weight matrices.}
#'   \item{`G_prior`}{`list[M]` of `list[S_m]` ash-shaped prior
#'     records (mutated by the IBSS loop).}
#'   \item{`sigma2`}{`list[M]`. Element m is either a scalar
#'     (when `residual_variance_scope = "per_outcome"`, the
#'     default) or a length-`S_m` vector of per-scale residual
#'     variances (when `residual_variance_scope = "per_scale"`).
#'     The richer per-scale shape is mfsusieR-specific;
#'     downstream code reading `sigma2` should handle both.}
#'   \item{`elbo`}{numeric vector of ELBO values per iteration.}
#'   \item{`niter`, `converged`}{IBSS termination metadata.}
#'   \item{`pip`}{length-p posterior inclusion probabilities.}
#'   \item{`sets`}{credible sets via `susie_get_cs`.}
#'   \item{`fitted`}{`list[M]` of running per-outcome fits in the
#'     wavelet domain. Used by `model_init` to warm-start a
#'     follow-up call.}
#'   \item{`Y_grid`}{`list[M]`, the post-remap position-space
#'     response on the padded grid. Attached when
#'     `attach_smoothing_inputs = TRUE` (default).}
#'   \item{`X_eff`}{`list[L]` of per-effect alpha-weighted X
#'     aggregates. Attached when `attach_smoothing_inputs = TRUE`.}
#'   \item{`lbf_variable_outcome`}{`L x p x M` array of per-(effect, variant,
#'     outcome) log Bayes factors. Parallels `lbf_variable` (which is
#'     `L x p`, the joint composite summed across scales and outcomes);
#'     `lbf_variable_outcome` keeps the M axis intact. Always attached;
#'     consumed by `susie_post_outcome_configuration(fit, by = "outcome")`.}
#' }
#'
#' @references
#' Manuscript: methods/online_method.tex (mfsusieR algorithm).
#' @export
mfsusie <- function(X, Y,
                    pos                       = NULL,
                    L                         = 20,
                    prior_variance_grid       = NULL,
                    prior_variance_scope      = c("per_outcome",
                                                  "per_scale"),
                    null_prior_init           = 0,
                    cross_outcome_prior       = NULL,
                    prior_weights             = NULL,
                    residual_variance         = NULL,
                    residual_variance_scope   = c("per_outcome",
                                                  "per_scale"),
                    standardize               = TRUE,
                    intercept                 = TRUE,
                    max_iter                  = 50,
                    tol                       = 1e-4,
                    coverage                  = 0.95,
                    min_abs_corr              = 0.5,
                    L_greedy                  = 5,
                    greedy_lbf_cutoff          = 0.1,
                    estimate_prior_variance   = TRUE,
                    convergence_method        = c("pip", "elbo"),
                    pip_stall_window          = 5L,
                    estimate_residual_variance = TRUE,
                    verbose                   = TRUE,
                    track_fit                 = FALSE,
                    max_padded_log2           = 10,
                    wavelet_basis_order       = 10,
                    wavelet_family            = "DaubLeAsymm",
                    wavelet_magnitude_cutoff  = 0,
                    wavelet_qnorm             = TRUE,
                    control_mixsqp            = NULL,
                    mixture_null_weight               = 0.05,
                    mixsqp_alpha_eps          = 5e-5,
                    model_init                = NULL,
                    small_sample_correction   = FALSE,
                    attach_smoothing_inputs   = TRUE) {
  if (!is.logical(small_sample_correction) ||
      length(small_sample_correction) != 1L ||
      is.na(small_sample_correction)) {
    stop("`small_sample_correction` must be `TRUE` or `FALSE`.")
  }
  prior_variance_scope    <- match.arg(prior_variance_scope)
  residual_variance_scope <- match.arg(residual_variance_scope)
  convergence_method      <- match.arg(convergence_method)

  # Translate `estimate_prior_variance` to susieR's internal
  # vocabulary. susieR's `single_effect_regression.default` skips
  # the per-effect prior update entirely when
  # `estimate_prior_method == "none"`; `"optim"` routes through
  # our mfsusieR override `optimize_prior_variance.mf_individual`,
  # which runs the mixsqp M step on `pi_V` per (outcome, scale).
  estimate_prior_method <- if (isTRUE(estimate_prior_variance)) "optim" else "none"

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
    X                        = X,
    Y                        = Y,
    pos                      = pos,
    max_padded_log2          = max_padded_log2,
    wavelet_basis_order      = wavelet_basis_order,
    wavelet_family           = wavelet_family,
    standardize              = standardize,
    intercept                = intercept,
    wavelet_magnitude_cutoff = wavelet_magnitude_cutoff,
    wavelet_qnorm            = wavelet_qnorm,
    verbose                  = verbose
  )

  # 2. Build the prior (scale-mixture-of-normals + cross-outcome
  #    combiner). The G_prior carried on the fit is mutated per
  #    effect by `optimize_prior_variance.mf_individual`.
  prior <- mf_prior_scale_mixture(
    data,
    prior_variance_grid  = prior_variance_grid,
    prior_variance_scope = prior_variance_scope,
    null_prior_init    = null_prior_init
  )

  # 3. Assemble params for `susie_workhorse`.
  params <- list(
    L                          = L,
    prior_weights              = prior_weights,
    prior                      = prior,
    cross_outcome_prior        = cross_outcome_prior,
    residual_variance          = residual_variance,
    residual_variance_scope    = residual_variance_scope,
    estimate_residual_variance = isTRUE(estimate_residual_variance),
    estimate_prior_variance    = (estimate_prior_method != "none"),
    estimate_prior_method      = estimate_prior_method,   # forwarded to single_effect_regression scaffolding
    convergence_method         = convergence_method,
    pip_stall_window           = pip_stall_window,
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
    mixture_null_weight        = mixture_null_weight,
    mixsqp_alpha_eps           = mixsqp_alpha_eps,
    control_mixsqp             = control_mixsqp,
    L_greedy                   = L_greedy,
    greedy_lbf_cutoff                    = greedy_lbf_cutoff,
    refine                     = FALSE,
    unmappable_effects         = "none",
    # susieR's `get_objective.default` short-circuits on `params$use_NIG`
    # before reaching the standard ELBO branch; we don't use the
    # NIG / shrinkage prior in mfsusieR, so wire it explicitly to FALSE
    # rather than leaving it NULL (NULL would trip an `&& nrow(...) == 1`
    # short-circuit when use_NIG is undefined).
    use_NIG                    = FALSE,
    residual_variance_lowerbound = 0,
    residual_variance_upperbound = Inf,
    model_init                 = model_init,
    small_sample_correction    = small_sample_correction,
    small_sample_df            = if (small_sample_correction) data$n - 1L
                                 else NULL
  )

  # 4. Run the susieR workhorse. All per-effect and per-iteration
  #    work dispatches to the .mf_individual / .mfsusie S3 methods
  #    registered by `.onLoad`.
  fit <- susie_workhorse(data, params)

  # Mask the iter-1 ELBO. mfsusie initialises sigma2 = var(Y) and
  # does the first closed-form M-step at the end of iter 1; the
  # iter-1 ELBO is computed against the inflated initial sigma2 and
  # therefore sits well above iter-2's value. Marking it `NA_real_`
  # keeps it out of `diff(fit$elbo)` plots / downstream consumers
  # that don't already drop the leading element. `tail(fit$elbo, 1)`
  # (used by `susie_get_objective`, refinement loops, etc.) is
  # unaffected.
  if (length(fit$elbo) >= 1L) fit$elbo[1] <- NA_real_

  # 5. Attach the smoothing inputs unless the caller opted out.
  #    `Y_grid[[m]]` is the post-remap, position-space Y on the
  #    padded n x T_basis[m] grid (raw, no centring or scaling).
  #    `X_eff[[l]]` is the alpha-weighted aggregate raw-X column
  #    `X_raw %*% alpha[l, ]` and serves as the per-effect
  #    regression predictor for `mf_post_smooth()`.
  if (attach_smoothing_inputs) {
    fit$Y_grid <- vector("list", data$M)
    for (m in seq_len(data$M)) {
      fit$Y_grid[[m]] <- data$Y[[m]]
    }
    # Use the post-greedy effect count from the fit; in the
    # L-greedy outer loop the workhorse may return fewer effects
    # than the requested `L`.
    L_eff <- nrow(fit$alpha)
    fit$X_eff <- vector("list", L_eff)
    for (l in seq_len(L_eff)) {
      # data$X is centered + standardized; multiply by data$csd
      # (X_scale) to get the raw-X column equivalent of the
      # alpha-weighted aggregate.
      fit$X_eff[[l]] <- as.numeric(
        data$X %*% (fit$alpha[l, ] * data$csd))
    }
  }

  # 6. Stash the wavelet pipeline metadata + X scaling on the fit so
  #    `predict.mfsusie`, `coef.mfsusie`, `fitted.mfsusie` can
  #    project posterior coefficients back to the original Y scale
  #    without the user passing `data` again.
  fit$dwt_meta <- list(
    M               = data$M,
    pos             = data$pos,
    T_basis        = data$T_basis,
    scale_index     = data$scale_index,
    csd           = data$csd,
    X_center        = data$wavelet_meta$X_center,
    X_scale         = data$wavelet_meta$X_scale,
    column_center   = data$wavelet_meta$column_center,
    column_scale    = data$wavelet_meta$column_scale,
    wavelet_filter  = data$wavelet_meta$filter_number,
    wavelet_family  = data$wavelet_meta$family,
    outcome_names   = names(data$D) %||% names(Y)
  )

  fit
}

# ---------------------------------------------------------------
# Single-outcome wrapper. Canonicalises (Y, X, pos) into
# (X, list(Y), list(pos)) and forwards to `mfsusie()`. No numerics.

#' Single-outcome functional SuSiE
#'
#' Convenience wrapper around `mfsusie()` for the single-outcome
#' case (`M = 1`). The response `Y` may be scalar (`T = 1`,
#' susieR-degenerate path) or functional (`T > 1`, wavelet-basis
#' path). Internally `fsusie(Y, X, pos, ...)` is equivalent to
#' `mfsusie(X, list(Y), list(pos), ...)`.
#'
#' Multi-outcome-only arguments (e.g. `cross_outcome_prior`) are
#' rejected with an explicit error to keep the wrapper honest about
#' its `M = 1` scope.
#'
#' @param Y numeric matrix `n x T` of functional responses on a
#'   regular position grid (or numeric vector for `T = 1`).
#' @param X numeric matrix `n x p` of covariates.
#' @param pos optional numeric vector of length `T` recording
#'   sampling positions. Defaults to `seq_len(ncol(Y))`.
#' @param ... forwarded to `mfsusie()`. See `?mfsusie` for the full
#'   parameter set; arguments that only make sense in the
#'   multi-outcome case (currently just `cross_outcome_prior`)
#'   error out.
#'
#' @return a list of class `c("mfsusie", "susie")`. See `?mfsusie`
#'   for the documented field set.
#' @export
fsusie <- function(Y, X, pos = NULL, ...) {
  args <- list(...)

  # Reject arguments that only make sense for M >= 2.
  forbidden_mv <- c("cross_outcome_prior")
  bad <- intersect(names(args), forbidden_mv)
  if (length(bad) > 0L) {
    stop(sprintf(
      "`fsusie()` is the single-outcome wrapper. The following arguments are only meaningful for multi-outcome fits and SHALL not be passed here: %s. Use `mfsusie()` directly.",
      paste(bad, collapse = ", ")))
  }

  if (!is.matrix(Y)) {
    if (is.numeric(Y)) Y <- matrix(Y, ncol = 1)
    else stop("`Y` must be a numeric matrix or numeric vector.")
  }

  pos_arg <- if (is.null(pos)) NULL else list(pos)
  do.call(mfsusie, c(list(X = X, Y = list(Y), pos = pos_arg), args))
}
