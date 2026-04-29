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
# (post-processors, separate API). The IBSS loop itself
# is susieR's, end-to-end.

#' Multi-functional, multi-outcome SuSiE
#'
#' Fit the multi-functional, multi-outcome Sum of Single Effects
#' (mfSuSiE) regression model. Each outcome `m` carries a per-sample
#' functional response `Y_m` of length `T_m`; `T_m` may differ across
#' outcomes, and either `T_m = 1` (univariate response) or `T_m > 1`
#' (functional response on a regular grid) is supported. mfsusieR runs
#' a per-outcome wavelet decomposition, then dispatches into susieR's
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
#' @param null_prior_weight numeric, weight on the null prior
#'   component. Default 2.
#' @param cross_outcome_prior optional cross-outcome combiner
#'   object. Defaults to the trivial independence combiner
#'   (`cross_outcome_prior_independent()`).
#' @param prior_weights optional length-p numeric vector, the
#'   variable-selection prior. Defaults to uniform `1/p`.
#' @param residual_variance optional list of length `M`, initial
#'   residual variance per outcome. Defaults to the per-outcome
#'   sample variance.
#' @param residual_variance_scope `"per_scale"` (default)
#'   or `"per_outcome"`. Controls the sigma2 update shape.
#' @param standardize logical, scale `X` columns to unit variance.
#' @param intercept logical, center `X` columns to mean zero.
#' @param max_iter integer, maximum IBSS iterations.
#' @param tol numeric, ELBO change tolerance for convergence.
#' @param coverage numeric in (0, 1), credible-set coverage.
#' @param min_abs_corr numeric, minimum variable-to-variable correlation
#'   inside a credible set (CS purity threshold).
#' @param L_greedy integer or `NULL`. When non-`NULL`, run susieR's
#'   greedy outer loop, growing the effect count from `L_greedy` up
#'   to `L` in linear steps until the fit saturates
#'   (`min(lbf) < greedy_lbf_cutoff`).
#' @param greedy_lbf_cutoff numeric saturation threshold for the
#'   greedy outer loop. Default 0.1.
#' @param estimate_prior_variance logical. When `TRUE` (default),
#'   run the per-effect empirical-Bayes update of the mixture
#'   weights `pi_V[[m]]` per (outcome, scale) using the `mixsqp`
#'   solver. When `FALSE`, hold the prior fixed at the initial
#'   `prior_variance_grid` / `null_prior_weight`; useful when
#'   collapsing `mfsusie()` to `susieR::susie()` for sanity
#'   checks.
#' @param convergence_method one of `"pip"` (default) or
#'   `"elbo"`. Selects the IBSS convergence criterion.
#'   `"pip"` stops when `max(abs(prev_alpha - alpha))` falls
#'   below `tol` (PIP-difference convergence) with stall
#'   detection -- the default because alpha is robust to the
#'   small per-iteration ELBO oscillations that arise from
#'   mixsqp's approximate M-step (Generalized EM residual:
#'   `ELBO(t+1) >= ELBO(t) - O(eps_mixsqp * L * M * S_m)`,
#'   Neal & Hinton 1998). `"elbo"` stops when the change in
#'   ELBO falls below `tol`; the per-iteration ELBO is a
#'   coherent variational free energy thanks to the
#'   `get_objective.mfsusie` post-iteration KL refresh
#'   (`refresh_lbf_kl.mf_individual`) which evaluates per-
#'   effect KL[l] against the iter-final pi_V rather than the
#'   pi_V state at the moment effect l was updated.
#' @param pip_stall_window integer. Number of consecutive
#'   iterations without PIP-difference improvement after which
#'   the PIP-based convergence path declares convergence even
#'   when `tol` has not been reached. Default 5. Only consulted
#'   when `convergence_method = "pip"` (or when ELBO produces
#'   non-finite values and the PIP fallback fires).
#' @param estimate_residual_variance logical. When `TRUE`
#'   (default), update `sigma2` per IBSS iteration via
#'   `update_variance_components`. When `FALSE`, hold `sigma2`
#'   at the initial value across iterations; useful for
#'   sensitivity analysis.
#' @param verbose logical.
#' @param track_fit logical, retain a per-iteration tracking list on
#'   the fit. Default `FALSE`.
#' @param max_padded_log2 integer, log2 cap on the post-remap grid
#'   length per outcome. Default 10.
#' @param wavelet_filter_number integer; see
#'   `filter.select`.
#' @param wavelet_family character; see `wd`.
#' @param low_count_filter non-negative numeric. Wavelet-domain
#'   columns with `median(|column|) <= low_count_filter` are
#'   flagged as uninformative and treated as `Bhat = 0`,
#'   `Shat = 1` at every IBSS iteration. Useful for
#'   sparse-coverage assays where many wavelet columns carry
#'   negligible signal. Default `0`; with the default the set
#'   is non-empty only when the response has at least one
#'   wavelet column whose absolute-value median is exactly zero.
#' @param quantile_norm logical. When `TRUE`, applies a
#'   column-wise rank-based normal quantile transform to the
#'   wavelet-domain response before the IBSS loop. Useful for
#'   non-Gaussian wavelet coefficients arising from heavy-
#'   tailed assays. Default `FALSE`.
#' @param control_mixsqp optional named list of `mixsqp` control
#'   arguments forwarded to the per-(outcome, scale) M-step.
#' @param mixsqp_null_penalty numeric, per-coefficient pseudo-count
#'   on the null component of the mixsqp M-step. The total null
#'   pseudo-count fed to mixsqp is `mixsqp_null_penalty * idx_size`,
#'   versus a total data weight of `idx_size` from the `zeta`-
#'   weighted SNP rows; the ratio `mixsqp_null_penalty : 1` is
#'   therefore the null:data weight ratio at the M-step. Default
#'   `0.1` corresponds to a 10% null:data ratio, anchored to
#'   `ashr::ash.workhorse(nullweight = 10)` over a typical
#'   `idx_size = 128` (ratio ~0.08). For multi-outcome fits the
#'   penalty is scaled by the number of outcomes `M` so that the
#'   null:data balance stays invariant as the per-effect SNP
#'   posterior concentrates linearly in `M`; single-outcome fits
#'   are unchanged.
#'   Default 0.7.
#' @param mixsqp_alpha_eps numeric, threshold below which a SNP's
#'   per-effect posterior `alpha[l, j]` is dropped from the
#'   mixsqp M-step input. The truncation error on the M-step
#'   gradient is bounded by `sum_{j outside} alpha_j *
#'   max_k(L_jk)`, well under floating-point precision for
#'   typical concentrated SuSiE posteriors. Set to `0` to use
#'   every SNP regardless of `alpha[l, j]`. Default `1e-6`.
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
#' @param attach_lbf_variable_outcome logical. When `TRUE` (default), the
#'   fit carries `lbf_variable_outcome`, an `L x p x M` array of per-(effect,
#'   variant, outcome) log Bayes factors, populated from the IBSS
#'   sweep. Consumed by `mf_post_outcome_configuration(fit,
#'   by = "outcome")` (and by `susieR::susie_post_outcome_configuration`
#'   directly). Set `FALSE` to drop these — extra storage is
#'   `L * p * M` doubles. When `FALSE`, recover the array by
#'   calling `mf_post_outcome_configuration(fit, X = X, Y = Y, ...)`
#'   with the original data.
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
#'     (when `residual_variance_scope = "per_outcome"`) or a
#'     length-`S_m` vector of per-scale residual variances
#'     (when `residual_variance_scope = "per_scale"`, the
#'     default). The richer per-scale shape is mfsusieR-
#'     specific; downstream code reading `sigma2` should
#'     handle both.}
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
#'     `lbf_variable_outcome` keeps the M axis intact. Attached when
#'     `attach_lbf_variable_outcome = TRUE` (default). Consumed by
#'     `susie_post_outcome_configuration(fit, by = "outcome")`.}
#' }
#'
#' @references
#' Manuscript: methods/online_method.tex (mfsusieR algorithm).
#' @export
mfsusie <- function(X, Y,
                    pos                       = NULL,
                    L                         = 10,
                    prior_variance_grid       = NULL,
                    prior_variance_scope      = c("per_outcome",
                                                  "per_scale"),
                    null_prior_weight         = 2,
                    cross_outcome_prior       = NULL,
                    prior_weights             = NULL,
                    residual_variance         = NULL,
                    residual_variance_scope   = c("per_outcome",
                                                  "per_scale"),
                    standardize               = TRUE,
                    intercept                 = TRUE,
                    max_iter                  = 100,
                    tol                       = 1e-4,
                    coverage                  = 0.95,
                    min_abs_corr              = 0.5,
                    L_greedy                  = NULL,
                    greedy_lbf_cutoff          = 0.1,
                    estimate_prior_variance   = TRUE,
                    convergence_method        = c("pip", "elbo"),
                    pip_stall_window          = 5L,
                    estimate_residual_variance = TRUE,
                    verbose                   = FALSE,
                    track_fit                 = FALSE,
                    max_padded_log2           = 10,
                    wavelet_filter_number     = 10,
                    wavelet_family            = "DaubLeAsymm",
                    low_count_filter          = 0,
                    quantile_norm             = FALSE,
                    control_mixsqp            = NULL,
                    mixsqp_null_penalty       = 0.1,
                    mixsqp_alpha_eps          = 1e-6,
                    model_init                = NULL,
                    small_sample_correction   = FALSE,
                    attach_smoothing_inputs   = TRUE,
                    attach_lbf_variable_outcome        = TRUE) {
  if (!is.logical(small_sample_correction) ||
      length(small_sample_correction) != 1L ||
      is.na(small_sample_correction)) {
    stop("`small_sample_correction` must be `TRUE` or `FALSE`.")
  }
  if (!is.logical(attach_lbf_variable_outcome) ||
      length(attach_lbf_variable_outcome) != 1L ||
      is.na(attach_lbf_variable_outcome)) {
    stop("`attach_lbf_variable_outcome` must be `TRUE` or `FALSE`.")
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
    X                     = X,
    Y                     = Y,
    pos                   = pos,
    max_padded_log2       = max_padded_log2,
    wavelet_filter_number = wavelet_filter_number,
    wavelet_family        = wavelet_family,
    standardize           = standardize,
    intercept             = intercept,
    low_count_filter      = low_count_filter,
    quantile_norm         = quantile_norm,
    verbose               = verbose
  )

  # 2. Build the prior (scale-mixture-of-normals + cross-outcome
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
    mixsqp_null_penalty        = mixsqp_null_penalty,
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
                                 else NULL,
    attach_lbf_variable_outcome         = isTRUE(attach_lbf_variable_outcome)
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
    L <- nrow(fit$alpha)
    fit$X_eff <- vector("list", L)
    for (l in seq_len(L)) {
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
