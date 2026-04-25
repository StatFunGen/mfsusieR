# susieR backbone

This note maps `susieR` as it stands on `master`. mfsusieR will sit on top of this
package, so any field, function, or S3 dispatch point we plan to override has to
be located here first.

Source root: `~/GIT/susieR`. All citations are file paths relative to that root.

## 1. Core data structures

### Public fit object (`class(fit) == "susie"`)

The class is set at `R/iterative_bayesian_stepwise_selection.R:96`. The fit object
is a list. The fields below are the ones a downstream package needs to know
about; the matrices are sized by `L` (effects) and `p` (variables).

| Field | Shape | Source | Meaning |
|---|---|---|---|
| `alpha` | L by p | `R/model_methods.R:81` (init), updated each iter | Posterior inclusion probability per (effect, variable) |
| `mu` | L by p | same | Posterior mean conditional on inclusion |
| `mu2` | L by p | same | Posterior second moment conditional on inclusion |
| `V` | length L | same | Per-effect prior variance (scalar in default mode) |
| `sigma2` | scalar | same | Residual variance |
| `lbf` | length L | per-effect SER | Log Bayes factor of the effect |
| `lbf_variable` | L by p | per-effect SER | Variable-level log BF used to form alpha |
| `KL` | length L | per-effect SER | KL divergence of q(b_l) from prior |
| `pi` | length p | input or default | Prior weights on the p variables |
| `null_weight` | scalar | input | Mass on the null variable column |
| `elbo` | length niter+1 | per iteration | ELBO trace (objective at each iter) |
| `niter`, `converged` | scalar | termination | Iteration count and stop reason |
| `pip` | length p | `susie_get_pip()` | `1 - prod_l (1 - alpha[l, ])` |
| `sets` | list | `susie_get_cs()` | Credible-set indices, coverage, purity, cs_index |
| `intercept`, `fitted`, `Xr` | scalar / n / n | finalize | Fitted values on the original scale |
| `X_column_scale_factors` | length p | finalize | For unstandardizing coefficients |

Optional fields that appear conditionally: `theta`, `tau2`, `X_theta` for
unmappable-effects mode; `predictor_weights` for some SER variants;
`slot_weights` for slot-prior runs; `trace` if `track_fit = TRUE`.

Evidence: `R/model_methods.R:81-98` (`initialize_matrices.default`),
`R/iterative_bayesian_stepwise_selection.R:329-382` (`ibss_finalize`).

### Internal state passed through the loop

The same list object that becomes the fit return is the internal state
(`model` in the workhorse, sometimes `s` in older code paths). It carries
extra fields that the IBSS loop reads and writes between iterations:

- `runtime$prev_elbo`, `runtime$prev_alpha`, `runtime$best_pip_diff`,
  `runtime$stall_count` — convergence bookkeeping (`R/model_methods.R:135-230`).
- `residuals`, `fitted_without_l`, `raw_residuals` — per-effect partial residual
  intermediates set by `compute_residuals()` (`R/individual_data_methods.R:104`)
  and consumed by `single_effect_regression()`.
- `c_hat_state` — slot activity (only when a `slot_prior` is supplied).

### Credible-set object

`susie_get_cs()` (`R/susie_get_functions.R:277-415`) returns
`list(cs, coverage, purity, cs_index, requested_coverage)`. `cs` is a list of
integer index vectors. `purity` is added if a correlation matrix or X is
provided. `cs_index` records which of the original L effects produced each
returned CS (drops null effects with low V or failed purity).

## 2. IBSS inner-loop shape

The exposed orchestration is `susie_workhorse()` (`R/susie_workhorse.R:14-80`).
The relevant skeleton (`R/susie_workhorse.R:31-55`):

```
for (iter in seq_len(params$max_iter)) {
  tracking <- track_ibss_fit(data, params, model, tracking, iter, elbo)
  model    <- ibss_fit(data, params, model)              # one sweep over L
  elbo[iter + 1] <- get_objective(data, params, model)
  model    <- check_convergence(data, params, model, elbo, iter)
  model$runtime$prev_elbo  <- elbo[iter + 1]
  model$runtime$prev_alpha <- model$alpha
  if (model$converged) break
  model    <- update_model_variance(data, params, model)  # sigma2 update
}
```

`ibss_fit()` (`R/iterative_bayesian_stepwise_selection.R:171-220`) sweeps `l = 1..L`,
calling `single_effect_update()` for each effect (and updating slot weights when
a slot prior is active).

`single_effect_update()` (`R/single_effect_regression.R:196-208`) is three steps:

1. `compute_residuals(data, params, model, l)` builds the residual y minus the
   sum of the other L-1 effects, then forms the SER sufficient statistics.
2. `single_effect_regression(data, params, model, l)` runs the SER, optimizes
   `V[l]` if requested, computes `lbf`, `lbf_variable`, `alpha[l, ]`, `mu[l, ]`,
   `mu2[l, ]`, and `KL[l]`.
3. `update_fitted_values(data, params, model, l)` writes the new `Xr`
   contribution of effect l back into the running fit.

ELBO is the sum of `Eloglik(data, model)` and `-sum(model$KL)`
(`R/generic_methods.R:171-228`). Convergence (`R/model_methods.R:135-230`) is
either `0 <= ELBO_diff < tol` (default `convergence_method = "elbo"`,
`tol = 1e-4`) or `max|alpha - alpha_prev| < tol` for `convergence_method = "pip"`.

Order of operations within ONE outer iteration:

1. Snapshot trace (if `track_fit = TRUE`).
2. Sweep through l = 1..L: residual, SER, fitted update.
3. Validate prior variance after the sweep
   (`R/iterative_bayesian_stepwise_selection.R:217`).
4. Compute objective, store it.
5. Check convergence; on stop, break before variance update.
6. Update `sigma2` (and other variance components for downstream classes).

## 3. User-facing parameters

The signature at `R/susie.R:304-338`:

- Data: `X`, `y`, `intercept = TRUE`, `na.rm = FALSE`.
- Model size: `L = min(10, ncol(X))`.
- Prior on effect sizes: `scaled_prior_variance = 0.2`, `prior_weights = NULL`,
  `null_weight = 0`, `prior_variance_grid = NULL`, `mixture_weights = NULL`.
- Residual variance: `residual_variance = NULL`, `estimate_residual_variance = TRUE`,
  `estimate_residual_method = "MoM"`, `residual_variance_lowerbound`,
  `residual_variance_upperbound = Inf`.
- Prior variance estimation: `estimate_prior_variance = TRUE`,
  `estimate_prior_method = "optim"`, `check_null_threshold = 0`,
  `prior_tol = 1e-9`.
- Convergence: `max_iter = 100`, `tol = 1e-4`, `convergence_method = "elbo"`,
  `verbose = FALSE`.
- Standardization and post-processing: `standardize = TRUE`,
  `compute_univariate_zscore = FALSE`, `coverage = 0.95`, `min_abs_corr = 0.5`,
  `n_purity = 100`, `refine = FALSE`, `track_fit = FALSE`,
  `unmappable_effects = "none"`, `model_init = NULL`, `init_only = FALSE`,
  `slot_prior = NULL`, `alpha0`, `beta0`.

Arguments documented as advanced or experimental in roxygen on this branch:
`unmappable_effects`, `slot_prior`, `refine`, `compute_univariate_zscore`,
`alpha0`/`beta0` (NIG path), `convergence_method = "pip"` (off by default).

## 4. Exports vs internals

Selected exports (full list in `NAMESPACE:1-108`):

- Entry points: `susie`, `susie_ss`, `susie_rss`, `susie_auto`, `susie_workhorse`,
  `susie_trendfilter`, `susie_init_coef`.
- S3 methods: `predict.susie`, `coef.susie`, `summary.susie`,
  `print.summary.susie`.
- Inference accessors: `susie_get_pip`, `susie_get_cs`, `susie_get_objective`,
  `susie_get_prior_variance`, `susie_get_residual_variance`,
  `susie_get_posterior_mean`, `susie_get_posterior_sd`, `susie_get_lfsr`,
  `susie_get_niter`, `susie_get_posterior_samples`.
- Plot helpers: `susie_plot`, `susie_plot_iteration`, `susie_plot_changepoint`.
- Pluggable bits: `slot_prior_betabinom`, `slot_prior_poisson`,
  `block_coordinate_ascent`, `ibss_initialize`, `ibss_finalize`, `get_objective`.
- Utilities: `compute_suff_stat`, `univariate_regression`, `calc_z`,
  `get_cs_correlation`, `estimate_s_rss`, `is_symmetric_matrix`.
- Mr.ASH path: `mr.ash`, `mr.ash.rss`, `coef.mr.ash`, `predict.mr.ash`,
  `get.full.posterior`.

Internal helpers a downstream package will want to call (most are reachable via
`susieR:::`):

- Constructors: `individual_data_constructor`, `sufficient_stats_constructor`,
  `summary_stats_constructor` (`R/susie_constructors.R`).
- Initialization: `initialize_susie_model`, `initialize_matrices`,
  `initialize_fitted` (`R/individual_data_methods.R:42`,
  `R/model_methods.R:81`).
- SER and partial residuals: `single_effect_regression`,
  `single_effect_update`, `compute_residuals`, `compute_ser_statistics`,
  `optimize_prior_variance`, `calculate_posterior_moments`, `compute_kl`,
  `loglik`, `neg_loglik` (`R/single_effect_regression.R`,
  `R/individual_data_methods.R`, `R/generic_methods.R`).
- Variance updates and ELBO: `Eloglik`, `get_objective`,
  `update_model_variance`, `update_variance_components`, `check_convergence`
  (`R/model_methods.R:109-230`, `R/generic_methods.R`).
- Finalization: `ibss_finalize` (`R/iterative_bayesian_stepwise_selection.R:329`).
- CS / PIP: `susie_get_pip`, `susie_get_cs`, `in_CS`, `get_purity`
  (`R/susie_get_functions.R`, `R/susie_utils.R`).

The S3 generics that take `data, params, model, ...` are the extension points.
`mvsusieR/refactor-s3` and `susieAnn` plug in by registering new methods on
their own data classes; mfsusieR will do the same.

## 5. Evidence pointers

- `R/susie.R:304-338` - public signature.
- `R/susie_workhorse.R:31-55` - outer IBSS loop.
- `R/iterative_bayesian_stepwise_selection.R:171-220` - per-effect sweep.
- `R/iterative_bayesian_stepwise_selection.R:329-382` - finalize.
- `R/single_effect_regression.R:17-69` - SER body.
- `R/single_effect_regression.R:196-208` - residual / SER / fitted-update wrapper.
- `R/model_methods.R:81-98` - initialize_matrices.
- `R/model_methods.R:135-230` - convergence check.
- `R/model_methods.R:109` - update_model_variance dispatcher.
- `R/individual_data_methods.R:104-160` - residual + SER stats + posterior moments.
- `R/generic_methods.R:171-228` - KL, Eloglik, loglik.
- `R/susie_get_functions.R:277-535` - CS extraction and PIP.
- `NAMESPACE:1-108` - full export list.

## Implications for mfsusieR

The S3 dispatch on `(data, params, model)` is the seam mfsusieR will hook into.
A multi-functional data class can implement `compute_residuals`,
`compute_ser_statistics`, `calculate_posterior_moments`, `compute_kl`,
`loglik`, `get_objective`, `update_variance_components`, and
`update_model_variance` for wavelet-coefficient sufficient statistics, and
inherit the entire IBSS scaffolding. The fit object should keep the standard
`alpha`, `mu`, `mu2`, `KL`, `lbf`, `lbf_variable`, `pi`, `sigma2`, `elbo`,
`niter`, `converged`, `sets`, and `pip` so that `predict.susie`, `coef.susie`,
and the `susie_get_*` accessors continue to work where their semantics still
hold. Quantities specific to functional traits (per-scale variance, per-trait
intercepts, reconstructed curves, per-trait residual covariance) become
additional fields without breaking the base contract. CS construction and PIP
formulas are the points where the FDR investigation in Phase 5 will likely
need to override `susie_get_cs` and `susie_get_pip` rather than reuse them
unmodified.
