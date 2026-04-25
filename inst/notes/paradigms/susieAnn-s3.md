# susieAnn S3 paradigm

This is paradigm reference #2 for mfsusieR: how to extend `susieR` by adding a
new prior structure (functional annotations driving inclusion through a
Gamma-Poisson activity indicator) and a pluggable predictor interface, without
extending the susie fit-object class hierarchy.

Source root: `~/GIT/susieAnn` on branch `main`. All citations are relative to
that root.

## 1. Package layout

- `susie_ann.R` - public `susie_ann()` entry and the outer EM loop.
- `estep.R` - per-block IBSS step, parallelized via `future_lapply`. Builds
  the slot-prior and `prior_weights` and calls `susieR::susie_rss` for each
  block.
- `mstep.R` - mixture-weight update (mixsqp), predictor refit, nu update,
  intercept correction.
- `annotation_predictor.R` - the predictor S3 system: generic constructors,
  generics, four built-in implementations (logistic, glmnet, lightgbm,
  catboost), and a custom escape hatch.
- `posterior.R` - result-object construction, CS collection, `print.susie_ann`.
- `utils.R` - prior-variance grid (via ashr), block indexing, omega
  initialization, polygenicity init, nu init, global PIP aggregation, ELBO.
- `ld_loader.R` - LD provider wrapper around `pecotmr::ld_loader`.
- `zzz.R` - caches `susieR:::warning_message` and
  `susieR:::slot_prior_poisson` into the package namespace.

No file mirrors a susieR file. Public entry: `susie_ann()`
(`R/susie_ann.R:86-316`).

## 2. Class hierarchy (S3)

The fit object carries one new class: `susie_ann`
(`R/posterior.R:23-64`, set at `R/posterior.R:62`). The susie objects produced
per block are kept inside `fits` unmodified, so per-block postprocessing with
`susieR::susie_get_pip` and friends still works.

The interesting hierarchy is on the predictor side
(`R/annotation_predictor.R`):

- Base class: `susie_ann_predictor`.
- Subclasses (also as S3): `predictor_auto`, `predictor_logistic`,
  `predictor_glmnet`, `predictor_lightgbm`, `predictor_catboost`,
  `predictor_custom`. Each is a list with `type`, `calibrate`, `screen_k`,
  `temperature`, `params`; `predictor_custom` also stores `fit_fn` and
  `predict_fn`.

Generics on predictors (`R/annotation_predictor.R:149-176`):

- `predictor_fit(predictor, X, y, ...)` - returns a fitted model.
- `predictor_predict(predictor, model, X, ...)` - returns probabilities.
- `predictor_shap(predictor, model, X, ...)` - returns importance/effects.

Methods exist for each subclass: logistic
(`:456-484`), glmnet (`:493-530`), lightgbm (`:541-578`), catboost
(`:588-623`), custom (`:632-645`), auto (`:185-207`).

The shared LOCO cross-validation pipeline (`predictor_loco`,
`R/annotation_predictor.R:230-294`) is not generic - all predictors run the
same outer pipeline and only the per-fold fit and predict dispatch on class.

`print.susie_ann` is at `R/posterior.R:129-161`. There are no `predict.susie_ann`
or `coef.susie_ann` methods; downstream code is expected to look at the
per-block `fits` or the aggregated `pip` and `omega`.

## 3. Delegation to susieR

susieAnn does not reimplement IBSS. The inner loop is `susieR::susie_rss`,
called twice per E-step pass (`R/estep.R:73`, `R/estep.R:99`) using
`do.call(susie_rss, rss_args)`. Two modes:

- `prior_variance_method = "local_estimate"` (`estep.R:59-78`): annotated
  inclusion probabilities are passed in `prior_weights`, susieR estimates the
  effect-size prior variance internally.
- `prior_variance_method = "mixture"` (`estep.R:80-108`): a global prior-
  variance grid plus mixture weights is passed; the Gamma-Poisson slot prior
  produced by `susieR:::slot_prior_poisson(C = lambda_g, nu = nu, ...)` is
  also supplied. susieR's IBSS handles slot c_hat updates internally.

Caching in `zzz.R:21-32` is small and explicit:
`warning_message` and `slot_prior_poisson` from susieR. No `mashr` use.

There is no reimplementation of SER, residual variance estimation, CS
construction, or LD routing - these all happen inside `susie_rss`.

## 4. Annotation-driven prior

The activity indicator is the Gamma-Poisson slot prior provided by susieR
(`susieR:::slot_prior_poisson`). susieAnn supplies its parameters:

- `lambda_g` per block: sum of `omega` across the block's variants
  (`estep.R:50`). This is the prior expected number of causal variants.
- `nu`: scalar overdispersion, estimated each EM iteration.
- `c_hat_init`: warm-started from the previous EM iteration's posterior
  (`estep.R:88`).
- `skip_threshold_multiplier`: prunes slots whose lambda-weighted lbf is
  small (`estep.R:88`).

The slot prior is rebuilt fresh each E-step, not reused as an object.

The annotation predictor, by contrast, is reused across iterations. Each
M-step refits it on (annotation, current PIP) pairs
(`mstep.R:92-95`) by calling `predictor_loco` (`annotation_predictor.R:230-294`):

1. LOCO folds: leave-one-chromosome-out, falling back to random 5-fold if
   chromosomes are not informative (`:334-350`).
2. Temperature sharpening: PIPs raised to power `1/T` to act as soft labels
   (`:255`).
3. Threshold and weighting: training restricted to variants with sharpened
   PIP > 1e-4; observation weights set to `sharpened_pip + 1/P`
   (`:260-267`).
4. Optional annotation screening by marginal correlation with y
   (`:281-283`).
5. Per-fold `predictor_fit` then `predictor_predict` (`:284-285`).
6. Optional isotonic calibration of out-of-fold omega against y
   (`:289-290`).

`omega` is the predicted causal probability per variant. `lambda_g` for the
next E-step is the per-block sum of `omega`.

Overdispersion is estimated by `update_nu` (`mstep.R:155-172`), a Newton
optimization of the marginal log-density of the Gamma-Poisson with sufficient
statistics `a_g` and `b_g` (Gamma posterior shape and rate from the slot
prior's posterior) and `lambda_g`. Called from `susie_ann.R:257`.

Annotation effects (the interpretable enrichment side) come out through
`predictor_shap`. For logistic this is `coef(model)` (`:481-483`), for
glmnet it is the regularized coefficients (`:527-530`), for tree methods
SHAP via the `shapviz` package. The M-step does not use these directly; they
are kept on `predictor$model` for downstream inspection.

## 5. IBSS loop differences vs susieR

susieAnn does not modify the IBSS loop body; it modifies the inputs.
Per E-step iteration, for block g:

1. `lambda_g <- sum(omega[block_idx[[g]]])`.
2. `pi_g <- omega[block_idx[[g]]] / lambda_g`.
3. Pass `prior_weights = pi_g` to `susie_rss`. In mixture mode, also pass the
   prior-variance grid, mixture weights, and a fresh `slot_prior_poisson`.
4. susieR runs IBSS, returns a `susie` fit including alpha, mu, c_hat
   posterior (`a_g`, `b_g` if slot_prior was used).
5. Per-block PIP and posterior c_hat are passed back to the M-step.

State that lives outside the susie objects, at the `susie_ann` level:
`omega` (P-vector), `w` (K-vector mixture weights), `nu` (scalar), `c_hat`
(L by G), the predictor object, `lambda_g`, `Lcg_hat`, `E_mu_g`, `a_g`,
`b_g`. Warm-starting passes the previous block's susie fit through
`model_init` (`estep.R:108-ish`) so IBSS starts from the previous posterior.

Outer convergence (`R/susie_ann.R:300-303`) is on the global maximum PIP
change across all blocks and variants, default tol 1e-3. Per-block IBSS uses
its own `tol_ibss` (default from susie_rss).

## 6. Pluggable predictor interface

The predictor abstraction is an S3 object. To plug in a new predictor a user
supplies (a) `predictor_fit` and (b) `predictor_predict`; `predictor_shap` is
optional and only used downstream. The escape hatch is
`custom_annotation_predictor(fit_fn, predict_fn, ...)`
(`R/annotation_predictor.R:112-132`), which constructs a `predictor_custom`
whose `predictor_fit.predictor_custom` simply calls `fit_fn(X, y, weights)`.

The shared LOCO pipeline lives outside the dispatch (`predictor_loco` is not
a generic), so all predictors share label sharpening, weighting, screening,
and isotonic calibration. Only the per-fold model-fitting and prediction
dispatch on class.

`predictor_auto` (`:185-207`) selects logistic, glmnet, or catboost/lightgbm
by annotation dimension `D` (`:88-100`). This is implemented as an S3 method
on `predictor_auto` that resolves the concrete class and then dispatches
back into the same generic.

## 7. Public API signature

```
susie_ann(z_list, R_loader, n,
          annotations = NULL, chrom = NULL,
          L = 30,
          pi_causal = "ashr",
          prior_variance_method = "mixture",
          nu_init = "auto",
          predictor = annotation_predictor("auto"),
          large_predictor = c("catboost", "lightgbm"),
          max_iter = 100, max_iter_local = 1,
          recheck_interval = 5L, skip_threshold_multiplier = 1.5,
          tol = 1e-3, warm_start = TRUE, verbose = TRUE, ...)
```

(`R/susie_ann.R:86-104`.) Diff vs `susieR::susie()`: takes summary
statistics across many blocks rather than X and Y. New arguments:
`annotations`, `chrom`, `pi_causal`, `prior_variance_method`, `nu_init`,
`predictor`, `large_predictor`, `max_iter_local`, `recheck_interval`,
`skip_threshold_multiplier`. Removes individual-data inputs in favor of a
list of z-vectors and an LD loader.

## 8. Return object

`susie_ann` (`R/posterior.R:23-64`):

- `fits` - list of G susie fits, each a normal `susie` object.
- `pip`, `cs` - global posterior inclusion probabilities and credible sets.
- `omega` - P-vector of annotation-derived causal probabilities.
- `w` - K-vector mixture weights (NULL in local-estimate mode).
- `nu` - scalar overdispersion (NULL in local-estimate mode).
- `grid` - K-vector of prior variance grid (NULL in local-estimate mode).
- `predictor` - fitted predictor (NULL if no annotations).
- `c_hat` - L by G slot-activity posteriors.
- `block_idx` - mapping from blocks to global variant indices.
- `lambda_g`, `Lcg_hat`, `E_mu_g`, `a_g`, `b_g` - per-block diagnostics.
- `polygenicity` = mean(omega).
- `elbo_trace`, `pip_trace` - convergence traces.

The per-block susie objects are unchanged. There is no class extension at
the fit level.

## 9. Contrast with mvsusieR/refactor-s3

Both packages extend susieR by registering S3 methods, but at different
points:

- mvsusieR extends the *response* structure: new data classes
  (`mv_individual`, `mv_ss`), new model class (`mvsusie`), full method
  override surface (~25 generics on the data class, 6 on the model). The
  IBSS loop is the same, but every per-iteration computation is overridden.
- susieAnn extends the *prior* structure without extending the fit class:
  no `mv_*` style data class, no `susie_ann` model methods, no override of
  the IBSS internals. Annotation enrichment is wired in by passing
  `prior_weights` and a `slot_prior_poisson` object into the public
  `susie_rss` interface, and by holding global state at the EM level.

The prior reuse pattern also differs. mvsusieR's `mash_prior` is a long-lived
object that the SER reads each iteration. susieAnn's slot_prior_poisson is
rebuilt every E-step. Both reuse the predictor-side state across iterations
(mvsusieR's `pi_V` mixture weights, susieAnn's `predictor` model and
`omega`).

## 10. Evidence pointers

- `R/susie_ann.R:86-316` - public entry and EM loop.
- `R/estep.R:38-128` - per-block IBSS with annotation prior.
- `R/estep.R:84-89` - slot prior assembly.
- `R/annotation_predictor.R:149-176` - predictor S3 generics.
- `R/annotation_predictor.R:230-294` - shared LOCO pipeline.
- `R/annotation_predictor.R:456-484` - logistic implementation.
- `R/annotation_predictor.R:88-100` - auto-predictor selection.
- `R/mstep.R:18-172` - mixture weights, predictor refit, nu update.
- `R/mstep.R:111-136` - intercept correction.
- `R/posterior.R:23-64` - `build_susie_ann_result` and class assignment.
- `R/utils.R:17-204` - prior grid, polygenicity, omega init, global PIP.
- `R/zzz.R:18-32` - cached susieR internals.

## Implications for mfsusieR

susieAnn shows that an extension is not obligated to extend the fit class.
For mfsusieR this matters in two places. First, anything in the
multi-functional pipeline that is naturally a *post-processing* step on an
otherwise-vanilla susie run (per-trait CS aggregation, smoothing,
transcript-level summaries) can stay as outer-loop code that wraps susie
fits, the way susieAnn wraps per-block fits. Second, the predictor S3
pattern is the right shape for any pluggable component mfsusieR ends up
needing (e.g., a wavelet-basis adapter or an alternative scale-prior
estimator), since it gives users a clean way to swap implementations
without subclassing the fit object. Where mfsusieR does need to override
per-iteration math (residual covariance across positions, scale-specific
prior variance), the mvsusieR pattern is the better fit; the two paradigms
are not in conflict and mfsusieR can use both.
