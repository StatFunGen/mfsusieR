# Expose `estimate_residual_method = "NIG"` for small-sample correction

## Why

The default residual-variance estimator (Method-of-Moments)
underestimates `sigma^2` when the sample size is small,
inflating per-variable Bayes factors and posterior inclusion
probabilities. `susieR` already implements a Normal-Inverse-
Gamma posterior on `sigma^2` selected by
`estimate_residual_method = "NIG"`; the SER step integrates
over `sigma^2` instead of conditioning on a point estimate.
The NIG path is the recommended small-sample correction in
the SuSiE backbone.

`mfsusie()` does not currently expose `estimate_residual_method`.
The legacy functional fine-mapping pathway has a different
small-sample correction (`cor_small = TRUE`) that switches
the per-variable Bayes factor from a Normal marginal
likelihood to a Student's t marginal likelihood. The two
mechanisms are not equivalent; we adopt the susieR NIG
because (i) it is the actively maintained correction in our
backbone, (ii) it integrates `sigma^2` at the SER stage
rather than per-variable, and (iii) it is what
`susie_workhorse` already supports.

## What changes

### 1. Investigation step (task 0)

Audit whether the mfsusieR S3 overrides for the SER step
correctly thread `model$rv` (the NIG residual-variance
posterior parameters) when `params$use_NIG = TRUE`. Read:

- `R/individual_data_methods.R::single_effect_regression.mf_individual`
- `R/individual_data_methods.R::update_model_variance.mf_individual`
- `R/individual_data_methods.R::Eloglik.mf_individual`

against the corresponding susieR defaults. The audit
outcome shapes the rest of the change:

- If the overrides correctly forward `model$rv` and
  `(alpha0, beta0, tau)`, the change is a plumb-through:
  add the public argument and map to `params$use_NIG`.
- If the overrides ignore `model$rv` (likely, given the
  per-(outcome, scale) prior structure adds a dimension
  not present in susieR's scalar SER), extend the
  overrides to honor the NIG path on each
  (outcome, scale) pair. The extension mirrors
  susieR's pattern.

### 2. Public argument

```
mfsusie(X, Y, pos = NULL, L = ...,
        estimate_residual_method = c("MoM", "MLE", "NIG"),
        ...)
```

The default matches susieR (`"MoM"`). Forwarded through
`fsusie()`. Mapped to `params$use_NIG <-
(estimate_residual_method == "NIG")` plus
`params$alpha0`, `params$beta0` per susieR's convention at
the workhorse boundary.

### 3. Tests

`tests/testthat/test_nig_residual_variance.R` (new):

- Bit-identity vs `susieR::susie(..., estimate_residual_method
  = "NIG")` for the M=1, T_1=1 case at tolerance `<= 1e-12`.
- Structural test on n = 80 multi-outcome data: NIG yields
  a tighter `sigma2` and smaller per-variable PIPs at the
  null variants compared to the MoM default. The test
  asserts the direction of the difference (NIG `sigma2 >=
  MoM sigma2` and NIG PIP at known-null variants `<= MoM
  PIP at the same variants`).

### 4. Refactor-exceptions ledger

Entry for `mvf.susie.alpha::multfsusie(... cor_small = TRUE)`:
"replaced-by-`susieR::estimate_residual_method = \"NIG\"`
exposed via `mfsusie(..., estimate_residual_method = ...)`.
The legacy `cor_small` switch implements a Wakefield-vs-
Johnson Bayes-factor toggle (Student's t marginal likelihood
per variable). The susieR NIG path implements a Normal-
Inverse-Gamma posterior on the residual variance integrated
at the SER stage. Both target small-sample correction; we
adopt NIG because it is the actively maintained backbone
mechanism."

## Impact

- Investigation step has uncertain duration; outcome shapes
  the implementation surface (plumb-only vs override fix).
- New: `tests/testthat/test_nig_residual_variance.R`,
  refactor-exceptions entry.
- Changed: `R/mfsusie.R`, `R/fsusie.R` signatures; optionally
  `R/individual_data_methods.R` SER-step overrides.
- DESCRIPTION: no changes (susieR already imported).
- Specs: new
  `inst/openspec/specs/mf-residual-variance-method/spec.md`.

## Out of scope

- Porting `cor_small` (Student's t Bayes-factor variant).
  Documented as the alternative we chose not to adopt.
- Vignette demonstration; that lives in the separate
  `add-practical-data-applications-vignette` change which
  showcases NIG on real ATAC-seq data.
