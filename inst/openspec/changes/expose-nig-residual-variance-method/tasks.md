# Tasks

## 0. Investigation (gates the rest)

- [ ] 0.1 Read `R/individual_data_methods.R::single_effect_regression.mf_individual`,
      `update_model_variance.mf_individual`,
      `Eloglik.mf_individual` against the corresponding
      defaults in susieR.
- [ ] 0.2 Determine whether `model$rv` and the
      `(alpha0, beta0, tau)` tuple are forwarded by the
      mfsusieR overrides when `params$use_NIG = TRUE`.
- [ ] 0.3 Decide implementation surface based on 0.2:
        outcome A (overrides forward correctly):
          plumb-through only.
        outcome B (overrides ignore the NIG path):
          extend the overrides per task 2 below.
- [ ] 0.4 Write the audit findings into a short report at
      `inst/notes/sessions/<date>-nig-audit.md`.

## 1. Public API

- [ ] 1.1 Add
      `estimate_residual_method = c("MoM", "MLE", "NIG")`
      to `mfsusie()` signature with roxygen.
- [ ] 1.2 Forward through `fsusie()`.
- [ ] 1.3 In `mfsusie()` body, set
      `params$use_NIG <- (estimate_residual_method == "NIG")`
      plus `params$alpha0` and `params$beta0` per susieR's
      convention.

## 2. Override extensions (only if outcome B)

- [ ] 2.1 Update `single_effect_regression.mf_individual` to
      thread `model$rv` per (outcome, scale).
- [ ] 2.2 Update `update_model_variance.mf_individual` to
      maintain the per-(outcome, scale) NIG posterior
      `(alpha0_t, beta0_t)`.
- [ ] 2.3 Update `Eloglik.mf_individual` to integrate over
      `sigma^2` under the NIG posterior.

## 3. Tests

- [ ] 3.1 `tests/testthat/test_nig_residual_variance.R`
      (new). Bit-identity vs
      `susieR::susie(..., estimate_residual_method = "NIG")`
      for the M=1, T_1=1 case at tolerance `<= 1e-12`.
- [ ] 3.2 Structural test on n = 80 multi-outcome data:
      NIG `sigma2` >= MoM `sigma2`; NIG PIP at known-null
      variants <= MoM PIP at the same variants.

## 4. Refactor-exceptions ledger

- [ ] 4.1 Add the entry mapping `cor_small = TRUE` ->
      `estimate_residual_method = "NIG"` with the
      Wakefield-vs-Johnson rationale.

## 5. Spec delta

- [ ] 5.1 Create
      `inst/openspec/changes/expose-nig-residual-variance-method/specs/mf-residual-variance-method/spec.md`.

## 6. Build + archive

- [ ] 6.1 `devtools::test()` passes.
- [ ] 6.2 `devtools::document()`.
- [ ] 6.3 Push, CI green.
- [ ] 6.4 `openspec archive expose-nig-residual-variance-method`.
