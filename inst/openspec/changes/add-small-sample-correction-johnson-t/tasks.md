# Tasks

## 0. Audit (Stage 4a investigation, complete)

- [x] 0.1 Audit the mfsusieR S3 SER overrides against the susieR
      NIG hooks. Outcome: every NIG hook in susieR is shadowed by
      an `mf_individual` override; plumb-through is impossible.
- [x] 0.2 Identify the legacy `cor_small` Johnson-t kernel in
      `mvf.susie.alpha/R/computational_routine.R::log_BFu` and
      `fsusieR/R/computational_functions.R::log_BF.mixture_normal_per_scale`.
- [x] 0.3 Decide redirect to Johnson-t (this proposal).

## 1. Dependency

- [ ] 1.1 Add `LaplacesDemon` to `DESCRIPTION` Imports.
- [ ] 1.2 Add `r-laplacesdemon = "*"` to `pixi.toml` core
      dependencies.
- [ ] 1.3 Add `@importFrom LaplacesDemon dstp` to
      `R/mfsusieR-package.R`.

## 2. Kernel

- [ ] 2.1 Add `mixture_log_bf_per_scale_johnson` in
      `R/posterior_mixture.R`. Mirror the Normal kernel
      structure but use `LaplacesDemon::dstp` per mixture
      component.

## 3. Public API

- [ ] 3.1 Add `small_sample_correction = FALSE` to `mfsusie()`
      signature with roxygen describing the BF-correction
      contract.
- [ ] 3.2 Forward through `fsusie()` via `...`.
- [ ] 3.3 Validate input (length-1 logical).
- [ ] 3.4 Pass `params$small_sample_correction` and
      `params$small_sample_df = data$n - 1` when the flag is
      set.

## 4. Loglik branch

- [ ] 4.1 In `loglik.mf_individual`, branch on
      `params$small_sample_correction` and dispatch to
      `mixture_log_bf_per_scale_johnson` when set.

## 5. Tests

- [ ] 5.1 `tests/testthat/test_small_sample_correction.R`
      (new). Default-FALSE bit-identity vs default call.
- [ ] 5.2 Johnson-t runs to convergence on small `n`.
- [ ] 5.3 Structural correction: signal recovered, aggregate
      null PIP not larger than Wakefield.
- [ ] 5.4 Argument validation.
- [ ] 5.5 Fidelity test: per-variable LBFs from the Johnson
      kernel summed across scales match `fsusieR::log_BF` with
      `df = n - 1` at `tolerance <= 1e-12`.

## 6. Refactor-exceptions ledger

- [ ] 6.1 Add the entry mapping `cor_small = TRUE` ->
      `small_sample_correction = TRUE` with the rationale for
      the rename and the documented decision not to port susieR
      NIG.

## 7. Spec delta

- [ ] 7.1 Update
      `inst/openspec/changes/add-small-sample-correction-johnson-t/specs/mf-small-sample-correction/spec.md`
      with the new requirements.

## 8. Build + archive

- [ ] 8.1 `devtools::test()` passes.
- [ ] 8.2 `devtools::document()`.
- [ ] 8.3 Push, CI green.
- [ ] 8.4 `openspec archive add-small-sample-correction-johnson-t`.
