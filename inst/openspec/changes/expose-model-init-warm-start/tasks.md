# Tasks

## 1. Public API

- [ ] 1.1 Add `model_init = NULL` to `mfsusie()` signature
      with roxygen describing the warm-start contract.
- [ ] 1.2 Forward `model_init` through `fsusie()`.
- [ ] 1.3 Pass `params$model_init <- model_init` in the
      `mfsusie()` body's params assembly.

## 2. Override fix

- [ ] 2.1 Update `ibss_initialize.mf_individual` in
      `R/ibss_methods.R`: when `params$model_init` is non-
      `NULL`, copy `alpha`, `mu`, `mu2`, `KL`, and `sigma2`
      from the fit into the new model. Mirror
      `ibss_initialize.default` from susieR.
- [ ] 2.2 Validate L compatibility: if
      `nrow(model_init$alpha) != params$L`, error with
      clear message identifying both values.

## 3. Tests

- [ ] 3.1 `tests/testthat/test_warm_start.R` (new).
      Convergence-shortcut test: assert `niter_warm <= 1`
      vs `niter_cold >= 4` from the same fixture.
- [ ] 3.2 Bit-identity test against `susieR::susie(...,
      s_init = m_prev)` for the M=1, T_1=1 case.
- [ ] 3.3 Reference test against
      `mvf.susie.alpha::multfsusie(...,
      multfsusie.obj = m_prev)`. Tolerance per the existing
      C3 contract.
- [ ] 3.4 Error-path test: `model_init` with mismatched L
      errors cleanly.

## 4. Vignette

- [ ] 4.1 Drop the "planned extension" caveat in
      `vignettes/mfsusie_long_running_fits.Rmd`.
- [ ] 4.2 Add a worked example: cap `max_iter = 4`, save the
      partial fit, resume via `model_init`, demonstrate
      <= 2 additional iterations to convergence.

## 5. Spec delta

- [ ] 5.1 Create
      `inst/openspec/changes/expose-model-init-warm-start/specs/mf-warm-start/spec.md`
      with requirements and scenarios.

## 6. Refactor-exceptions ledger

- [ ] 6.1 Add entry mapping `multfsusie(... multfsusie.obj
      = m_prev)` -> `mfsusie(..., model_init = m_prev)`
      with the rationale for retiring the legacy
      `max_step` workaround.

## 7. Build + archive

- [ ] 7.1 `devtools::test()` passes including the new tests.
- [ ] 7.2 `devtools::document()` regenerates man + NAMESPACE.
- [ ] 7.3 Push, CI green.
- [ ] 7.4 `openspec archive expose-model-init-warm-start`.
