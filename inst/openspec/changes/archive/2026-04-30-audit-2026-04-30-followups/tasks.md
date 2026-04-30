# Tasks

Severity: LOW. Five independent polish items from the 2026-04-30
cross-package audit. Each task is a standalone PR; ordering does
not matter.

## 1. A-2. Warm-start validator

- [x] 1.1 Add `validate_init.mf_individual` (or extend
      `expand_model_init_to_L`) to NA/Inf-check `alpha`, `mu`,
      `mu2`, `pi_V`, `fitted_g_per_effect`, `sigma2`, `V` and
      enforce list-of-list shape sanity.
- [x] 1.2 New unit test exercising the validator on a corrupted
      `model_init` fixture.

## 2. A-3. Post-hook S3 dispatch cleanup

- [x] 2.1 In `post_loglik_prior_hook.mf_individual`
      (`R/individual_data_methods.R:701-714`) replace the direct
      `optimize_prior_variance.mf_individual(...)` call with the
      dispatched `optimize_prior_variance(...)` form.
- [x] 2.2 Drop the unused `moments = get_post(model, l)`
      argument.
- [x] 2.3 Confirm the existing test suite still passes.

## 3. B-5. CS purity / coverage parity smoke test

- [x] 3.1 New `tests/testthat/test_cs_parity_fsusier.R` that
      fits both `mfsusieR::fsusie` and `fsusieR::susiF` on the
      same toy fixture, normalises the per-CS structures, and
      asserts equality at the appropriate tolerance.

## 4. C-7. Combiner default method

- [x] 4.1 Add `combine_outcome_lbfs.default` in
      `R/prior_cross_outcome.R` that errors with a message naming
      the registered combiner classes (introspect via
      `methods("combine_outcome_lbfs")`).
- [x] 4.2 New unit test asserting the helpful error message.

## 5. C-8. L_greedy expansion regression test

- [x] 5.1 New `tests/testthat/test_l_greedy_expansion.R` fitting
      with `L = 5, L_greedy = 5` then `L = 10, L_greedy = 0` on
      the same seed and asserting the alpha/mu/mu2/pip agree at
      `1e-8`.

## 6. Validate

- [x] 6.1 `openspec validate audit-2026-04-30-followups` from
      `inst/`.
- [x] 6.2 Full `devtools::test()` clean.
