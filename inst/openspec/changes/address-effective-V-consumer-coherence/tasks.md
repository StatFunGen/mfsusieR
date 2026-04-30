# Tasks

Severity: MEDIUM. Should land before the next user-facing
release; not blocking unrelated feature work.

## 1. Code

- [x] 1.1 Decide between (a) deleting
      `trim_null_effects.mf_individual` so the susieR `.default`
      fires, or (b) overriding with mfsusieR-shaped zeroing of
      `alpha[l, ]`, `mu[[l]]`, `mu2[[l]]`, `lbf[l]`, `KL[l]`
      for `V[l] < params$prior_tol`. Read
      `susieR/R/model_methods.R:339-352` and confirm dispatch
      shape compatibility before choosing (a).
- [x] 1.2 Add the chosen implementation. Refresh the docstring at
      `R/ibss_methods.R:228-240` to remove the "V[l] = 1 always"
      claim.
- [x] 1.3 Add a one-line comment to
      `pre_loglik_prior_hook.mf_individual`
      (`R/individual_data_methods.R:677-689`) noting the V_init
      discard and the V=1 semantic for `loglik.mf_individual`.

## 2. Documentation

- [x] 2.1 Add a paragraph to `mfsusie()` return-value roxygen for
      `pip`, naming the V-based filter and the `prior_tol`
      argument.
- [x] 2.2 NEWS.md entry under the next release header.

## 3. Test

- [x] 3.1 New regression test (`tests/testthat/test_v_filter_pip.R`)
      on a null-only fixture asserting the V-filter behavior is
      deterministic and documented.

## 4. Validate

- [x] 4.1 `openspec validate address-effective-V-consumer-coherence`
      from `inst/`.
- [x] 4.2 Full `devtools::test()` clean.
