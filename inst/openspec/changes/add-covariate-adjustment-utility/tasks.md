# Tasks

## 1. Helper refactor (R/em_helpers.R)

- [ ] 1.1 Add `.compute_per_scale_likelihood_matrix(Bhat, Shat,
      prior, indx_lst, is_ebmvfr = FALSE)`. Default
      (`is_ebmvfr = FALSE`) prepends the SuSiE null-component
      penalty row (existing behavior). `is_ebmvfr = TRUE`
      omits it. Refactor the existing SuSiE call site to
      route through this helper.
- [ ] 1.2 Add `.m_step_per_scale(L, w, indx_lst, init_pi0_w,
      control_mixsqp, null_prior_weight, is_ebmvfr = FALSE)`.
      Default behavior preserves the SuSiE M-step's
      `null_prior_weight` weight prepend; `is_ebmvfr = TRUE`
      omits it. Refactor the existing SuSiE call site to
      route through this helper.
- [ ] 1.3 Smoke tests in `tests/testthat/test_em_helpers.R`
      exercising both `is_ebmvfr = TRUE` and
      `is_ebmvfr = FALSE` and asserting the structural
      differences (penalty-row count; nullweight prepend).

## 2. New file R/adjust_covariates.R

- [ ] 2.1 `mf_residualize_ols(Y, Z, X = NULL)`. Closed-form
      Frisch-Waugh-Lovell. Returns
      `list(Y_adjusted, X_adjusted, fitted_func)`. ~10 lines
      base R.
- [ ] 2.2 Internal `.mf_get_er2_ebmvfr_adjust(Y_wd, Z, fitted_wc,
      fitted_wc2)`. Implements the upstream-buggy formula
      `sum((Y - Z*postF)^2) - sum(postF^2) + sum(postF2)` and
      normalizes by `1/(n*T_padded)`. Pattern B kept-upstream-
      bug. ~5 lines.
- [ ] 2.3 Internal `.mf_residualize_wavelet_eb_inner(Y_wd, Z,
      prior, control)`. The Gauss-Seidel outer loop honoring
      every gotcha from `inst/notes/refactor-exceptions.md`'s
      EBmvFR-adjust entry: per-covariate partial residual,
      `compute_marginal_bhat_shat`, `pmax(Shat, 1e-32)`,
      `ashr::postmean`/`postsd`, `.mf_get_er2_ebmvfr_adjust`,
      `.compute_per_scale_likelihood_matrix(is_ebmvfr = TRUE)`,
      `.m_step_per_scale(is_ebmvfr = TRUE)`,
      convergence check `||pi - pi_prev||/log(K) < 1e-3`.
      ~80 lines.
- [ ] 2.4 `mf_residualize_wavelet_eb(Y, Z, X = NULL,
      wavelet_filter_number, wavelet_family, max_iter, tol,
      null_prior_weight, init_pi0_w, control_mixsqp,
      grid_mult)`. Wraps the inner: centers Y, scales Z,
      runs DWT, calls `init_scale_mixture_prior_default`,
      runs the Gauss-Seidel inner, then per-covariate inverse
      DWT to produce `fitted_func`, then
      `Y_adjusted = Y - Z %*% fitted_func`. Optionally
      `X_adjusted = (I - H_Z) X`. ~50 lines.
- [ ] 2.5 Public `mf_adjust_for_covariates(Y, Z, X = NULL,
      method = c("wavelet_eb", "ols"), ...)`. Dispatches.
      Errors on `thresh_lowcount > 0` and
      `quantile_trans = TRUE`. ~20 lines + roxygen.
- [ ] 2.6 Roxygen + manuscript citation in
      `R/adjust_covariates.R` per writing-style policy. NO
      port-source mentions in user-facing roxygen.

## 3. Refactor-exceptions ledger

- [ ] 3.1 Add carve-out entry: EBmvFR adjust path is in scope
      as of round 2; full EBmvFR model remains out of scope.
      Cite the helper-reuse map.
- [ ] 3.2 Add Pattern B entry for the kept-upstream `get_ER2`
      bug. Cite `operation_on_EBmvFR_obj.R:124-132` and
      Denault's TODO at line 74. Note that `sigma2` is dead-
      end in the EBmvFR loop so the bug only affects the
      stored sigma2 field.

## 4. Vignette work

- [ ] 4.1 Rename `vignettes/fsusie_covariates_and_coloc.Rmd` to
      `vignettes/fsusie_colocalization.Rmd`. Drop the
      covariates section. Update the intro paragraph to
      reflect the narrowed scope. Update the
      `\VignetteIndexEntry` line.
- [ ] 4.2 New `vignettes/fsusie_covariates_adjustment.Rmd`.
      Author: William Denault. Sections per design.md D7.
      Uses `mf_adjust_for_covariates()` and `mfsusie_plot()`.
      No port-source mentions.

## 5. Tests

- [ ] 5.1 `tests/testthat/test_adjust_covariates.R`: apple-to-
      apple bit-identity vs `fsusieR:::EBmvFR(adjust = TRUE)`
      across `n in {100, 200}`, `T in {64, 128}`,
      `K in {2, 3, 5}` at tolerance `<= 1e-12`. Gated
      `skip_if_not_installed("fsusieR")`.
- [ ] 5.2 OLS closed-form sanity: residual-orthogonality
      (`max(abs(t(Z) %*% Y_adj)) < 1e-12`), idempotency
      (`adjust(adjust(Y))$Y_adjusted == adjust(Y)$Y_adjusted`).
- [ ] 5.3 Error-path tests: `thresh_lowcount > 0` errors;
      `quantile_trans = TRUE` errors; bad shapes error;
      method dispatch.

## 6. Public API + Suggests

- [ ] 6.1 `NAMESPACE`: export `mf_adjust_for_covariates`,
      `mf_residualize_ols`, `mf_residualize_wavelet_eb`.
- [ ] 6.2 `DESCRIPTION` Suggests: add `fsusieR` (test-only).
- [ ] 6.3 No `pixi.toml` change (fsusieR installed via the
      existing `r-remotes` machinery in CI).

## 7. Spec delta

- [ ] 7.1 Create `inst/openspec/changes/add-covariate-adjustment-utility/specs/mf-covariate-adjustment/spec.md`.
      Requirements + scenarios for:
      - `mf_adjust_for_covariates` signature, return shape.
      - `wavelet_eb` bit-identity contract at upstream defaults.
      - `ols` closed-form contract.
      - Error contract for unsupported preprocessing flags.
      - X-argument FWL extension (returns `X_adjusted` when X
        is supplied).

## 8. Excalidraw design diagram

- [ ] 8.1 `design/covariate-adjustment-flow.excalidraw`:
      data-flow diagram showing the pipeline from `(Y, Z)`
      through DWT, EB regression per scale, mixture-weights
      EM, inverse DWT, to `Y_adjusted`. Reference helpers and
      pre-existing primitives by file:function.

## 9. Build + tests + archive

- [ ] 9.1 `devtools::document()` regenerates man + NAMESPACE.
- [ ] 9.2 `devtools::test()` passes (smoke + apple-to-apple).
- [ ] 9.3 Vignettes render cleanly under pixi env.
- [ ] 9.4 Push, let CI run.
- [ ] 9.5 `openspec archive add-covariate-adjustment-utility`
      once green.
