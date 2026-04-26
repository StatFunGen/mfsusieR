# Tasks

## 1. Public API

- [x] 1.1 Rename `cross_modality_prior` argument in `mfsusie()` to
      `cross_outcome_prior`.
- [x] 1.2 Rename `prior_variance_scope` value vocabulary to
      `c("per_outcome", "per_scale")`; flip default to first value.
- [x] 1.3 Rename `residual_variance_method` -> `residual_variance_scope`;
      values `c("per_outcome", "per_scale")`; flip default.
- [x] 1.4 Update `fsusie()` docstring (passes `...` to `mfsusie()`).
- [x] 1.5 Rename `nullweight` -> `mixsqp_null_penalty` in
      `mfsusie()` argument and `params` forwarding.

## 2. Internal data + model fields

- [x] 2.1 `data$predictor_weights` -> `data$xtx_diag` everywhere.
- [x] 2.2 `data$csd_X` -> `data$csd`.
- [x] 2.3 `data$T_padded` -> `data$T_basis`.
- [x] 2.4 `model$cross_modality_combiner` -> `model$cross_outcome_combiner`.
- [x] 2.5 `fit$mf_meta` -> `fit$dwt_meta`.

## 3. Internal helpers + S3 generics

- [x] 3.1 `combine_modality_lbfs.*` -> `combine_outcome_lbfs.*` (S3 generic).
- [x] 3.2 `cross_modality_prior_independent()` -> `cross_outcome_prior_independent()`.
- [x] 3.3 `mf_per_modality_bhat_shat` -> `mf_per_outcome_bhat_shat`.
- [x] 3.4 `mf_invert_per_modality` -> `mf_invert_per_outcome`.
- [x] 3.5 Update `R/zzz.R::.onLoad` S3-method-registration list.
- [x] 3.6 Rename file `R/prior_cross_modality.R` -> `R/prior_cross_outcome.R`.

## 4. Prose audit (roxygen, comments)

- [x] 4.1 Replace "SNP" / "SNPs" -> "variable" / "variables" in R/.
- [x] 4.2 Replace "modality" / "modalities" -> "outcome" / "outcomes"
      in R/.
- [x] 4.3 Replace "trait" / "traits" -> "outcome" / "outcomes" in R/
      and vignettes.
- [x] 4.4 Apply same renames to `src/` (cpp11 file headers).

## 5. Tests

- [x] 5.1 Update tests that set legacy value names.
- [x] 5.2 Update tests that read renamed data / model fields.
- [x] 5.3 Update tests that call renamed internal helpers.
- [x] 5.4 Verify default-flip implications on C1/C2/C3 contracts.
- [x] 5.5 Update `tests/testthat/scripts/make_fixtures.R` and
      regenerate `scenario_minimal.rds`.
- [x] 5.6 Skip the first ELBO diff in the smoke monotonicity test:
      the initial sigma2 (var(Y) guess) gets replaced by its first
      closed-form update, which is a one-shot non-monotone step.
      Mirrors mvsusieR pattern.

## 6. Vignettes

- [x] 6.1 Apply prose rename across all vignettes.
- [x] 6.2 Update explicit `prior_variance_scope` /
      `residual_variance_scope` calls (only one site:
      getting_started.Rmd's C1 degeneracy block).

## 7. Modularity audit guard

- [x] 7.1 D12 forbidden-string list unchanged; modularity guard
      passes after the rename.

## 8. Build

- [x] 8.1 `devtools::document()` to refresh man / NAMESPACE.
- [x] 8.2 `devtools::test()` passes.
- [x] 8.3 `R CMD check --no-manual` returns Status: OK.
