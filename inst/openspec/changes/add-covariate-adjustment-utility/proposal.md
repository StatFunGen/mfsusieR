# Add `mf_adjust_for_covariates()` and split the covariates vignette

## Why

Covariate adjustment for functional fine-mapping requires more
than position-by-position OLS residualization. fsusieR exposes
`EBmvFR(Y, X = Z, adjust = TRUE)` for this purpose: a wavelet-
domain empirical-Bayes regression that fits per-coefficient
mixture-of-normals priors on the regression coefficients of Y
on Z, then returns `Y_adjusted = Y - Z %*% fitted_func`.

The current `vignettes/fsusie_covariates_and_coloc.Rmd`
implements covariate adjustment inline as a position-space OLS
hat-matrix, which (i) is not what `fsusieR::Adjusting_covariate.Rmd`
demonstrates, (ii) discards the wavelet-domain shrinkage that
fsusieR's recipe gives, and (iii) has no public utility a user
can call.

mfsusieR needs a public `mf_adjust_for_covariates()` that
reproduces the upstream wavelet-EB result bit-identically at
upstream defaults. The new function lives alongside
`mf_post_smooth()` as a standard pre-fit utility for users.

## What changes

### 1. Public utility `mf_adjust_for_covariates()`

```
mf_adjust_for_covariates(Y, Z, X = NULL,
                         method = c("wavelet_eb", "ols"),
                         wavelet_filter_number = 1L,
                         wavelet_family        = "DaubExPhase",
                         max_iter              = 100L,
                         tol                   = 1e-3,
                         null_prior_weight     = 10,
                         init_pi0_w            = 1,
                         control_mixsqp        = list(verbose = FALSE,
                                                      eps = 1e-6,
                                                      numiter.em = 4L),
                         grid_mult             = sqrt(2))
```

Returns `list(Y_adjusted, X_adjusted, fitted_func, sigma2,
niter, converged, method)`.

- `method = "wavelet_eb"` (default): wavelet-domain empirical-
  Bayes regression per coefficient. Bit-matches
  `fsusieR::EBmvFR(Y, X = Z, adjust = TRUE)` at upstream
  defaults at tolerance `<= 1e-12`.
- `method = "ols"`: closed-form Frisch-Waugh-Lovell
  residualization. `Y_adjusted = (I - Z(Z'Z)^{-1}Z') Y` and
  optionally `X_adjusted = (I - Z(Z'Z)^{-1}Z') X`. ~10 lines of
  base R.

The optional `X` argument lets users residualize the genotype
matrix as well (FWL), addressing covariates correlated with
genotype. Upstream's `EBmvFR(adjust = TRUE)` does not do this;
we expose it as an option since FWL is the statistically
correct treatment when Z and X are not orthogonal.

### 2. Internal helpers (in `R/em_helpers.R`, refactored)

The wavelet-EB path needs two operations our existing IBSS
machinery already does in slightly different form:

- `.compute_per_scale_likelihood_matrix(Bhat, Shat, prior,
                                        indx_lst,
                                        is_ebmvfr = FALSE)`
  Wraps the L_mixsq logic. The `is_ebmvfr` flag suppresses the
  null-component penalty row that the SuSiE M-step prepends.
- `.m_step_per_scale(L, w, indx_lst, init_pi0_w,
                     control_mixsqp, null_prior_weight,
                     is_ebmvfr = FALSE)`
  Wraps the m_step logic. The `is_ebmvfr` flag suppresses the
  `null_prior_weight` weight prepend used in the SuSiE path.

Both are added under `R/em_helpers.R` so both the existing
SuSiE prior-update path AND the new EBmvFR-adjust outer loop
exercise the same code.

### 3. Reused primitives (no edits)

- `init_scale_mixture_prior_default()` (R/prior_scale_mixture.R)
- `compute_marginal_bhat_shat()` (susieR import)
- `ashr::postmean()`, `ashr::postsd()`
- `mixsqp::mixsqp()`
- `col_scale()`, `mf_dwt_matrix()`, `mf_invert_dwt()`
  (R/utils_wavelet.R)

### 4. Refactor-exceptions ledger update

Two new entries in `inst/notes/refactor-exceptions.md`:

- **In-scope carve-out for the EBmvFR adjust path.** Hard rule
  #1 declared all of `EBmvFR.R`, `EBmvFR_workhorse.R`, and
  `operation_on_EBmvFR_obj.R` out of scope. The carve-out
  scopes the per-coefficient EB regression core needed for
  `adjust = TRUE` only; the full EB-multivariate-functional-
  regression model as a SuSiE alternative remains out of scope.
  We do not port `EBmvFR.R` verbatim; we assemble the algorithm
  from existing primitives plus the two helpers in §2 above.
- **Kept-upstream-bug entry for `get_ER2.EBmvFR`.** fsusieR's
  `operation_on_EBmvFR_obj.R:124-132` formula
  `sum((Y - X*postF)^2) - sum(postF^2) + sum(postF2)` is
  missing the `predictor_weights` factor. Denault flagged this
  himself in the TODO at line 74. We replicate the buggy
  formula in `.mf_get_er2_ebmvfr_adjust()` (Pattern B per
  refactor-discipline) so the apple-to-apple bit-identity
  contract holds. Note: `sigma2` is dead-end in EBmvFR's outer
  loop (`cal_Bhat_Shat` ignores its `resid_var` argument), so
  the bug only affects the stored `sigma2` field, not the
  returned `Y_adjusted`. We do NOT replicate this bug in the
  SuSiE / mfsusie path (`mf_get_ER2_per_position` keeps its
  corrected formula).

### 5. New vignette

`vignettes/fsusie_covariates_adjustment.Rmd` (split out from
the existing combined vignette):

- Author: William Denault (the upstream vignette's author).
- Reproduces the structure of `fsusieR/vignettes/Adjusting_covariate.Rmd`
  with our `mf_adjust_for_covariates()` in place of `EBmvFR`.
- Same simulation recipe and seed as upstream so the
  vignette's numerical output matches upstream's: same
  `Y`, `Z`, `X`, same shape of adjusted curves, same
  fine-mapping result on the adjusted response.
- After computing `Y_adjusted`, the vignette runs
  `fsusie(Y_adjusted, X, pos)` and reports PIPs / credible
  sets / effect curves. These SHALL match the equivalent
  output of upstream's vignette under the same simulation
  seed, with the only acceptable deviation being numerical
  drift driven by the `get_ER2.EBmvFR` bug we keep at
  Pattern B (which affects `sigma2` only, and is dead-end in
  the EBmvFR loop, so the downstream `fsusie()` fit on
  `Y_adjusted` is bit-identical).
- No mention of any port-source package in vignette prose
  per the standing rule.

### 6. Renamed vignette

`vignettes/fsusie_covariates_and_coloc.Rmd` ->
`vignettes/fsusie_colocalization.Rmd`. The covariates section
is removed (it lived in the combined vignette only because we
had no proper utility); the remaining content is the
colocalization workflow already in place.

### 7. Unit tests

`tests/testthat/test_adjust_covariates.R`:

- Apple-to-apple bit-identity vs `fsusieR:::EBmvFR(Y, X = Z,
  adjust = TRUE)$Y_adjusted` and `$fitted_func` at tolerance
  `<= 1e-8`. Coverage: `n in {100, 200}`, `T in {64, 128}`,
  `K in {2, 3, 5}`. Gated `skip_if_not_installed("fsusieR")`.
- Math-correctness sanity for `method = "ols"`: residuals
  orthogonal to Z, idempotency (`adjust(adjust(Y)) == adjust(Y)`).
- Error paths: bad shapes, unsupported flags
  (`thresh_lowcount > 0`, `quantile_trans = TRUE`), method
  dispatch.
- Refactored EM helpers: smoke tests in
  `tests/testthat/test_em_helpers.R`, exercising
  `is_ebmvfr = TRUE` and `is_ebmvfr = FALSE` against expected
  shapes.

## Impact

- New: `R/adjust_covariates.R` (~150 lines),
  `tests/testthat/test_adjust_covariates.R`, new vignette,
  small additions to `R/em_helpers.R`.
- Renamed: `vignettes/fsusie_covariates_and_coloc.Rmd` ->
  `vignettes/fsusie_colocalization.Rmd`.
- DESCRIPTION: add `fsusieR` to Suggests (test-only;
  `skip_if_not_installed`).
- pixi.toml: no change (`fsusieR` is installed in CI for the
  apple-to-apple tests via the existing `r-remotes` machinery).
- No public-API breakage. Pure additive plus a vignette
  rename.
- Specs: new `inst/openspec/specs/mf-covariate-adjustment/spec.md`.

## Out of scope

- `thresh_lowcount > 0` (low-count wavelet-coefficient filtering).
  Errors with "not supported in v1". Tracked for follow-up.
- `quantile_trans = TRUE` (normal-quantile transform on wavelet
  coefficients). Errors with "not supported in v1". Tracked
  for follow-up.
- ELBO computation in EBmvFR. Upstream `stop()`s on
  `cal_obj = TRUE` so this is non-issue at defaults.
- Full `EBmvFR()` as a SuSiE alternative (the original
  out-of-scope EBmvFR ruling stays for that purpose).
