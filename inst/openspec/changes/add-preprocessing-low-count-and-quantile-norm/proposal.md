# Add `low_count_filter` and `quantile_norm` preprocessing options

## Why

Two preprocessing options control how mfsusieR treats wavelet-
domain data before fine-mapping: a low-count filter that masks
sparse columns and a column-wise quantile transform that
maps wavelet coefficients onto the standard normal scale. Both
are present in the upstream functional fine-mapping pathways
that mfsusieR replaces; mfsusieR currently exposes neither and
its `mf_adjust_for_covariates` rejects them outright. Adding
both options closes a real gap in the public interface and
strengthens the apple-to-apple contract on data with sparse
or zero-median wavelet columns.

The current state has two consequences worth flagging:

1. On data with at least one zero-median wavelet column, the
   fine-mapping result deviates from the upstream functional
   fine-mapping result even when both packages are called with
   the upstream default `thresh_lowcount = 0`. The upstream
   default still triggers the median-zero mask via
   `which_lowcount(Y_f, 0)`; mfsusieR does not run that step
   and therefore differs.
2. Cell-type ATAC-seq and other sparse-coverage applications
   benefit from quantile normalization on the wavelet
   coefficients to suppress heavy-tailed residuals; mfsusieR
   has no entry point for this step.

## What changes

### 1. Public arguments

Add to `mfsusie()`, `fsusie()`, and
`mf_adjust_for_covariates()`:

```
low_count_filter = 0
quantile_norm    = FALSE
```

`low_count_filter` is a non-negative threshold; columns of
the wavelet matrix with `median(|column|) <= low_count_filter`
are masked. `quantile_norm = TRUE` applies a column-wise
rank-based normal quantile transform to the wavelet matrix
before the IBSS loop.

### 2. New helpers in `R/utils_wavelet.R`

```
mf_low_count_indices(Y_wd, threshold = 0)
  -> integer vector of column indices with
     median(|.|) <= threshold
mf_quantile_normalize(Y_wd)
  -> column-wise INT(rank(., ties.method = "random")) under
     `set.seed(1)` to bit-match upstream
```

Twelve and six lines respectively. Both ported as semantic
features, not as algorithmic surgery on susieR.

### 3. Marginal regression stays in susieR

The IBSS loop calls `susieR::compute_marginal_bhat_shat()`
unchanged. The low-count mask is applied via a thin wrapper:

```
mf_compute_bhat_shat_with_low_count_mask(
  X, Y_wd, lowc_idx)
  bs <- susieR::compute_marginal_bhat_shat(X, Y_wd)
  if (length(lowc_idx) > 0L) {
    bs$Bhat[, lowc_idx] <- 0
    bs$Shat[, lowc_idx] <- 1
  }
  bs
```

This matches the upstream pipeline's post-compute mask
(Bhat -> 0, Shat -> 1 on flagged columns) without
re-implementing marginal regression in mfsusieR.

### 4. Pipeline placement

After the existing per-outcome DWT in
`R/individual_data_class.R::create_mf_individual`:

```
1. col_scale + DWT on Y -> W$D, W$C
2. Y_wd <- cbind(W$D, W$C)
3. lowc_idx <- mf_low_count_indices(Y_wd, low_count_filter)
4. if (quantile_norm) Y_wd <- mf_quantile_normalize(Y_wd)
5. data$lowc_idx[[m]] <- lowc_idx
6. data$D[[m]] <- Y_wd      (replaces the unmasked matrix)
```

The IBSS path threads `data$lowc_idx[[m]]` into both
`init_scale_mixture_prior_default` (already accepts a
low-count argument) and the per-iteration
`mf_compute_bhat_shat_with_low_count_mask`. The mask is also
applied during prior weight EM updates so the masked columns
do not contribute to the mixsqp likelihood matrix.

### 5. `mf_adjust_for_covariates` updates

The v1 reject in `R/adjust_covariates.R::mf_adjust_for_covariates`
that errors on `thresh_lowcount > 0` and `quantile_trans =
TRUE` is lifted. The wavelet-EB residualization path threads
`lowc_idx` and `Y_wd_normalized` through the same helper
flow as the main IBSS pipeline. When `quantile_norm = TRUE`
the function returns `Y_adjusted` on the original (un-INT'd)
position scale by inverting the transform after the
residualization step so the downstream `fsusie(Y_adjusted, X)`
call works without unit confusion.

### 6. Roxygen and vignettes

- Public-API roxygen for `mfsusie()`, `fsusie()`, and
  `mf_adjust_for_covariates()` documents both arguments
  with their semantics and references the upstream-equivalent
  default values.
- `vignettes/mfsusie_intro.Rmd` and
  `vignettes/fsusie_intro.Rmd` get short paragraphs noting
  when each option is appropriate. The full demonstration
  lives in the `practical_data_applications` vignette
  (separate change `add-practical-data-applications-vignette`).

### 7. Tests

`tests/testthat/test_preprocessing.R` (new):

- Bit-identity vs upstream functional fine-mapping with both
  flags toggled across `n in {84, 200}`,
  `T in {64, 128, 256}`, and `low_count_filter in {0, 0.1, 0.5}`,
  `quantile_norm in {FALSE, TRUE}`. Tolerance `<= 1e-12` per
  the tightened C2 contract.
- A targeted test on synthetic data with at least one
  zero-median wavelet column at `low_count_filter = 0` to
  verify the upstream mask is applied (this is the off-spec
  case the current mfsusieR misses).
- Unit tests for `mf_low_count_indices` and
  `mf_quantile_normalize` against `which_lowcount` and
  `Quantile_transform` on randomized inputs.

## Impact

- New: `R/utils_wavelet.R` extensions (~30 lines),
  `R/preprocessing.R` (small wrapper file, ~30 lines), tests.
- Changed: `R/individual_data_class.R` (DWT pipeline tail),
  `R/mfsusie.R` and `R/fsusie.R` public signatures,
  `R/adjust_covariates.R` (lift v1 rejects), three vignettes
  (intro mentions only).
- C2 contract: tightened. The test suite now exercises the
  zero-median-column case where pre-port mfsusieR was
  off-spec.
- DESCRIPTION: no changes (susieR / wavethresh / mixsqp / ashr
  already imported).
- Specs: new `inst/openspec/specs/mf-preprocessing/spec.md`.

## Out of scope

- No reimplementation of marginal regression in mfsusieR;
  `susieR::compute_marginal_bhat_shat` stays the engine.
- No quantile-normalization-aware variant of the post-fit
  smoother (out of scope for round 3; the smoothers operate on
  fitted curves, not on the IBSS data path).
