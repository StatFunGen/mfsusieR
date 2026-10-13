## Why

The permutation FDR study (meeting 2026-09-29, Section 4) found that a
wavelet-level TensorQTL pre-filter at p < 1e-5 removes 191/200 null
regions while passing all 200 real regions, reducing wavelet FDR to 0%
for all CIP=OFF per-outcome configs without cutting real HP-CS yield.

The meeting of 2026-10-06 directed that this screen be integrated into
the mfSuSiE R package as a unified pipeline function: users supply one
or more regions (each with X, Y, pos), the package runs the association
screen internally, and only passing regions enter the full mfsusie() fit.
Users should not have to manage a separate upstream TensorQTL workflow.
The same pipeline should accept genome-wide candidate regions and handle
the parallel split (e.g. by chromosome) that the production run requires.

## What Changes

**New public functions in R/mfsusie_pipeline.R:**

1. `mfsusie_screen(X, Y, pos, wavelet_qnorm = TRUE, p_cutoff = 1e-5, ...)`
   - Screens a single region via wavelet-level marginal association tests.
   - Applies `mf_dwt()` per outcome, optionally `mf_quantile_normalize()`,
     then computes Pearson correlations between each (wavelet column, SNP)
     pair and converts to a two-sided t-test p-value (df = n - 2).
   - Returns a list with `min_pval`, `passed` (logical), `n_tested` (integer
     count of (column, outcome, SNP) triples), and `region_id` (passed through).
   - Does not call mfsusie(); intended as a standalone QC utility.

2. `mfsusie_pipeline(regions, use_wavelet_filter = TRUE, wavelet_qnorm = TRUE,
                     p_cutoff = 1e-5, n_cores = 1L, verbose = TRUE, ...)`
   - Accepts `regions`: a named or unnamed list where each element has fields
     `X`, `Y`, `pos`, and optionally `region_id` (character scalar).
   - For each region: (a) screen with `mfsusie_screen()` if
     `use_wavelet_filter = TRUE`; (b) if passed, call `mfsusie(X, Y, pos, ...)`
     with `wavelet_qnorm = wavelet_qnorm` forwarded; (c) if failed, return a
     structured empty result.
   - Supports `n_cores > 1` via `parallel::mclapply` for multi-core execution.
   - Returns `list(results = ..., funnel = ...)`:
     - `results`: named list (by region_id or index) of per-region outputs,
       each an `mfsusie` fit object (passed regions) or a structured empty
       object (failed/no-discovery regions) with consistent class and fields.
     - `funnel`: data.frame with one row per region and columns
       `region_id`, `status` ("screen_rejected" | "fit_completed" |
       "fit_failed"), `min_pval`, `n_cs`, `n_hp_cs`, `cell_types_entered`,
       `cell_types_active` (post M-reduction, if applicable).

3. `mfsusie_empty_result(region_id = NULL, status = "screen_rejected",
                          min_pval = NA_real_)`
   - Returns an `mfsusie` object with class `c("mfsusie", "susie")` and the
     same top-level fields as a fitted object, but with empty/NULL contents.
   - Carries a `.$screen_result` sub-list with `region_id`, `status`,
     and `min_pval`.
   - Downstream code (e.g. HP-CS extraction, predict.mfsusie) must not error
     on these objects.

**No changes to existing public API.** `mfsusie()`, `fsusie()`, and all
existing S3 methods are unchanged. The `wavelet_qnorm` parameter already
accepted by `mfsusie()` is reused consistently.

## Non-Goals

- This change does not modify the IBSS fitting loop or any method registered
  in R/zzz.R.
- Parallel via `future` or `BiocParallel` is out of scope; `parallel::mclapply`
  (forking, Unix/macOS only) is sufficient for the chromosome-split use case.
  Windows users fall back to `n_cores = 1L`.
- The screen does not implement Bonferroni correction or other FWER controls;
  the p_cutoff is an empirical operating threshold.
- No new DESCRIPTION Imports beyond `parallel` (already available in base R).

## Backward Compatibility

None of the existing functions change. `mfsusie_pipeline()` is additive.
