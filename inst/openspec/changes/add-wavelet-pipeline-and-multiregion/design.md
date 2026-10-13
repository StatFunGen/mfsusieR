## Context

The permutation FDR benchmark (200 real regions, 200 permuted regions)
showed that a wavelet-level marginal association screen at p < 1e-5
passes all 200 real regions and only 9/200 permuted regions, reducing
wavelet FDR to 0% for CIP=OFF per-outcome configs. The meeting of
2026-10-06 requested this screen be packaged as a usable function
shipped with mfsusieR, with a multi-region pipeline entry that handles
the screen + fit + empty-result logic in one call.

The code in `R/dwt.R` already provides `mf_dwt()` for the wavelet
transform and `R/utils_wavelet.R` provides `mf_quantile_normalize()`.
The new pipeline reuses these internals without duplicating them.

## Goals / Non-Goals

**Goals:**

- Ship `mfsusie_screen()` as a standalone region-level QC function.
- Ship `mfsusie_pipeline()` as a multi-region entry that screens then fits.
- Provide `mfsusie_empty_result()` for consistent return type on failed
  or screened-out regions.
- Record an analysis funnel (input -> screen pass -> fit complete -> CS/HP-CS)
  for every `mfsusie_pipeline()` call.
- Flexible parameters: `use_wavelet_filter`, `wavelet_qnorm`, `p_cutoff`,
  `n_cores`.
- `wavelet_qnorm` in the screen and in fitting use the same default (TRUE)
  so the representation tested and fitted are consistent.

**Non-Goals:**

- No changes to existing functions or S3 dispatch.
- No Bonferroni correction; `p_cutoff` is an empirical threshold.
- No `future` / `BiocParallel` parallel backend.

## API Contract

### D1. Region input format

Each element of `regions` is a named list:

```r
list(
  X         = <numeric matrix n x p>,
  Y         = <list of M numeric matrices n x T_m>,
  pos       = <list of M numeric vectors length T_m>,
  region_id = <character scalar, optional>
)
```

When `region_id` is absent, the pipeline uses the element's name in the
`regions` list, or `"region_<index>"` as a fallback.

### D2. mfsusie_screen() return value

```r
list(
  region_id = character(1),
  min_pval  = numeric(1),   # min over all (wavelet col, outcome, SNP) pairs
  n_tested  = integer(1),   # total (wavelet col x outcome x SNP) tested
  passed    = logical(1)    # min_pval < p_cutoff
)
```

The screen skips outcomes with `T_m = 1` (scalar phenotypes have no
wavelet expansion; they are included in mfsusie() fitting unchanged).
For constant-variance wavelet columns (all values identical after QN),
the correlation is 0 and the p-value is 1; they do not affect `min_pval`.

### D3. mfsusie_empty_result() fields

```r
structure(
  list(
    alpha        = NULL,
    mu           = NULL,
    mu2          = NULL,
    lbf          = NULL,
    lbf_variable = NULL,
    V            = NULL,
    sigma2       = NULL,
    pi           = NULL,
    pip          = numeric(0),
    sets         = list(cs = list(), purity = list(), cs_index = integer(0)),
    converged    = NA,
    niter        = 0L,
    screen_result = list(
      region_id = character(1),
      status    = character(1),  # "screen_rejected" | "fit_failed"
      min_pval  = numeric(1)
    )
  ),
  class = c("mfsusie", "susie")
)
```

`predict.mfsusie`, `coef.mfsusie`, and `fitted.mfsusie` must not error
on this object; they should return zero-length or NA results.

### D4. mfsusie_pipeline() funnel data.frame

Columns: `region_id` (character), `status` (character), `min_pval`
(numeric), `n_cs` (integer), `n_hp_cs` (integer),
`cell_types_entered` (integer), `cell_types_active` (integer or NA
when screen rejected).

`status` values: "screen_rejected", "fit_completed", "fit_failed".
"fit_completed" includes regions with zero CS (no discovery).

### D5. Parallelism

`n_cores = 1L` uses `lapply` (no forking). `n_cores > 1L` uses
`parallel::mclapply(..., mc.cores = n_cores)`. On Windows,
`parallel::mclapply` falls back to `lapply` automatically; the function
emits a warning when `n_cores > 1L` is requested on `.Platform$OS.type
== "windows"`.

### D6. Consistency of wavelet_qnorm

The `wavelet_qnorm` argument is forwarded to both `mfsusie_screen()`
(used in the internal DWT before the association test) and to
`mfsusie()` (used during IBSS fitting). Defaults are both TRUE. A
mismatch between screen and fit representations is a bug; the pipeline
enforces identical values.

## Implementation Notes

The correlation-based p-value computation in `mfsusie_screen()` is:

```r
r   <- cor(D_m, X_scaled)               # T_basis x p matrix
t   <- r * sqrt((n - 2) / (1 - r^2))
pv  <- 2 * pt(-abs(t), df = n - 2)
min_pval <- min(pv, na.rm = TRUE)
```

This gives the same p-value as `lm(D_m[, j] ~ X[, k])` but avoids
the overhead of p separate lm calls per wavelet column.

`X` should be centered (mean 0) before computing correlations to match
the centering that `create_mf_individual()` applies before IBSS. The
screen centers X internally; it does not modify the caller's X.

## File Layout

- `R/mfsusie_pipeline.R` — all three new functions plus roxygen docs.
- No other files change.

## Test Plan

A new test file `tests/testthat/test_mfsusie_pipeline.R` covers:

- `mfsusie_screen()` with a functional outcome (T_m > 1, wavelet_qnorm
  TRUE and FALSE) returns a valid list with correct field names.
- `mfsusie_screen()` returns `passed = FALSE` for a shuffled X (null
  region); smoke test, not a probabilistic guarantee.
- `mfsusie_pipeline()` on a list of two regions (one real, one null)
  returns a funnel with both statuses represented.
- `mfsusie_empty_result()` has class `c("mfsusie", "susie")` and does not
  error under `predict.mfsusie`.
- `mfsusie_pipeline()` with `use_wavelet_filter = FALSE` fits all regions
  regardless of association evidence.
