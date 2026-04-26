# Add the `practical_data_applications` vignette

## Why

The three preprocessing and small-sample-correction features
landing in round 3 (`low_count_filter`, `quantile_norm`, and
`small_sample_correction`) need a vignette that demonstrates
all three on a real-world dataset and provides the brief
mathematical statement of the Johnson 2005 scaled Student-t
marginal Bayes factor used by the small-sample correction.

The dataset at `~/Documents/fsusie_test/` provides 169
multi-region multi-cell-type ATAC-seq example datasets (n = 84
individuals; 6 cell-type phenotypes per region; T = 10240
bins per cell-type; ~3500 SNPs per region after MAF filter).
The small sample size (n = 84) is exactly the regime where
the Johnson-t correction matters; the sparse-coverage
structure is the regime where `low_count_filter` and
`quantile_norm` matter.

The original ATAC-seq counts and individual-level data are
not redistributable. The vignette ships only fitted-model
output from a small set of selected regions; raw inputs are
never displayed.

## What changes

### 1. New vignette

`vignettes/practical_data_applications.Rmd`. Authors:
Anjing Liu, Gao Wang, William Denault.

Sections:

1. **The small-sample problem.** Frame the issue: at small
   `n`, the Wakefield Normal marginal Bayes factor in the SER
   step under-propagates residual-variance uncertainty,
   inflating per-variable PIPs at null variants.
2. **Johnson-t correction.** Brief statement of the Johnson
   2005 scaled Student-t marginal:
   - Wakefield: per-component density
     `dnorm(Bhat; 0, sqrt(sigma_k^2 + Shat^2))`.
   - Johnson-t: per-component density
     `dstp(Bhat; tau = 1 / (sigma_k^2 + Shat^2), nu = df)`
     with `df = n - 1`.
   - Demonstration:
     `mfsusie(..., small_sample_correction = TRUE)` on a
     representative region. Compare PIPs against the
     Wakefield default; show that the correction tightens
     PIPs at the null variants.
3. **Low-count filtering for sparse-coverage cell types.**
   Some cell types in the dataset have median coverage
   close to zero in many wavelet columns. Without
   `low_count_filter > 0`, these columns contribute
   noise to the prior init and the SER step. Demonstrate
   the effect:
   `mfsusie(..., low_count_filter = 0.5)` vs the default,
   PIP comparison.
4. **Quantile normalization for non-Gaussian wavelet
   coefficients.** Sparse-count assays produce heavy-
   tailed wavelet coefficients. Demonstrate
   `quantile_norm = TRUE` rendering the coefficients
   approximately Gaussian and the downstream effect on
   the fit.

### 2. Data flow

- `data-raw/make_practical_dataset.R` (new): loads the rds
  and the analysis script from `~/Documents/fsusie_test/`,
  selects a representative region, runs `mfsusie()` under
  each setting (Wakefield default, small-sample correction,
  low-count filter, quantile normalization), saves the
  fitted-model objects to `data/practical_fits.rda`. Raw
  inputs are NOT saved with the package.
- The vignette loads `data(practical_fits)` and renders
  plots and summary tables. The original `X`, `Y`, and
  individual-level data never appear in the vignette
  HTML.

### 3. Document the data

`R/data.R` gets a `practical_fits` roxygen entry describing
the fields present and a clear note that the underlying
ATAC-seq data is not redistributable; the dataset ships
fitted-model summaries only.

### 4. pkgdown wiring

`_pkgdown.yml` articles list adds
`practical_data_applications` under the "Single-outcome
fine-mapping" or "Multi-outcome fine-mapping" section
(it covers multi-outcome data, so the latter).

## Impact

- New: `vignettes/practical_data_applications.Rmd`,
  `data-raw/make_practical_dataset.R`,
  `data/practical_fits.rda`, roxygen entry in `R/data.R`.
- Changed: `_pkgdown.yml` articles list.
- Specs: new
  `inst/openspec/specs/mf-vignette-practical/spec.md`.

## Out of scope

- The actual implementation of `low_count_filter`,
  `quantile_norm`, and `small_sample_correction` — those
  land in the prior changes
  (`add-preprocessing-low-count-and-quantile-norm`,
  `add-small-sample-correction-johnson-t`). This change
  depends on those landing first.
- Demonstrating warm-start. That belongs in the
  `mfsusie_long_running_fits` vignette which is updated
  in `expose-model-init-warm-start`.
