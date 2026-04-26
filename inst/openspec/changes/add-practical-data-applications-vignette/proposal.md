# Add the `practical_data_applications` vignette

## Why

The three preprocessing and small-sample-correction features
landing in round 3 (`low_count_filter`, `quantile_norm`, and
`estimate_residual_method = "NIG"`) need a vignette that
demonstrates all three on a real-world dataset and provides
the brief mathematical derivation for the NIG correction.

The dataset at `~/Documents/fsusie_test/` provides 169
multi-region multi-cell-type ATAC-seq fixtures (n = 84
individuals; 6 cell-type phenotypes per region; T = 10240
bins per cell-type; ~3500 SNPs per region after MAF filter).
The small sample size (n = 84) is exactly the regime where
the NIG correction matters; the sparse-coverage structure is
the regime where `low_count_filter` and `quantile_norm` matter.

The original ATAC-seq counts and individual-level data are
not redistributable. The vignette ships only fitted-model
output from a small set of selected regions; raw inputs are
never displayed.

## What changes

### 1. New vignette

`vignettes/practical_data_applications.Rmd`. Authors:
Anjing Liu, Gao Wang, William Denault.

Sections:

1. **The small-sample problem.** Frame the issue: when
   n is comparable to or smaller than the number of
   wavelet-domain observations per outcome, the Method-of-
   Moments residual variance estimator underestimates
   `sigma^2`, inflating per-variable Bayes factors and
   PIPs.
2. **NIG correction.** Brief mathematical derivation:
   - Prior: `sigma^2 ~ Inverse-Gamma(alpha_0, beta_0)`,
     posterior `sigma^2 | data ~ Inverse-Gamma(alpha_0', beta_0')`.
   - SER Bayes factor integrates over `sigma^2` under
     the posterior instead of conditioning on a point
     estimate; the resulting per-variable BF takes the
     form of a Student's t marginal likelihood.
   - Notation matches the SuSiE manuscript and the
     `susieR` source.
   - Demonstration:
     `mfsusie(..., estimate_residual_method = "NIG")` on
     a representative region. Compare PIPs against the
     MoM default; show that NIG produces tighter PIPs at
     the null variants.
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
  selects 2-3 representative regions, runs `mfsusie()`
  under each setting (MoM default, NIG, low-count-filter,
  quantile-norm), saves the fitted-model objects to
  `data/practical_fits.rda`. Raw inputs are NOT saved
  with the package.
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
  `quantile_norm`, and NIG — those land in the prior
  changes
  (`add-preprocessing-low-count-and-quantile-norm`,
  `expose-nig-residual-variance-method`). This change
  depends on those landing first.
- Demonstrating warm-start. That belongs in the
  `mfsusie_long_running_fits` vignette which is updated
  in `expose-model-init-warm-start`.
