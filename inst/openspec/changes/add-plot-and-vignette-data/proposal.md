# Add `mfsusie_plot()`, `mf_post_smooth()` (scalewise / TI / HMM),
# packaged example data, and refreshed vignettes

## Why

Vignettes used inline simulations and ad-hoc base-R plotting. The
SuSiE family has paired plot helpers (`susieR::susie_plot`,
`mvsusieR::mvsusie_plot`, `fsusieR::plot_susiF`); mfsusieR was
missing the equivalent. Effect curves were the raw inverse-DWT of
the wavelet posterior, which is correct but visually noisy and
without credible bands.

fSuSiE provides three post-processing smoothers
(`smash_regression`, `TI_regression`, `HMM_regression`). mfsusieR
needs the same three so users get publication-quality figures
matching fSuSiE / mfSuSiE paper figures.

## What changes

### 1. Plot API (single function)

- Re-export `susieR::susie_plot()` for the standard PIP plot.
- New `mfsusie_plot(fit, m = NULL, ...)` (base R only, no
  `ggplot2` / `cowplot` dependency). Auto-detects `M`:
  - `M = 1` -> 2-panel column (PIP top, effect bottom)
  - `M > 1` -> dense grid (PIP top-left, M effect panels)
  - `m = i` -> single-outcome effect panel only
- S3 method `plot.mfsusie()` dispatches to `mfsusie_plot()`.
- LFSR overlay: when `fit$lfsr_curves[[m]][[l]]` is populated
  (HMM smoother), the plot draws it as a secondary black trace
  with a dashed `0.05` reference line, mirroring
  `fsusieR::plot_susiF_effect`'s style.
- Color palette: Okabe-Ito (modern, colorblind-friendly).

### 2. `mf_post_smooth(fit, method = c("scalewise", "TI", "HMM"))`

Single user-facing entry point with method dispatch. Operates on
`fit` only (no `X`, no `Y`). The `mfsusie()` / `fsusie()` wrapper
saves what's needed:

- `fit$residuals[[m]]` -- wavelet-domain residual `D - fitted`.
- `fit$lead_X[[l]]` -- per-effect lead-variable column of the
  centred + scaled `X`. Length-`n` numeric vectors, one per
  effect.

Both are always populated. The `save_residuals` argument is
removed from the public API (cost is < 1 MB per fit).

Three method implementations:

- **`"scalewise"`** -- per-scale soft-threshold of the lead
  variable's wavelet posterior mean at `factor * sigma_s *
  sqrt(2 log T)`. Pointwise band from `mu2 - mu^2` Parseval-
  inverted. Fast, uses no residuals.
- **`"TI"`** -- faithful port of `fsusieR::univariate_TI_regression`.
  For each effect: subtract other effects from the
  position-space response, regress on the lead variable,
  stationary-wavelet-transform, scalewise `ashr::ash` shrinkage
  on detail coefficients, cycle-spinning average via
  `wavethresh::av.basis`. Pointwise variance via the ported
  `wd.var` / `AvBasis.var` / `convert.var` helpers (~300 lines
  added to `R/utils_wavelet.R`).
- **`"HMM"`** -- faithful port of
  `fsusieR::univariate_HMM_regression`. Per-position regression on
  lead variable, then `fit_hmm` posterior smoothing. Populates
  `fit$lfsr_curves[[m]][[l]]` in addition to
  `effect_curves` / `credible_bands`.

### 3. cpp11 acceleration (deferred, but planned)

The TI and HMM smoothers loop over wavelet shifts and per-position
regressions; both can be slow on T_basis = 1024. cpp11 ports of
the inner loops are tracked as a follow-up under
`add-cpp11-smoothers` (a separate OpenSpec change). Pure-R
versions ship now and are correct.

### 4. Packaged example data

Two `.rda` files under `data/` so vignettes load via `data(<name>)`
without external paths:

- `data(fsusie_methyl_example)` -- DNAm fixture from
  `~/GIT/fsusie-experiments/plot/CR1_CR2/CR1_CR2_obj.RData`,
  trimmed to ~250 KB. Sample IDs replaced with anonymous
  `S001..S<n>`; PHI fields stripped; only public reference
  identifiers (chrom, position, alleles, CpG positions) retained.
- `data(mfsusie_joint_example)` -- multi-outcome fixture, similarly
  trimmed and de-identified. Source: trimmed from CASS4 region.

Build scripts under `data-raw/`; `.Rbuildignore`d.

### 5. Vignette refresh

Eight vignettes:

- `getting_started.Rmd` -- public-API tour; uses
  `data(fsusie_methyl_example)` for an `mfsusie_plot()` preview.
- `fsusie_intro.Rmd` -- adapts `fsusieR::fsusie_intro`.
- `fsusie_covariates_and_coloc.Rmd` -- adapts
  `fsusieR::Adjusting_covariate` + `Coloc_fsusie`.
- `fsusie_dnam_case_study.Rmd` -- uses
  `data(fsusie_methyl_example)` (CR1/CR2 region).
- `fsusie_why_functional.Rmd` -- adapts
  `fsusieR::Limitation_SuSiE` + `susie_top_PC_fails`.
- `mfsusie_intro.Rmd` -- uses `data(mfsusie_joint_example)` (CASS4).
- `mfsusie_long_running_fits.Rmd` -- keep current.
- `post_processing.Rmd` -- demonstrates all three smoothers with
  before / after plots, plus an LFSR curve panel for HMM.

Every vignette uses `mfsusie_plot()` (or `susie_plot()` for PIPs)
instead of ad-hoc `matplot()`. Titles describe topics, not entry
points (e.g., "DNA methylation QTL fine-mapping" not
"DNA methylation QTL fine-mapping with `fsusie()`").

### 6. NAMESPACE / Imports

- `export(mfsusie_plot, mf_post_smooth, susie_plot)` (last is a
  re-export from susieR).
- `S3method(plot, mfsusie)`.
- New `@importFrom wavethresh av.basis convert nlevelsWT` (already
  partially imported).
- `ashr::ash` is already imported.
- No new package Imports.

## Impact

- Affected specs: new `inst/openspec/specs/mf-plot/spec.md` and
  `inst/openspec/specs/mf-post-processing/spec.md`.
- Affected code: new `R/utils_wavelet.R` extensions (variance
  helpers); `R/mfsusie_methods.R` gains the three smoother
  implementations and the unified `mfsusie_plot()` dispatch path
  (already merged into `R/mfsusie_methods.R` in the current
  session); `R/mfsusie.R` always populates residuals + lead_X.
- Affected tests: new `test_post_smooth.R` covering smoke +
  bit-identity vs `fsusieR` (skip-when-not-installed) for TI / HMM
  point estimates at `<= 1e-8`; smoke for scalewise.
- Affected vignettes: all 8 refreshed.
- No behaviour change in `fsusie()` / `mfsusie()`. Pure additive.
- Pre-1.0; no deprecation.

## Out of scope (own future change)

- cpp11 acceleration of TI/HMM inner loops --
  `add-cpp11-smoothers`.
- Smash smoother (`fsusieR::smash_regression`) -- adds `smashr`
  Imports; track in `add-smash-smoother`.
- Gviz integration -- `add-gviz-plotting`.
