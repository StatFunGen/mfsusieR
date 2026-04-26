# Add `fsusie_plot()` / `mfsusie_plot()`, packaged example data, and refreshed vignettes

## Why

Vignettes currently render with simulated data inline, so figures are
unstable across knit runs and not comparable to the figures shown in
the fSuSiE paper / `fsusieR` and `mvf.susie.alpha` original vignettes.
Plotting in vignettes is ad-hoc (`matplot`, `plot(fit$pip)`, base
graphics) without a unified API.

The SuSiE family already has paired plotting helpers:

- `susieR::susie_plot()` — PIP plot for scalar SuSiE
- `mvsusieR::mvsusie_plot()` — PIP plot + per-outcome effect / lfsr
  plot for multi-trait scalar SuSiE
- `fsusieR::plot_susiF()` — PIP plot + per-effect functional curve
  with credible bands (uses ggplot2 + cowplot)

mfsusieR is missing the equivalent. Vignettes can't show the actual
methylation / chromatin-track figures readers expect from this
package family.

## What changes

### 1. Plot API

Three plotting functions, all base R (no ggplot, no cowplot):

- **PIP plot for either fsusie or mfsusie fit**: re-export
  `susieR::susie_plot()` from mfsusieR. The fit object inherits from
  `c("mfsusie", "susie")` so `susie_plot(fit, y = "PIP")` already
  works. We just expose it.

- **`fsusie_plot(fit, ...)`**: single-outcome functional effect plot.
  Takes a fit from `fsusie()` (M=1, T_m > 1). Plots the per-effect
  curves against the post-remap position grid, colored by credible
  set, with optional credible bands (when post-processing is
  available; section 2 below) and optional affected-region shading.

- **`mfsusie_plot(fit, m = NULL, ...)`**: multi-outcome version.
  When `m` is `NULL` (default), tiles one panel per outcome via
  `par(mfrow = ...)`; the user can also pass a single `m` index to
  see one outcome. Each panel shows the per-effect curves on the
  outcome's position grid colored by CS, plus a top-right legend
  with the variable index of the sentinel for each CS. Scalar
  outcomes (`T_m = 1`) get a dot-plot with horizontal CS bars
  instead of curves.

- **`plot.mfsusie(x, ...)`** S3 method: dispatch on the fit class so
  `plot(fit)` Just Works. Defaults to `mfsusie_plot()`; degenerates
  to `fsusie_plot()` when `M == 1`.

### 2. Post-processed effect storage on fit objects

For credible-band shading on `fsusie_plot()` / `mfsusie_plot()`, the
fit needs:

- `fit$effect_curves[[m]]` — a length-`L` list of length-`T_m`
  vectors holding the post-remap-grid-to-original-grid-back-projected
  effect for each effect `l`, optionally smoothed.
- `fit$credible_bands[[m]]` (optional) — a length-`L` list of
  `T_m x 2` matrices `[lower, upper]` holding the credible-band
  endpoints. When absent, the plot draws curves only, no bands.
- `fit$lfsr_curves[[m]]` (optional) — local false sign rate curves
  per effect. When present, the plot shows them as a secondary
  black curve with a 0.05 threshold dashed line.

These slots are populated by the post-processor (`mf_post_smooth()`,
PR group 6b) when run. The plot functions handle absence gracefully:
they fall back to plotting the wavelet-domain effect curves
projected through `mf_invert_dwt()` (which `coef.mfsusie()` already
does).

### 3. Packaged example data

Add two datasets under `data/` so vignettes are reproducible.
`susieR::N3finemapping` is already available via the `susieR`
Imports and is used directly without re-packaging:

- `data/fsusie_methyl_example.rda` — small DNAm fixture (n ~= 200,
  p ~= 100, T_m ~= 64) for `fsusie_dnam_case_study.Rmd`. Either a
  trimmed version of `fsusie-experiments/fig_1_data` or a fixed-seed
  simulation saved as RDA.
- `data/mfsusie_joint_example.rda` — joint-fine-mapping fixture
  (M = 5: 2 functional + 3 scalar; n ~= 200; p ~= 100). Either ported
  from `mvf.susie.alpha` simulation outputs (via
  `mvf.susie.alpha::simu_effect_multfsusie` with fixed seed, run
  once and saved) or built locally in
  `data-raw/make_mfsusie_data.R`.

Each rda file is added to the `data/` directory and roxygen
documentation added in `R/data.R`. Total data size kept under
~500KB to stay within CRAN's 5MB package size budget.

### 4. Vignette refresh

All seven vignettes:

- Replace inline simulations with `data(<example>)` loads from the
  packaged datasets above. Existing inline simulations stay as a
  one-paragraph "or simulate your own" callout.
- Replace ad-hoc `matplot()` / `plot()` calls with `susie_plot()`,
  `fsusie_plot()`, or `mfsusie_plot()` as appropriate.
- Add reference figures: each vignette ends with a "Reference"
  section that names the original fsusieR / mvf.susie.alpha vignette
  whose figure ours reproduces.

Specific vignette-by-vignette mapping:

- `getting_started.Rmd` — small fsusie_plot() + mfsusie_plot()
  preview using the packaged data. Keep the `susieR` C1 numerical
  identity demo.
- `fsusie_intro.Rmd` — adapts `fsusieR::fsusie_intro` with our
  packaged data + `fsusie_plot()`.
- `fsusie_covariates_and_coloc.Rmd` — adapts
  `fsusieR::Adjusting_covariate` + `fsusieR::Coloc_fsusie`.
- `fsusie_dnam_case_study.Rmd` — adapts `fsusieR::methyl_demo`
  using `data(fsusie_methyl_example)`.
- `fsusie_why_functional.Rmd` — adapts `fsusieR::Limitation_SuSiE`
  + `susie_top_PC_fails`.
- `mfsusie_intro.Rmd` — adapts
  `mvf.susie.alpha::joint_functional_and_univariate_fine_mapping`
  using `data(mfsusie_joint_example)`.
- `mfsusie_long_running_fits.Rmd` — keep current content; add an
  `mfsusie_plot()` preview at the end.

### 5. NAMESPACE / Imports

- Export `fsusie_plot`, `mfsusie_plot`, `plot.mfsusie`.
- Re-export `susieR::susie_plot` so users do `mfsusieR::susie_plot()`
  without loading susieR explicitly.
- No new package Imports. Plot uses `graphics::*` (base) and
  `grDevices::*` (already loaded by R).

## Impact

- Affected specs: new `mf-plot/spec.md` (added under
  `inst/openspec/specs/`).
- Affected code: ~3 new R files (`R/fsusie_plot.R`,
  `R/mfsusie_plot.R`, `R/data.R`); `R/mfsusie_methods.R` gains
  `plot.mfsusie`; `data/*.rda` added; all vignettes rewritten.
- Tests: smoke tests that the plot functions don't error
  (`expect_silent(plot(fit))` style) plus snapshot tests for the
  data-loading paths.
- No behaviour change in `fsusie()` / `mfsusie()`. Pure additive.
- Pre-1.0; no deprecation needed.
