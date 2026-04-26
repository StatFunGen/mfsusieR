# Design

## Decisions locked

- **Smoother**: ship a minimal TIWT post-processor now,
  `mf_post_smooth(fit, method = "tiwt")`. Adds
  `$effect_curves` (smoothed position-space curves) and
  `$credible_bands` (`mean +/- 2*sd` pointwise) to the fit.
  Full smash / HMM ports stay deferred. New vignette
  `post_processing.Rmd` shows before / after.
- **Color palette**: Okabe-Ito (colorblind-friendly, modern).
  Vector: `c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
  "#D55E00", "#CC79A7", "#999999")`. Scales by recycling.
- **Multi-outcome layout**: dense grid via
  `par(mfrow = c(rows, cols))` with `rows * cols >= M + 1`,
  `cols = ceiling(sqrt(M + 1))`. Top-left panel is the shared PIP
  plot; remaining panels are per-outcome effect plots.
- **Real-data examples**: trim CASS4 (`data/fig_4_data/CASS4/`) and
  CR1/CR2 (`data/fig_4_data/CR1_CR2_obj.RData`) from
  `~/GIT/fsusie-experiments`. Each trimmed to ~250 KB by
  subsampling variables + individuals + positions to a
  vignette-scale example dataset; total package size budget ~500 KB.
- **De-identification (binding)**: the build scripts under
  `data-raw/` SHALL strip every personally identifying field from
  the genotype matrix `X` and the molecular phenotype matrix `Y`
  before saving. Concretely: `rownames(X)` and `rownames(Y)` are
  replaced with anonymous indices `S001, S002, ...`; any
  `colnames` carrying participant IDs (e.g., subject IDs in
  pheno matrices) are stripped or replaced; covariates with PHI
  (sex, age, batch, cohort) are dropped. Only the SNP IDs (chr,
  position, alleles) and CpG positions are retained as
  identifiers, since those are public reference information.
  A `de_identified = TRUE` attribute is set on each saved object
  for traceability.

## Plot function structure

Three layers, base R only:

1. **Low-level helpers** in `R/fsusie_plot.R`:
   - `mf_cs_colors(L, palette = "fsusie")` — color vector
     keyed by CS index. fsusie palette mirrors `fsusieR`'s for
     visual continuity with the paper figures; alternative palettes
     are easy to add.
   - `mf_pip_colors(fit)` — per-variable color (CS color if in a
     CS, grey otherwise).
   - `mf_effect_curve(fit, l, m)` — returns the length-T_m position-
     space curve for effect l, outcome m. Reads
     `fit$effect_curves[[m]][[l]]` if present, else inverts the
     wavelet representation via the existing `coef.mfsusie()` path
     (so callers don't need to run the post-processor first).
   - `mf_credible_band(fit, l, m)` — returns NULL if absent, else a
     T_m x 2 matrix `[lower, upper]`.

2. **Single-outcome plot** `fsusie_plot(fit, ...)`:
   - Two panels stacked via `par(mfrow = c(2, 1))`: PIP (top) +
     functional effect (bottom).
   - PIP panel: identical layout to `susie_plot(fit, y = "PIP")`.
   - Effect panel: one curve per CS, colored by `mf_cs_colors`.
     Credible band drawn as a transparent ribbon via
     `polygon()` with `col = adjustcolor(cs_color, alpha.f = 0.3)`.
     Affected-region shading via `rect()` along the x-axis.

3. **Multi-outcome plot** `mfsusie_plot(fit, m = NULL, ...)`:
   - Tile via `par(mfrow = ...)`. Default to `c(M+1, 1)`: top
     panel is the shared PIP plot (across-outcome PIPs), bottom
     M panels are per-outcome effect curves.
   - When `m` is a single index, plot just that outcome's effect
     panel.
   - Scalar outcomes (T_m = 1) render as the mvsusie-style point
     plot: variables on x-axis, effect on y-axis, dots colored by
     CS, with a horizontal CS bar at the top of the panel.

## Data files

- `data/fsusie_methyl_example.rda`: a list with elements
  `X` (n x p genotype matrix), `Y` (n x T_m methylation matrix),
  `pos` (length-T_m CpG positions), `signal_idx` (the true causal
  variable indices), `metadata` (`data.frame` with chr, gene name,
  source citation).
- `data/mfsusie_joint_example.rda`: a list with elements
  `X` (n x p), `Y` (length-M list of n x T_m matrices),
  `pos` (length-M list of position vectors), `signal_idx`
  (length-2 vector of true causal variable indices),
  `outcome_names` (character vector of length M).

Build scripts under `data-raw/`:

- `data-raw/make_fsusie_methyl_example.R` — pulls a real region
  from `~/GIT/fsusie-experiments/fig_1_data` if available; falls
  back to a fixed-seed simulation otherwise. Trims to ~500 KB.
- `data-raw/make_mfsusie_joint_example.R` — uses
  `mvf.susie.alpha::simu_effect_multfsusie` with seed 1, trimmed
  to ~500 KB.

Both scripts are committed but not run by R CMD check (they live
in `data-raw/` which is `.Rbuildignore`d).

## Why no ggplot

- mvsusie_plot uses ggplot + cowplot + ggrepel. That's heavy for a
  package that already has cpp11 + susieR + wavethresh. Base R
  graphics suffice for our two-panel layouts.
- Removes a runtime dependency; smaller install footprint; faster
  vignette knit.
- mvsusie_plot can stay ggplot-based; mfsusie_plot is independently
  designed.

## Plot API contract (lives in mf-plot/spec.md)

The new `inst/openspec/specs/mf-plot/spec.md` will document:
- `fsusie_plot(fit, ...)`: SHALL return invisibly, SHALL handle
  fits with or without post-processed `$effect_curves` /
  `$credible_bands` / `$lfsr_curves` slots, SHALL color CSes via
  `mf_cs_colors`.
- `mfsusie_plot(fit, m = NULL, ...)`: SHALL tile per-outcome
  panels when `m` is `NULL`, SHALL handle scalar outcomes
  (`T_m = 1`) and functional outcomes (`T_m > 1`) in the same
  call, SHALL accept either an `mfsusie` or an `mfsusie`
  derived-class fit.
- Both functions SHALL be base-R only (no `ggplot2`, no
  `cowplot`).

## What stays out of scope

- Post-processed smoothers (smash, TI, HMM) and TIWT credible
  bands. These are PR group 6b in the
  `add-mfsusier-s3-architecture` change. The plot functions
  draft draws curves only when those slots are absent; bands
  light up automatically when the smoother runs.
- Interactive / Plotly plots. Outside scope.
- `Gviz` integration. Outside scope (fsusieR has it; we defer).
