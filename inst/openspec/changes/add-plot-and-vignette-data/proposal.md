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

### 1. Plot API

`mfsusie_plot()` covers PIP + per-CS effect curves. A separate
`mfsusie_plot_lfsr()` covers the per-CS local-false-sign-rate
bubble grid (HMM-smoother output). Both are base R only, no
`ggplot2` / `cowplot` dependency in the package.

#### `mfsusie_plot(fit, ...)`

Auto-detects `M`:
- `M = 1` -> 2-panel column (PIP top, effect bottom)
- `M > 1` -> dense grid (PIP top-left, M effect panels)
- `m = i` -> single-outcome effect panel only (override)

Effect panel shape is parametrized:
- `effect_style = c("band", "errorbar")`. `"band"` (default):
  solid mean curve plus translucent ribbon for the credible
  band (current behavior). `"errorbar"`: per-position dot at
  the mean with a vertical bar to the credible band, mirroring
  the per-CS errorbar plots used in published methylation
  case-study figures.
- `facet_cs = c("auto", "stack", "overlay")`. `"overlay"`
  (current behavior): all CSes on one panel. `"stack"`: one
  row per CS within the effect panel. `"auto"` (default):
  stack when there are >= 3 CSes or when affected-region
  masks are pairwise disjoint; overlay otherwise.
- `show_lfsr_curve = TRUE` keeps the secondary-axis lfsr
  overlay for HMM-smoothed fits.
- `lfsr_threshold = 0.01` (changed from 0.05); applies to the
  secondary-axis dashed reference line.

The `m` argument is retained for single-outcome focus on
multi-outcome fits.

#### `mfsusie_plot_lfsr(fit, ...)`

Bubble grid: rows are CSes (or rows-per-outcome when M > 1),
columns are positions, dot size = `-log10(lfsr)`, color
indicates `lfsr <= lfsr_threshold` by default. The function
handles M = 1 and M > 1 transparently (no `m` argument); for
M > 1 it tiles one bubble panel per outcome with the same
layout policy as `mfsusie_plot()`.

The optional `truth` argument accepts a length-`T_m` boolean
vector (`M = 1`) or a length-`M` list of length-`T_m` boolean
vectors (`M > 1`) and recolors dots by ground truth. The
argument is intended for simulated data where the truly
affected positions are known by construction; in real-data
use `truth = NULL` and color follows the threshold rule.
`lfsr_threshold` defaults to `0.01`.

S3 method `plot.mfsusie()` dispatches to `mfsusie_plot()`.
Color palette: Okabe-Ito (modern, colorblind-friendly).

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

Existing example datasets stay; two are enhanced in place. No
new example datasets.

- `dnam_example` -- enhanced to `n = 100`, `p = 12`, `T = 32`
  with three causal SNPs producing two CpG clusters. The
  smaller, three-causal example is what the methylation
  case-study vignette uses for the per-CpG-susie versus
  fSuSiE comparison.
- `rnaseq_example` -- bumped from one to two causal SNPs at
  positions 25 and 75, with smooth random per-position
  effects generated by wavelet-coefficient shrinkage. Two CSes
  give the intro vignette enough structure to demonstrate
  multi-CS fits, post-smoothing, and prior comparison without
  adding a second example dataset.
- The intro vignette's uneven-sampling section subsamples
  positions inline (`pos <- sort(sample(...))`) on
  `rnaseq_example`; no second example dataset is needed.

Build scripts under `data-raw/`; `.Rbuildignore`d.

### 5. Vignette refresh

`fsusie_intro.Rmd` and `fsusie_dnam_case_study.Rmd` are
substantively merged with the corresponding fSuSiE-paper
case-study content. The vignettes never name or compare to
any port-source package; reproducibility versus port sources
is enforced by unit tests, not by visible vignette text.

- `fsusie_intro.Rmd` (single author: William Denault) merges
  the existing intro with the math / theory framing,
  uneven-sampling section, and prior-comparison section.
  Uses `data(rnaseq_example)` (two causals; see Section 4).
  Prior comparison: `prior_variance_scope = "per_scale"` vs
  `"per_outcome"` with a scatter of PIPs and a runtime
  comparison; no upstream-package mention.
- `fsusie_dnam_case_study.Rmd` (authors: William Denault and
  Gao Wang) merges the existing dnam case study with the
  small-scale methylation comparison: data inspection, GWAS by
  SNP / by CpG, per-CpG susie fits showing the failure mode,
  a single fSuSiE fit recovering all CSes, the lfsr bubble
  panel (`mfsusie_plot_lfsr()` with the `truth` mask passed
  in from the simulation), and the TI-smoother errorbar plot
  (`mfsusie_plot(effect_style = "errorbar", facet_cs = "stack")`).
  Uses `data(dnam_example)` (enhanced; see Section 4).
- Other vignettes keep their current structure with only the
  new plot-API options pulled through where useful.

ggplot2 / cowplot / reshape2 are vignette-only render
dependencies and are NOT in `DESCRIPTION` Suggests. They are
added to `pixi.toml` so the CI vignette renderer has them.
Vignette chunks that use them gate with
`eval = requireNamespace("ggplot2", quietly = TRUE)`.

Every vignette uses `mfsusie_plot()` and
`mfsusie_plot_lfsr()` for fit visualization. `susieR::susie_plot`
is not used in vignettes. Titles describe topics, not entry
points.

### 6. NAMESPACE / Imports

- `export(mfsusie_plot, mfsusie_plot_lfsr, mf_post_smooth)`.
- `S3method(plot, mfsusie)`.
- `@importFrom wavethresh av.basis convert nlevelsWT` (already
  imported).
- `ashr::ash` is already imported.
- No new package Imports. ggplot2 / cowplot / reshape2 are
  pixi-only (CI vignette render); not in DESCRIPTION.

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
