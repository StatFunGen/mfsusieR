# mf-plot capability

## ADDED Requirements

### Requirement: `mfsusie_plot()` is the unified PIP + effect-curve plot

`mfsusie_plot()` SHALL be the public entry point for visualising
PIPs and per-CS effect curves on an `mfsusie` / `fsusie()` fit.
It SHALL handle single-outcome (`M = 1`) and multi-outcome
(`M > 1`) fits transparently and SHALL NOT depend on `ggplot2`,
`cowplot`, or any non-base graphics package.

#### Scenario: signature

```
mfsusie_plot(fit, m = NULL, pos = NULL,
             effect_style   = c("band", "errorbar"),
             facet_cs       = c("auto", "stack", "overlay"),
             show_grid_dots = FALSE,
             show_lfsr_curve = TRUE,
             show_affected_region = TRUE,
             lfsr_threshold = 0.01,
             lwd = 1.5, add_legend = TRUE, ...)
```

The function SHALL return `invisible(NULL)`.

#### Scenario: layout

When `m` is `NULL`:
- `M == 1` -> 2-panel column (PIP top, effect bottom).
- `M > 1`  -> dense grid (PIP top-left, M effect panels).

When `m` is a positive integer in `1..M` -> single effect panel
for outcome `m`, no PIP panel.

#### Scenario: effect_style

`effect_style = "band"` (default) SHALL draw the per-CS mean
effect curve as a solid line and the credible band as a
translucent ribbon when `fit$credible_bands[[m]][[l]]` is
populated; without credible bands it draws the mean curve only.

`effect_style = "errorbar"` SHALL draw a per-position dot at
the mean and a vertical bar to the credible band lower / upper
endpoints. It SHALL emit a clear error when
`fit$credible_bands` is unpopulated.

#### Scenario: facet_cs

`facet_cs = "overlay"` SHALL draw all CSes on a single effect
panel.

`facet_cs = "stack"` SHALL draw one row per CS within the effect
panel using `layout()` so the rows share the position axis.

`facet_cs = "auto"` (default) SHALL select `"stack"` when
`length(fit$sets$cs) >= 3` OR when the affected-region masks
(positions where each CS's credible band excludes zero) are
pairwise disjoint; otherwise `"overlay"`.

#### Scenario: lfsr secondary-axis overlay

When `fit$lfsr_curves[[m]][[l]]` is populated and
`show_lfsr_curve = TRUE`, the effect panel SHALL overlay each
CS's lfsr curve on a secondary y-axis with a dashed reference
line at `lfsr_threshold` (default `0.01`).

### Requirement: `mfsusie_plot_lfsr()` is the per-CS lfsr bubble grid

`mfsusie_plot_lfsr()` SHALL be a separate public function for
the per-CS lfsr bubble grid produced by HMM post-processing.
It SHALL handle `M = 1` and `M > 1` fits transparently and
SHALL NOT take an `m` argument.

#### Scenario: signature

```
mfsusie_plot_lfsr(fit,
                  lfsr_threshold = 0.01,
                  truth          = NULL,
                  cex_max        = 6, ...)
```

The function SHALL return `invisible(NULL)`.

#### Scenario: layout

`M == 1` -> single bubble panel: rows = CSes, columns =
positions, dot size = `-log10(fit$lfsr_curves[[1]][[l]])`,
clamped to `[0, cex_max]`.

`M > 1` -> one bubble panel per outcome, tiled with the same
layout policy as `mfsusie_plot()` (dense grid).

#### Scenario: color rule

When `truth = NULL` (default), dots SHALL be colored by
`lfsr <= lfsr_threshold` (two-color: above-threshold vs
at-or-below).

When `truth` is supplied, color SHALL reflect the ground-truth
mask:
- `M == 1`: `truth` is a length-`T_1` boolean vector.
- `M > 1`: `truth` is a length-`M` list, each entry a
  length-`T_m` boolean vector.

The `truth` argument is intended for vignette code where
ground truth is known by simulation; in real-data use it is
`NULL`.

#### Scenario: missing lfsr

The function SHALL emit a clear error when
`fit$lfsr_curves` is `NULL` or empty (i.e., the fit was not
post-processed with `mf_post_smooth(method = "HMM")`).

### Requirement: S3 `plot.mfsusie` method dispatches to `mfsusie_plot()`

`plot(fit)` on an `mfsusie` / `fsusie()` fit SHALL call
`mfsusie_plot(fit, ...)` and forward `...`.

#### Scenario: dispatch and arg forwarding

`plot(fit, effect_style = "errorbar")` SHALL produce the same
output as `mfsusie_plot(fit, effect_style = "errorbar")` and
SHALL return `invisible(NULL)`.

### Requirement: vignettes use `mfsusie_plot()` and `mfsusie_plot_lfsr()` exclusively

Package vignettes SHALL use `mfsusie_plot()` and
`mfsusie_plot_lfsr()` for fit visualization. They SHALL NOT
call `susieR::susie_plot()` or any non-package plot
function on a fit object. They SHALL NOT mention port-source
packages by name in any plot caption, axis label, or prose.

#### Scenario: vignette grep audit

A grep over `vignettes/*.Rmd` for `susie_plot`, `fsusieR`,
`mvf.susie.alpha`, or `plot_susiF` SHALL return zero matches.
Apple-to-apple comparisons against port sources SHALL appear
only under `tests/testthat/` and `inst/notes/`.
