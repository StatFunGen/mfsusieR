# Rewrite `fsusie_why_functional` vignette using the
# Limitation_SuSiE simulation and add an mvsusieR comparison

## Why

The current `vignettes/fsusie_why_functional.Rmd` was written
from scratch with two Gaussian-shaped causal effects. The
upstream functional fine-mapping motivation vignette used a
single causal SNP with a cosine effect on positions 41-128 of
T = 128, which is the simulation that mvsusieR's diagnostic
prototype (`mvsusieR/inst/prototypes/diagnose_mvsusie_issue_61.R`)
is also written against. Reusing the same simulation lets the
vignette compare per-CpG SuSiE, top-PC SuSiE, mvSuSiE (default
prior + true-prior settings per issue 61), and fSuSiE on
identical data, with effect-curve plots overlaying every
method's recovered effect onto the truth.

mvSuSiE is a SuSiE-suite peer; its inclusion is natural in a
"why functional" vignette and shows the user when each tool in
the suite is the right choice.

## What changes

### 1. Rewrite `vignettes/fsusie_why_functional.Rmd`

- Authors: Gao Wang, Anjing Liu, William Denault (existing
  author list retained).
- Adopt the upstream simulation verbatim (same recipe, same
  seed): `n = 100` from `susieR::N3finemapping$X[1:100, ]`,
  `T = 128`, single causal SNP at position 700, cosine effect
  on positions 41-128. The numerical output for the methods
  that overlap with upstream (per-CpG SuSiE, fSuSiE) SHALL
  match upstream at tolerance `<= 1e-12`.
- Methods compared, each with PIPs and recovered effect curve:
  1. **Per-CpG SuSiE.** Run `susie(X, y[, t])` for each
     `t in 1..T`, aggregate PIPs by per-variable max. Show
     fragmented CSes.
  2. **Top-PC SuSiE for `pc in 1..5`.** Project Y onto top-`pc`
     principal components, run `susie(X, Y_proj)` on the
     stacked columns. Five separate fits. Most fail because
     the principal components do not align with the variable
     effect direction.
  3. **mvSuSiE (default prior + estimate residual variance).**
     The configuration that gives all-zero PIPs per issue 61.
     Show the failure mode.
  4. **mvSuSiE (true outer-product prior + fixed identity
     residual variance).** The configuration that recovers
     the correct CS per issue 61. Plot recovered effect curve
     via `coef(fit_mv)[lead_snp + 1, ]` overlaid on truth.
  5. **fSuSiE.** Single `fsusie()` call. Plot via
     `mfsusie_plot(fit_s)` after `mf_post_smooth()`.

### 2. Add mvsusieR dependency

- `pixi.toml [dependencies]`: add `r-mvsusier = ">=0.2.0"`
  (cross-platform; dnachun ships both linux-64 and osx-arm64
  builds; no Bioconductor deps).
- `DESCRIPTION` Suggests: add `mvsusieR (>= 0.2.0)`.
- Vignette chunks that use mvsusieR gate with
  `eval = requireNamespace("mvsusieR", quietly = TRUE)`.

### 3. No new R code

The vignette uses public APIs from `mfsusieR`, `susieR`, and
`mvsusieR`. No package code changes; no test changes.

The mvsusie effect-curve plot is an inline base-R `plot()` call
in the vignette using `coef(fit_mv)`; we do not extend our
public plot API for cross-package fits.

### 4. Excalidraw design diagram

`design/methods-comparison-flow.excalidraw` shows the five
methods and where each succeeds / fails on the cosine
simulation. Renders via the official Excalidraw MCP in Claude
Desktop.

## Impact

- Rewritten: `vignettes/fsusie_why_functional.Rmd` (one file).
- DESCRIPTION: `Suggests: mvsusieR (>= 0.2.0)`.
- pixi.toml: `r-mvsusier = ">=0.2.0"` under cross-platform
  `[dependencies]`.
- No R code changes, no test changes, no public API changes.
- Specs: new `inst/openspec/specs/mf-vignette-why-functional/spec.md`
  documenting the vignette content contract.

## Out of scope

- Any change to `mfsusie_plot()` or `mfsusie_plot_lfsr()`. The
  mvsusie effect-curve overlay is inline `plot()`.
- Any port-source mention in the vignette. mvsusieR is a
  SuSiE-suite peer, NOT a port source; we name it directly.
  fsusieR is a port source and stays unmentioned.
- Issue 61 follow-up code in mvsusieR. The vignette uses
  whatever mvsusieR ≥ 0.2.0 ships.
