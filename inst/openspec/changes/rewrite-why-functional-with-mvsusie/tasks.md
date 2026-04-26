# Tasks

## 1. Dependencies

- [ ] 1.1 `pixi.toml [dependencies]`: add `r-mvsusier = ">=0.2.0"`
      under cross-platform deps (dnachun ships both linux-64
      and osx-arm64; no Bioconductor deps).
- [ ] 1.2 `DESCRIPTION` Suggests: add `mvsusieR (>= 0.2.0)`.
- [ ] 1.3 Vignette chunks gate with
      `eval = requireNamespace("mvsusieR", quietly = TRUE)`.

## 2. Vignette rewrite

- [ ] 2.1 Open `vignettes/fsusie_why_functional.Rmd` and
      replace the body. Authors stay: Gao Wang, Anjing Liu,
      William Denault.
- [ ] 2.2 Section: simulation. Cosine effect on positions
      41..128 of T = 128, single causal SNP at position 700.
      Use `susieR::N3finemapping` with `n = 100`.
- [ ] 2.3 Section: per-CpG SuSiE. `susie(X, y[, t])` for each
      `t`, aggregate PIPs by `apply(pip_per_pos, 1, max)`.
      Plot fragmented PIP pattern.
- [ ] 2.4 Section: top-PC SuSiE for `pc in 1..5`. Five
      separate fits; show PIP for each.
- [ ] 2.5 Section: mvsusie default. `mvsusie(X, y, L = 1)`;
      print PIPs and CS; confirm the failure mode.
- [ ] 2.6 Section: mvsusie with true outer-product prior +
      fixed identity residual variance. The recovery
      configuration. Print PIPs and CS; plot effect curve via
      `coef(fit_mv)[lead + 1, ]` overlaid on truth.
- [ ] 2.7 Section: fsusie. Single `fsusie(y, X, L = 1)` plus
      `mf_post_smooth(fit, method = "TI")`. Plot via
      `mfsusie_plot(fit_s)`.
- [ ] 2.8 Closing section: when to use which tool in the
      SuSiE suite. Brief, mechanical, no marketing.

## 3. Spec delta

- [ ] 3.1 Create `inst/openspec/changes/rewrite-why-functional-with-mvsusie/specs/mf-vignette-why-functional/spec.md`.
      Document the vignette content contract: simulation
      recipe, methods compared, no port-source mentions.

## 4. Excalidraw design diagram

- [ ] 4.1 `design/methods-comparison-flow.excalidraw`:
      five-method comparison fan from `(X, y, truth)` to
      success/failure outcomes.

## 5. Build + render + archive

- [ ] 5.1 Vignette renders cleanly under pixi env on linux-64
      (mvsusieR available).
- [ ] 5.2 Vignette renders cleanly under pixi env on
      osx-arm64 (mvsusieR available; alternatively gates
      cleanly if absent).
- [ ] 5.3 Push, let CI run.
- [ ] 5.4 `openspec archive rewrite-why-functional-with-mvsusie`
      once green.
