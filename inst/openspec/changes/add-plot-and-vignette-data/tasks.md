# Tasks

## 1. Plot API

- [x] 1.1 Implement `R/mfsusie_plot.R` with
      `mfsusie_plot(fit, m = NULL, ...)`. Base R only. Tiles
      `par(mfrow = ...)` over outcomes when `m` is `NULL`.
- [x] 1.2 Add S3 method `plot.mfsusie(x, ...)` in
      `R/mfsusie_methods.R` dispatching to `mfsusie_plot()`.
- [x] 1.3 Color palette helper `mf_cs_colors(L)` (Okabe-Ito).
- [x] 1.4 Smoke tests for `mfsusie_plot()`.
- [ ] 1.5 Add `effect_style = c("band", "errorbar")` to
      `mfsusie_plot()`. `"band"` keeps the current ribbon;
      `"errorbar"` draws a per-position dot with vertical
      bar to the credible band. Default `"band"`. Fails fast
      when `fit$credible_bands` is unpopulated and
      `effect_style = "errorbar"`.
- [ ] 1.6 Add `facet_cs = c("auto", "stack", "overlay")` to
      `mfsusie_plot()`. `"overlay"` keeps the current single
      panel. `"stack"` draws one row per CS within the effect
      panel using `layout()` rather than nested `par(mfrow)`.
      `"auto"` defaults to `"stack"` when `length(cs) >= 3`
      OR when affected-region masks (band excludes zero) are
      pairwise disjoint; `"overlay"` otherwise.
- [ ] 1.7 New `R/mfsusie_plot_lfsr.R` with
      `mfsusie_plot_lfsr(fit, lfsr_threshold = 0.01,
      truth = NULL, cex_max = 6, ...)`. No `m` argument.
      Bubble grid: rows = CSes, columns = positions, dot
      size = `-log10(lfsr)`, color = `lfsr <= threshold`
      (or `truth` when supplied). Handles M = 1 (single
      panel) and M > 1 (one panel per outcome, layout policy
      mirrors `mfsusie_plot()`). Errors when
      `fit$lfsr_curves` is unpopulated.
- [ ] 1.8 Change `mfsusie_plot()` default `lfsr_threshold`
      from `0.05` to `0.01` (and propagate to docs).
- [ ] 1.9 Smoke tests in `tests/testthat/test_plot.R`:
      `effect_style` x `facet_cs` x (post-smoothed / not),
      `mfsusie_plot_lfsr()` for M = 1 and M > 1, with and
      without `truth`. Each call wrapped in
      `expect_silent()` against a tiny fit.

## 2. Post-processed effect storage

- [x] 2.1 `fit$effect_curves`, `fit$credible_bands`,
      `fit$lfsr_curves` slots populated by `mf_post_smooth()`.
- [x] 2.2 Document the slots in `mfsusie()` `@return`.

## 3. Packaged example data (enhance, no new datasets)

- [ ] 3.1 Update `data-raw/make_data.R`: enhance
      `dnam_example` to `n = 100`, `p = 12`, `T = 32` with
      three causal SNPs producing two CpG clusters. Include
      a length-T boolean `truth_mask` per CS so the
      methylation vignette can pass `truth = ...` to
      `mfsusie_plot_lfsr()`.
- [ ] 3.2 Update `data-raw/make_data.R`: bump `rnaseq_example`
      from one to two causal SNPs at positions 25 and 75
      with smooth gene-body shapes via the same
      `simu_IBSS_per_level()`-style construction the intro
      vignette demonstrates.
- [ ] 3.3 Regenerate `.rda` files; update roxygen blocks in
      `R/data.R` with the new sizes / causal counts.
- [ ] 3.4 Verify package size remains reasonable.

## 4. Vignette merges

- [ ] 4.1 Merge `vignettes/fsusie_intro.Rmd`. Single author
      William Denault. Sections: overview / math model;
      example dataset (`data(rnaseq_example)`); fitting;
      inspecting the fit (`mfsusie_plot()`); smoothed effect
      curves (`mf_post_smooth(method = "TI")`); uneven
      sampling positions (subsample inline); prior choices
      (`prior_variance_scope = "per_scale"` vs
      `"per_outcome"`, scatter of PIPs, runtime comparison).
      No mention of any port-source package. No
      `susie_plot()`.
- [ ] 4.2 Merge `vignettes/fsusie_dnam_case_study.Rmd`.
      Authors: William Denault and Gao Wang. Sections:
      overview; example dataset (`data(dnam_example)`);
      data inspection (ggplot2-gated); per-CpG association
      tests + GWAS-by-SNP / GWAS-by-CpG plots
      (ggplot2-gated); per-CpG `susie()` failure-mode demo
      (ggplot2-gated); single fSuSiE fit; `mfsusie_plot()`;
      `mfsusie_plot_lfsr(truth = ...)`;
      `mf_post_smooth(method = "TI")` followed by
      `mfsusie_plot(effect_style = "errorbar", facet_cs =
      "stack")`; prediction (retain). No port-source
      mention. No `susie_plot()`.
- [ ] 4.3 Audit other vignettes
      (`fsusie_covariates_and_coloc`, `fsusie_why_functional`,
      `mfsusie_intro`, `mfsusie_long_running_fits`,
      `post_processing`, `getting_started`, `fsusie_gtex_case_study`)
      and remove any `susie_plot()` calls and any port-source
      references. Pull the new plot-API options through where
      they help.

## 5. Spec deltas

- [ ] 5.1 Add `inst/openspec/specs/mf-plot/spec.md` documenting
      the plot API: `mfsusie_plot()` requirements (PIP panel,
      effect panel, `effect_style`, `facet_cs`,
      `show_lfsr_curve`, `lfsr_threshold` default 0.01) and
      `mfsusie_plot_lfsr()` requirements (no `m`, M = 1 / M > 1
      handling, `truth` semantics).

## 6. CI / pixi

- [ ] 6.1 Add `r-ggplot2`, `r-cowplot`, `r-reshape2` to
      `pixi.toml` only (NOT to DESCRIPTION Suggests).
      Vignette chunks gate with
      `eval = requireNamespace("ggplot2", quietly = TRUE)`.

## 7. Verification (developer-only, never in vignettes)

- [ ] 7.1 Confirm `prior_variance_scope = "per_scale"` and
      `"per_outcome"` compute identical priors to
      `mixture_normal_per_scale` and `mixture_normal` in the
      port source (`R/prior_scale_mixture.R::init_scale_mixture_prior_default`
      already verified). Document the mapping in
      `inst/notes/refactor-exceptions.md` if not already
      there. No vignette text mentions either source.
- [ ] 7.2 Add unit-test coverage for the smooth-shape
      simulation used by `rnaseq_example` (apple-to-apple
      PIP / alpha / mu / lbf at the documented contract
      tolerance) if not already covered. Same for the
      enhanced `dnam_example` (`n = 100`, `T = 32`, three
      causal SNPs).

## 8. Build + archive

- [ ] 8.1 `devtools::test()` passes (smoke + apple-to-apple).
- [ ] 8.2 `devtools::document()` regenerates man + NAMESPACE.
- [ ] 8.3 Vignettes render cleanly under pixi env.
- [ ] 8.4 Push, let CI run.
- [ ] 8.5 Archive once everything is green.
