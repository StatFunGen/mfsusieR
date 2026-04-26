# Tasks

## 1. Plot API

- [ ] 1.1 Re-export `susieR::susie_plot` from mfsusieR (add to
      `R/mfsusieR-package.R` `@importFrom` + `@export susie_plot`).
- [ ] 1.2 Implement `R/fsusie_plot.R` with `fsusie_plot(fit, ...)`.
      Base R only. Args: `fit`, `effect = "all"`, `cred_band = TRUE`,
      `show_affected_region = TRUE`, `lfsr_curve = TRUE`,
      `show_grid_dots = FALSE`, `pos_label = NULL`, `legend = TRUE`.
- [ ] 1.3 Implement `R/mfsusie_plot.R` with
      `mfsusie_plot(fit, m = NULL, ...)`. Base R only. Tiles
      `par(mfrow = ...)` over outcomes when `m` is `NULL`.
      Functional outcomes get curves; scalar outcomes get the
      mvsusie-style point + CS-bar plot.
- [ ] 1.4 Add S3 method `plot.mfsusie(x, ...)` in
      `R/mfsusie_methods.R`. Dispatches: `M == 1` -> `fsusie_plot`,
      else -> `mfsusie_plot`.
- [ ] 1.5 Color palette helper `mf_cs_colors(L)` returning a
      length-L vector matching `fsusieR`'s palette
      (`c("black", "dodgerblue2", "green4", "#6A3D9A", "#FF7F00",
        "darkorchid", "khaki2", "khaki4", "deeppink1")`).
- [ ] 1.6 Tests: smoke tests in
      `tests/testthat/test_plot.R`. Each public plot call wrapped
      in `expect_silent()` against a tiny fit. Verify the function
      returns invisibly (`invisible(NULL)` per base-R convention).

## 2. Post-processed effect storage

- [ ] 2.1 Add (optional) `fit$effect_curves`, `fit$credible_bands`,
      `fit$lfsr_curves` slots to the documented mfsusie fit
      contract. Slots are populated by `mf_post_smooth()` (PR
      group 6b) and read by the plot functions.
- [ ] 2.2 Document the slots in the `mfsusie()` `@return` block.
      When absent, plot functions fall back to the wavelet-domain
      effect through `coef.mfsusie()`.

## 3. Packaged example data

- [ ] 3.1 Create `data-raw/make_fsusie_methyl_example.R`.
      Builds a small DNAm fixture from
      `fsusie-experiments/fig_1_data` (or fixed-seed simulation).
      Saves to `data/fsusie_methyl_example.rda`.
- [ ] 3.2 Create `data-raw/make_mfsusie_joint_example.R`. Builds a
      multi-outcome fixture (M=5: 2 functional + 3 scalar) using
      `mvf.susie.alpha::simu_effect_multfsusie` with fixed seed,
      saved to `data/mfsusie_joint_example.rda`.
- [ ] 3.3 Document data sets in `R/data.R` with `@docType data`
      blocks; one per dataset.
- [ ] 3.4 Add `data-raw/` to `.Rbuildignore`.
- [ ] 3.5 Verify package size remains reasonable
      (`tools::checkBuiltPackage()` size summary).

## 4. Vignette refresh

- [ ] 4.1 `getting_started.Rmd`: load `data(N3finemapping)` from
      susieR, plus `data(fsusie_methyl_example)` for an
      `fsusie_plot()` preview, plus `data(mfsusie_joint_example)`
      for an `mfsusie_plot()` preview. Keep the `susieR` C1
      numerical identity demo.
- [ ] 4.2 `fsusie_intro.Rmd`: rewrite using
      `data(fsusie_methyl_example)` + `fsusie_plot()`. Mirror the
      structure of `fsusieR::fsusie_intro.Rmd`.
- [ ] 4.3 `fsusie_covariates_and_coloc.Rmd`: rewrite using
      packaged data + `fsusie_plot()`. Mirror
      `fsusieR::Adjusting_covariate.Rmd` + `Coloc_fsusie.Rmd`.
- [ ] 4.4 `fsusie_dnam_case_study.Rmd`: rewrite using
      `data(fsusie_methyl_example)` + `fsusie_plot()`. Mirror
      `fsusieR::methyl_demo.Rmd`.
- [ ] 4.5 `fsusie_why_functional.Rmd`: keep the simulation-based
      structure (it demonstrates a contrast); add `fsusie_plot()`
      and `susie_plot()` for visual comparison.
- [ ] 4.6 `mfsusie_intro.Rmd`: rewrite using
      `data(mfsusie_joint_example)` + `mfsusie_plot()`. Mirror
      `mvf.susie.alpha::joint_functional_and_univariate_fine_mapping.Rmd`.
- [ ] 4.7 `mfsusie_long_running_fits.Rmd`: keep current; add an
      `mfsusie_plot()` preview at the end.

## 5. Spec + NAMESPACE

- [ ] 5.1 Add `inst/openspec/specs/mf-plot/spec.md` documenting the
      plot API contract.
- [ ] 5.2 NAMESPACE: `export(fsusie_plot)`,
      `export(mfsusie_plot)`, `S3method(plot, mfsusie)`,
      `export(susie_plot)` (re-exported from susieR).
- [ ] 5.3 R CMD check (CI) passes; vignettes knit cleanly.

## 6. Build + archive

- [ ] 6.1 `devtools::test()` passes (smoke tests for plot
      functions).
- [ ] 6.2 `devtools::document()` regenerates man + NAMESPACE.
- [ ] 6.3 Push, let CI run on the change branch.
- [ ] 6.4 Archive the OpenSpec change once everything is green.
