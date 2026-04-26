# Tasks

Depends on: `add-preprocessing-low-count-and-quantile-norm`
and `expose-nig-residual-variance-method` archived first.

## 1. Data preparation

- [ ] 1.1 `data-raw/make_practical_dataset.R` (new). Load the
      rds at `~/Documents/fsusie_test/` and the analysis
      function. Select 2-3 representative regions exhibiting
      a clean small-sample signal.
- [ ] 1.2 For each selected region, run `mfsusie()` four
      times (default MoM, NIG, low-count-filter, quantile-
      norm). Save fitted-model objects to
      `data/practical_fits.rda`.
- [ ] 1.3 Confirm package-size impact remains reasonable
      (target: bundled rda <= 5 MB; if larger, downsample
      regions or store summary fields only).
- [ ] 1.4 Add `R/data.R` roxygen entry for `practical_fits`
      with field documentation and the not-redistributable
      note.

## 2. Vignette

- [ ] 2.1 `vignettes/practical_data_applications.Rmd` (new).
      Authors: Anjing Liu, Gao Wang, William Denault.
- [ ] 2.2 Section: small-sample problem framing.
- [ ] 2.3 Section: NIG correction with brief math
      derivation. Demonstration block loads
      `data(practical_fits)`, plots PIPs MoM vs NIG.
- [ ] 2.4 Section: `low_count_filter` demo. PIP comparison
      before vs after.
- [ ] 2.5 Section: `quantile_norm` demo. Before vs after
      comparison.
- [ ] 2.6 Audit: vignette never displays raw `X`, `Y`, or
      individual-level data; only fitted-model summaries
      and plots.

## 3. pkgdown

- [ ] 3.1 Add `practical_data_applications` to the articles
      list in `_pkgdown.yml`.

## 4. Spec delta

- [ ] 4.1 Create
      `inst/openspec/changes/add-practical-data-applications-vignette/specs/mf-vignette-practical/spec.md`
      with the no-raw-data and section-coverage
      requirements.

## 5. Build + archive

- [ ] 5.1 Vignette renders cleanly under pixi env on
      linux-64 (smashr / coloc / Bioconductor deps
      available).
- [ ] 5.2 Vignette renders on osx-arm64 with appropriate
      requireNamespace gates (where needed).
- [ ] 5.3 Push, CI green.
- [ ] 5.4 `openspec archive add-practical-data-applications-vignette`.
