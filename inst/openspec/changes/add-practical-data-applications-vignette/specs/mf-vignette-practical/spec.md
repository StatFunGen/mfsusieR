# mf-vignette-practical capability

## ADDED Requirements

### Requirement: vignette demonstrates the three round-3 features on real data

`vignettes/practical_data_applications.Rmd` SHALL contain sections demonstrating `mfsusie(..., small_sample_correction = TRUE)`, `mfsusie(..., low_count_filter = c)`, and `mfsusie(..., quantile_norm = TRUE)` on a real-world multi-outcome ATAC-seq example dataset. Each section SHALL provide a PIP comparison between the option enabled and the default.

#### Scenario: section coverage

A section-header grep over the vignette returns at least one header for each of: small-sample correction, low-count filtering, quantile normalization. Each section SHALL contain at least one code chunk that loads `data(practical_fits)` and at least one plot or table comparing the option enabled vs the default.

### Requirement: no raw data shipped

The vignette SHALL NOT display individual-level genotype, phenotype, cell-count, or read-count data. Only fitted-model summaries (PIPs, credible sets, smoothed effect curves) and plots derived from those summaries SHALL appear in the rendered output.

#### Scenario: raw-data audit

A grep over the vignette source for `head(X)`, `head(Y)`, `print(X)`, `print(Y)`, or `str(Y_f)` returns zero matches. The fitted-model object loaded via `data(practical_fits)` SHALL contain only summary fields and SHALL NOT contain `X` or `Y` matrices from the source dataset.

### Requirement: small-sample correction section explains the Johnson-t marginal

The small-sample correction section of the vignette SHALL include a brief mathematical statement of the Johnson 2005 scaled Student-t marginal Bayes factor and the role of the degrees-of-freedom parameter `df = n - 1`. The statement SHALL contrast the Wakefield Normal marginal with the Student-t marginal and cite the source of the correction (`mvf.susie.alpha::multfsusie(..., cor_small = TRUE)` and Johnson 2005 JRSSB).

#### Scenario: derivation present

The vignette source contains at least one display-math block (e.g., `$$ ... $$`) under the small-sample correction section that references either the Student-t marginal density or the contrast with the Normal marginal.
