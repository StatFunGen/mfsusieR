# mf-vignette-practical capability

## ADDED Requirements

### Requirement: vignette demonstrates the three round-3 features on real data

`vignettes/practical_data_applications.Rmd` SHALL contain
sections demonstrating
`mfsusie(..., estimate_residual_method = "NIG")`,
`mfsusie(..., low_count_filter = c)`, and
`mfsusie(..., quantile_norm = TRUE)` on a real-world
multi-outcome ATAC-seq example dataset. Each section SHALL provide a
PIP comparison between the option enabled and the default.

#### Scenario: section coverage

A section-header grep over the vignette returns at least one
header for each of: NIG correction, low-count filtering,
quantile normalization. Each section SHALL contain at least
one code chunk that loads `data(practical_fits)` and at
least one plot or table comparing the option enabled vs the
default.

### Requirement: no raw data shipped

The vignette SHALL NOT display individual-level genotype,
phenotype, cell-count, or read-count data. Only fitted-model
summaries (PIPs, credible sets, smoothed effect curves) and
plots derived from those summaries SHALL appear in the
rendered output.

#### Scenario: raw-data audit

A grep over the vignette source for `head(X)`, `head(Y)`,
`print(X)`, `print(Y)`, or `str(Y_f)` returns zero matches.
The fitted-model object loaded via `data(practical_fits)`
SHALL contain only summary fields and SHALL NOT contain `X`
or `Y` matrices from the source dataset.

### Requirement: NIG section includes a brief mathematical derivation

The NIG section of the vignette SHALL include a brief
mathematical derivation of the Normal-Inverse-Gamma posterior
on `sigma^2` and the resulting form of the SER Bayes factor
when `sigma^2` is integrated out. The derivation SHALL match
the SuSiE manuscript notation and cite the manuscript by file
path.

#### Scenario: derivation present

The vignette source contains at least one display-math block
(e.g., `$$ ... $$`) under the NIG section that references
either the prior `Inverse-Gamma(alpha_0, beta_0)` form or the
posterior conjugate update.
