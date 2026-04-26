# mfsusieR

[![CI](https://github.com/StatFunGen/mfsusieR/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/StatFunGen/mfsusieR/actions/workflows/ci.yml)
[![Website](https://github.com/StatFunGen/mfsusieR/actions/workflows/dispatch_pkgdown_build.yml/badge.svg?branch=main)](https://statfungen.github.io/mfsusieR)
[![Coverage](https://img.shields.io/endpoint?url=https://statfungen.github.io/mfsusieR/coverage.json)](https://github.com/StatFunGen/mfsusieR/actions/workflows/ci.yml)
[![License](https://img.shields.io/badge/License-BSD_3--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)
[![GitHub release](https://img.shields.io/github/v/release/StatFunGen/mfsusieR?include_prereleases&sort=semver)](https://github.com/StatFunGen/mfsusieR/releases)

The `mfsusieR` package implements multi-outcome functional regression
using the Sum of Single Effects (mfSuSiE) model. The method is general
for Bayesian variable selection in functional-regression problems, but
was motivated and developed for fine-mapping spatially correlated
genomic features (e.g., DNA methylation profiles, chromatin
accessibility tracks) where the outcome is a curve along the genome
and a sparse set of variables affects several correlated positions
jointly. Built on the `susieR` backbone via S3 dispatch.

The package exposes two main functions:

- `fsusie(Y, X, pos = NULL, ...)` for fine-mapping a **single**
  response. `Y` may be scalar or functional (a matrix sampled at a
  grid of positions).
- `mfsusie(X, Y, pos = NULL, ...)` for fine-mapping **multiple**
  responses jointly. `Y` is a list of length `M` outcomes; each
  element is a matrix `n x T_m` (with `T_m = 1` for scalar outcomes
  and `T_m > 1` for functional outcomes; `T_m` may differ across
  outcomes).

When `Y` is scalar, `fsusie()` reduces to a version of the SuSiE
model. We expose this case here for two reasons: (i) it provides a
sanity-check path against `susieR` (the C1 contract test suite locks
exact element-wise equivalence), and (ii) it allows scalar outcomes
to be fit jointly with functional outcomes through `mfsusie()`. Users
analyzing a single scalar response on its own should use the
[`susieR`](https://github.com/stephenslab/susieR) package directly.
Users with **multiple correlated scalar outcomes** without a spatial
structure (e.g., several QTL traits across tissues) should use
[`mvsusieR`](https://github.com/stephenslab/mvsusieR) instead.
mfsusieR is the right choice when one or more outcomes are
functional, or when functional and scalar outcomes are jointly
modelled.

## Installation

```R
# install.packages("remotes")
remotes::install_github("stephenslab/susieR")
remotes::install_github("StatFunGen/mfsusieR")
```

`mfsusieR` currently depends on the GitHub master of `susieR` for the
per-iteration S3 generics. CRAN and conda releases of both packages
are planned; once those land the dependency on the GitHub master will
be dropped and a single `install.packages()` (or conda install) will
be sufficient.

## Quick start

See the [pkgdown website](https://statfungen.github.io/mfsusieR/),
in particular the
[Getting Started](https://statfungen.github.io/mfsusieR/articles/getting_started.html)
vignette, for worked examples in both single- and multi-outcome
contexts.

## Citing this work

If you use `mfsusieR::fsusie()` (single-outcome functional
fine-mapping) in your work, please cite:

> Denault, W.R.P., Sun, H., Carbonetto, P., Liu, A., De Jager, L.P.,
> Bennett, D., The Alzheimer's Disease Functional Genomics
> Consortium, Wang, G. & Stephens, M. (2025). fSuSiE enables
> fine-mapping of QTLs from genome-scale molecular profiles.
> *bioRxiv* [DOI:
> 10.1101/2025.08.17.670732](https://doi.org/10.1101/2025.08.17.670732)

If you use `mfsusieR::mfsusie()` (multi-outcome joint fine-mapping)
in your work, please cite:

> Liu, A., Sun, H., De Jager, L.P., Bennett, D., The Alzheimer's
> Disease Functional Genomics Consortium, Wang, G. & Denault, W.R.P.
> (2025). mfSuSiE enables multi-cell-type fine-mapping and multi-omic
> integration of chromatin accessibility QTLs in aging brain.
> *bioRxiv*. *(in preparation; DOI will be added on submission)*

For the underlying SuSiE backbone and engineering improvements that
mfsusieR builds on, please also cite:

> McCreight, A., Cho, Y., Li, R., Nachun, D., Gan, H-Y., Carbonetto,
> P., Stephens, M., Denault, W.R.P. & Wang, G. (2025). SuSiE 2.0:
> improved methods and implementations for genetic fine-mapping and
> phenotype prediction. *Submitting to Genome Biology*.

## License

BSD 3-Clause. See `LICENSE`.

## Issues and contributions

Please [file an issue](https://github.com/StatFunGen/mfsusieR/issues)
for bug reports, questions, or suggestions.
