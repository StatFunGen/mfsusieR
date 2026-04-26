# mfsusieR

[![CI](https://github.com/StatFunGen/mfsusieR/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/StatFunGen/mfsusieR/actions/workflows/ci.yml)
[![Website](https://github.com/StatFunGen/mfsusieR/actions/workflows/dispatch_pkgdown_build.yml/badge.svg?branch=main)](https://statfungen.github.io/mfsusieR)
[![License](https://img.shields.io/badge/License-BSD_3--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)
[![GitHub release](https://img.shields.io/github/v/release/StatFunGen/mfsusieR?include_prereleases&sort=semver)](https://github.com/StatFunGen/mfsusieR/releases)

The `mfsusieR` package implements the multi-modality Sum of Single
Functions (mfSuSiE) model for joint fine-mapping of multiple
functional and scalar molecular phenotypes. It generalises the Sum
of Single Functions (fSuSiE) and the scalar Sum of Single Effects
(SuSiE) frameworks, sharing wavelet-based functional priors across
modalities so that a SNP affecting several spatially structured
traits is prioritised over a SNP affecting just one. Built on the
`susieR` backbone via S3 dispatch.

The package exposes two public entry points:

- `fsusie(Y, X, pos = NULL, ...)` for fine-mapping a **single**
  response (scalar or functional).
- `mfsusie(X, Y, pos = NULL, ...)` for fine-mapping **multiple**
  responses jointly. `Y` is a list of length `M` modalities; each
  element is a matrix `n x T_m` (with `T_m = 1` for scalar
  modalities).

## System Requirements

### Software Dependencies

- R (>= 4.0.0)
- Required R packages: `susieR` (>= 0.16.1), `wavethresh`, `mixsqp`,
  `ashr`, `matrixStats`
- Compiled code via `cpp11`

### Operating Systems

- macOS 12.0+
- Ubuntu 20.04+
- Windows 10+

### Hardware

No special hardware required. Standard desktop with 8 GB RAM is
sufficient for typical analyses.

## Installation

```R
# install.packages("remotes")
remotes::install_github("stephenslab/susieR")
remotes::install_github("StatFunGen/mfsusieR")
```

`mfsusieR` depends on the GitHub master of `susieR` for the per-
iteration S3 generics; the conda or CRAN release is older.

## Quick start

```R
library(mfsusieR)

# Single-modality functional fine-mapping
fit_f <- fsusie(Y, X, pos = positions, L = 5)

# Multi-modality joint fine-mapping (Y is a length-M list of matrices)
fit_m <- mfsusie(X, Y_list, L = 5,
                 prior_variance_scope = "per_modality")

# Standard SuSiE accessors and S3 methods
fit_m$pip
fit_m$sets$cs
predict(fit_m, newx = X_test)
coef(fit_m)
summary(fit_m)
```

See the [pkgdown website](https://statfungen.github.io/mfsusieR/) for
worked examples in both univariate and multivariate contexts.

## Citing this work

If you use `mfsusieR::fsusie()` (single-modality functional) in your
work, please cite:

> Denault, W.R.P., Sun, H., Carbonetto, P., Li, A., De Jager, L.P.,
> Bennett, D., The Alzheimer's Disease Functional Genomics
> Consortium, Wang, G. & Stephens, M. (2025). fSuSiE enables
> fine-mapping of QTLs from genome-scale molecular profiles.
> *bioRxiv* [DOI:
> 10.1101/2025.08.17.670732](https://doi.org/10.1101/2025.08.17.670732)

If you use `mfsusieR::mfsusie()` (multi-modality joint fine-mapping)
in your work, please cite:

> Liu, A., Sun, H., De Jager, L.P., Bennett, D., The Alzheimer's
> Disease Functional Genomics Consortium, Wang, G. & Denault, W.R.P.
> (2025). mfSuSiE enables multi-cell-type fine-mapping and multi-omic
> integration of chromatin accessibility QTLs in aging brain.
> *bioRxiv*. *(forthcoming; DOI to be assigned)*

For the underlying SuSiE backbone, please also cite:

> Wang, G., Sarkar, A., Carbonetto, P. & Stephens, M. (2020). A simple
> new approach to variable selection in regression, with application
> to genetic fine mapping. *Journal of the Royal Statistical Society,
> Series B* 82, 1273-1300.
> [DOI: 10.1111/rssb.12388](https://doi.org/10.1111/rssb.12388)

## License

BSD 3-Clause. See `LICENSE`.

## Issues and contributions

Please [file an issue](https://github.com/StatFunGen/mfsusieR/issues)
for bug reports, questions, or suggestions.
