# mf-coverage capability

## ADDED Requirements

### Requirement: package coverage targets

`covr::package_coverage()` SHALL report:

- Overall coverage `>= 95%`.
- Per-file coverage `>= 90%` for every file in `R/` except
  `R/zzz.R` (exempt; `.onLoad` hooks are untestable through
  covr).

#### Scenario: coverage measurement

A local run of
`R -e 'cat(covr::percent_coverage(covr::package_coverage(quiet = TRUE)))'`
returns a value `>= 95`. The same run with
`covr::file_coverage` reports per-file percentages
`>= 90` for all files except `R/zzz.R`.

### Requirement: predict / coef fidelity reference test exists

`tests/testthat/test_predict_coef_fidelity.R` SHALL assert
bit-identity at tolerance `<= 1e-12` between
`mfsusieR::predict.mfsusie` / `coef.mfsusie` (M=1 case) and
the corresponding upstream functional fine-mapping methods,
gated on `skip_if_not_installed("fsusieR")`.

#### Scenario: file present and exercised

The file exists and is exercised in
`devtools::test(filter = "predict_coef_fidelity")` without
errors when fsusieR is installed.
