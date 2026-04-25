## Why

mfsusieR is a from-scratch refactor that replaces both
`mvf.susie.alpha::multfsusie` (multi-modality functional SuSiE) and
`fsusieR::susiF` (single-modality functional SuSiE) with one S3 package
on top of the `susieR` backbone. The single-modality case
(`M = 1, T_1 > 1`) is a special case of the multi-modality case
exposed through a thin wrapper. The univariate case
(`M = 1, T_1 = 1`) is a deeper degeneracy that reduces to
`susieR::susie`. Phase 1 paradigm notes (`inst/notes/paradigms/`)
catalogued the port sources and the two architecture references. The
math anchor is the manuscript at
`/home/gw/GIT/MultifSuSiE_Manuscript`. This change locks the
architecture so Phase 3 (port and implementation under the
Claude-native review loop) and Phase 4 (reference unit tests) have a
stable target.

The package is held to three apple-to-apple equivalence contracts.
Each contract has its own test file and its own tolerance. When a
port source has a known bug, mfsusieR fixes it, the corresponding
test asserts the deviation explicitly, and the OpenSpec change that
authorized the deviation is cited inline.

| Contract | Inputs | Reference | Tolerance |
|---|---|---|---|
| C1: scalar SuSiE | `M = 1`, `T_1 = 1`, single-component prior | `susieR::susie` | `<= 1e-10` |
| C2: single-modality functional SuSiE | `M = 1`, `T_1 > 1` | `fsusieR::susiF` | `<= 1e-8` |
| C3: multi-modality functional SuSiE | `M >= 1`, legacy variance mode | `mvf.susie.alpha::multfsusie` | `<= 1e-8` |

## What Changes

- **Two public entry points.** `mfsusieR::mfsusie(X, Y, pos, ...)` is the
  multi-modality entry; `mfsusieR::fsusie(Y, X, pos, ...)` is a thin
  wrapper that accepts a single-phenotype `(Y, X)` input in the same
  argument order as `fsusieR::susiF` and calls `mfsusie()` with `M = 1`.
  The wrapper is the user-facing migration path off fsusieR.
- **Single `mf_individual` data class.** Univariate (`T_m = 1`),
  single-modality (`M = 1`), and multi-modality with ragged `T_m` all
  flow through one class with shape branches, mirroring how
  `mvsusieR/refactor-s3` collapses `R = 1` vs `R > 1` inside one
  `mv_individual`.
- **Single `mfsusie` model class.** `class(fit) = c("mfsusie", "susie")`
  so existing `susieR` accessors with well-defined semantics for our
  outputs (e.g., `susieR::susie_get_pip`) reuse without override.
- **IBSS plug-in via S3 method registration on `mf_individual` and
  `mfsusie`,** following the `mvsusieR/refactor-s3` paradigm.
- **Two prior layers.** A per-(scale, modality) wavelet-coefficient
  scale-mixture-of-normals prior (manuscript `methods/derivation.tex`
  eq:additive_model with `pi_{k,s,m}`, `sigma_{k,s,m}^2`), and an
  optional cross-modality covariance prior (off by default). The
  susieAnn pluggable-predictor pattern (paradigm reference #2) is the
  template for the cross-modality plug-in seam.
- **Wavelet pathway ported into `mfsusieR/R/`.** Routines we previously
  intended to call as `fsusieR::*` (`remap_data`, `colScale`,
  `gen_wavelet_indx`, `init_prior.default`) are ported into
  `R/utils_wavelet.R` and `R/prior_scale_mixture.R` so mfsusieR can
  replace fsusieR as a runtime dependency. `wavethresh::wd` and
  `wavethresh::wr` remain as math primitives. The DWT is computed once
  at data-class construction and cached.
- **Decoupled post-processing.** `mfsusie()` does NOT take a
  `post_processing` argument. Smoothing and credible bands are
  separate functions on the fit object: `mf_post_smooth(fit, method =
  c("smash", "TI", "HMM"))` and `mf_credible_bands(fit, method =
  "TIWT", level = 0.95)`. The fit stores per-modality residuals by
  default (`save_residuals = TRUE`) so post-processors take only the
  fit; for very-large data the user opts out with `save_residuals =
  FALSE` and passes `(X, Y)` to the post-processors directly. Roxygen
  on each post-processor explains the residual contract in
  statgen-writing-style prose.
- **CS-then-PIP ordering invariant.** PIPs are a function of the
  *final* alpha state after any CS filtering, not before. This is the
  most-likely FDR-miscalibration source identified in
  `inst/notes/paradigms/mvf-original.md` section 8.
- **Per-(scale, modality) residual variance is the v1 default.** The
  legacy `residual_variance_method = "shared_per_modality"` is
  available and is the mode used by Phase 4 contract C3 fidelity
  tests.
- **Forbidden public-API names.** `nullweight`, `gridmult`,
  `max_scale`, `max_SNP_EM`, `cov_lev`, `min_purity`, `filter_cs`,
  `filter.number`, `family` (without scope), and the rest of the
  legacy non-snake_case names from `mvf.susie.alpha::multfsusie` and
  `fsusieR::susiF` SHALL NOT appear in any public mfsusieR signature.
  A naming test asserts this programmatically.
- **EXPLICITLY OUT OF SCOPE.** Sufficient-statistics path
  (`mfsusie_ss`); concrete cross-modality covariance prior
  implementation (only the plug-in seam ships in v1, the mash-style
  default comes in a follow-up); FDR investigation (Phase 5) and the
  fixes that follow it (Phase 6); Rcpp ports (Phase 7); fsusieR's
  `EBmvFR` algorithm (a different model, no SuSiE structure, lives in
  a follow-up if at all).

## Capabilities

### New Capabilities

- `mf-data-class`. The `mf_individual` class, its constructor,
  wavelet-pathway state (DWT cache, scale-index lists, column-scale
  cache), per-modality residual storage, and degenerate-case handling
  (`M = 1`, `T_m = 1`).
- `mf-prior`. Default per-(scale, modality)
  scale-mixture-of-normals prior; pluggable interface for the optional
  cross-modality covariance prior; constructors; fit-time
  empirical-Bayes update via mixsqp.
- `mf-ibss`. IBSS integration via S3 method registration on
  `mf_individual` and `mfsusie`. Per-iteration flow (residual, SER
  stats, posterior moments, KL, fitted update, variance components,
  ELBO). Per-(scale, modality) residual variance update. The three
  apple-to-apple equivalence contracts (C1, C2, C3) are required of
  the IBSS path.
- `mf-public-api`. The `mfsusie()` and `fsusie()` entry functions,
  their argument sets, defaults, return-object structure, and
  standard S3 methods (`coef`, `predict`, `fitted`, `summary`,
  `print`).
- `mf-credible-sets`. PIP and credible-set computation. Strict
  ordering rule that PIPs are computed from the *post-filter* alpha
  state.
- `mf-post-processing`. Decoupled smoothing and credible-band
  functions on the fit object. Residual contract (post-processors
  read `fit$residuals` when `save_residuals = TRUE`; otherwise the
  user passes `X` and `Y`). Apple-to-apple equivalence with
  `fsusieR::susiF(..., post_processing = "...")` followed by the
  matching post-processor.

### Modified Capabilities

None. mfsusieR has no prior `inst/openspec/specs/` entries.

## Impact

- **DESCRIPTION.** Imports: `susieR`, `wavethresh`, `mixsqp`, `ashr`,
  `methods`, `stats`. fsusieR is removed; routines we need are ported
  into `mfsusieR/R/`. Any smashr / ebnm / Rfast dependencies that
  fsusieR pulled in for smash/TI/HMM smoothers are evaluated and
  added only when the smoother in question is ported (PR group 6b).
- **New code under `R/`.** Data class, IBSS S3 methods, public API
  (`mfsusie`, `fsusie`), credible-set helpers, prior constructors,
  ported wavelet helpers (`R/utils_wavelet.R`), ported smoothers
  (`R/post_processing.R`), output formatter, S3 method registrations
  in `R/zzz.R`.
- **New tests under `tests/testthat/`.**
  - `test_susier_degeneracy.R` (contract C1, tolerance 1e-10).
  - `test_fsusier_degeneracy.R` (contract C2, tolerance 1e-8).
  - `test_mvf_alpha_fidelity_*.R` (contract C3, tolerance 1e-8) with
    one file per ported routine.
  - `test_post_processing_*.R` smoothers and credible-bands fidelity
    against `fsusieR::susiF(..., post_processing = "...")`.
  - `test_public_api_naming.R` (forbidden-name check).
  - `test_modularity.R` (no hand-rolled DWT / mixsqp / CS construction
    inside `R/`).
- **Refactor-exceptions ledger.** `inst/notes/refactor-exceptions.md`
  records every line of `mvf.susie.alpha/R/*.R` AND `fsusieR/R/*.R`
  (susiF path only; EBmvFR is out of scope) that is intentionally not
  ported, with reason. Created at the start of Phase 3, populated as
  PR groups land.
- **Roxygen on every function.**
  - `@manuscript_ref methods/<file>.tex eq:<label>` for every formula
    in main code.
  - Original-code references (`mvf.susie.alpha`, `multfsusie`,
    `fsusieR`, `susiF`, "original implementation",
    `@references_original`) FORBIDDEN under `R/`. Provenance lives in
    test-file headers, `inst/notes/refactor-exceptions.md`, and
    free-form prose in dev notes.
- **No change to read-only reference packages.** `susieR`, `mvsusieR`,
  `susieAnn`, `mvf.susie.alpha`, and `fsusieR` are read-only. The
  single authorized exception is the `feature/L_greedy` branch on
  `~/GIT/susieR` for the `L_greedy` argument addition (Migration Plan
  in `design.md`).
