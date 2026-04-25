## Why

mfsusieR is a from-scratch refactor of `mvf.susie.alpha`, a multi-functional
SuSiE implementation with a known FDR miscalibration and a layout that does
not align with the S3 conventions adopted by sister packages
(`mvsusieR/refactor-s3` and `susieAnn`). Phase 1 paradigm notes
(`inst/notes/paradigms/`) catalogued the port source and the two reference
paradigms. This change locks the architecture so Phase 3 (port and
implementation under the Claude-native review loop) and Phase 4 (reference unit tests)
have a stable target. The package's math is anchored to the manuscript at
`/home/gw/GIT/MultifSuSiE_Manuscript`.

## What Changes

- Define `mfsusie()` as a public entry that fits a multi-modality functional
  regression model on top of `susieR::susie_workhorse`, with per-modality
  per-scale prior variance and mixture weights (manuscript
  `methods/derivation.tex` line 47, eq:additive_model).
- Introduce a single `mf_individual` data class (no separate
  multi-functional vs univariate path; univariate trait is a degenerate
  case of `T_m = 1` and single-modality is `M = 1`, mirroring how
  `mvsusieR/refactor-s3` collapses R = 1 vs R > 1 inside one
  `mv_individual` class).
- Introduce a fit-object class `c("mfsusie", "susie")` so existing
  `susieR` accessors that have well-defined semantics for our outputs
  (e.g., `susie_get_pip` shape) can be reused where appropriate.
- Plug into the `susieR` IBSS loop via S3 method registrations on
  `mf_individual` and on the `mfsusie` model class, following the
  `mvsusieR/refactor-s3` paradigm (paradigm reference #1, see
  `inst/notes/paradigms/mvsusieR-s3.md`).
- Compose two prior structures: a per-(scale, modality) wavelet-coefficient
  scale-mixture-of-normals prior (manuscript eq:additive_model with
  `pi_{k,s,m}`, `sigma_{k,s,m}^2`), and an optional cross-modality covariance
  prior (mash-style, off by default) for users who need shared effects across
  modalities. The susieAnn pluggable-predictor pattern (paradigm reference
  #2) is the design template for the cross-modality prior plug-in.
- Port the wavelet pathway from `mvf.susie.alpha`: forward DWT via
  `wavethresh`, scale indexing via `fsusieR::gen_wavelet_indx`, inverse DWT
  for reconstruction, optional TIWT step for credible bands. The
  reorganization respects the manuscript's per-(scale, modality)
  decomposition (manuscript `online_method.tex` lines 218-272).
- **BREAKING vs `mvf.susie.alpha`**: the function is `mfsusie`, not
  `multfsusie`. Argument names follow CLAUDE.md naming rules (snake_case,
  no abbreviations outside the SuSiE lineage). Default residual variance
  is per-(scale, modality), not global per modality — this is the
  manuscript-aligned choice and the most plausible structural fix for the
  FDR miscalibration Anjing observed.
- **EXPLICITLY OUT OF SCOPE**: sufficient-statistics path
  (`mfsusie_ss`); cross-modality covariance prior implementation (only the
  plug-in seam is added in v1, the default `mash`-style implementation
  comes in a follow-up change); FDR-bug investigation (Phase 5) and the
  fixes that follow it (Phase 6); Rcpp ports (Phase 7).

## Capabilities

### New Capabilities

- `mf-data-class`: the `mf_individual` data class, its constructors,
  wavelet pathway state (DWT cache, scale index lists), and degenerate-case
  handling (M = 1, T_m = 1).
- `mf-prior`: prior object hierarchy. Default per-(scale, modality)
  scale-mixture-of-normals prior; pluggable interface for an optional
  cross-modality covariance prior; constructor functions; fit-time
  empirical-Bayes update via mixsqp.
- `mf-ibss`: IBSS integration via S3 method registration. Per-iteration
  flow (residual, SER stats, posterior moments, KL, fitted update,
  variance components, ELBO). Per-(scale, modality) residual variance
  update.
- `mf-public-api`: the `mfsusie()` entry function, its argument set,
  defaults, and return-object structure (class, fields).
- `mf-credible-sets`: PIP and credible-set computation. Strict ordering
  rule that PIPs are computed from the *final* alpha state after any CS
  filtering, not before, fixing the most-likely FDR-miscalibration source
  identified in `inst/notes/paradigms/mvf-original.md` section 8.

### Modified Capabilities

None. mfsusieR has no prior `openspec/specs/` entries.

## Impact

- New code under `R/` (data class, IBSS S3 methods, public API,
  credible-set helpers, prior constructors, output formatter).
- New tests under `tests/testthat/` comparing all numeric outputs against
  `mvf.susie.alpha::multfsusie` per CLAUDE.md Phase 4 rules and the
  test-fidelity memory entry. Tolerances: 1e-6 for posterior summaries,
  1e-10 for deterministic intermediates.
- Roxygen on every function: `@references_original
  mvf.susie.alpha/R/...` for ported logic and `@manuscript_ref
  methods/<file>.tex eq:label` for any formula.
- DESCRIPTION imports: `susieR`, `fsusieR` (delegated wavelet helpers),
  `wavethresh`, `mixsqp`, `ashr` (for the scale-mixture prior, matching
  the original code's choice).
- No change to the four read-only reference packages (`susieR`,
  `mvsusieR`, `susieAnn`, `mvf.susie.alpha`). They remain reference only.
