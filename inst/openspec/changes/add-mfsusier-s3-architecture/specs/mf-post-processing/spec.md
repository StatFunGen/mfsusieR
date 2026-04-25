## ADDED Requirements

### Requirement: post-processing functions are decoupled from mfsusie()

`mfsusie()` and `fsusie()` SHALL NOT accept a `post_processing`
argument. Smoothing of effect curves and credible-band computation
SHALL be exposed as separate functions on the fit object:

- `mf_post_smooth(fit, method = c("smash", "TI", "HMM"), ..., X = NULL, Y = NULL)`
- `mf_credible_bands(fit, method = "TIWT", level = 0.95, ..., X = NULL, Y = NULL)`

The decoupling rationale (mechanism-first, statgen-writing-style
prose) is documented in roxygen on each function and in
`design.md` D8c: separating fitting (a variational estimate of
effect locations and posterior moments) from summarizing (smoothed
point estimates and credible bands) lets users add or change
smoothers without rerunning the IBSS loop and matches the
susieR/mvsusieR convention where `coef`, `predict`, `summary` are
post-fit accessors.

#### Scenario: mfsusie() does not accept post_processing

- **WHEN** `formalArgs(mfsusie)` and `formalArgs(fsusie)` are
  inspected
- **THEN** neither SHALL contain `post_processing`

#### Scenario: post-processors are exported

- **WHEN** the package is loaded
- **THEN** `mf_post_smooth` and `mf_credible_bands` SHALL be
  exported and discoverable via `?mf_post_smooth`, `?mf_credible_bands`

### Requirement: residual contract

Each post-processor SHALL read per-modality residuals from
`fit$residuals` when present (the `save_residuals = TRUE` default
case, per `mf-data-class/spec.md`). When `fit$residuals` is NULL
(opt-out path), the post-processor SHALL require the user to pass
`X` and `Y` explicitly and SHALL recompute the residuals on the
fly.

#### Scenario: default path reads from fit$residuals

- **WHEN** `mfsusie(X, Y, pos)` is called with default
  `save_residuals = TRUE`, then `mf_post_smooth(fit, method = m)`
  is invoked without passing X or Y
- **THEN** the post-processor SHALL succeed, read residuals from
  `fit$residuals`, and return the augmented fit; if X or Y are
  also passed and disagree with `fit`'s identity, the
  post-processor SHALL warn

#### Scenario: opt-out path requires X and Y

- **WHEN** `mfsusie(X, Y, pos, save_residuals = FALSE)` is called,
  then `mf_post_smooth(fit, method = m)` is invoked without passing
  X or Y
- **THEN** the post-processor SHALL stop with an informative error
  message naming the missing arguments

#### Scenario: opt-out path with X and Y produces identical output

- **WHEN** the same fit is processed two ways: once with
  `save_residuals = TRUE` and a `(fit, method)` call, once with
  `save_residuals = FALSE` and a `(fit, method, X, Y)` call
- **THEN** the two augmented fits SHALL be element-wise identical
  at tolerance `1e-12` (the only difference between the two paths
  is where the residuals come from)

### Requirement: smash / TI / HMM dispatch

`mf_post_smooth(fit, method)` SHALL dispatch on `method` to one of
three smoother implementations ported from
`fsusieR/R/operation_on_susiF_obj.R` (`smash_regression`,
`TI_regression`, `HMM_regression`) into `R/post_processing.R`.

#### Scenario: each method is callable

- **WHEN** `mf_post_smooth(fit, method = m)` is invoked for each
  `m in c("smash", "TI", "HMM")` on the standard fixture
- **THEN** each call SHALL return an augmented fit object of class
  `c("mfsusie", "susie")` with the smoother's output stored in a
  documented field (e.g., `fit$fitted_func` or the per-modality
  smoothed-curve cache); the call SHALL NOT mutate the original
  fit

#### Scenario: method dispatch is the only signature variant

- **WHEN** `args(mf_post_smooth)` is inspected
- **THEN** the public signature SHALL be a single function with
  `method` as the dispatch argument, not three separate functions
  (`mf_smash`, `mf_TI`, `mf_HMM`); per `design.md` D8c, single
  dispatch was chosen for namespace cleanliness and to make adding
  a fourth smoother a non-breaking change

### Requirement: TIWT credible bands

`mf_credible_bands(fit, method = "TIWT", level = 0.95)` SHALL
compute per-effect translation-invariant wavelet-transform credible
bands per the manuscript's `online_method.tex` step 6, ported from
`fsusieR/R/operation_on_susiF_obj.R` into `R/post_processing.R`.

#### Scenario: TIWT bands attach per-effect

- **WHEN** `mf_credible_bands(fit, method = "TIWT", level = 0.95)`
  is called on a fit with L effects and M modalities
- **THEN** the returned fit SHALL contain a field `cred_band` that
  is a list of length L, each entry a list of length M, each inner
  entry a `2 x T_m` matrix of upper / lower band per position

### Requirement: contract C2-with-post-processing equivalence with fsusieR::susiF

For each `method`, the apple-to-apple contract MUST hold:

```
mfsusie(X, Y, pos) |> mf_post_smooth(method = m)
                                 ≡   to <= 1e-8
fsusieR::susiF(Y, X, pos, post_processing = m)
```

per `design.md` D11d. The same contract SHALL hold for
`mf_credible_bands(fit, method = "TIWT")` against the band output
that `fsusieR::susiF` produces internally for the TI path.

#### Scenario: smash equivalence

- **WHEN** `mfsusie(X, Y, pos) |> mf_post_smooth(method = "smash")`
  is run on a single-modality (`M = 1, T_1 > 1`) fixture and seed,
  and `fsusieR::susiF(Y, X, pos, post_processing = "smash")` is
  run on the same fixture and seed with arguments aligned per the
  test-file header
- **THEN** the two augmented fits SHALL agree element-wise on
  every numeric output (`alpha`, `mu`, `mu2`, `lbf`, `lbf_variable`,
  `KL`, `sigma2`, `elbo`, `niter`, `pip`, plus the smoother-specific
  output `fitted_func`) at tolerance `<= 1e-8`

#### Scenario: TI equivalence

- **WHEN** the same comparison is run with `method = "TI"`
- **THEN** the same element-wise tolerance `<= 1e-8` SHALL hold,
  including on the TIWT credible bands when those are part of
  fsusieR's TI output

#### Scenario: HMM equivalence

- **WHEN** the same comparison is run with `method = "HMM"`
- **THEN** the same element-wise tolerance `<= 1e-8` SHALL hold,
  except that fsusieR's HMM path sets `cred_band <- NULL`; the
  mfsusieR side SHALL match this

#### Scenario: intentional fsusieR-bug fixes are asserted explicitly

- **WHEN** mfsusieR fixes a smoother bug present in fsusieR (e.g.,
  a missing edge-case in the HMM forward pass, found during the
  D13 audit of the smoother helpers)
- **THEN** the C2-with-post-processing test SHALL include an
  explicit `expect_*` assertion that records the deviation, with
  a comment citing the OpenSpec change that authorized the fix;
  tolerance SHALL NOT be loosened to mask the deviation

### Requirement: post-processor port-quality audit

The smoother and credible-band routines SHALL be ported from
`fsusieR/R/operation_on_susiF_obj.R` and supporting helpers
(`fsusieR::TI_regression`, `fsusieR::smash_regression`,
`fsusieR::HMM_regression`) under the port-quality audit in
`design.md` D13. Style improvements logged in the commit message;
behaviour-preserving only.

#### Scenario: audit recorded in commit message

- **WHEN** the PR landing the post-processors is reviewed
- **THEN** the commit message SHALL state that the D13 audit was
  run (simplify skill or Explore subagent) and SHALL list any
  intentionally omitted lines with entries added to
  `inst/notes/refactor-exceptions.md`

#### Scenario: numerical drift caught by C2 test

- **WHEN** the audit edit unintentionally changes numerics
- **THEN** the contract C2-with-post-processing test SHALL fail at
  tolerance `<= 1e-8`, the audit edit SHALL be reverted, and the
  PR SHALL re-run the audit with the failing change identified
