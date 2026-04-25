# Refactor exceptions ledger

This file records every line of the two port sources that mfsusieR
intentionally does NOT carry over. Doctrine lives in
`inst/openspec/changes/add-mfsusier-s3-architecture/design.md`
D11e. Reviewer pass on every Phase 3 numerical PR confirms the
ledger is up to date; PRs that omit lines without a corresponding
entry SHALL be blocked.

## Port sources in scope

- `mvf.susie.alpha/R/multfsusie.R` and supporting
  `operation_on_multfsusie_*.R`, `EM.R`, `ELBO_mutlfsusie.R`,
  `computational_routine.R`, `utils_wavelet_transform.R`, and the
  parts of `utils.R` and `utils_formatting.R` that the IBSS path
  touches.
- `fsusieR/R/susiF.R` and `susiF_workhorse.R` and supporting
  `operation_on_susiF_obj.R`, `wavelet_utils.R`, the parts of
  `operation_on_prior.R` that init the scale-mixture prior, the
  per-method smoothers in `operation_on_susiF_obj.R`
  (`TI_regression`, `smash_regression`, `HMM_regression`), and
  the parts of `utils.R` needed by the susiF path.

## Out of scope (file-level entries; no per-line walk required)

- fsusieR/R/EBmvFR.R
  Behavior: EB multivariate functional regression, a different
    model with no SuSiE structure.
  Decision: out-of-scope-EBmvFR
  Reason: design.md D5 and CLAUDE.md hard rule #2 record EBmvFR as
    an entirely separate algorithm. Not in v1; possibly addressed
    in a follow-up OpenSpec change after the SuSiE-track ports
    stabilize.
- fsusieR/R/EBmvFR_workhorse.R
  Behavior: Inner loop for EBmvFR.
  Decision: out-of-scope-EBmvFR
  Reason: As above.
- fsusieR/R/operation_on_EBmvFR_obj.R
  Behavior: Operations on the EBmvFR fit object.
  Decision: out-of-scope-EBmvFR
  Reason: As above.

## Entry schema for in-scope omissions

```
- <port_source>/R/<file>.R:<lo>-<hi>
  Behavior: <one-line description of what the original lines do>
  Decision: omit | replaced-by-<upstream> | deferred-to-<phase> |
            out-of-scope-EBmvFR
  Reason: <one-paragraph justification, citing OpenSpec change or
          manuscript section>
```

## In-scope omissions

(Empty until PR groups 2-7 land. PRs SHALL append entries here as
they go, per the doctrine in design.md D11e and the reviewer
checklist in `inst/notes/review-loop-methodology.md`.)
