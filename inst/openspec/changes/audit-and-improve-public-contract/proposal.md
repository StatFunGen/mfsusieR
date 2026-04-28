# Audit and improve mfsusieR's public contract

## Why

The core mfsusieR algorithm is in place and the perf work has
landed. Before broader use, the package needs a cleanup pass:
the fit object has accumulated fields, several S3 methods may
duplicate logic that `susieR`'s default already provides,
verbose output is sparse, the HMM credible band is suppressed
on a hunch that it is too wide, and a representative scan
case (n=84, p≈3500, M=6, T=128) reportedly times out at 5
minutes. The recent perf work (cache + subsetting + cpp11)
likely changed the timeout picture; we have not re-measured.

This change captures the audit-and-improve work as one
coordinated effort. Implementation will split into focused
commits per section.

## Scope (seven sections)

### Section 1: Trim and document the fit object

Walk every field on a typical `mfsusie` fit. For each field,
decide:
- Keep as is (with one-line description in the return-value
  roxygen).
- Remove (redundant with `susieR`'s base shape, derivable from
  other fields, or never read by any public method or vignette).
- Move to `summary(fit)` (informational, not a fit-state field).
- Rename (current name diverges from `susieR`'s convention).

Result: a smaller, susieR-shape-compatible fit object with a
documented field set. Fields that are mfsusieR-specific are
clearly flagged as extensions.

### Section 2: Feature gap with `fsusieR` and `mvf.susie.alpha`

Read every exported function and parameter of `fsusieR` and
`mvf.susie.alpha`. List unported features in a session note.
For each, classify as `would-port-if-asked`, `out-of-scope`,
or `port-now`. Pure investigation; no code changes here.

### Section 3: Maximize susieR backbone usage

Find S3 methods mfsusieR overrides. For each, decide whether
the body meaningfully diverges from `susieR`'s default or is
trivially redundant. Drop trivial overrides; let `susieR`'s
default fire. Document the survivors with a one-line reason
each. Flag anything found in `susieR` that looks redundant or
miscalibrated for upstream discussion (no edits to `susieR`
without explicit approval).

### Section 4: Useful verbose output and convergence diagnostics

Currently `verbose = TRUE` prints essentially nothing
per-iteration. `susieR` prints ELBO and max alpha change.
Match that, optionally adding peak memory if cheap. Offer
`convergence_metric = c("pip_diff", "elbo")` and add a unit
test that the two stopping criteria agree on standard
fixtures within an iteration or two of each other.

### Section 5: HMM credible band

Re-derive the per-position credible band for HMM smoothing
from the manuscript's posterior-mixture formula. Compare
against `fsusieR::fit_hmm`'s band shape. If our band is
genuinely wider than the manuscript implies, it is a bug;
fix the formula. If both packages produce a wide band, the
manuscript band is correct as derived but the visualization
is misleading — in that case, document the band semantics
clearly and stop suppressing it.

### Section 6: Performance and convergence on the heavy fixture

Re-measure `n=84, p≈3500, M=6, T=128` after the recent perf
landing. If it still times out, profile and identify whether
the bottleneck is per-effect work that scales in `p` or
something else (initial wavelet basis cost? sigma2 update?
fitted-value reconstruction?). Apply targeted fixes (Rfast
where it cleanly wins, more cpp11 only if profile says so,
loop-invariant lifting). Also investigate the existing
50-iteration non-convergence warning on the test fixture:
identify whether it is a fixture issue (too tight tolerance,
too low PVE) or an algorithmic issue (init quality, Eloglik
disagreement with PIP-diff).

### Section 7: Comment polish and helper promotion

Two clean-up jobs:
- Audit and remove British-English spellings; align on
  American English everywhere (R, comments, docs, vignettes).
- Audit comments and roxygen for any reference to upstream
  packages or "port-fidelity" language. Comments must be
  self-contained: a reader who has never heard of fsusieR or
  mvf.susie.alpha should fully understand them. References
  to other packages live in `inst/notes/` or test file
  headers, not in `R/`.
- Audit internal helpers prefixed with `.`; promote (drop
  the `.`) those whose generality justifies a stable internal
  name. Pattern set by `univariate_smash_regression`. Do not
  export; just rename internally so future contributors can
  find them.

## Impact

- Fit object: smaller and more documented. Existing user code
  reading well-known fields (`pip`, `sets`, `alpha`, `mu`,
  `coef`) is unaffected. Code reading mfsusieR-specific
  internals may break; documented at the cleanup.
- S3 surface: smaller. No behavior change.
- Verbose: existing users see new per-iteration output when
  `verbose = TRUE`. Default unchanged.
- HMM: band displayed when populated; visualization may show
  uncertainty more honestly than before.
- Perf: targeted only at fixtures that show real headroom; no
  regressions tolerated on existing fixtures.
- Comments / naming: pure hygiene, no behavior change.

## Out of scope

- New algorithmic features (mfsusieR's algorithm is locked).
- Renaming user-facing arguments of `mfsusie()` or `fsusie()`
  (Phase 8 polish, separate change).
- Cross-package work in `susieR` (we list candidates; the
  upstream PRs are separate).
- Any C++ port beyond what section 6 measurably justifies.
