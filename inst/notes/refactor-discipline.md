# Refactor discipline for the mfsusieR port

Binding policy for all Phase 3+ implementation work that ports
routines from `mvf.susie.alpha` or `fsusieR/susiF`. The companion
file `inst/notes/test-fidelity-policy.md` covers the unit-test
contract that this discipline relies on as a safety net.

This file is the canonical source for the refactor discipline.
Copies in the per-Claude memory system
(`~/.claude/projects/.../memory/`) are short pointers at this file
so the policy is reproducible across sessions, machines, and
contributors.

## 0. Upstream-first principle (BINDING for all port work)

Before porting any routine from `mvf.susie.alpha` or `fsusieR`
into mfsusieR, FIRST check whether `susieR` already has a
similar helper or whether one can be added or generalised in
susieR for shared use. Code duplication from susieR is forbidden;
when a routine could live upstream and benefit other downstream
packages (mvsusieR, mfsusieR, future packages), it SHALL be added
to susieR rather than ported into each downstream individually.

**How to apply:**

- During the D13 author-pass audit, before opening a new file
  under `R/`, search susieR for related helpers (`grep`,
  `Glob`, search the relevant susieR source files). Look for
  functions that compute the same quantity, even if the call
  shape differs.
- If susieR has a close helper:
  - Use it directly when the call shape matches.
  - Generalise upstream (with author authorisation) if the
    shape needs adjustment (matrix-Y, additional kwargs,
    Rfast acceleration, etc.). Open a separate susieR commit;
    refactor mvsusieR to use the helper too where applicable.
    Document the upstream change in the mfsusieR change's
    Migration Plan section.
- If no susieR helper exists and the routine is shareable across
  packages, propose adding it to susieR. Examples landed under
  this principle:
  - `susieR::susie_workhorse(L_greedy, lbf_min)` (greedy outer
    loop, used by susieR / mvsusieR / mfsusieR; landed
    2026-04-25).
  - `susieR::compute_marginal_bhat_shat(X, Y, ...)` (per-position
    marginal OLS regression, used by susieR / mvsusieR /
    mfsusieR; landed 2026-04-25).
- Only port to mfsusieR locally when the routine is genuinely
  mfsusieR-specific (e.g., wavelet-domain wrappers, the
  per-modality prior plumbing). Document the local-port decision
  in the refactor-exceptions ledger.

**Why:** every port-source routine (cal_Bhat_Shat, init_prior,
greedy logic, etc.) that gets duplicated into mfsusieR creates
maintenance burden AND drifts from the canonical susieR /
mvsusieR implementation. The upstream-first principle keeps the
SuSiE family consistent and lets all downstream packages benefit
from improvements made in any one of them.

## 1. Mathematical fidelity first

Every line of the original is accounted for. Lines intentionally
omitted are logged in `inst/notes/refactor-exceptions.md` per
design.md D11e (entry schema: file/range, behavior, decision,
reason). When porting, translate the math faithfully; do NOT
"improve" the algorithm during the port. The OpenSpec change that
covers algorithmic changes ships separately.

## 2. Verify math independently of the original code

The manuscript at `~/GIT/MultifSuSiE_Manuscript` is the ground
truth. The port sources (`mvf.susie.alpha`, `fsusieR/susiF`) are
known to have bugs (the FDR miscalibration, the
PIP-after-CS-filter ordering, the `gen_wavelet_indx` `lev_res = 1`
edge case, the `col_scale` `d`-attribute redundant formula, the
`interpol_mat`/`remap_data` redundant min/max chain). When
reference and refactored disagree, be critical of BOTH and read
the manuscript before deciding which side is correct.

## 3. Binary apple/orange test philosophy

Apple-to-apple comparisons (same algorithm, same code path) use
`tolerance <= 1e-8`. Apple-to-orange comparisons (genuinely
different algorithms) are smoke tests only, no numerical
comparison. Tolerances of `1e-6`, `1e-2`, `5e-2` for
apple-to-apple comparisons are FORBIDDEN: they hide bugs.

A failing apple-to-apple test at `1e-8` is a bug to investigate,
not a tolerance to relax.

### Pattern B exception (NOT a violation of the binary rule)

When a port-source bug is fixed in mfsusieR (see section 4 below),
the unit test for the cleaned-up quantity uses
`tolerance = 1e-12` (or similarly tight). This is four orders of
magnitude tighter than the `1e-8` contract floor; it accepts the
port source's machine-epsilon rounding noise from the buggy
formula we corrected, while still catching real drift. The
forbidden graduated tolerances (`1e-6`, `1e-2`, `5e-2`) are
forbidden because they mask drift; `1e-12` does not.

## 4. mfsusieR ships clean code: port-source bugs are FIXED

Wasteful, buggy-in-spirit, or numerically noisy patterns from
`fsusieR` / `mvf.susie.alpha` SHALL be cleaned up during the port
rather than carried over for cosmetic bit-identity. Two canonical
patterns for handling such cleanups:

**Pattern A: behavioral deviation.** When the original code is
functionally wrong (e.g., the PIP-after-CS-filter ordering in
fsusieR), mfsusieR fixes the behaviour. The unit test DIRECTLY
ASSERTS the corrected behaviour (does NOT "expect" the old buggy
output) and includes a comment citing the OpenSpec change that
authorized the fix.

**Pattern B: ULP-level cleanup of a buggy formula.** When the
original code computes the right quantity via a wasteful or
rounding-noisy expression (e.g., fsusieR's `d` attribute via
`(n*cm^2 + (n-1)*csd^2 - n*cm^2)/csd^2` instead of the direct
`(n-1)`), mfsusieR computes the simplified form directly. The
unit test KEEPS asserting equivalence with the port source but
uses a relaxed tolerance (e.g., `1e-12`) well below the contract
floor, with an inline comment saying "this is an fsusieR-bug fix,
not numerical drift."

Every Pattern-A or Pattern-B cleanup gets a refactor-exceptions
ledger entry per design.md D11e with decision
`replaced-by-simplified-math`, `replaced-by-cleaner-message-text`,
`replaced-by-seq_len`, `replaced-by-single-step`,
`replaced-by-parameter-unification`, etc.

## 5. Refactor cleanups during a port land WITH unit tests in the same edit pass

When the D13 port-quality audit substep flags a cleanup beyond
cosmetic style (function-extraction, signature change, removal of
redundant intermediate steps, parameter unification across
helpers), the unit test that locks the new behaviour SHALL land
in the same author-pass edit, not as a follow-up commit.

For each cleanup: write a regression-net test that would fail if
someone reintroduces the redundancy. Example:

- `interpol_mat` was simplified to fold position-rescaling into
  one step; the regression test asserts
  `remap_data(...)$outing_grid == interpol_mat(...)$outing_grid`
  exactly, so reintroducing a separate rescale step in
  `remap_data` breaks the test.

The C2 / C3 / C1 contracts at the IBSS / fit level are the
ultimate safety net. Per-routine tests catch drift earlier.

## 6. No quick fixes

When tests fail, investigate root cause. Do NOT paper over with
`as.character()` casts, swallowing `tryCatch`, defensive `NA`
coercions, or graduated-tolerance loosening. If a fix feels
hacky, it is hacky. Pattern B's `1e-12` tolerance is not a quick
fix: it is a documented response to a port-source bug, paired
with a refactor-exceptions ledger entry and an inline comment in
the test.

## 7. No reference to original code in main package roxygen

`@references_original` tags, `mvf.susie.alpha` / `fsusieR` /
`susiF` / `multfsusie` mentions, and "original implementation"
comments are FORBIDDEN under `R/`. Provenance lives in:
- `inst/notes/refactor-exceptions.md` (line-level omission ledger)
- Test file headers (apple-to-apple comparison metadata)
- Free-form prose in session notes and paradigm notes

Manuscript citations are required in main roxygen and stay. The
canonical form is the built-in `@references` tag with a
`Manuscript:` prefix:
```r
#' @references
#' Manuscript: methods/<file>.tex eq:<label>
```
The built-in tag is used because roxygen2 rejects unknown tags;
the `Manuscript:` prefix makes the citation greppable for the
reviewer-checklist scan.

## 8. Structural cleanliness AFTER correctness

Group runtime state in `model$runtime`, distinguish permanent
output fields from runtime, remove dead code only after
correctness is established. Refactor-then-clean, not
clean-while-refactoring.

## 9. Optimization is a separate phase

Phase 7 only. Keep the pure-R reference next to any Rcpp version.
Reference tests live in `tests/testthat/test_rcpp_vs_r.R`.

## 10. Show plan before large changes

Don't surprise Gao with a big refactor. When told to stop, stop.

## 11. Review style: Gao reviews summaries, not files

After completing a unit of work (a refactor pass, an
audit-then-fix cycle, a spec rewrite, etc.) Claude SHALL produce
a self-contained summary in the chat covering: what was changed,
what was decided, what was skipped (with reasons), and what
tests / validation steps ran. Gao reviews the summary and replies
with inline quotes from it plus comments per quote. Each
quote-and-comment is a discrete piece of feedback that Claude
SHALL address in the next response (fix the code, justify the
skip, or escalate the question).

Summaries should answer: what changed, what decisions were made,
what was deferred (and why), what verification was run. For
audits: list each finding with severity, decision (fix / skip /
defer), and one-line justification per skip. For refactors: state
what was preserved, what was simplified, what numerical
contracts were re-verified. Be concise but complete; the summary
is the artifact Gao reads.

## Concrete examples (from PR group 2, 2026-04-25)

| Cleanup | Pattern | Test tolerance | Ledger entry |
|---|---|---|---|
| `col_scale` `d` attribute simplified to `rep(n - 1, ncol(x))` | B | `1e-12` on `attr(., "d")` | `replaced-by-simplified-math` |
| `gen_wavelet_indx` `lev_res = 1` edge fixed via `seq_len` | A (edge case below contract floor) | bit-identical for `lev_res >= 2` | `replaced-by-seq_len` |
| `interpol_mat` / `remap_data` redundant `min/max` chain folded into one step | B | `1e-12` on `outing_grid` | `replaced-by-single-step` |
| `dwt_matrix`'s `min.scale = 10` lifted to `max_scale` parameter | (cosmetic, not a bug) | bit-identical at default | `replaced-by-parameter-unification` |
| fsusieR's verbose-message multi-line text replaced with semicolon | (cosmetic, not numerical) | n/a | `replaced-by-cleaner-message-text` |
