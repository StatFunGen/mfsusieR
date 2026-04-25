# Claude-native review loop for Phase 3+ implementation

This is the implementation methodology for Phases 3, 4, 6, and 7. It runs
entirely inside Claude Code with no external review CLI.

## Origin

The structure is borrowed from external write-review-revise loops (e.g., the
Humanize plugin's RLCR design): an authoring agent writes code, a separate
reviewing agent looks at it with fresh context, findings come back, the
author addresses, repeat. The novelty here is that the reviewing agent is a
Claude Agent (Explore subagent) running in the same Claude Code session, not
an external CLI.

We do not depend on any external review tool. Earlier drafts of the project
brief referenced a `Humanize` plugin and a `codex` CLI; both have been
removed. References to either in older notes or commit messages are
historical.

## What problem the loop solves

A single authoring pass writes code, sees its own assumptions baked in, and
misses the bugs that those assumptions hide. A fresh-context reviewer reading
the same code without the authoring trace catches those bugs. The two
agents have to be different *contexts*, not different models, for this to
work; running the same Claude with reset context as the reviewer is
sufficient.

## The loop

For each task in `inst/openspec/changes/<name>/tasks.md`:

1. **Author pass.** Claude (the main agent in the session) writes code per
   the spec.md contract for the relevant capability. Tests land in the same
   commit. Author writes:
   - the code,
   - the tests,
   - a short PR description in the commit message,
   - any new manuscript or `mvf.susie.alpha` / `fsusieR::susiF` cross-reference
     roxygen tags.

   **Port-quality audit substep.** When the author pass copies a routine
   from a port source (`mvf.susie.alpha/R/*.R` or `fsusieR/R/*.R`), the
   author SHALL run a code-quality audit on the copied routine before
   advancing to the reviewer pass. Two acceptable mechanisms:
   (a) invoke the `simplify` skill on the changed files, or
   (b) spawn an `Agent` (Explore subagent) with a focused "audit ported
   routine for naming, dead code, redundant branches, and idiomatic R"
   prompt and address its findings.
   Style improvements logged in the commit message; behaviour-preserving
   only (numerical changes are out of scope and are caught by the
   apple-to-apple equivalence test). The fidelity test against the port
   source is the safety net: any post-audit numerical drift fails the
   test and reverts the audit edit.

2. **Reviewer pass.** Spawn an `Agent` (Explore subagent) with a focused
   prompt that contains only:
   - the task ID and one-line description,
   - the relevant capability's `spec.md`,
   - the file paths the author touched,
   - explicit checklist items (below).
   The reviewer reads the diff with no awareness of the author's reasoning
   and reports findings as structured output (file:line + issue).

3. **Address pass.** The author reads the reviewer's findings. For each:
   - fix the code if the finding is correct,
   - reply in the commit message or in a follow-up note if the finding is
     wrong (and explain why),
   - escalate to Gao if the finding identifies a contract issue (spec
     wrong, manuscript wrong, fixture wrong).

4. **Re-review or land.** If the address pass changed substantive code, run
   step 2 again with a fresh reviewer agent. Two consecutive clean reviewer
   passes mean the task can be marked done in tasks.md and the PR can be
   committed.

A "clean review" is one where every finding is either fixed or has a written
justification accepted by Gao.

## Reviewer checklist

The Agent reviewer prompt embeds this checklist. Every numerical PR walks
through it. The checklist absorbs principles from the code-refactor skill
at `~/Documents/obsidian/AI/general/agents/AGENT-code-refactor.md`.

1. **Spec adherence.** Does the code satisfy every requirement in the
   relevant `spec.md`? Cite the specific requirement if not.
2. **Modularity audit (design.md D10).** Does any function in the diff
   duplicate logic already provided by `susieR`, `fsusieR`, `mixsqp`,
   `ashr`, or `wavethresh`? Any reimplemented IBSS loop, DWT, or CS
   construction is a finding.
3. **Manuscript cross-reference (in main code).** Does every formula in
   the diff carry a manuscript citation in roxygen? The canonical
   form is `@references` followed by `Manuscript: methods/<file>.tex
   eq:<label>` (built-in roxygen tag, greppable via the
   `Manuscript:` prefix). Open
   `~/GIT/MultifSuSiE_Manuscript/methods/<file>.tex` and verify the
   cited equation says what the code computes. Verify the math
   independently of the original code; the original may have bugs.
4. **Original-code references are NOT in main code.** Per design.md
   D12, scan the `R/` portion of the diff for `@references_original`,
   `mvf.susie.alpha`, `multfsusie`, `fsusieR`, `susiF`,
   `original implementation`. Zero matches expected (except inside
   user-facing error/warning string literals). For test-file portions
   of the diff, verify the test file has a header comment block listing
   the `mvf.susie.alpha/R/<file>.R#L<lo>-L<hi>` or
   `fsusieR/R/<file>.R#L<lo>-L<hi>` ranges being compared.
5. **Refactor-exceptions completeness.** If the diff omits or replaces
   any line of `mvf.susie.alpha` or `fsusieR/susiF`, is there a
   corresponding entry in `inst/notes/refactor-exceptions.md`? Walk the
   cited original line range; for any line not implemented in the diff,
   the omission must be logged.
6. **Naming rules.** Snake_case throughout (functions, arguments, AND
   filenames in `R/` / `tests/testthat/`; hyphenated source filenames
   are forbidden). No abbreviations beyond `pip`, `cs`, `lbf`, `elbo`,
   `kl`. None of the names on the `mf-public-api/spec.md` forbidden list.
7. **CS-then-PIP ordering.** If the diff touches finalize logic, verify
   the five-step order in `mf-credible-sets/spec.md` is preserved.
8. **Binary tolerance classification.** For each test in the diff,
   classify as apple-to-apple (same algorithm, same code path) or
   apple-to-orange (genuinely different algorithms). Apple-to-apple
   tests SHALL use `tolerance <= 1e-8`. Apple-to-orange tests SHALL
   be smoke tests only (return shape, no NaN, ELBO monotone) with no
   numerical comparison. Tolerances of `1e-6`, `1e-2`, `5e-2` are
   forbidden.
9. **Test fidelity.** Per the test-fidelity memory entry, do the
   apple-to-apple tests compare *every* numeric output (alpha, mu,
   mu2, lbf, lbf_variable, KL, sigma2, elbo, niter, pip, cs
   membership) against `mvf.susie.alpha`, not just `length(cs)` and
   `pip`?
10. **Degenerate-case test.** If the diff touches IBSS or finalize,
    does the susieR-degeneracy test (`mf-ibss/spec.md` D11) still pass
    at `tolerance <= 1e-10`?
11. **No quick fixes.** Are there any `as.character()`,
    `tryCatch({...}, error = function(e) NA)`, hardcoded fallbacks,
    or "this works for now" comments that paper over a bug instead of
    investigating it? Each such pattern is a blocker until the root
    cause is identified.
12. **Roxygen cleanliness.** Does `roxygen2::roxygenise()` run without
    errors on the diff? Are documented defaults consistent with the
    function signatures? Are the manuscript citations rendered?

A finding can be a *blocker* (must fix) or a *note* (should fix in a
follow-up). The default is blocker; the reviewer downgrades to note only
when explicitly justified. When a reviewer and author disagree on
whether a finding is correct, the user (Gao) decides.

## Reviewer prompt template

```
You are a code reviewer. Read the diff at <commit-hash>.

Capability spec to verify against:
<paste spec.md contents>

Files touched:
<list>

Walk through the 10-item checklist in
inst/notes/review-loop-methodology.md exactly. For each item:
- "PASS" with one-line justification, OR
- "FINDING (blocker|note): file:line  one-line description  one-line
  recommended fix"

Do not paraphrase the spec. Do not propose new requirements. Do not
write code. Output structured findings only.
```

The author then reads the findings and decides per item.

## When NOT to use the loop

- **Pure boilerplate** (DESCRIPTION, NAMESPACE skeleton, empty
  `zzz.R::.onLoad`, test setup files). The reviewer's checklist has
  nothing to verify, and a fresh-context reviewer adds nothing. PR group
  1 in `tasks.md` is the canonical example: skip the reviewer pass.
- **Documentation-only changes** (vignettes, roxygen polish, README
  edits). Author + Gao review is sufficient.
- **Pure rename PRs** that preserve behavior. Tests catch behavior
  preservation; the loop adds no value.

For *every* PR that touches numerical routines or CS / PIP logic, the
loop runs. No exceptions.

## Reviewer model selection

The Agent tool's `model` parameter selects the model used by the
reviewer subagent. Default inheritance gives the reviewer the same
model the main session is using. For numerical PRs (Phase 3
implementation work, Phase 6 fixes, Phase 7 Rcpp) the reviewer
SHOULD be invoked with `model: "opus"` to maximize correctness
review quality on the binary tolerance / port-fidelity contracts.
For lighter passes (the parallel /simplify code-quality / efficiency
agents, documentation-only PRs) `sonnet` is sufficient. The
author-pass agent stays on the session model; only the reviewer
pass dials up.

If the user runs into rate limits on opus, fall back to sonnet for
the reviewer pass and note the substitution in the commit message.
Two consecutive clean reviewer passes are still required regardless
of model.

## Cost and runtime

The reviewer pass costs roughly one Agent invocation per task. Most tasks
in `tasks.md` need one or two cycles before clean. Phase 3's 51 tasks (12
PR groups) are likely 50-100 reviewer-passes worth of work over multiple
sessions.

There is no hard token budget on the loop. If a task is taking more than
three reviewer cycles to clear, that is a signal to escalate to Gao
(the spec is probably wrong or the requirement is ambiguous), not to
brute-force a fourth cycle.

## Coordination with OpenSpec

The loop and OpenSpec are orthogonal. OpenSpec defines *what* (the spec
contract); the loop is *how* (the authoring + review process). When the
loop finds a contract issue (spec.md is wrong, ambiguous, or contradicts
the manuscript), the right response is to open or update an OpenSpec
proposal, not to patch the spec inside the implementation PR.

## What changes for Gao

You review the same artifacts you would have reviewed under any
implementation methodology: the diff, the tests, the commit message. The
only added artifact is the reviewer's findings list, which is in the
session's chat log (or attached to the PR as a comment block) so you
can see the review trail without reading every line of code.
