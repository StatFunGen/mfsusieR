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
   - any new manuscript or `mvf.susie.alpha` cross-reference roxygen tags.

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
through it.

1. **Spec adherence.** Does the code satisfy every requirement in the
   relevant `spec.md`? Cite the specific requirement if not.
2. **Modularity audit (design.md D10).** Does any function in the diff
   duplicate logic already provided by `susieR`, `fsusieR`, `mixsqp`,
   `ashr`, or `wavethresh`? Any reimplemented IBSS loop, DWT, or CS
   construction is a finding.
3. **Manuscript cross-reference.** Does every formula in the diff carry an
   `@manuscript_ref methods/<file>.tex eq:<label>` tag? Open
   `~/GIT/MultifSuSiE_Manuscript/methods/<file>.tex` and verify the cited
   equation says what the code computes.
4. **Original-code cross-reference.** Does every routine ported from
   `mvf.susie.alpha` carry an
   `@references_original mvf.susie.alpha/R/<file>.R#L<lo>-L<hi>` tag? Open
   the cited lines and verify the port preserves the numerical behavior
   (or, if not, that the deviation is documented and authorized by an
   OpenSpec change).
5. **Naming rules.** Snake_case throughout. No abbreviations beyond `pip`,
   `cs`, `lbf`, `elbo`, `kl`. None of the names on the
   `mf-public-api/spec.md` forbidden list. Function-level cross-check.
6. **CS-then-PIP ordering.** If the diff touches finalize logic, verify
   the five-step order in `mf-credible-sets/spec.md` D9 is preserved.
7. **Test fidelity.** Per the test-fidelity memory entry, do the tests
   compare *every* numeric output (alpha, mu, mu2, lbf, KL, sigma2, elbo,
   pip, cs membership) against `mvf.susie.alpha` at the documented
   tolerance, not just `length(cs)` and `pip`?
8. **Degenerate-case test.** If the diff touches IBSS or finalize, does
   the susieR-degeneracy test (D11, `mf-ibss/spec.md`) still pass at
   tolerance 1e-10?
9. **Modularity test.** If the diff adds a new R file, did the author
   update or run `tests/testthat/test-modularity.R`?
10. **Roxygen cleanliness.** Does `roxygen2::roxygenise()` run without
    errors on the diff? Are documented defaults consistent with the
    function signatures?

A finding can be a *blocker* (must fix) or a *note* (should fix in a
follow-up). The default is blocker; the reviewer downgrades to note only
when explicitly justified.

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
