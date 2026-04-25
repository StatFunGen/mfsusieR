# mfsusieR — agent brief

You are Claude Code working on **mfsusieR**, the refactored multi-functional SuSiE package.
The mission is to port and reimplement William Denault's `mvf.susie.alpha` on top of the
`susieR` backbone, matching the S3 paradigm in `mvsusieR/refactor-s3`, with reference-grade
unit tests, a clean path to Rcpp optimization, and an interface worth living with.

All of these packages are under ~/GIT folder. `mvf.susie.alpha` is a very messy package. 

The methodology for mfsusie is under /home/gw/GIT/MultifSuSiE_Manuscript . It is an extension of
the fSuSiE method which is under ~/GIT/fsusieR 

We use OpenSpec for phase planning and a Claude-native multi-round review loop for
implementation work (see `inst/notes/review-loop-methodology.md`). As you can see we would
like to keep the mfsusieR code base as short as possible with most of the structure sharing
susieR. You can use zzz.R to load internal functions without having to expose from susieR.

## How to start every session

1. Run `openspec list --changes` and `openspec list --specs` to see state.
2. Decide which phase you are in (see roadmap below). A phase is complete only when its
   OpenSpec change is archived.
3. Never start phase N+1 while phase N is unarchived.
4. If no active change and Phase 1 notes are missing, start Phase 1.
5. If no active change and all Phase 1 notes exist but no architecture proposal has been
   archived, the user should start Phase 2 (do NOT open that proposal autonomously;
   Gao reviews the Phase 1 notes first).

## Repo map (all under `~/GIT/`)

- **`susieR/`** — the backbone. SuSiE + IBSS loop + utilities we depend on. Read-only.
- **`mvsusieR/`** (branch `refactor-s3`) — S3 paradigm reference #1: how to extend
  `susieR` with multi-trait structure (MASH prior, cross-condition sharing).
  Read-only.
- **`susieAnn/`** — S3 paradigm reference #2: how to extend `susieR` with an
  additional prior structure (functional annotations, Poisson-Gamma activity
  indicator, pluggable predictor interface). Read-only.
- **`mvf.susie.alpha/`** — William's original multi-functional SuSiE. This is the
  **port source**: every functional/wavelet pathway in mfsusieR comes from here,
  reorganized using the paradigms from `mvsusieR/refactor-s3` and `susieAnn`.
  Every numerical result in our tests must match this at appropriate tolerance
  unless we document a reason to deviate. Read-only.
- **`mfsusieR/`** — the target. This is where you write.

`fsusieR` is a separate package unrelated to this refactor. Do not consult it.

## Tools and when to use them

| Tool | Use for |
|------|---------|
| Regular Claude Code | Phase 1 exploration, reading code, drafting notes |
| OpenSpec (`openspec new change / validate / archive`, run from `inst/`) | Phase 2, 6, 7, 8 planning |
| Claude-native review loop (see `inst/notes/review-loop-methodology.md`) | Phase 3, 4, 6, 7 implementation |
| `targets` pipeline + `bench` + `profvis` | Phase 5 calibration, Phase 7 profiling |

The review loop borrows the structure of an external code-review-bot RLCR loop
(write -> review -> revise) but runs entirely inside Claude Code by spawning
fresh-context Agent (Explore subagent) reviewers between Claude's authoring
passes. There is no external review CLI to install.

For all notes please save to inst/notes under this folder ~/GIT/mfsusieR

## Phases

### Phase 1 — paradigm study (notes only, no code)

Produce four Markdown files under `inst/notes/paradigms/`. No OpenSpec proposal.
Each note answers: (1) core data structures, (2) IBSS inner-loop shape,
(3) user-facing parameters, (4) what mfsusieR can inherit or mirror vs. must add.

Required files:

- `inst/notes/paradigms/susieR-backbone.md` — what `susieR::susie()` exposes. The fit
  object (`alpha`, `mu`, `mu2`, `KL`, `lbf`, `pip`, `sets`). The IBSS inner loop.
  Which helpers are exported vs. internal. Which functions a wrapper package is
  expected to call directly.

- `inst/notes/paradigms/mvsusieR-s3.md` — how `mvsusieR/refactor-s3` extends `susieR`
  with multi-trait structure. Class hierarchy, which methods are overridden,
  dispatch points, how the MASH prior plugs into the IBSS loop. This is paradigm
  reference #1 for mfsusieR.

- `inst/notes/paradigms/susieAnn-s3.md` — how `susieAnn` extends `susieR` with a new
  prior structure. Class hierarchy, where the activity-indicator / annotation
  prior enters the IBSS update, how the predictor interface is abstracted, how
  overdispersion enters. This is paradigm reference #2 for mfsusieR. Note which
  patterns differ from `mvsusieR/refactor-s3` and why.

- `inst/notes/paradigms/mvf-original.md` — William's `mvf.susie.alpha`. Public API,
  key numerical pathways including the wavelet basis construction and
  scale-specific prior handling, the FDR miscalibration Anjing observed (point
  to specific functions if you can identify them), what's elegant vs. what's
  wart. This is the **port source**: all functional machinery comes from here.
  The paradigms above tell us how to reorganize it; this note tells us what
  needs reorganizing.

Each note ends with a short section titled "Implications for mfsusieR" that
states in 3-5 sentences how this paradigm or port source informs the
mfsusieR architecture. Keep these sections factual, not prescriptive; the
architecture proposal is Phase 2.

Tone: mechanism-first prose. No filler ("striking," "notably," "importantly,"
"crucial," "robust," "pervasive," "leverage," "delve," "seamless"). No em dashes.
Cite file paths and line numbers using the form `pkg/R/file.R:123`.

Exit: commit the four notes with message `paradigm study: <one-line summary>`.
Then stop. The user will review before Phase 2.

### Phase 2 — architecture proposal (OpenSpec, no code)

Once the user says Phase 1 is approved, run from `inst/`:

```
cd inst
openspec new change add-mfsusier-s3-architecture
```

The proposal must cover:
- Class hierarchy. Which classes inherit from `susie` fit objects, which are new.
  How functional signals are represented (list column? separate S3 class? method
  dispatch on a `functional_response` class?).
- Public API. Signature of `mfsusie()` with argument names chosen per the naming
  rules below. Sketches of `predict.mfsusie()`, `summary.mfsusie()`,
  `coef.mfsusie()`, `plot.mfsusie()`.
- Delegation map. For each behavior in `mvf.susie.alpha::multfsusie`, state whether
  it will be (a) delegated to `susieR`, (b) reimplemented in `mfsusieR` following
  the `mvsusieR/refactor-s3` pattern, (c) reimplemented in `mfsusieR` following the
  `susieAnn` pattern, (d) reimplemented bespoke because neither paradigm fits, or
  (e) dropped. Each (d) or (e) entry needs a reason.
- Prior handling. How multi-trait residual covariance enters (paradigm from
  `mvsusieR`); how wavelet-scale prior variance enters (ported from
  `mvf.susie.alpha`, organized per `susieAnn`'s pattern for plugging an extra
  prior structure into the IBSS update); how the two compose when both are active.
- Task list. Small, reviewable PRs. No task should take more than a day of work to implement.

Run `openspec validate` from `inst/`. Human review required. Note that OpenSpec
1.3.1 has no separate `apply` step: the proposal is the contract for Phase 3
once validated, and `openspec archive` is the final merge step that runs after
Phase 4 tests pass.

### Phase 3 — core port (Claude-native review loop)

After Phase 2 proposal is approved by Gao, start Phase 3 implementation under
the Claude-native review loop documented in
`inst/notes/review-loop-methodology.md`. The loop borrows the structure of an
external write-review-revise bot but runs entirely inside Claude Code: an
authoring pass writes code per the proposal's spec.md contract, a fresh-context
Agent (Explore subagent) reviews, findings come back, the author addresses,
the cycle repeats until clean.

Rules for Phase 3:
- Every function that has an analog in `mvf.susie.alpha` cites its origin in roxygen:
  ```r
  #' @references_original mvf.susie.alpha/R/multfsusie.R#L123-L187
  ```
- Every formula cites the manuscript:
  ```r
  #' @manuscript_ref methods/derivation.tex eq:post_f_mix
  ```
- Before marking a task done, a fresh-context Agent reviewer must diff behavior
  against `mvf.susie.alpha` on one scripted example and report match/mismatch
  at the documented tolerance, AND audit modularity against design.md D10.
- No new algorithmic ideas in this phase. Faithful port first; innovations
  come in Phase 6.
- When porting a numerical routine, preserve the original parameter names
  inside the internal function. Renaming happens only at the public-API
  boundary per Phase 2.

Exit: Phase 4 reference tests pass against the port at the documented
tolerances, then `openspec archive add-mfsusier-s3-architecture` (run from
`inst/`) merges the deltas into `inst/openspec/specs/`.

### Phase 4 — reference unit tests (interleaved with Phase 3 under the same review loop)

For every ported function, write a testthat test locking its behavior against
`mvf.susie.alpha` output on a fixed seed. The canonical test shape:

```r
test_that("mfsusie matches mvf.susie.alpha on toy scenario A", {
  skip_if_not_installed("mvf.susie.alpha")
  data <- readRDS(test_path("fixtures/scenario_A.rds"))

  set.seed(1)
  fit_new <- mfsusie(data$X, data$Y, data$pos, L = 10,
                     residual_variance_method = "shared_per_modality")

  set.seed(1)
  fit_old <- mvf.susie.alpha::multfsusie(data$Y, data$X, data$pos, L = 10)

  # Per memory test_fidelity_and_manuscript_xref: every numeric output, not
  # just aggregates.
  expect_equal(fit_new$alpha,        fit_old$alpha,        tolerance = 1e-6)
  expect_equal(fit_new$mu,           fit_old$mu,           tolerance = 1e-6)
  expect_equal(fit_new$mu2,          fit_old$mu2,          tolerance = 1e-6)
  expect_equal(fit_new$lbf,          fit_old$lbf,          tolerance = 1e-6)
  expect_equal(fit_new$lbf_variable, fit_old$lbf_variable, tolerance = 1e-6)
  expect_equal(fit_new$KL,           fit_old$KL,           tolerance = 1e-6)
  expect_equal(fit_new$sigma2,       fit_old$sigma2,       tolerance = 1e-6)
  expect_equal(fit_new$elbo,         fit_old$elbo,         tolerance = 1e-6)
  expect_equal(fit_new$pip,          fit_old$pip,          tolerance = 1e-6)
  expect_equal(length(fit_new$cs),   length(fit_old$cs))
})
```

Fixture rules:
- `tests/testthat/fixtures/*.rds` are generated by `tests/testthat/scripts/make-fixtures.R`
  with a recorded seed. The script is committed; fixtures are committed if small (<1MB each),
  otherwise generated in the test setup.
- Scenarios must include: single-trait functional, multi-trait scalar, multi-trait
  functional, high-LD, low-PVE, and at least one large-null case (no real signal).

Exit: `R CMD check --as-cran` returns 0 errors, 0 warnings, and only notes that are not
about our code (e.g., known susieR notes). `covr::package_coverage()` ≥ 80% on exported
functions. Both runs saved under `bench/ci-logs/` with a commit hash.

### Phase 5 — FDR calibration investigation (targets pipeline, notes only)

Goal: diagnose the FDR miscalibration Anjing observed. This is investigation, not
refactoring. No code changes; no OpenSpec change during the investigation.

Set up `bench/_targets.R`:
- Scenarios: single-trait functional, multi-trait functional, with/without correlated
  residuals, across sample sizes {100, 400, 1000}, PVE {0.01, 0.05, 0.1},
  `max_scale` {6, 8, 10}.
- Per combination: observed FDR at nominal 0.05 PIP threshold, PIP calibration curve
  (binned observed vs. claimed), credible-set coverage.
- Reproducible seeds. Output: one `html` per scenario family in `inst/notes/investigations/fdr-plots/`.

For each hypothesis, write `inst/notes/investigations/fdr-<short-name>.md` with sections:
*Hypothesis, Evidence, Test, Result, Conclusion*. One hypothesis per file. Commit as
you go. Do not collapse multiple hypotheses into one note.

Exit: a summary note `inst/notes/investigations/fdr-summary.md` naming the specific
components responsible (residual variance estimation? scale-specific prior variance?
CS construction? screening threshold?). Each identified issue is a candidate OpenSpec
proposal in Phase 6.

### Phase 6 — targeted fixes (OpenSpec + Claude-native review loop, one per issue)

For each issue from Phase 5, open a focused OpenSpec change with a measurable
success criterion. Example:

> `change: fix-residual-variance-under-correlated-noise`
> Success: observed FDR at 0.05 PIP threshold matches nominal within ±0.01 across
> all scenarios defined in bench/scenarios/fdr/*.

Implement under the same Claude-native review loop as Phase 3 (see
`inst/notes/review-loop-methodology.md`). No changes outside the scope declared
in the proposal's delta spec.

### Phase 7 — Rcpp optimization (profiling first, proposal second, implementation last)

Only after Phase 3-6 land and behavior is locked. Workflow:

1. Profile. Run `bench/profiling/profile-<scenario>.R` with `profvis`, save flamegraphs
   under `bench/profiling/flamegraphs/` with a commit hash.
2. Identify the top 3 hot paths from the flamegraph.
3. For each hot path, open an OpenSpec proposal: `perf-<component>-rcpp`. Proposal
   must state the current wall-clock on the bench harness, the target wall-clock,
   and the C++ sketch (headers, main loop, Armadillo types).
4. Implement under the Claude-native review loop. Every Rcpp function has a
   pure-R reference implementation that the test suite compares against at
   tolerance 1e-10.

Do not start Rcpp work before Phase 6 is archived.

### Phase 8 — interface polish

Parameter and function names are frozen after Phase 2. If during Phase 7 you see a
name that reads badly (`max_SNP_EM`, `null_weight`, `max_scale`), open a single
OpenSpec change `rename-public-api-v1` with a `lifecycle::deprecate_warn()` path
and keep the old names working for one minor version.

## Naming rules (binding for Phases 2 and 8)

- Snake_case for function and argument names.
- No abbreviations unless they are established in the SuSiE lineage (`pip`, `cs`,
  `lbf`, `elbo`, `kl`). `max_SNP_EM` is not established — rename.
- Parameter names describe what the parameter IS, not what it does. `prior_variance`
  not `V`. `max_iter` not `max_SNP_EM`. `null_prior_weight` not `null_weight`.
- Public S3 methods: `mfsusie()`, `predict.mfsusie()`, `summary.mfsusie()`,
  `coef.mfsusie()`, `plot.mfsusie()`, `print.mfsusie()`.
- Internal helpers prefixed with `.` or in `R/utils-*.R`.
- No boolean flags named in the negative (`no_cache = FALSE` becomes `cache = TRUE`).

## Writing rules (all markdown, roxygen, and commit messages)

The lab's `statgen-writing-style` conventions apply. In short:

- No em dashes. Use commas, periods, or parentheses.
- No filler: "striking," "notably," "importantly," "crucial," "robust,"
  "pervasive," "leverage," "delve," "seamless," "powerful," "unique," "significant"
  when used as decoration rather than statistical significance.
- Mechanism-first. Describe what a thing does before why it matters.
- No review-style framing in Results or in commit bodies.
- Quantitative claims cite their source. Commit hash, fixture filename, or table ID.
- Avoid colons in subheadings.

## Commits and PRs

- Conventional-commit prefixes: `feat:`, `fix:`, `test:`, `docs:`, `refactor:`,
  `perf:`, `chore:`.
- One logical change per commit. One feature per PR.
- Every PR description links to the OpenSpec change it implements. If there is no
  OpenSpec change, there should be no PR (except for typo-level `docs:` fixes).
- Tests and code land in the same PR.

## Hard rules

1. `mvf.susie.alpha/`, `susieR/`, `mvsusieR/`, `susieAnn/` are read-only reference.
   Do not edit. The single authorized exception is the `feature/L_greedy` branch on
   `~/GIT/susieR`, used to add the `L_greedy` argument to `susie_workhorse` per the
   external coordination plan in
   `inst/openspec/changes/add-mfsusier-s3-architecture/design.md` (Migration Plan).
2. `fsusieR` is not a paradigm source. Do NOT model mfsusieR architecture decisions
   on fsusieR's structure. Do NOT cite fsusieR alongside mvsusieR/refactor-s3 or
   susieAnn as a paradigm reference. Specific *runtime utility* delegations to
   fsusieR (e.g., `fsusieR::remap_data`, `fsusieR::colScale`,
   `fsusieR::gen_wavelet_indx`, `fsusieR::init_prior.default`, smash/TI/HMM
   smoothers) ARE allowed and are documented in design.md D5 and D8. The rule is
   about design influence, not function calls.
3. Do not run expensive simulations without declaring an estimated runtime at the
   top of the script. If projected wall-clock exceeds 30 minutes, ask the user first.
4. Do not start phase N+1 while phase N is not archived in OpenSpec. Note that
   Phase 1 has no OpenSpec change (notes only), so the Phase 1 -> Phase 2 transition
   only requires Gao's review of the paradigm notes; for all other transitions, the
   Phase N OpenSpec change must be archived first.
5. If a design decision comes up that is not covered by an OpenSpec proposal, stop
   and open one. Do not decide silently.
6. If you find yourself reframing a numerical discrepancy against `mvf.susie.alpha`
   as "acceptable," that is the signal to ask the user, not to proceed.

## Reporting at end of session

Before ending a session, write `inst/notes/sessions/YYYY-MM-DD-HHMM.md` with:
- What was done.
- What phase it was.
- OpenSpec state after the session (what's active, what's archived).
- What the next session should pick up.
- Token usage estimate if available.

Then stop. The next session reads that note first.
