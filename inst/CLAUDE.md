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
- **`mvf.susie.alpha/`** — William's original multi-functional SuSiE. **Port
  source #1**: every multi-modality functional pathway in mfsusieR comes
  from here, reorganized using the paradigms from `mvsusieR/refactor-s3`
  and `susieAnn`. Apple-to-apple equivalence contract #3 (M ≥ 1, legacy
  variance mode) compares against `mvf.susie.alpha::multfsusie`.
  Read-only.
- **`fsusieR/`** — single-modality functional SuSiE. **Port source #2**:
  the routines we need (`remap_data`, `colScale`, `gen_wavelet_indx`,
  `init_prior.default`, smash / TI / HMM smoothers) are ported into
  `mfsusieR/R/`, not delegated. Apple-to-apple equivalence contract #2
  (M = 1, T_1 > 1) compares against `fsusieR::susiF`. fsusieR is a port
  source ONLY; it is not a paradigm source for architecture decisions
  (mvsusieR and susieAnn fill that role) and not a runtime dependency.
  fsusieR's `EBmvFR` is out of scope. Read-only.
- **`mfsusieR/`** — the target. This is where you write.

## Tools and when to use them

| Tool | Use for |
|------|---------|
| Regular Claude Code | Phase 1 exploration, reading code, drafting notes |
| OpenSpec skills (`/openspec-explore`, `/openspec-propose`, `/openspec-apply-change`, `/openspec-archive-change`) | Phase 2, 6, 7, 8 design + implementation. Skills are the **preferred entry points** for OpenSpec work. |
| OpenSpec CLI (`openspec validate`, `openspec status`, `openspec list`, run from `inst/`) | In-place refinement of already-scaffolded changes; safety-net validation after manual edits |
| Claude-native review loop (see `inst/notes/review-loop-methodology.md`) | Phase 3, 4, 6, 7 implementation |
| `targets` pipeline + `bench` + `profvis` | Phase 5 calibration, Phase 7 profiling |

### OpenSpec skill flow (binding for all future design work)

When a design task requires opening or substantially revising an
OpenSpec change, invoke the OpenSpec skills rather than hand-editing
artifacts from a blank page. The mapping:

- **Start of a design phase, requirements are unclear or in flux:**
  invoke `/openspec-explore` (or the experimental
  `/opsx:explore`). The skill puts the session in explore mode, a
  thinking partner for clarifying requirements with Gao before any
  artifact is drafted.
- **Scaffolding a new change:** invoke `/openspec-propose` (or
  `/opsx:propose`). The skill generates `proposal.md`, `design.md`,
  `tasks.md`, and the spec capabilities directory in one step. Refine
  manually after.
- **Implementing tasks from an existing change:** invoke
  `/openspec-apply-change` (or `/opsx:apply`). The skill walks the
  task list and helps drive the Claude-native review loop on each
  PR group.
- **Archiving a completed change:** invoke
  `/openspec-archive-change` (or `/opsx:archive`) after Phase 4
  reference tests pass. The skill performs the merge of spec deltas
  into `inst/openspec/specs/`.
- **In-place refinement of an already-scaffolded change** (the case
  during the 2026-04-25 architectural pivot): hand-editing +
  `openspec validate` is acceptable. Skills add little when the
  artifacts already exist and only need section-level rewrites.

For every step where the skill flow applies, prefer the skill. If
the skill does not match the work being done (e.g., a tiny typo fix
in a spec, or a one-line rename), hand-editing + `openspec validate`
is the right choice. When in doubt, ask Gao before opening an
OpenSpec change without a skill.

The review loop borrows the structure of an external code-review-bot
RLCR loop (write -> review -> revise) but runs entirely inside
Claude Code by spawning fresh-context Agent (Explore subagent)
reviewers between Claude's authoring passes. There is no external
review CLI to install.

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

Rules for Phase 3 (binding; mirrors the code-refactor skill at
`~/Documents/obsidian/AI/general/agents/AGENT-code-refactor.md`):

- **Mathematical fidelity first.** Every line of the original is
  accounted for. Lines intentionally omitted are logged in
  `inst/notes/refactor-exceptions.md` with file/range, behavior,
  decision, reason. Verify math independently of code (read the
  manuscript and supplementary material; the original may have bugs).
- **Manuscript citations** are required in main code roxygen.
  Use the built-in `@references` tag with a `Manuscript:` prefix
  so roxygen2 accepts the tag and the citation is greppable:
  ```r
  #' @references
  #' Manuscript: methods/derivation.tex eq:post_f_mix
  ```
- **Original-code references are FORBIDDEN in main code roxygen.**
  No `@references_original`, no `# from mvf.susie.alpha/...`, no
  `# original implementation` comments under `R/`. Provenance lives
  in (a) test file headers (apple-to-apple comparison metadata),
  (b) `inst/notes/refactor-exceptions.md`, (c) free-form prose in
  session notes and paradigm notes.
- **Binary tolerance philosophy.** Apple-to-apple comparisons (same
  algorithm, same code path) use tolerance `<= 1e-8`. Apple-to-orange
  comparisons (genuinely different algorithms) are smoke tests only,
  not numerical. Graduated tolerances (`1e-2`, `5e-2`, `1e-6`) are
  forbidden. A failing apple-to-apple test at `1e-8` is a bug to
  investigate, not a tolerance to relax.

  **The ONLY exception is a documented port-source bug fix.** When
  mfsusieR fixes a buggy or numerically noisy formula in
  `mvf.susie.alpha` / `fsusieR`, the affected unit test EITHER
  (a) directly asserts the new corrected behaviour and does NOT
  expect the old buggy output (Pattern A), OR (b) keeps asserting
  equivalence with the port source but uses a relaxed tolerance
  (e.g., `1e-12`) well below the contract floor, with an inline
  comment saying "this is a port-source-bug fix, not numerical
  drift" and a refactor-exceptions ledger entry per design.md
  D11e (Pattern B). The relaxed `1e-12` tolerance is NOT a
  graduated tolerance: it is four orders of magnitude tighter
  than the contract floor and catches real drift while accepting
  the port source's machine-epsilon rounding noise. The full
  policy lives in `inst/notes/refactor-discipline.md` sections 3
  and 4.
- **Reviewer pass on every numerical PR**: fresh-context Agent
  reviewer runs the 10-item checklist in
  `inst/notes/review-loop-methodology.md`, which includes the
  modularity audit (design.md D10), the manuscript cross-check, the
  refactor-exceptions completeness check, and the binary-tolerance
  test classification.
- **No new algorithmic ideas in this phase.** Faithful port first;
  innovations come in Phase 6 after the FDR investigation.
- **No quick fixes.** When tests fail, investigate root cause. Do not
  paper over with `as.character()`, `tryCatch` swallowers, or
  loosened tolerances. If a fix feels hacky, it is hacky.
- **Public-API parameter names** follow the harmonization with susieR
  / mvsusieR documented in design.md D7. Internal function
  parameter names within ported routines preserve the original to
  ease line-by-line review; renaming happens only at the public
  boundary.

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
  # Apple-to-apple: same algorithm, same code path -> tolerance <= 1e-8.
  # See design.md D11b (binary tolerance philosophy).
  expect_equal(fit_new$alpha,        fit_old$alpha,        tolerance = 1e-8)
  expect_equal(fit_new$mu,           fit_old$mu,           tolerance = 1e-8)
  expect_equal(fit_new$mu2,          fit_old$mu2,          tolerance = 1e-8)
  expect_equal(fit_new$lbf,          fit_old$lbf,          tolerance = 1e-8)
  expect_equal(fit_new$lbf_variable, fit_old$lbf_variable, tolerance = 1e-8)
  expect_equal(fit_new$KL,           fit_old$KL,           tolerance = 1e-8)
  expect_equal(fit_new$sigma2,       fit_old$sigma2,       tolerance = 1e-8)
  expect_equal(fit_new$elbo,         fit_old$elbo,         tolerance = 1e-8)
  expect_equal(fit_new$pip,          fit_old$pip,          tolerance = 1e-8)
  expect_equal(length(fit_new$cs),   length(fit_old$cs))
})
```

Fixture rules:
- `tests/testthat/fixtures/*.rds` are generated by `tests/testthat/scripts/make_fixtures.R`
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
- Snake_case for filenames in `R/`, `tests/testthat/`, and
  `tests/testthat/scripts/`. Hyphens are not allowed in source filenames.
  ISO date filenames in `inst/notes/sessions/` (`YYYY-MM-DD-HHMM.md`) keep
  their hyphens by convention.
- No abbreviations unless they are established in the SuSiE lineage (`pip`, `cs`,
  `lbf`, `elbo`, `kl`). `max_SNP_EM` is not established — rename.
- Parameter names describe what the parameter IS, not what it does. `prior_variance`
  not `V`. `max_iter` not `max_SNP_EM`. `null_prior_weight` not `null_weight`.
- Public S3 methods: `mfsusie()`, `predict.mfsusie()`, `summary.mfsusie()`,
  `coef.mfsusie()`, `plot.mfsusie()`, `print.mfsusie()`.
- Internal helpers prefixed with `.` or in `R/utils_*.R`.
- Class-related code may be split by lifecycle phase, not by
  arbitrary "class vs methods" axis. The split mfsusieR uses for
  `mf_individual`:
  - `R/individual_data_class.R` -- constructor (`create_mf_individual`)
  - `R/individual_data_methods.R` -- per-effect SER-step methods
    (`compute_residuals`, `compute_ser_statistics`, `loglik`,
    `calculate_posterior_moments`, `optimize_prior_variance`,
    `update_fitted_values`, etc.)
  - `R/ibss_methods.R` -- per-iteration IBSS-loop methods
    (`update_variance_components`, `update_model_variance`,
    `Eloglik`, `get_objective`, `ibss_initialize`, `cleanup_model`,
    finalize accessors)
  Same pattern for `mfsusie` model class:
  - `R/model_class.R` -- constructor + per-effect accessors
    (`get_alpha_l`, `get_posterior_mean_l`, etc.)
  - `R/mfsusie_methods.R` -- user-facing S3 methods
    (`predict.mfsusie`, `coef.mfsusie`, `summary.mfsusie`,
    `print.mfsusie`)
  The lifecycle split keeps file sizes manageable (none > 600
  lines) and groups code by when it runs rather than which class
  it belongs to. Each file has a one-paragraph header that names
  the lifecycle phase it covers. Section divider format inside a
  file:
  ```r
  # =====================================================================
  # Per-effect getter / setter accessors on the `<class>` model
  # =====================================================================
  ```
- No boolean flags named in the negative (`no_cache = FALSE` becomes `cache = TRUE`).
- The full public-API harmonization with susieR and mvsusieR (which choices
  match, which diverge, and why) lives in design.md D7. Phase 8 polish follows
  that table.

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

1. `mvf.susie.alpha/`, `fsusieR/`, `susieR/`, `mvsusieR/`, `susieAnn/` are
   read-only reference. Do not edit. The single authorized exception is the
   `feature/L_greedy` branch on `~/GIT/susieR`, used to add the `L_greedy`
   argument to `susie_workhorse` per the external coordination plan in
   `inst/openspec/changes/add-mfsusier-s3-architecture/design.md` (Migration
   Plan).
2. `fsusieR` is a port source, not a runtime dependency. mfsusieR replaces
   fsusieR by making the single-modality functional case (`M = 1`, `T_1 > 1`)
   a special case of `mfsusie()`, surfaced through the thin wrapper
   `mfsusieR::fsusie()`. mfsusieR therefore SHALL NOT have `fsusieR` in
   `DESCRIPTION` Imports, SHALL NOT call `fsusieR::*` at runtime, and SHALL
   NOT model architecture decisions on fsusieR's structure (mvsusieR and
   susieAnn remain the paradigm references). Routines we need from fsusieR
   (`remap_data`, `colScale`, `gen_wavelet_indx`, `init_prior.default`,
   smash / TI / HMM smoothers) are ported into `R/utils_wavelet.R` and
   `R/post_processing.R`. The fsusieR/susiF M=1 path is one of the three
   apple-to-apple equivalence contracts (alongside susieR scalar and
   mvf.susie.alpha legacy). When fsusieR has a bug, we fix it; the
   equivalence test asserts the deviation explicitly with a citation to
   the OpenSpec change that authorized it.
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
