# OpenSpec primer for mfsusieR

A short orientation. Read this before approving Phase 2.

## What OpenSpec is

OpenSpec is a planning layer that sits on top of git. It tracks two trees in
the repo, both kept under `inst/openspec/`:

- `inst/openspec/specs/` - the source of truth for what the package currently
  delivers. One folder per *capability*. Each contains a `spec.md` with
  testable requirements written as `WHEN ... THE SYSTEM SHALL ...`.
- `inst/openspec/changes/` - one folder per *open proposal*. Each holds:
  - `proposal.md` - why this change, what shifts.
  - `design.md` - architecture, trade-offs, S3 dispatch decisions, etc.
  - `tasks.md` - small, reviewable implementation steps.
  - delta spec fragments (`specs/<capability>/spec.md`) showing what gets
    `ADDED`, `MODIFIED`, or `REMOVED` in `inst/openspec/specs/`.

Once a change is applied, its deltas merge into `inst/openspec/specs/`, and
the proposal moves to `inst/openspec/changes/archive/`. The archive plus the
current `inst/openspec/specs/` together describe the package's full design
history.

The `inst/openspec/` directory is excluded from package builds via
`.Rbuildignore` so it does not ship with the installed package.

## Working directory for openspec commands

The openspec CLI looks for `./openspec/` in the current working directory
and has no flag to relocate it. Because we keep it under `inst/`, all
openspec commands must be run from `inst/`:

```bash
cd ~/GIT/mfsusieR/inst
openspec list --changes
openspec new change <name>
openspec validate <name>
openspec apply <name>
openspec archive <name>
```

The Claude Code slash commands (`/opsx:propose`, etc.) wrap the CLI; they
likewise need to run with `inst/` as the working directory. If a command
fails with "No OpenSpec changes directory found", check cwd first.

## The four-step lifecycle

For each phase or each focused fix, this is the sequence:

1. `openspec new change <name>` (or `/opsx:propose <description>` inside
   Claude Code). Scaffolds `openspec/changes/<name>/` with placeholder
   artifacts.
2. Iterate on `proposal.md`, `design.md`, `tasks.md`, and the delta specs
   until they read well. Gao reviews on this step.
3. `openspec validate <name>` lints structure, required sections, and delta
   format. No state change.
4. `openspec apply <name>` merges the deltas into `openspec/specs/` (the
   change is now considered "in flight"). After implementation lands and
   tests pass, `openspec archive <name>` moves the change to
   `openspec/changes/archive/`.

`openspec list --changes` and `openspec list --specs` show current state.

## How mfsusieR's phases map to OpenSpec changes

Per the project CLAUDE.md, the package goes through eight phases. OpenSpec is
involved from Phase 2 onward:

| Phase | OpenSpec change name | Notes |
|---|---|---|
| 1 | none | paradigm notes only, already committed |
| 2 | `add-mfsusieR-s3-architecture` | one big proposal: classes, public API, delegation map, prior composition, task list |
| 3 | (same as Phase 2) | Phase 3 is the implementation that *applies* the Phase 2 change |
| 4 | (same as Phase 2) | reference tests interleaved; lands under the same change |
| 5 | none | FDR investigation, notes only |
| 6 | one per identified bug, e.g. `fix-residual-variance-under-correlated-noise` | each with a measurable success criterion |
| 7 | one per Rcpp hot path, e.g. `perf-ser-stats-rcpp` | profiling first, then proposal |
| 8 | `rename-public-api-v1` | with `lifecycle::deprecate_warn()` shims |

Phases 1 and 5 are notes-only because they are research, not commitments.
Everything else gets a proposal.

## What a good proposal contains

For Phase 2 specifically (just so it's concrete):

- **proposal.md** - "we extend susieR with an S3 hierarchy that combines the
  mvsusieR pattern (data-class method overrides) and the susieAnn pattern
  (predictor-style pluggability)". Names the FDR-bug investigation as a
  later change, not as in-scope here.
- **design.md** - class hierarchy diagram, public API signature with
  argument names and defaults, delegation map (per behavior in
  mvf.susie.alpha::multfsusie, state inherit / mvsusie-style override /
  susieAnn-style plug-in / bespoke / drop), prior composition story,
  alternatives considered.
- **tasks.md** - small PRs. None should take more than a day. Each PR has a
  clear deliverable and lists which tests gate its merge.
- **delta specs** - one spec.md per new capability under
  `inst/openspec/changes/add-mfsusieR-s3-architecture/specs/<capability>/`.
  Capabilities I expect for mfsusieR: `multifunctional-data-class`,
  `multifunctional-prior`, `multifunctional-ibss`, `mfsusie-public-api`,
  `mfsusie-credible-sets`. Each spec.md uses `## ADDED Requirements`
  sections.

## Slash commands installed by `openspec init`

The init step also installs four Claude Code slash commands and four
matching skills under `.claude/`:

- `/opsx:propose <name-or-description>` - scaffolds a new change and walks
  through filling in the artifacts.
- `/opsx:apply <name>` - merges deltas into specs.
- `/opsx:archive <name>` - moves an applied change to archive.
- `/opsx:explore <name>` - read-only inspection.

These are convenience wrappers; the underlying `openspec` CLI does the same
thing. Either is fine.

## Failure modes to watch for

- *Skipping validate before apply*: validate is cheap, run it.
- *Cramming multiple capabilities into one change*: each change should map
  to one architectural decision or one bug fix. Phase 2 is the exception
  (it sets up the whole architecture in one go), and it is large because
  the package is starting from scratch.
- *Editing applied specs directly*: once `openspec/specs/` has content, do
  not hand-edit it; produce a new change with `MODIFIED` deltas instead.
  This keeps history readable.
- *Phases out of order*: CLAUDE.md hard rule #4. Do not start Phase N+1
  while Phase N is unarchived. Treat the OpenSpec archive as the gate.

## Where to find more

- `openspec --help` for the CLI.
- `https://github.com/Fission-AI/OpenSpec` for the upstream.
- `.claude/skills/openspec-propose/SKILL.md` and the others for the agent-
  side procedures.
