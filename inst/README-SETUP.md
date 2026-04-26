# mfsusieR refactor kit — how to use

Three files in this kit:

1. `bootstrap-mfsusieR.sh` — one-time setup. Run on your workstation.
2. `CLAUDE.md` — the agent brief. Goes into `~/GIT/mfsusieR/CLAUDE.md`.
3. This file — instructions for you, not the agent.

## Prereqs (check these first)

- `git`, `R --version` (≥ 4.3)
- `pixi` (https://pixi.sh) for tool installs. nodejs and openspec are
  installed via pixi by the bootstrap script.
- Claude Code installed: https://docs.claude.com/en/docs/claude-code

If any of these are missing, install before running bootstrap. The script will fail fast
and tell you which one is missing.

## One-time setup

```bash
# 1. Drop the two files somewhere
cp bootstrap-mfsusieR.sh ~/
cp CLAUDE.md /tmp/mfsusieR-CLAUDE.md    # we'll move it in after clone

# 2. Run the bootstrap
chmod +x ~/bootstrap-mfsusieR.sh
~/bootstrap-mfsusieR.sh

# 3. Install the CLAUDE.md into the repo
cp /tmp/mfsusieR-CLAUDE.md ~/GIT/mfsusieR/CLAUDE.md
cd ~/GIT/mfsusieR
git add CLAUDE.md
git commit -m "chore: add agent brief for refactor workflow"

# 4. Open Claude Code in the repo
claude
```

## Kicking off Phase 1 (the first session)

Tell Claude Code:

> Read CLAUDE.md. Execute Phase 1 only. Stop when the four paradigm notes are
> committed. Do not start Phase 2.

This session is regular Claude Code (no OpenSpec proposal yet). It reads
`susieR`, `mvsusieR/refactor-s3`, `susieAnn`, and `mvf.susie.alpha`, then writes
four Markdown notes under `inst/notes/paradigms/`. Expected: 30-90 minutes of
wall time.

When it finishes, you read the four notes. This is the most important human review
in the whole project. Everything downstream assumes the paradigm analysis is right.
Things to check:
- Is the `mvsusieR/refactor-s3` S3 pattern described accurately? (paradigm reference #1)
- Is the `susieAnn` S3 pattern described accurately, with clear contrast vs. mvsusieR?
  (paradigm reference #2)
- Does the `mvf-original.md` note identify plausible candidates for the FDR bug?
- Does the `susieR-backbone.md` note correctly describe what's exported vs. internal?
- Does each note's "Implications for mfsusieR" section make sense?

If anything is off, correct the notes by hand, commit the corrections, and tell the
next session to re-read the notes before Phase 2.

## Kicking off Phase 2 (next session)

After Phase 1 notes are reviewed and committed, a new session:

> Phase 1 notes are approved. Execute Phase 2 only. Produce the OpenSpec proposal
> `add-mfsusier-s3-architecture` (kebab-case, lowercase) under `inst/openspec/changes/`
> and run `openspec validate` from `inst/`. Do not archive.

Output: one OpenSpec change folder with `proposal.md`, `design.md`, delta specs,
and `tasks.md`. You review, edit in place, then approve. Note that OpenSpec
1.3.1 has no separate `apply` command: validated specs are the contract for
Phase 3, and `openspec archive` is the final step that runs after Phase 4
tests pass.

## Kicking off Phase 3 (the implementation session)

After you've approved the architecture proposal:

> Phase 2 architecture is approved. Start Phase 3 implementation under the
> Claude-native review loop documented in `inst/notes/review-loop-methodology.md`.
> Begin with PR group 1 (package skeleton).

Phase 3 uses a multi-round review loop entirely inside Claude Code: an
authoring pass writes code per the proposal spec contract, a fresh-context
Agent (Explore subagent) reviews, findings come back, the author addresses,
the cycle repeats until clean. Multiple sessions will be needed to complete
Phase 3.

## What the agent WILL NOT do on its own

- Open any OpenSpec proposal without explicit approval.
- Start Phase N+1 before Phase N's OpenSpec change is archived (with the
  exception of Phase 1 → Phase 2, since Phase 1 produces no OpenSpec change).
- Run simulations projected to exceed 30 min wall-clock without asking.
- Edit `mvf.susie.alpha/`, `susieR/` (except the authorized `feature/L_greedy`
  branch), `mvsusieR/`, or `susieAnn/`. These are read-only reference.
- Start Rcpp porting before Phase 6 is done.

## Handy commands after bootstrap

```bash
# see OpenSpec state (must be run from inst/)
cd ~/GIT/mfsusieR/inst && openspec list --changes && openspec list --specs

# run the test suite
R -e 'devtools::test()'

# run R CMD check
R -e 'devtools::check()'

# run the bench harness (Phase 5+)
R -e 'targets::tar_make(script = "bench/_targets.R")'

# pick up where last session left off
ls -t inst/notes/sessions/ | head -1
```

## Honest scope note

The roadmap covers eight phases. Phases 1-4 alone are probably two to four weeks of
calendar time if you're reviewing attentively. Don't expect one weekend to get through
the port. The right mental model is: you're setting up a pipeline where each session
moves the state forward one phase (or partial phase), and OpenSpec + the session notes
carry state between them. The value isn't any one session, it's that no session loses
what the previous one learned.
