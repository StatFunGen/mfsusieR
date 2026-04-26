#!/usr/bin/env bash
# bootstrap-mfsusieR.sh
#
# One-time setup for the mfsusieR refactor workflow.
# Safe to re-run; idempotent where possible.
#
# Prereqs you must have installed BEFORE running this:
#   - git
#   - pixi (https://pixi.sh)  -- used to install nodejs (which provides npm
#                                for the openspec install)
#   - Claude Code  (https://docs.claude.com/en/docs/claude-code)
#   - R 4.3+ with devtools, testthat, roxygen2, covr (can be installed via pixi)
#
# Usage:
#   chmod +x bootstrap-mfsusieR.sh
#   ./bootstrap-mfsusieR.sh
#
# Env overrides:
#   GIT_ROOT (default: $HOME/GIT)

set -euo pipefail

GIT_ROOT="${GIT_ROOT:-$HOME/GIT}"
mkdir -p "$GIT_ROOT"
cd "$GIT_ROOT"

say() { printf "\n\033[1;34m==>\033[0m %s\n" "$*"; }
warn() { printf "\n\033[1;33mWARN:\033[0m %s\n" "$*"; }

say "Git root: $GIT_ROOT"

# ---------------------------------------------------------------------
# 1. Clone reference and target repositories
# ---------------------------------------------------------------------
# Format:  name|url|branch
#
# Paradigm references:
#   - susieR             : backbone we build on
#   - mvsusieR refactor-s3 : S3 pattern for multi-trait extension of susieR
#   - susieAnn           : S3 pattern for annotation-aware extension of susieR
#
# Manuscript:
#   - MultifSuSiE_Manuscript : LaTeX source, the math anchor for mfsusieR
#
# Port source:
#   - mvf.susie.alpha    : William's original multi-functional SuSiE.
#                          All functional/wavelet machinery comes from here;
#                          we reorganize it using the S3 patterns above.
#
# Target:
#   - mfsusieR           : the repo we write
#
# EDIT the susieAnn URL below to match where the repo actually lives.

REPOS=(
  "mfsusieR|https://github.com/StatFunGen/mfsusieR.git|main"
  "mvf.susie.alpha|https://github.com/william-denault/mvf.susie.alpha.git|main"
  "mvsusieR|https://github.com/stephenslab/mvsusieR.git|refactor-s3"
  "susieR|https://github.com/stephenslab/susieR.git|master"
  "susieAnn|https://github.com/StatFunGen/susieAnn.git|main"
  "MultifSuSiE_Manuscript|https://github.com/StatFunGen/MultifSuSiE_Manuscript.git|main"
)

for entry in "${REPOS[@]}"; do
  IFS='|' read -r name url branch <<< "$entry"
  if [[ -d "$name/.git" ]]; then
    say "$name: already cloned, fetching"
    (cd "$name" && git fetch --all --prune --quiet)
  else
    say "Cloning $name"
    git clone --quiet "$url" "$name" || warn "clone failed for $url -- continuing"
  fi
  if [[ -d "$name/.git" ]]; then
    (cd "$name" && git checkout --quiet "$branch" 2>/dev/null \
      || warn "branch '$branch' not found for $name; staying on default")
  fi
done

# ---------------------------------------------------------------------
# 2. nodejs (via pixi) and OpenSpec (via npm) -- planning layer
# ---------------------------------------------------------------------
say "Checking pixi"
command -v pixi >/dev/null 2>&1 \
  || { echo "Install pixi first: https://pixi.sh"; exit 1; }

if ! command -v node >/dev/null 2>&1; then
  say "Installing nodejs via pixi global"
  pixi global install nodejs
fi

if ! command -v openspec >/dev/null 2>&1; then
  say "Installing OpenSpec via npm"
  npm install -g @fission-ai/openspec

  # pixi global expose only sees conda-installed binaries, so the
  # npm-installed openspec binary needs a manual symlink onto PATH.
  if [[ -x "$HOME/.pixi/envs/nodejs/bin/openspec" ]]; then
    mkdir -p "$HOME/.pixi/bin"
    ln -sf "$HOME/.pixi/envs/nodejs/bin/openspec" "$HOME/.pixi/bin/openspec"
  fi
fi

say "OpenSpec present: $(openspec --version 2>/dev/null || echo unknown)"

# ---------------------------------------------------------------------
# 3. Claude Code sanity check
# ---------------------------------------------------------------------
command -v claude >/dev/null 2>&1 \
  || { echo "Install Claude Code first: https://docs.claude.com/en/docs/claude-code"; exit 1; }

# ---------------------------------------------------------------------
# 4. Initialize OpenSpec inside mfsusieR/inst
# ---------------------------------------------------------------------
# OpenSpec lives under inst/openspec/ (excluded from package builds via
# .Rbuildignore). The CLI resolves ./openspec/ from the cwd and has no
# flag to relocate it, so all openspec commands are run from inst/.

cd "$GIT_ROOT/mfsusieR"

if [[ ! -d "inst/openspec" ]]; then
  mkdir -p inst
  (cd inst && say "Initializing OpenSpec in mfsusieR/inst" && \
   openspec init --tools claude || warn "openspec init returned non-zero -- check the repo state")
fi

# ---------------------------------------------------------------------
# 5. Skeleton directories the CLAUDE.md references
# ---------------------------------------------------------------------
mkdir -p inst/notes/paradigms inst/notes/investigations inst/notes/sessions \
         bench/scenarios bench/profiling \
         tests/testthat/fixtures R

# ---------------------------------------------------------------------
# 6. Done
# ---------------------------------------------------------------------
cat <<'EOF'

==========================================================
Bootstrap complete.

1. Open a Claude Code session in the mfsusieR repo:

   cd ~/GIT/mfsusieR
   claude

2. Inside Claude Code, kick off Phase 1 by telling Claude:

   "Read CLAUDE.md and execute Phase 1 only. Stop when the four
    paradigm notes are committed. Do not start Phase 2."

   Phase 1 produces inst/notes/paradigms/*.md only. No code work.

3. After Phase 1 review, run a new session and ask for Phase 2.
   Phase 2 produces an OpenSpec proposal at
   inst/openspec/changes/add-mfsusier-s3-architecture/.

4. Phase 3 (implementation) uses a Claude-native multi-round review
   loop. See inst/notes/review-loop-methodology.md after Phase 1
   lands.

==========================================================
EOF
