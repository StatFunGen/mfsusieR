# SERIOUS Cross-package audit after the SER hook + per-effect prior refactor

## Why

Production has hit several silent-numerical and convergence
regressions over the last week, all rooted in mfsusieR drifting
from `fsusieR` / `mvf.susie.alpha` behavior on the SER hot path:

- A pre-loglik M-step that collapses pi to all-null on hard
  workloads (cascade bug, fixed by switching to a post-loglik
  M-step in the new hook design).
- A shared-G_prior-across-L design that let one effect's M-step
  poison the BF for every other effect (fixed by adding the
  `model$pi_V[[l]]` and `model$fitted_g_per_effect[[l]][[m]][[s]]`
  per-effect sidecars and the pre-loglik G_prior swap).
- A V scalar that arrived in `loglik` as a length-1 list instead
  of a numeric scalar (fixed by hardcoding V = 1 in both hooks
  for the mfsusie path).
- A 50-iter convergence stall driven by spurious-effect drift on
  noise correlations (partial fix via the inner-EM loop in
  `post_loglik_prior_hook.mf_individual` with
  `max_inner_em_steps`; outer-iter count drops from 50 to 6 on
  the getting_started multi case but the algorithm still takes
  ~2x more outer iters than fsusieR on the same fixture).

Each of these surfaced AT THE USER, not in CI. Each cost
multiple hours to track down. The accumulating delta from
upstream is a serious threat to release quality.

The prior multi-agent audit (2026-04-26) covered round 3
(preprocessing, model_init, scalewise SD fix). Since then ~12
commits landed substantially modifying the SER scaffolding,
prior storage, M-step dispatch, and convergence path:

```
9af590e Improve inner em steps count
daf47c3 Add back an inner em to pre-hook to sync pip and alpha
a8397e6 fix documentation
2e8e11f Fix V deminsions bug
6a54da0 Add hook to prior
78aa3f7 Minor fixes
c22f6ce Change qnorm back to FALSE by default; add standardize to TRUE by default
4274418 Update vignettes
6f3ad9d Use same variants kept for mixsqp as ebnm
2528620 Add EBNM priors
e266733 Fix customized EM bugs
9e103f3 Draft the per_scale_normal prior
```

A FRESH cold-start audit is required. The prior round's audit
prompts are not reusable as-is because the SER scaffolding and
prior storage shape have both changed.

## Known divergences (do NOT re-flag)

The following are documented and intentional. Audit agents must
treat them as accepted and not flag them.

### From `~/Documents/mfsusier-mvf-divergences/divergences.tex`

- **Divergence 1 — Expected residual sum of squares.** mfsusieR
  has the correct `ER2_{m,t}` decomposition; mvf.susie.alpha
  misses the `predictor_weights` factor in the bias correction
  and aggregates effects via sum-then-square. Pattern A
  regression test at
  `tests/testthat/test_mvf_alpha_degeneracy.R:114-129`.
- **Divergence 2 — ELBO sign on the SER posterior expected
  log-likelihood.** mfsusieR follows the correct susieR
  decomposition `Eloglik - sum(KL)`; mvf has a sign error in
  `cal_KL_l` (`R/ELBO_mutlfsusie.R:269-282`) that the entropy
  term `o = sum_l[log(p) + H(alpha_l)]` cannot offset.
- **Divergence 3 — Per-variant SE estimator.** mfsusieR uses
  the SuSiE backbone `Shat^2 = sigma2 / colSums(X^2)` (joint
  sigma2 per outcome / per scale); fsusieR / mvf use a
  per-(variable, outcome) marginal residual variance from a
  separate univariate regression. Different model formulation,
  not a bug.

### From `inst/notes/refactor-exceptions.md`

Every entry under "In-scope omissions" (PR groups 2, 4, 6 plus
the 2026-04-26 cross-package audit fix-now bundle plus the
2026-04-26 covariate-adjustment / em_helpers entries). Agents
read this file in full before flagging.

## What changes

### 1. Three Opus agents in parallel (Phase A)

Each agent runs from a fresh context. Each receives:

1. The known-divergence ledger above (verbatim).
2. The current `inst/notes/refactor-exceptions.md`.
3. A pointer to the recent commit list above so the agent
   knows what changed since 2026-04-26.
4. A specific scope (which upstream package to compare
   against and which mfsusieR pathway to trace).

#### Agent A — mfsusieR vs susieR

Focus: SER scaffolding integration. The recently-landed
hook generics (`pre_loglik_prior_hook`, `post_loglik_prior_hook`)
are exported from susieR and dispatched in mfsusieR. Audit
the integration end-to-end:

- Each hook's signature, return contract, and side-effect
  surface match what susieR's `single_effect_regression`
  expects.
- mfsusieR's `optimize_prior_variance.mf_individual`,
  `compute_residuals.mf_individual`,
  `compute_ser_statistics.mf_individual`,
  `loglik.mf_individual`,
  `calculate_posterior_moments.mf_individual`,
  `compute_kl.mf_individual`,
  `update_fitted_values.mf_individual`,
  `Eloglik.mfsusie`,
  `get_objective.mfsusie`,
  `check_convergence.mfsusie` (or default),
  `trim_null_effects.mf_individual` (or default),
  `cleanup_extra_fields.mf_individual`,
  `cleanup_model.mf_individual` are consistent with the
  susieR generic contract.
- The per-effect persistent state (`model$pi_V[[l]]`,
  `model$fitted_g_per_effect[[l]][[m]][[s]]`) is correctly
  preserved across IBSS iters, L_greedy expansions, and
  warm-start `model_init` paths.
- The hardcoded `V = 1` returned by both mfsusieR hooks does
  not break any susieR consumer that reads `model$V[l]`
  (e.g., `susie_get_pip` filtering, `trim_null_effects`).
- No mfsusieR helper duplicates a susieR helper that should
  be delegated.

Report at
`inst/notes/sessions/<date>-cross-package-audit-a-susier-posthooks.md`.

#### Agent B — mfsusieR vs fsusieR

Focus: single-modality functional case (M=1, T_1>1). Trace
each numerical pathway and report any residual divergence
not in the known-ledger.

- Data prep: `mf_dwt` / `mf_quantile_normalize` /
  `mf_low_count_indices` vs `DWT2 + colScale + (optional)
  qnorm`. Account for our `wavelet_standardize = TRUE`
  default vs fsusieR's implicit double-`colScale`.
- Prior init: `init_scale_mixture_prior_default` (deterministic
  10th-percentile-Shat anchor) vs `init_prior.default`
  (ash-driven). Differences in grid construction and
  cold-start pi already noted; flag any further drift.
- SER step: mfsusieR's
  `compute_residuals + compute_ser_statistics + loglik` vs
  fsusieR's `cal_partial_resid + cal_Bhat_Shat + log_BF`.
- M-step: `mf_em_m_step_per_scale` vs `m_step.lik_mixture_normal`,
  including the `nullweight * idx_size` weighting and the
  mixsqp control args.
- **Inner-EM control surface**: `max_inner_em_steps` and the
  early-exit on `|d(model$lbf[l])| < params$tol` vs fsusieR's
  `max_step_EM` and `abs(newloglik-oldloglik) >= espsilon`.
  Confirm semantics match (both treat the count as a CAP on
  inner cycles, both exit early on per-effect lbf
  convergence). Flag any subtle off-by-one.
- Variant thinning: mfsusieR's `alpha_thin_eps` vs fsusieR's
  `max_SNP_EM`. Both cap mixsqp's L-matrix row count, by
  threshold vs rank. Confirm neither drops the truly-signal
  variants.
- Residual variance update; ELBO; smoothers (TI, scalewise,
  smash, HMM); CS construction; PIP filter; refinement.

Report at `<date>-cross-package-audit-b-fsusier-posthooks.md`.

#### Agent C — mfsusieR vs mvf.susie.alpha

Focus: multi-modality functional case (M>1). Same pathways
as B with the cross-modality combiner added:

- `combine_outcome_lbfs.cross_outcome_prior_independent`
  vs the `f_logBF + u_logBF` aggregation in `multfsusie`'s
  EM loop.
- The per-(outcome, scale) sigma2 aggregation
  (`update_variance_components.mf_individual`) vs mvf's
  scope.
- The per-effect variational posterior shape: mfsusieR keeps
  `model$mu[[l]][[m]]` as a `p x T_basis[m]` matrix per
  outcome m; mvf has a slightly different layout in
  `multfsusie.obj`. Flag any algorithmically meaningful
  difference (not just naming).
- The inner-EM loop interaction with the multi-outcome path
  (mfsusieR's loglik aggregates per-outcome lbf into
  `model$lbf[l]` via the cross-outcome combiner; the
  inner-EM check uses that aggregated value).

Report at `<date>-cross-package-audit-c-mvf-posthooks.md`.

### 2. Triage (Phase B)

Findings classified into:

- `accept-known`: matches a divergence already in the
  divergence ledger or refactor-exceptions. Cite the
  reference and move on.
- `add-to-divergence-ledger`: real and intentional but NOT
  yet in the ledger. Add an entry to
  `~/Documents/mfsusier-mvf-divergences/divergences.tex`
  and a refactor-exceptions stanza in the same session.
- `fix-now`: clear bug. Open a targeted OpenSpec change in
  the same session and link it from this proposal.
- `track-for-later`: real, bounded, not blocking release.
  Open a follow-up OpenSpec with a deferred task list.

### 3. Audit summary (Phase C)

`inst/notes/cross-package-audit-summary-posthooks.md`
collates the three reports plus the triage table plus the
list of follow-up OpenSpec changes opened.

## Impact

- **Severity**: HIGH. Halt feature work until Phase B
  triage completes.
- **Source-code changes**: zero from this proposal directly;
  drives zero or more follow-up changes.
- **Documentation**: three session-note reports + a summary +
  any new divergence-ledger / refactor-exceptions entries.
- **Process precedent**: every future round of substantial
  algorithmic change in `R/individual_data_methods.R`,
  `R/model_class.R`, `R/ibss_methods.R`, or
  `R/em_helpers.R` triggers a cross-package audit before
  the round is archived.

## Out of scope

- Vignette text comparison. Numeric / code only.
- Performance comparisons. Separate concern.
- Smoother-internal numerical comparison beyond a smoke
  check that the smoothers produce sane output.
