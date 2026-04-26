# Cross-package audit summary

Date: 2026-04-26
Inputs: three independent fresh-context opus-model audits, read-only.

- `sessions/2026-04-26-cross-package-audit-a.md` — vs susieR (delegation + divergence)
- `sessions/2026-04-26-cross-package-audit-b.md` — vs fsusieR (M = 1 functional)
- `sessions/2026-04-26-cross-package-audit-c.md` — vs mvf.susie.alpha (multi-modality functional)

## Headline

**No undocumented numerical correctness divergences.** Every arithmetic difference is either a shape generalization (per-(scale, outcome) instead of scalar), a different algorithm (mixture-of-normals vs single-Gaussian) inherent to mfsusieR's design, or a Pattern A / Pattern B port-source bug fix already in the `inst/notes/refactor-exceptions.md` ledger.

The findings split into:

- **5 `fix-now`** items: 3 ledger gaps + 2 small code cleanups. None blocks Phase 4 fidelity.
- **~15 `track-for-later`** items: documentation gaps, cosmetic delegation misses, missing end-to-end coverage tests.
- **~25 `accept-as-documented`** items.

## `fix-now` items (5)

| ID | Source | Description | Type |
|---|---|---|---|
| B-1.1 | Audit B | `mf_invert_variance_curve` is a Pattern A correction of upstream's `(invert_dwt(sqrt(var_w)))^2` formula but lacks a `refactor-exceptions.md` entry. | Ledger entry |
| B-4.3 | Audit B | `mf_residualize_wavelet_eb` defaults (`wavelet_filter_number = 1L, wavelet_family = "DaubExPhase"`) disagree with `mf_adjust_for_covariates` wrapper defaults (10, "DaubLeAsymm"). Wrapper passes its values through, so impact is only on direct calls to the inner. | Default alignment |
| C-1.3 | Audit C | mvf.susie.alpha scales `nullweight` by `K_f + K_u` in the M step (`R/EM.R:58-65, 102`); mfsusieR passes `mixsqp_null_penalty` unscaled. Upstream's inline comment claims this addresses FDR miscalibration under M > 1. Investigate whether mfsusieR's per-(outcome, scale) M step needs the same scaling. | FDR investigation |
| C-4.4 | Audit C | `get_objective.mf_individual` omits the `o = sum(alpha_l * log(p / pmax(alpha_l, 1e-6)))` entropy term that mvf.susie.alpha adds in `ELBO_mutlfsusie.R:269-282`. mfsusieR follows susieR's clean ELBO decomposition; needs ledger entry documenting the deliberate omission. | Ledger entry |
| C-5.3 | Audit C | `R/post_smooth_smash.R:108` uses local `mf_per_position_bhat_shat` instead of `susieR::compute_marginal_bhat_shat`. The HMM branch (`R/post_smooth_hmm.R:286`) correctly uses the susieR helper. | Delegation cleanup |

## `track-for-later` items (high-signal subset)

Architecture / behaviour gaps:

- **A-F2**: `check_convergence.mf_individual` silently disables susieR's PIP-convergence and verbose tabular output (medium severity).
- **A-F5**: `ibss_initialize.mf_individual` skips `slot_prior` init silently.
- **A-F13**: `compute_Xb` / `compute_Xty` not cached in `zzz.R`, so `compute_residuals` and `update_fitted_values` use plain `crossprod` / `%*%` instead of susieR's sparse-aware helpers (forecloses sparse-X support).
- **B-3.1, B-3.2**: `.post_smooth_smash` and `.post_smooth_ti` lack the upstream `n_iter = 5` outer Gauss-Seidel loop. Affects L > 1 fits only; current tests use L = 1 and pass at `1e-12`. Worth flagging before users invoke smash / TI on multi-effect fits.
- **C-1.3, C-1.4**: documentation gaps for `max_SNP_EM` truncation drop and inner-EM-loop collapse.

Cosmetic delegation misses:

- **A-F1, A-F4**: `validate_prior.mf_individual` and `configure_data.mf_individual` are byte-equivalent to defaults.
- **A-F3**: `cleanup_model.mf_individual` re-implements field-stripping rather than `NextMethod()`.
- **A-F8 / C-5.3 / C-5.4**: smash kernel uses local `mf_per_position_bhat_shat` and `mf_col_scale` duplicates of helpers that already exist in the package (and in susieR for the Bhat/Shat case).

Stale ledger entries:

- **A-F11**: `lowc_wc filtering deferred` is stale; `low_count_filter` is now a public argument and the prior init excludes flagged columns.
- **A-F12**: `fit_hmm` HMM-mu-subset entry's "Resolution (2026-04-26)" suggests upstream caught up but the entry-status header is ambiguous.

Documentation:

- **A-F7 / C-4.5**: `update_model_variance.mf_individual` roxygen truncated.
- **A-F10**: `summary.mfsusie` schema diverges from susieR (Phase 8 polish channel).

Test coverage:

- **B-5.1**: `test_adjust_covariates.R` uses graduated tolerances `< 0.1`, `< 0.05` that conflict with the binary apple/orange philosophy in `refactor-discipline.md` section 3. Justified inline by Pattern A bug-magnitude framing; rewrite as structural assertions.
- **B-5.2 / C-6.2**: no top-level end-to-end fidelity test against `fsusieR::susiF` or `mvf.susie.alpha::multfsusie` confirmed (per-component tests cover the kernels; the full IBSS-step bit-identity is blocked by the documented sigma2 divergence).

## Pattern A / B exceptions verified load-bearing

All 24 ledger entries cross-checked. Two are stale (A-F11, A-F12). The remaining 22 are correct and current:

- `EBmvFR-MLE_wc2`, `EBmvFR-ER2` (Pattern A; bug fixes against fsusieR EBmvFR)
- `multfsusie-ER2` (Pattern A; bug fix against mvf.susie.alpha)
- `HMM-mu-subset` (resolved 2026-04-26 by fsusieR PR #31; ledger header should be updated)
- All wavelet pipeline simplifications (Pattern B, `1e-12` ULP cleanup)
- All renames (`gridmult` -> `grid_multiplier`, etc.)
- All scope decisions (out-of-scope-EBmvFR, etc.)

## Triage rationale

No `fix-now` item changes numerical correctness on the supported paths. Three are documentation-only (B-1.1, C-1.3 investigation, C-4.4); two are small code cleanups (B-4.3 default alignment, C-5.3 delegation). All five are bounded edits suitable for a single follow-up OpenSpec change.

The proposed follow-up: open one consolidated change `apply-cross-package-audit-fixnow` covering the 5 items, plus the two ledger updates (A-F11 stale entry, A-F12 status header). Estimated scope: ~50 lines of code, ~3 new refactor-exceptions entries. Track-for-later items are queued for Phase 8 polish or specific follow-up changes (sparse-X support, slot_prior interop, smoother outer-loop n_iter, schema harmonization).

## Decisions

- **Audit run**: complete and documented.
- **Fix-now items**: open a single follow-up OpenSpec change `apply-cross-package-audit-fixnow` to land them. Out of scope for the audit proposal itself per its "Out of scope" section.
- **Track-for-later items**: triaged into the audit notes. To be picked up individually as their own OpenSpec changes when prioritized.
- **Refactor-exceptions ledger updates** (stale entries A-F11, A-F12) bundle into the same follow-up change.
