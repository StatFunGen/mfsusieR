# Cross-package audit Agent C — mfsusieR vs mvf.susie.alpha (2026-04-30)

## Summary
- Total findings: 9 (by severity: fix-now=1, track-for-later=4, needs-trace=2, accept-pre-existing=2)
- Pathways audited: data prep / prior init / SER step + cross-modality combiner / M-step / inner-EM / variant thinning / sigma2 aggregation / ELBO / posterior shape / cleanup / smoothers+CS / small-sample correction
- Agent confidence: data prep — high; prior init — high; SER + combiner — high; M-step — high; inner-EM — high; variant thinning — medium; sigma2 — high; ELBO — high; posterior shape — high; cleanup — medium; smoothers/CS — medium; small-sample — high

## Findings

### Finding C-1: refactor-exceptions ledger entry is stale on parameter name + default
- **Severity**: fix-now (ledger correctness)
- **Pathway**: M-step on multi-modality (M-fold scaling)
- **Location**: `mfsusieR/R/individual_data_methods.R:773-774`; `mfsusieR/inst/notes/refactor-exceptions.md:380-382`; `mfsusieR/R/mfsusie.R:155-164` (docstring); `mvf.susie.alpha/R/EM.R:64-65, 102` (mvf reference).
- **Verified by reading**: yes
- **Description**: The 2026-04-26 fix-now ledger entry codifies the M-fold scaling as `mixsqp_null_penalty <- (params$mixsqp_null_penalty %||% 0.7) * max(1L, data$M)`. The shipped code uses `mixture_null_weight <- (params$mixture_null_weight %||% 0.05) * max(1L, data$M)` — a different parameter name (`mixsqp_null_penalty` → `mixture_null_weight`) and a different default (0.7 → 0.05). The M-fold multiplier is intact. mvf's default is `nullweight = 0.7` (`R/multfsusie.R:184`); a 14× drop relative to upstream is a substantive numerical choice that warrants either (a) preserving the 0.7 default, or (b) updating the ledger to record the deliberate 0.05 baseline + scaling. Currently the ledger's reference value mismatches the code.
- **Triage hint**: add-to-ledger
- **Suggested action**: Update the ledger entry's parameter name to `mixture_null_weight` and the default value to `0.05`. If the 0.05 default is itself a deliberate departure from mvf's `nullweight = 0.7`, add a one-paragraph derivation pointer (otherwise the read of `0.05` vs `0.7` reads as an undocumented mfsusieR-specific tuning).

### Finding C-2: variant thinning mechanism is fundamentally different from mvf
- **Severity**: track-for-later
- **Pathway**: variant thinning (M-step input)
- **Location**: `mfsusieR/R/individual_data_methods.R:587-591` (alpha-threshold); `mvf.susie.alpha/R/EM.R:67-71` (top-K).
- **Verified by reading**: yes
- **Description**: mvf's `EM_pi_multsusie` truncates the M-step input to the top `max_SNP_EM = 100` SNPs by `lBF` (`R/EM.R:67-71`); under heavy LD, the cutoff varies per IBSS iter and modality. mfsusieR uses `keep_idx <- which(zeta_l > alpha_eps)` with `alpha_thin_eps = 5e-5` (default), which is alpha-driven rather than lBF-driven and is variable in count. For dense priors with diffuse alpha (early IBSS iterations), mfsusieR keeps far more variables than mvf's hard top-100; in the saturated-signal regime (alpha concentrated on one SNP), mfsusieR keeps fewer. The previous Audit C report (2026-04-26 §1.2) already track-flagged this; not yet ledger'd.
- **Triage hint**: add-to-ledger
- **Suggested action**: ledger entry naming the alpha-threshold vs top-K design choice and asserting the M-step output is robust to either thinning mechanism.

### Finding C-3: model-state scratchpads survive cleanup onto the returned fit
- **Severity**: track-for-later
- **Pathway**: cross-modality cleanup and finalization
- **Location**: `mfsusieR/R/individual_data_methods.R:1001-1003` (`cleanup_extra_fields.mf_individual` returns `character(0L)`); `mfsusieR/R/individual_data_methods.R:235-237` (`ser_cache`); `mfsusieR/R/individual_data_methods.R:499-538` (`refresh_iter_cache.mf_individual` writes `iter_cache`); `susieR/R/generic_methods.R:430-445` (`cleanup_model.default` strip list).
- **Verified by reading**: yes
- **Description**: The current `cleanup_extra_fields.mf_individual` returns `character(0L)` (commit 605f395 changed it from `c("V")` so model$V can carry the per-effect effective slab variance). Two large scratchpads survive cleanup onto the returned fit: `model$ser_cache` (a list with `betahat[[m]]` and `shat2[[m]]` of size `p × T_basis[m]` per outcome) and `model$iter_cache` (with `sdmat[[m]][[s]]` and `log_sdmat[[m]][[s]]` of size `(|keep_idx| × |idx|) × K` per (m, s) for mixsqp prior classes). Neither is consumed by post-fit accessors (`predict`, `coef`, `summary`, `fitted`), only by the IBSS hot path. On the heavy fixture (n=84, p=3500, M=6, T=128, see `2026-04-28-audit-findings.md`) these two slots can be tens of MB. The 2026-04-28 audit notes do not mention this specific cleanup gap.
- **Triage hint**: add-to-ledger
- **Suggested action**: extend `cleanup_extra_fields.mf_individual` to `c("iter_cache", "ser_cache")` (and verify `expand_model_init_to_L` does not need them on warm-restart).

### Finding C-4: model-state-only V cleanup vs susieR prior_tol semantics
- **Severity**: needs-trace
- **Pathway**: cross-modality cleanup and finalization
- **Location**: `mfsusieR/R/individual_data_methods.R:725, 734-748` (`.effective_V_l` writes mean-over-(m, s) of `sum_k pi[l,m,s,k] * var_k`); `mfsusieR/R/mfsusie.R:358` (`prior_tol = 1e-9`); `susieR/R/susie_get_functions.R:486-498` (`susie_get_pip` filters effects with `res$V <= prior_tol`).
- **Verified by reading**: yes
- **Description**: The post-loglik hook now writes `model$V[l] = .effective_V_l(...)` (commit 605f395). `susie_get_pip` in `ibss_finalize` uses `prior_tol = 1e-9` to drop effects with `model$V <= prior_tol`. The new effective slab variance can be very small in practice (mean over (m, s) of `sum_k pi_k * var_k` collapses toward zero when the per-effect prior concentrates on the null component). Under `prior_variance_scope = "per_outcome"` with a single `mixture_normal` group and an iterated null-collapsed prior, `pi_null * 0 + pi_nonnull[K] * sd[K]^2` can fall below `1e-9` for very weak signals, silently dropping the effect from `pip`. mvf has no analog: it filters CSes by purity / lbf_min in `out_prep`, not by V. I cannot verify whether this is empirically reachable without running fits; flagging for trace.
- **Triage hint**: needs-trace
- **Suggested action**: write a regression test that fits a synthetic null-only fixture (no causal signal) and asserts `model$pip` matches what `susie_get_pip(model, prior_tol = 0)` returns, or relaxes mfsusieR's `params$prior_tol` to 0 by default.

### Finding C-5: pre-hook reads model$V[l] but discards it; V channel dual-encodes
- **Severity**: track-for-later
- **Pathway**: per-effect SER step (V state contract)
- **Location**: `mfsusieR/R/individual_data_methods.R:677-689` (pre-hook returns `list(V = 1, model = model)`); `susieR/R/single_effect_regression.R:37, 40-43` (V flow into pre-hook); `mfsusieR/R/individual_data_methods.R:725` (post-hook returns effective V).
- **Verified by reading**: yes
- **Description**: After commit 605f395, model$V[l] dual-encodes two semantically different things at different points in the IBSS step: V=1 from the prior pre-hook (consumed by `loglik` as the V-scale multiplier on the mixture sd grid), and V=`.effective_V_l(...)` from the post-hook (consumed by `set_prior_variance_l` and lives until the next pre-hook overwrites it). The pre-hook's `list(V = 1, model = model)` (line 688) silently drops the V_init it received — correct, because V=1 is the load-bearing semantic for `loglik.mf_individual`. The dual-encoding works but is non-obvious; the pre-hook docstring would benefit from an explicit "V_init is intentionally discarded; effective slab variance lives only on the model object between IBSS iters" note. Not a numerical bug.
- **Triage hint**: track-for-later
- **Suggested action**: docstring note on `pre_loglik_prior_hook.mf_individual` explaining the V-discard.

### Finding C-6: inner-EM tol semantics differ from mvf's tol
- **Severity**: track-for-later
- **Pathway**: inner-EM control surface
- **Location**: `mfsusieR/R/individual_data_methods.R:701-723` (inner-EM loop); `mvf.susie.alpha/R/EM.R:33, 81-120` (mvf inner EM with `espsilon = 0.0001`, `max_step = 1` default).
- **Verified by reading**: yes
- **Description**: mfsusieR's inner-EM loop reads `inner_tol <- params$tol %||% 1e-4`, sharing the IBSS outer-convergence tol with the inner per-effect lbf-change cutoff. mvf's `EM_pi_multsusie` has its own `espsilon = 0.0001` (`R/EM.R:33`) and `max_step = 1` default for the inner EM (only one M-step pass per IBSS effect). The inner-EM loop is mfsusieR-only conceptually — mvf's `max_step = 1` makes it a single M-step. mfsusieR's `max_inner_em_steps = 5L` default plus the always-on first M-step gives up to 6 inner cycles. Higher fidelity to (alpha, pi) in lockstep, more compute. No numerical contract violation, but the semantics of `tol` differ from mvf (where `tol` is purely the outer ELBO-change cutoff).
- **Triage hint**: add-to-ledger
- **Suggested action**: ledger entry recording that mfsusieR's `tol` controls both inner-EM and outer convergence; if a separate `inner_em_tol` is desired, expose it on the public surface.

### Finding C-7: cross-outcome combiner's class hierarchy registration incomplete
- **Severity**: needs-trace
- **Pathway**: per-effect SER step (cross-modality combiner extension point)
- **Location**: `mfsusieR/R/prior_cross_outcome.R:50-56` (combine method has only one class registered).
- **Verified by reading**: yes
- **Description**: `cross_outcome_prior_independent()` constructs an object of class `c("mf_prior_cross_outcome_independent", "mf_prior_cross_outcome")`. `combine_outcome_lbfs` is a generic with one `.mf_prior_cross_outcome_independent` method registered via `@exportS3Method`. There is no `combine_outcome_lbfs.mf_prior_cross_outcome` (parent class) method or `combine_outcome_lbfs.default` method. If a user creates a class that inherits from `mf_prior_cross_outcome` but not from `_independent`, dispatch will fall through to `UseMethod` with no method and error. The seam exists for future combiners (per the file-header comment), but the protocol for "register a new combiner" is not defended by a default. Not a bug today (only the independent combiner is shipped); a tripping hazard for the seam's stated purpose. The 2026-04-28 audit confirms this is the intended extension surface.
- **Triage hint**: needs-trace
- **Suggested action**: add a `combine_outcome_lbfs.default` method that errors with a helpful message naming the registered combiner classes, or document in the seam header how downstream code should subclass.

### Finding C-8: mu cache rebuild on warm-start uses zero-state rather than uniform alpha
- **Severity**: track-for-later
- **Pathway**: per-effect variational posterior shape (warm-start rehydration)
- **Location**: `mfsusieR/R/ibss_methods.R:96-119` (`expand_model_init_to_L`).
- **Verified by reading**: yes
- **Description**: On L_greedy expansion, `expand_model_init_to_L` appends `L_diff` extra effects with `alpha = 1/p` uniform but `mu` and `mu2` set to zero matrices. The first SER step on the new effect therefore sees a partial residual `R_l = D_m - sum_{l' < L_prev} X * alpha_l' * mu_l'` (correct, because the new effects contribute nothing at zero state), but the `loglik.mf_individual` is then driven against the same per-effect mixture prior shape as the existing effects (the seed_pi from `mi$pi_V[[1L]]`). For a fixed-warm-start path (no L expansion) this is fine; for L_greedy (the default `L_greedy = 5`) the new effects' M-step inputs are influenced by `seed_pi` from effect 1 — which may have already concentrated on the null. The expansion seeds new effects' priors with the lowest-effect prior, not a generic init prior. Not algorithmically wrong (it is a defensible warm-start convention), but mvf has no analog (mvf's `multfsusie.obj` resume fully populates per-effect state from the prior call, never grows L mid-run — see refactor-exceptions PR group `R/multfsusie.R::multfsusie.obj` entry).
- **Triage hint**: needs-trace
- **Suggested action**: regression test on L_greedy growth from 5 to 10, asserting the expansion-warm-start fits within `1e-8` of a cold-start fit at L=10.

### Finding C-9: lbf_variable_outcome stores in a 3D array; mvf stores list-of-list
- **Severity**: accept-pre-existing
- **Pathway**: per-effect variational posterior shape
- **Location**: `mfsusieR/R/model_class.R:118` (allocates `array(NA_real_, dim = c(L, p, M))`); `mfsusieR/R/individual_data_methods.R:240-246` (writes); `mvf.susie.alpha/R/operation_on_multfsusie_obj.R` (`update_lBF.multfsusie` stores `multfsusie.obj$lBF_per_trait[[l]]` as a list with `f_logBF` and `u_logBF` matrices).
- **Verified by reading**: yes
- **Description**: mfsusieR collapses the f-trait/u-trait split of mvf into a unified M-axis via the 3D array `lbf_variable_outcome`. mvf's storage is `list[L]` of `list(f_logBF, u_logBF)`, where each entry is a `K_f × p` or `K_u × p` matrix. Algorithmically the layouts encode the same information (per-(effect, variable, outcome) log-BF), with mfsusieR using the unified-modality convention from CLAUDE.md hard rule #2. Not a divergence: the layout choice is dictated by the unified `mfsusie` design.
- **Triage hint**: accept-known
- **Suggested action**: none.

## Pathways with no findings

- **Multi-modality data prep (`mf_dwt` + `create_mf_individual`)**: silent because-clean. Ledger PR group 2 (utils_wavelet_transform.R DWT2, pack_dwt, multfsusie.R inline DWT pipeline) covers the structural reorg; the iteration order, slot naming, and step list match. `Y_remapped`, `D`, `scale_index`, `T_basis`, `pos_grid`, `column_center`, `column_scale`, `wavelet_center`, `wavelet_scale`, `lowc_idx`, `na_idx`, `xtx_diag_list` all match the per-modality counterparts in mvf's `multfsusie()` body.
- **Per-modality prior init (`mf_prior_scale_mixture`)**: silent because-clean. The per-(m, s) groups, the V_grid construction, and the pi initialization all match the equivalence asserted in the existing audit (`2026-04-26-cross-package-audit-c.md` §3). The new `per_scale_normal` / `per_scale_laplace` paths in `init_ebnm_prior_per_scale` are mfsusieR-only and slot cleanly into the per-modality loop (verified by reading the `groups_m` plumbing).
- **Per-effect SER step + cross-modality combiner (`loglik.mf_individual`)**: silent because-clean. Per-modality loop sums per-(m, s) BFs across scales; `combine_outcome_lbfs.mf_prior_cross_outcome_independent` reduces with `+`; the equivalence with mvf's `apply(lBF_per_trait$f_logBF, 2, sum) + apply(lBF_per_trait$u_logBF, 2, sum)` is asserted at `1e-15` in `tests/testthat/test_mvf_alpha_fidelity.R:80-92`. The inner-EM convergence check (line 723) reads `model$lbf[l]` which is the cross-modality aggregated `lbf_model` (line 226).
- **M-step on multi-modality (`.opv_mixsqp` and `.opv_ebnm_point`)**: silent because-clean (modulo C-1 ledger correction). Both mfsusieR's per-(m, s) M-step and mvf's `m_step_multsusie` -> `fsusieR::m_step.lik_mixture_normal_per_scale` -> `scale_m_step` end up running one mixsqp solve per (m, s) using the same `zeta` (the per-effect alpha). The M-fold scaling (`max(1L, data$M)`) is preserved.
- **Per-(outcome, scale) sigma2 aggregation (`update_variance_components.mf_individual`)**: silent because-clean. The per-scale aggregation matches the formula in the ledger (PR group 6, `R/individual_data_methods.R::update_variance_components.mf_individual`), and the per-outcome path collapses to mvf's legacy divisor `n * T_basis[m]`.
- **ELBO on multi-modality (`Eloglik.mf_individual`, `get_objective.mfsusie`)**: silent because-clean. Per-modality contributions are aggregated correctly; the `o` entropy term from mvf is omitted per the ledger entry (Divergence 2). The `refresh_lbf_kl.mf_individual` path makes the ELBO a coherent variational free energy.
- **Smoothers and CS in multi-modality**: silent because-not-checked end-to-end against mvf, but the credible-band formula divergence is already documented in the 2026-04-26 ledger entry (`mf_invert_variance_curve` replacing `abs(invert_dwt(sqrt(var_w)))`). The `.post_smooth_*` methods iterate per modality; CS construction goes through `susie_get_cs(model, X, ...)` which is shared with susieR (no per-modality CS construction in mvf either).
- **Small-sample correction (`mixture_log_bf_per_scale_johnson`)**: silent because-clean. The kernel matches `fsusieR::log_BF(..., df = n - 1)` at `1e-14` per `tests/testthat/test_small_sample_correction.R:115-116`, well below the ledger's stated `1e-12` floor. Posterior moments are intentionally unchanged (the correction acts on selection only).

## Limits of this audit

- Did not run any IBSS fits; all findings are static-analysis-derived. C-3, C-4, C-7, C-8 in particular would benefit from end-to-end runs to confirm the failure modes are reachable.
- Did not audit fsusieR's `EBmvFR` paths or mfsusieR's `R/adjust_covariates.R` (out of scope per refactor-exceptions PR group 4 and the EBmvFR file-level entries).
- Did not re-trace credible-band SD math (already covered by the 2026-04-26 audit fix-now bundle and audit C-1.1 / C-1.2).
- Did not exhaustively read `R/post_smooth_smash.R` and `R/post_smooth_hmm.R` — they were covered in the 2026-04-26 audit C §5 (HMM-mu-subset bit-identity and smash duplicate-helper findings).
- Did not run `tests/testthat/` to verify the assertions cited; relied on test-file reading.
- The `lbf_variable_outcome` field is L×p×M — for very large M this could be a large array on the fit, but that's a memory-cost call, not a numerical concern.
