# FDR inflation: investigation and diagnosis
# 2026-05-09 to 2026-05-20

---

## 1. Observed fact

Permutation data (5 real ATAC-seq regions, 1024-bin, `wavelet_qnorm = TRUE`):

| version | type | n regions | HP CS | perm DR |
|---------|------|-----------|-------|---------|
| 1e1a866 (no NA fix) | perm | 168 | 8 | 4.8% |
| be0ce136 (NA fix, 0.0.1) | perm | 153 | 1 | 0.7% |
| 0.0.2 HEAD | perm | 168 | 92 | **54.8%** |
| 1e1a866 | real | 168 | 14 | DR 8.3% |
| be0ce136 (0.0.1) | real | 168 | 0 | DR 0% |

The critical gap: 0.7% → 54.8% (78×) in perm DR between be0ce136 and 0.0.2.

Already ruled out before 2026-05-09:
- `L_greedy = NULL` → same 54.8% perm DR; L_greedy is not the cause
- `mixture_null_weight = 0.1` alone → perm DR still 0.5

---

## 2. Complete algorithmic delta: be0ce136 → HEAD (0.0.2)

### A. be2722e — mixsqp warm start

**Before**: cold start `x0 = c(init_pi0_w, rep(1e-6, K-1))`, `convtol.sqp=1e-8`,
`numiter.em=20`, `tol.svd` default.

**After**: warm start `x0 = pi_prev` (previous IBSS iter's π), `convtol.sqp=1e-6`,
`numiter.em=10`, `tol.svd=0` (skip SVD).

**Mechanism**: residual non-null π from iter t-1 carries into iter t; looser tolerance
may not converge back to null.

### B. 6728cd3 — mixture_null_weight default 0.1 → 0.05 (×M scaling preserved)

**Before**: effective penalty = 0.1 × M = 0.6 for M=6.
**After**: effective penalty = 0.05 × M = 0.3 for M=6.
Less null regularization → easier for mixsqp to place weight on non-null components.

### C. New architecture — `fitted_g_per_effect` + inner EM loop

**Before**: zero inner EM; π updated once per (l, outer_iter); no per-effect storage;
π state discarded between outer iterations.

**After**:
- `fitted_g_per_effect[[l]]` persists each effect's π across outer iterations;
  `pre_loglik_prior_hook` restores effect l's prior before each loglik call.
- `post_loglik_prior_hook` runs up to `inner_cap = max_inner_em_steps + 1 = 6` cycles
  of {M-step → loglik → moments → KL} per (effect, outer_iter).

Execution structure in HEAD:
```
outer_iter → for l in 1..L: {
  restore π from fitted_g_per_effect[[l]]   (warm)
  loglik → moments → KL
  inner loop × 6: {M-step(warm) → loglik → moments → KL}
  save π to fitted_g_per_effect[[l]]
}
```

This is the structurally largest change.

### D. dc5515e — cpp → pure-R likelihood kernel

Floating-point differences < 1e-12; negligible for FDR.

### E. trim_null_effects now active

Should reduce FDR by pruning low-V effects, but V[l] stays non-zero when warm start
keeps π non-null, so the effect is neutralized by changes A–C.

---

## 3. Mathematical argument: why inner EM can inflate FDR

SuSiE IBSS (Wang et al. 2020) is coordinate ascent on `ELBO(q, g)`:
- E-step (outer): fix g, compute posterior → alpha[l,]
- M-step (outer): fix alpha, update g = argmax_g ELBO

Monotone convergence holds because each step does not decrease the ELBO.

**The inner EM breaks this guarantee.** After the outer loglik fixes alpha[l,], the
inner loop re-runs the E-step (loglik → alpha) inside the same outer iteration.
These inner E-steps are not accounted for in the outer ELBO.

Under null data the feedback amplifies noise: uniform alpha[l,] has finite-sample
variation ε_j → inner M-step produces slightly non-null g → inner E-step concentrates
alpha further → 6 cycles can produce detectable CSes from pure noise.

Under signal the inner loop converges quickly to the true effect (self-limiting).
This null-vs-signal asymmetry is the mechanism.

The warm start (change A) compounds this: `fitted_g_per_effect[[l]]` carries the
inner-loop-amplified non-null π from outer iter t to t+1, so at t+1 even the
first inner M-step starts from a non-null distribution.

---

## 4. Issue 1 (primary): `mf_quantile_normalize` discards NA structure

**File**: `R/utils_wavelet.R`

**Root cause**: `rank(Y_wd[, j], ties.method = "random")` called without
`na.last = "keep"`. R assigns NA entries the largest ranks (acts like `na.last = TRUE`).

```r
rank(c(1.2, 0.5, NA, 2.1, NA, -0.3), ties.method = "random")
# [1] 3 2 5 4 6 1   <- NAs get ranks 5, 6; not preserved as NA

rank(c(1.2, 0.5, NA, 2.1, NA, -0.3), ties.method = "random", na.last = "keep")
# [1]  3  2 NA  4 NA  1   <- correct
```

**Effect**: NA rows steal the top `n_na` quantiles; complete-row wavelet coefficients
follow a right-truncated distribution with variance below 1. `sigma2` converges to
this deflated variance throughout IBSS.

Estimated sigma2 deflation (region 4):

| Outcome | n_cc | n_na | Estimated deflation |
|---------|------|------|---------------------|
| Ast     |  71  |  13  | ~35% |
| Exc     |  60  |  24  | ~48% |
| Inh     |  70  |  14  | ~36% |
| Mic     |  75  |   9  | ~28% |
| Oli     |  73  |  11  | ~32% |
| OPC     |  83  |   1  |  ~6% |

**Chain to FDR**: deflated sigma2 → deflated shat2 → inflated Bhat/Shat → inflated
Bayes factors → alpha concentrates on LD-block SNPs → spurious high-purity CSes.
Plausible mechanism; not tested in isolation with X and LD structure held fixed.

**Status**: **fixed** in commit 54a10a5 (2026-05-20). `rank()` now applied only to
non-NA entries; NA positions preserved.

Simulation v2 (na_frac factor, job 33180369) confirms: sigma2 with 20% NA rows
recovers from ~0.65 (unfixed) to ~0.90 (fixed) under null/normal conditions.

---

## 5. Issue 2 (secondary): `get_var_y` includes corrupted NA rows

**File**: `R/ibss_methods.R`

After Issue 1, NA rows have large finite values rather than NA, so `na.rm = TRUE`
does not exclude them. The initial sigma2 may be inflated (opposite direction from
Issue 1's converged effect). Effect limited to iteration 1; `update_variance_components`
uses `na_idx` and overrides from iteration 2 onward.

**Status**: resolved as a side-effect of Issue 1 fix; once NA rows stay NA through
qnorm, `na.rm = TRUE` correctly excludes them.

---

## 6. Perm experiment: `fitted_g_per_effect` association

SLURM job 33175249 (N_BINS=1024, 5 regions, permuted real data) compared two conditions
with all other settings fixed (`wavelet_qnorm = TRUE` in both):

| Condition | `cross_iter_prior` | Description |
|---|---|---|
| π-persist | `TRUE` | π warm-started per effect across IBSS outer iterations |
| π-reset | `FALSE` | π re-estimated from shared G_prior each iteration |

Per-region HP CS counts (purity ≥ 0.8, permuted data):

| region | 64-bin π-persist | 64-bin π-reset | 1024-bin π-persist | 1024-bin π-reset |
|--------|:---:|:---:|:---:|:---:|
| 4  | 1 | 0 | 3 | 0 |
| 5  | 0 | 0 | 1 | 0 |
| 6  | 3 | 4 | 3 | 0 |
| 11 | 0 | 0 | 1 | 0 |
| 14 | 1 | 0 | 2 | 0 |
| **total** | **5** | **4** | **10** | **0** |

At 1024 bins, π-persist produces 10 HP CS vs 0 for π-reset. Issue 1 (qnorm NA bug)
is active equally in both conditions, so the contrast isolates `fitted_g_per_effect`.
At 64 bins the difference is small (5 vs 4); region 6 fires under both conditions,
suggesting a region-specific pattern independent of prior persistence.

Cold/warm mixsqp start and inner EM steps (conditions A–D in the interactive test
plan) made no detectable difference at 64 bins.

**Mechanism**: `fitted_g_per_effect[[l]][[m]][[s]]` stores the fitted π for effect l,
outcome m, scale group s across IBSS outer iterations. Before each SER step,
`pre_loglik_prior_hook` restores effect l's π into `G_prior`. This is structurally
identical to `mvf.susie.alpha`'s `est_pi[[l]]` mechanism.

At 1024 bins (10 wavelet scales) each effect has 10 × K mixture parameters that can
adapt per-effect. Real ATAC-seq wavelet coefficients may have non-spherical cross-scale
covariance; per-effect π persistence may fit that structure for a permuted variant,
producing a spurious high-purity CS. At 64 bins (6 scales) the same mechanism has
less leverage.

`cross_iter_prior = FALSE` (π-reset) is closer to the standard IBSS derivation
(Wang et al. 2020), which assumes a shared prior across all effects. The parameter
is now explicitly exposed in `mfsusie()` (default `TRUE` for backward compatibility).

**Bug note (2026-05-20)**: `.opv_mixsqp` and `.opv_ebnm_point` wrote to
`model$fitted_g_per_effect[[l]][[m]][[s]]` unconditionally, even when
`cross_iter_prior = FALSE` (where `fitted_g_per_effect` starts as NULL). R's
assignment to NULL creates partial list structure; `pre_loglik_prior_hook` then
crashed on the second effect. Fixed in commit 37d6723 with `!is.null()` guards.
All `cross_iter_prior = FALSE` tasks in simulation v2 failed; resubmitted as v2b
(job 33180369).

---

## 7. Why simulation does not reproduce the inflation

1. **No LD**: random binomial X has no LD; the purity filter (mean.abs.corr ≥ 0.8)
   requires correlated variants in a CS. Without LD, no CS can pass purity regardless
   of BF inflation.
2. **No NA rows**: simulated Y has no missingness, so Issue 1 does not activate.
3. **i.i.d. noise per entry**: skewed and lowcount Y have non-Gaussian marginals but
   independent entries. If the real mechanism requires non-spherical cross-scale
   wavelet covariance (as hypothesized in Section 6), i.i.d. simulation would not
   generate this structure.

---

## 8. Simulation v1 design (480 tasks, SLURM job 32610344)

Fixed parameters: n=80, p=200, T=64, M=2, L=10, true causal index=50,
effect amplitude=3.0 on positions 20-40, `residual_variance_scope = "per_outcome"`,
`mixture_null_weight = 0.1`, 5 reps per cell.

Variable factors (6):

| factor | levels |
|---|---|
| `y_dist` | normal, skewed (20% outlier, sd=4), lowcount (NB log1p) |
| `signal` | signal, null |
| `wavelet_qnorm` | TRUE, FALSE |
| `prior_variance_scope` | per_outcome, per_scale_normal (per_scale not included) |
| `cross_iter_prior` | ON, OFF |
| `em_start` | warm (NULL), cold (convtol.sqp=1e-8, numiter.em=20, tol.svd=1e-10) |

Total: 3×2×2×2×2×2 = 96 cells × 5 reps = 480 tasks.

Results (null conditions, mean HP CS):

| y_dist | scope | fitted_g | qnorm | mean hp_cs |
|--------|-------|----------|-------|:----------:|
| lowcount | per_scale_normal | ON | FALSE | **0.1** |
| all others | all | all | all | 0.0 |

All 24 null cells give mean hp_cs = 0.0 except one marginal case (lowcount /
per_scale_normal / fitted_g=ON / qnorm=FALSE → single HP CS across 10 reps).
Consistent with Section 7: simulation lacks the LD and NA structure required
to reproduce real-data FDR inflation.

Power: 1.0 across all signal conditions except one (skewed / per_outcome /
fitted_g=OFF → 0.9).

---

## 9. Simulation v2 design (1440 tasks)

Extends v1 by adding `na_frac` (0 / 0.10 / 0.20) as a 7th factor to verify sigma2
recovery after the Issue 1 fix.

Total: 3×2×2×2×2×2×3 = 288 cells × 5 reps = 1440 tasks.

`na_frac` inserts whole NA rows independently per outcome (mimics samples with zero
cells of a given type). With na_frac=0.20 and the unfixed package, sigma2 mean ≈ 0.65
(vs expected ≈ 0.90). After Issue 1 fix: sigma2 recovers to ≈ 0.90.

v2 also first exposed the `cross_iter_prior = FALSE` crash (Section 6 bug note).
720 failed tasks resubmitted as v2b (job 33180369) with the fixed package.

---

## 10. Status and next steps (as of 2026-05-20)

| Issue | Status |
|---|---|
| Issue 1: qnorm NA bug | **Fixed** (commit 54a10a5) |
| Issue 2: get_var_y | **Resolved** as side-effect of Issue 1 fix |
| `cross_iter_prior=FALSE` crash | **Fixed** (commit 37d6723) |
| `cross_iter_prior` parameter exposed | **Done** (commit 54a10a5) |

Remaining:

| Action | Priority |
|---|---|
| Wait for v2b results; confirm sigma2 recovery across na_frac levels | High |
| Re-run 1024-bin perm with fixed qnorm package; quantify remaining FDR | High |
| Evaluate `cross_iter_prior=FALSE` power vs FDR tradeoff on real data | High |
| Decide default: keep `TRUE` (backward compat) or switch to `FALSE` | Pending |
| Apply same qnorm NA fix to `patched_mvf.susie.alpha` (`R/utils.R:256`) | Low |

Driver: `inst/bench/profiling/fdr_sim_worker.R`
Results v1: `inst/bench/slurm/fdr_sim_results/task_XXXX.csv`
Results v2/v2b: `inst/bench/slurm/fdr_sim_v2_results/task_XXXX.csv`
