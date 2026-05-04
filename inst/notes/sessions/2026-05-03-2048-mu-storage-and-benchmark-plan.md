# 2026-05-03 mfsusieR PR plan: per_scale_normal benchmark + mu/mu2 storage

Source documents
- Slack + meeting digest at `/home/anjing.liu/mydata/anjing.liu/project/mfsusie/notes` (consolidated 4-29 conversation, 5-03 update).
- Latest binding directive (5-03): build PR around `per_scale_normal`; run 6-combination benchmark; pull v0.0.2; install ebnm via pixi; update susieR from GitHub.

## Decisions confirmed with user

| Item | Decision |
|---|---|
| `save_mu_method` default | `"complete"` (no surprise breakage of existing fits / warm-start) |
| `save_mu_method` modes | three-way: `"complete"` (default, p×T, can warm-start), `"alpha_collapsed"` (1D = α %*% μ, lossless for coef/predict/post_smooth, no warm-start), `"lead"` (1D = μ_l,j*, biased cheap-coef, no warm-start). Plus `mf_thin(fit, method)` post-fit helper for users who want both a complete checkpoint and a thinned distribution copy. |
| Benchmark fixture | reuse `data-raw/make_practical_dataset.R` |
| Top-level priority | PR around `per_scale_normal` (Gao 5-03), 6-grid benchmark first |

## Fact-finding results

1. **`mixsqp_null_penalty` does not exist in current main.** `grep -rn mixsqp_null_penalty R/ man/` returns zero hits. All current usages are `mixture_null_weight` (`R/mfsusie.R:278/292/296/368`, `R/em_helpers.R:90/109/123/136`, `R/individual_data_methods.R:749/752/799`, `man/mfsusie.Rd:41/211`). The rename Gao mentioned in chat has not landed. Two possibilities:
   - (a) Gao plans the rename but it is not yet in code. Action: open a renaming PR with `lifecycle::deprecate_warn` alias.
   - (b) Gao's "penalty" refers to mixsqp's internal `weight` argument, not the public mfsusie() arg. Action: confirm with Gao.
2. **William's "top"-like storage is α-weighted, not lead-SNP only.** `mvf.susie.alpha/R/operation_on_multfsusie_obj.R:2080-2096` (`update_cal_fit_func`) builds `multfsusie.obj$fitted_func[[l]][[k]] = wr( α_l %*% fitted_wc[[l]][[k]] )`. That is the correct posterior mean, collapsed to a single 1D curve per (effect, outcome). William keeps both `fitted_wc` (full p×T) and `fitted_func` (1D). He does NOT trim to lead SNP.
   - Gao's Slack said "save only top per CS (refer to williams code) ... wrong but useful: mu_top".
   - Two possible readings: (A) lead-SNP `μ_l,j*` (user's instruction, "wrong but useful" = biased toward lead); (B) α-collapsed 1D (William's actual pattern, correct posterior mean).
   - Per user instruction, proceed with (A) lead-SNP. Flag this as Open Question OQ-1 below.

## Priority order (binding, user-confirmed 2026-05-03)

P0: pull latest packages (susieR github + ebnm via pixi).
P1: fix the mu/mu2 storage issue (the `save_mu_method` PR — was PR-2, now PR-1).
P2: `per_scale_normal` validation + 6-grid benchmark (was PR-1, now PR-2).
P3: `mixture_null_weight` → `mixsqp_null_penalty` rename (blocked on OQ-2).

The benchmark in P2 produces interpretable FDR/power numbers ONLY after P1 lands, because object-size and downstream coef/predict/post_smooth behavior depend on `save_mu_method`. Running benchmarks before P1 would just have to be redone.

## PR scope

**PR-1 (P1): `save_mu_method = c("complete", "top")` + downstream compat**

- (P0 prerequisite, done before this PR) susieR pulled from GitHub master, ebnm available in pixi env.
- New arg `save_mu_method` on `mfsusie()` and `fsusie()`, default `"complete"` (no behavior change for existing callers).
- `top` mode trims `fit$mu[[l]][[m]]` and `fit$mu2[[l]][[m]]` to the lead SNP `j* = argmax_j fit$alpha[l, ]`; store `fit$top_index[l]` and `attr(fit, "save_mu_method")`.
- Helper `get_effect_curve(fit, l, m)` unifies coef/predict/post_smooth read paths (auto-dispatch on stored shape).
- `coef.mfsusie`: complete → `Σ_j α_lj μ_lj`; top → `μ_l,j*` (mark cheap-coef in attribute, not the strict posterior mean).
- `predict.mfsusie`, `fitted.mfsusie`: route through helper.
- `mf_post_smooth`: dispatch on mu dimension (p-dim → α-weight; 1-dim → use directly).
- `model_init` guard: if `attr(model_init, "save_mu_method") == "top"`, stop with a clear message (warm-start checkpoint requires complete state).
- Tests: `tests/testthat/test_save_mu_method.R` covering coef shape/value, predict, post_smooth, object size reduction, model_init error path, complete-mode numerical equivalence with old fits at tolerance `1e-12`.
- Docs: roxygen for new arg in `R/mfsusie.R`, regenerated `man/mfsusie.Rd`, vignette section in `vignettes/post_processing.Rmd` on storage modes, NEWS.md entry.
- Block on OQ-1 only if Gao prefers α-collapsed 1D over lead-SNP; current direction is lead-SNP per user instruction.

**PR-2 (P2): `per_scale_normal` validation + 6-grid benchmark**

- Confirm `per_scale_normal` semantics match Gao's vignette claim: faster, ignores `mixture_null_weight`, π_0 fit by `ebnm::ebnm_point_normal`. Add a sanity test if missing.
- Benchmark script: `inst/bench/profiling/benchmark_per_scale_normal_6grid.R` (declare estimated runtime in header per CLAUDE.md hard rule 3).
- Grid:
  ```r
  bench_grid <- rbind(
    expand.grid(wavelet_qnorm        = c(FALSE, TRUE),
                prior_variance_scope = "per_scale",
                mixture_null_weight  = c(0.05, 0)),
    data.frame   (wavelet_qnorm        = c(FALSE, TRUE),
                  prior_variance_scope = "per_scale_normal",
                  mixture_null_weight  = NA_real_)
  )
  # 6 rows total
  ```
- Metrics per cell: empirical FDR @ 0.05 PIP threshold, power, #CS, CS purity, CS coverage, runtime, memory, niter, warnings.
- Fixture source: `data-raw/make_practical_dataset.R`. If wall-clock projects > 30 min, ask user before running.
- Output: results table + summary memo at `inst/notes/sessions/2026-05-XX-per-scale-normal-benchmark-results.md`.

**PR-3 (P3, blocked on OQ-2): `mixture_null_weight` → `mixsqp_null_penalty` rename**

- Add new arg `mixsqp_null_penalty`; if user passes `mixture_null_weight`, `lifecycle::deprecate_warn("0.0.3", ...)` and forward.
- Replace internal usages.
- Update man, vignettes, NEWS.md, tests.
- Hold until OQ-2 resolved.

## Default-value table (current vs. proposed)

| Arg | Current default (`R/mfsusie.R`) | Proposed | Notes |
|---|---|---|---|
| `prior_variance_scope` | `"per_outcome"` | unchanged for now | benchmark may motivate switching to `"per_scale_normal"` post-PR-1 |
| `wavelet_qnorm` | `FALSE` | unchanged | Gao 5-03 confirmed FALSE pending rebenchmark |
| `mixture_null_weight` | `NULL` (→ 0.05 internally) | unchanged | rename pending OQ-2 |
| `null_prior_init` | `0` | unchanged | only init; EM overrides |
| `small_sample_correction` | `FALSE` | unchanged | issue #8, sensitivity only |
| `L`, `L_greedy`, `greedy_lbf_cutoff` | `20`, `5`, `0.1` | unchanged | per issue #11 |
| `save_mu_method` | (not present) | add as `"complete"` in PR-1 | |

## Open questions for Gao (OQ-N)

- ~~**OQ-1** (resolved 2026-05-03)~~: support both via three-way `save_mu_method = c("complete", "alpha_collapsed", "lead")` plus `mf_thin()` helper. `alpha_collapsed` is the William-style 1D summary that smoothers already consume (`α %*% μ`), lossless for downstream consumers. `lead` is the biased lead-SNP cheap-coef. Neither 1D mode supports `model_init` warm-start; users who want both warm-start and storage savings keep a complete checkpoint and a thinned copy via `mf_thin()`.
- ~~**OQ-2** (resolved 2026-05-03)~~: `mixsqp_null_penalty` is a deprecated former name; current name is `mixture_null_weight`. No rename PR needed.
- **OQ-3**: Where should the benchmark fixtures live — `inst/bench/profiling/` (current convention) or `tests/testthat/fixtures/`? PR-1 will use `inst/bench/profiling/` unless told otherwise.

## Three-layer warm start (terminology pin, do not conflate)

1. mixsqp internal warm start — already on by default (commit `be2722e perf(ibss): mixsqp warm start + ser_cache`).
2. `L_greedy` ramp 5 → L=20 — already default.
3. operational warm start via `model_init` (cheap fit → expensive refit) — works only with `save_mu_method = "complete"` (PR-1 will guard against top).

## Process discipline

- Every code change ships docs (roxygen + man) and tests in the same commit, per CLAUDE.md.
- Default-value changes go through this memo first; do not edit `R/mfsusie.R` defaults silently.
- Track open work in TaskList; mark each item completed only when its docs + tests are in.
- This memo lives at `inst/notes/sessions/2026-05-03-2048-mu-storage-and-benchmark-plan.md` and is the canonical handoff for next session.

## Vignette / test parameter audit (2026-05-03)

Audit driver: confirm all vignette and test code uses current public-API parameter names (`mixture_null_weight`, `prior_variance_scope`, `wavelet_qnorm`, `wavelet_magnitude_cutoff`, `null_prior_init`, `alpha_thin_eps`, ...), not legacy names from `mvf.susie.alpha` / `fsusieR` (`nullweight` as null prior, `null_weight`, `null_prior_weight`, `max_SNP_EM`, `low_count_filter`, `cor_small`, `maxit`, `cov_lev`, `min_purity`, `init_pi0_w`, `posthoc`).

Findings:
- Vignettes: the only legacy-name-looking string is `nullweight = 300` at `vignettes/post_processing.Rmd:269-276`. This is **not** a mfsusieR null-prior knob; it is `ashr::ash()`'s own `nullweight` argument forwarded through `mf_post_smooth(..., method = "ash", nullweight = ...)`'s `...`. Action: keep the name (ash's own), but tighten the surrounding prose to disambiguate from `mfsusie()`'s `mixture_null_weight`.
- `vignettes/fsusie_intro.Rmd:256`: uses `mixture_null_weight` correctly (current name).
- Tests: legacy names in `test_fsusier_degeneracy.R`, `test_mvf_alpha_degeneracy.R`, `test_cs_parity_fsusier.R` (`backfit`, `maxit`) are arguments of the apple-to-apple comparison targets (`mvf.susie.alpha::multfsusie`, `fsusieR::susiF`). These must NOT be renamed; they are the upstream packages' parameter names and the equivalence contract requires calling those signatures literally.
- `test_public_api_naming.R:15-25` already lists legacy names (`max_SNP_EM`, `max_scale`, `init_pi0_w`, `min_purity`, ...) as **forbidden** on the mfsusieR public API. No change needed.
- `test_per_scale_normal_degeneracy.R:684-706` uses `mixture_null_weight` correctly and validates the warning emitted when it is passed under `prior_variance_scope = "per_scale_normal"` (where it is ignored).
- `test_ti_uniform.R`, `test_smash_lw.R`: `nullweight = 300` again forwarded to `ashr::ash` via the smoother kernel, not a mfsusie() param. Keep.

Conclusion: only `vignettes/post_processing.Rmd:269-276` needs a small prose tightening. No code changes to vignettes or tests for parameter names. PR-1 will land that prose tightening alongside the `save_mu_method` documentation block.

## Status (will be updated as work progresses)

- [x] Pull latest mfsusieR main (at `fa53ad2`)
- [x] Read Slack + meeting digest, write this plan memo
- [x] Verify `mixsqp_null_penalty` absence
- [x] Verify William's `fitted_func` storage pattern
- [x] Update susieR from GitHub master in local R lib (0.16.1, commit 220a191)
- [x] Confirm `r-ebnm` is in pixi env (1.0.55)
- [x] PR-1 branch `fix-mu-storage`, full implementation, full test suite green (0 fail / 0 regression)
- [ ] Notify Gao on issue #11: three-way `save_mu_method` design + `mf_thin()` helper
- [ ] PR-1 6-grid run + results memo
- [ ] PR-2 branch (after PR-1 merges)

## susieR 0.16.1 trace-format compat (2026-05-03, surfaced during PR-1 P0)

Symptom: with susieR upgraded from the previously-installed older version to GitHub master 0.16.1 (commit `220a191`), `tests/testthat/test_ibss_methods_branches.R:35:3` ("`track_ibss_fit.mf_individual` exercises the recording branch when `track_fit = TRUE`") errors with:

```
Error in `make_susie_track_history(model)`: fit$trace is not a compact SuSiE track
```

Root cause: susieR 0.16.1 added `is_compact_track_snapshots(fit$trace)` to `ibss_finalize -> make_susie_track_history`. The validator requires each tracking list element to be a `list(alpha = data.frame, effect = data.frame, iteration = data.frame)`. The old `track_ibss_fit.mf_individual` wrote `list(alpha = matrix, sigma2 = list, pi_V = list, elbo = scalar)`, which the new validator rejects.

Fix: delegate to susieR's own `make_track_snapshot(model, iteration)` helper (cached as a package-level binding via `R/zzz.R::.onLoad`'s internal-cache list), and pre-sanitise `model$sigma2 <- NA_real_` for the snapshot copy because mfsusieR's `sigma2` is `list[M]` whereas susieR's `track_scalar` assumes a length-1 numeric. Real per-iteration sigma2 stays on the fit; only the trace records NA. ELBO trace is unaffected (lives on `fit$elbo`).

Verified: pre-fix, `test_ibss_methods_branches.R` is FAIL on both `main` and `fix-mu-storage` (the failure pre-dates the PR; surfaced because P0 upgraded susieR). Post-fix on `fix-mu-storage`, the test is PASS and the full suite is PASS (0 fail).

Scope: this is a susieR upstream tracking compat fix, not part of the original `save_mu_method` design. Landed in the same PR because P0 upgraded susieR and surfaced it; documented separately in `NEWS.md` and `R/ibss_methods.R::track_ibss_fit.mf_individual` roxygen.
