# 2026-04-28 audit-and-improve findings (read-only batch)

Two read-only audit agents + a heavy-fixture re-measurement.
Source: the active OpenSpec change
`audit-and-improve-public-contract`. Findings inform the
implementation order; no code changed in this session.

## A. S3 override audit

### A1. mfsusieR overrides — three classes

**Tier A: delete-and-inherit** (no upstream patch needed).
Pure deletions; susieR's default fires.

| Override | File:line | Reason |
|---|---|---|
| `validate_prior.mf_individual` | `R/ibss_methods.R:191` | body is `invisible(TRUE)`, byte-identical to `.default` |
| `check_convergence.mf_individual` | `R/ibss_methods.R:227` | trivial ELBO-only body; `.default` does ELBO + PIP fallback + stall detection + verbose tabular output. **Prerequisite**: land the `format_sigma2_summary` patch first so the default's `model$sigma2` print does not crash on the list-of-vectors shape. |
| `configure_data.mf_individual` | `R/ibss_methods.R:296` | identity body; `.default` is also identity. mfsusieR sets `unmappable_effects = "none"` upstream of dispatch, so `.individual`'s ash branch is bypassed regardless. |
| `get_cs.mf_individual` | `R/ibss_methods.R:337` | `susie_get_cs(model, X, coverage, min_abs_corr, n_purity %||% 100)`. mfsusieR sets `params$n_purity = 100` in the workhorse call, so `%||% 100` is dead. Class hierarchy `c("mf_individual", "individual")` falls through cleanly. |
| `neg_loglik.mf_individual` | `R/individual_data_methods.R:336` | `-loglik(...)` after `V <- exp(V_param)`. `.individual` is identical. |
| `get_alpha_l.mfsusie` | `R/model_class.R:135` | `model$alpha[l, ]`, byte-identical to `.default` |
| `set_prior_variance_l.mfsusie` | `R/model_class.R:187` | `model$V[l] <- V; model`, byte-identical to `.default` |
| `get_prior_variance_l.mfsusie` | `R/model_class.R:193` | `model$V[l]`, byte-identical to `.default` |
| `get_variable_names.mf_individual` | `R/ibss_methods.R:349` | likely byte-identical to `.individual`'s `assign_names(...)` call; verify before deletion |

**~9 deletions**, no upstream patch. Total LOC removed: ~50.

**Tier B: patch-susieR-and-delete** (small upstream hook
unlocks deletion).

| Override | susieR patch sketch | Deletes after patch |
|---|---|---|
| `track_ibss_fit.mf_individual` | add `track_extra_fields(model)` generic returning extra fields to merge into the per-iter snapshot; default `list()`; mfsusieR returns `list(pi_V = …, elbo = elbo[iter])` | yes |
| `cleanup_model.mf_individual` | add `cleanup_extra_fields(data)` generic returning extra field names to strip; default `character(0)`; mfsusieR returns `c("raw_residuals")` | yes |
| `compute_kl.mf_individual` | allow `SER_posterior_e_loglik` to return either scalar or list; `compute_kl.default` reads `result$eloglik %||% result`. Plus `compute_log_null_density` generic for the per-outcome `L_null`. | yes |
| `get_objective.mf_individual` | change `sum(model$KL)` to `sum(model$KL, na.rm = TRUE)` in `.default`. **Cheapest patch**, also unlocks the two mvsusieR deletions below. | yes |
| `update_model_variance.mf_individual` | factor the bounds-clamp into a `clamp_sigma2(model, params)` generic; default clamps scalar; mfsusieR returns model unchanged | yes |
| `check_convergence.default` extra-columns | add `format_sigma2_summary(model)` and `format_extra_diag(model)` generics. **User-pre-approved.** | unlocks A1's `check_convergence.mf_individual` deletion |

**~5–6 patches**, each ≤10 LOC in susieR, each deletes one
override here.

**Tier C: keep, document divergence.**
The remaining ~25 overrides have real, irreducible
divergence: multi-outcome list shapes, per-(scale, outcome)
sigma2, mixture-of-normals posterior structure. Each gets
a one-line rationale comment in its body.

### A2. mvsusieR sidebar audit (no edits; report only)

mvsusieR has analogous cleanup opportunities:

**Tier A — delete-and-inherit (5 methods):**
- `neg_loglik.mv_individual`
- `get_cs.mv_individual` (with hierarchy `c("mv_individual","individual")`)
- `get_alpha_l.mvsusie`, `get_prior_variance_l.mvsusie`,
  `set_prior_variance_l.mvsusie` (byte-equivalent to defaults)

**Tier B — patch-susieR-and-delete (5–6 methods):**
- `compute_kl.mv_individual`, `get_objective.mv_individual`,
  `cleanup_model.mv_individual`, `get_objective.mv_ss`,
  `get_cs.mv_ss`, possibly `get_variable_names.*`
- The `compute_kl` and `get_objective` patches are **shared
  with mfsusieR** — one upstream PR fixes overrides in both
  packages.

mvsusieR's `check_convergence` is not overridden, so it
already gets susieR's verbose tabular output and PIP
fallback for free. mfsusieR is the one that masks this.

The mvsusieR cleanup is reported but not actioned in this
change. Forward to mvsusieR maintainers as a sidebar.

### A3. Patch priority

1. **Tier A deletions** in mfsusieR — no coordination, do
   immediately. Recovers ~50 LOC.
2. **`format_sigma2_summary` + `format_extra_diag` upstream
   patch** (P1) — user-pre-approved; unlocks `check_convergence`
   deletion which is the single biggest diagnostic-quality
   win.
3. **`na.rm = TRUE` upstream patch in `get_objective.default`**
   (P2) — cheapest patch (one-line); deletes
   `get_objective.mf_individual` plus two mvsusieR overrides.
4. **`cleanup_extra_fields` upstream patch** (P3) — modest
   effort; deletes overrides in both packages.
5. **`compute_kl` patch** (P4) — larger contract change
   (`SER_posterior_e_loglik` return shape); needs deprecation
   cycle.
6. **`track_extra_fields` and `clamp_sigma2`** (P5, P6) —
   lower priority; debug feature and shape-specific
   respectively.

## B. Feature parity audit

### B1. Function-level parity with fsusieR + mvf.susie.alpha

Almost everything is ported (S3 dispatch on `mf_individual` /
`mfsusie` covers the ~50 internal `.multfsusie` /
`.susiF` methods). Three meaningful gaps:

| Gap | Source | Status | Action |
|---|---|---|---|
| `affected_reg` (fsusieR), `affected_reg_effect.multfsusie` | both upstream | **port-now** | ~30 LOC: convenience accessor returning per-(CS, outcome) position ranges where the credible band excludes zero. Reads existing `credible_bands` field; one-stop diagnostic for "where is the effect non-null?" |
| Per-iter verbose output | both upstream | covered by §4 of OpenSpec change (delete `check_convergence.mf_individual`) | scheduled |
| `convergence_method` + `pip_stall_window` exposure | susieR (not mfsusieR) | covered by §4b of OpenSpec change | scheduled |
| `posthoc = TRUE` (Yuan 2024 causal configurations) | mvf.susie.alpha | **port-if-asked** | speculative research direction; not a current request |
| `simu_IBSS_ash_vanilla` | fsusieR | **port-if-asked** | only used in vignettes/tests upstream; not user-facing |
| `cs_relation`, `cal_cor_cs` (CS overlap helpers) | fsusieR | **port-if-asked** | rarely invoked |

### B2. Argument-level parity on the main entrypoint

mfsusieR exposes effectively all of fsusieR's `susiF` and
mvf.susie.alpha's `multfsusie` arguments, with renames to
align with susieR's vocabulary
(`maxit`→`max_iter`, `cov_lev`→`coverage`,
`min_purity`→`min_abs_corr`, etc.). Three pending renames
(authorized in the active OpenSpec change):

- `mixture_weight_method` → `estimate_prior_variance` (TRUE/FALSE)
- `lbf_min` → `greed_lbf_cutoff`
- `prior` parameter naming: kept as `prior_variance_scope` for now

Other flagged-but-deferred renames (Phase 8):

- `prior_variance_scope = c("per_outcome", "per_scale")`
  vs upstream `prior = c("mixture_normal", "mixture_normal_per_scale")`
- `low_count_filter` (mfsusieR) vs `thresh_lowcount` (upstream)
- `wavelet_filter_number`, `wavelet_family` (mfsusieR) vs
  `filter.number`, `family` (upstream)

### B3. HMM credible band — the alleged "wide band" diagnosis

**fsusieR's wavelet-domain band formula has a real bug** (in
fsusieR's `update_cal_credible_band.susiF`,
`operation_on_susiF_obj.R:1845-1864`):

1. Aggregates per-SNP posterior **standard deviations**
   linearly weighted by alpha (`alpha %*% Smat`), rather
   than computing the mixture variance via the law of total
   variance. Inflates wavelet-domain variance.
2. Uses `3*sqrt(tt)` (3-sigma, ~99.7%) instead of
   `qnorm(0.975)*sqrt(tt)` (1.96, 95%). Bands are ~50%
   wider than nominal at the requested coverage level.

The `HMM_regression.susiF` path **does not produce a
credible band at all** — `get_fitted_effect` warns
"credible band option not available for post processing =
HMM or none". So the "wide HMM band" complaint actually
maps to fsusieR's DWT / smash / TI band path, not HMM
proper. fsusieR's HMM has no band; mfsusieR populates one
(via `mf_invert_variance_curve`).

**mfsusieR's band formula is mathematically correct** and
matches the manuscript:

- `var_w = mu2_w - mean_w^2` where
  `mu2_w = alpha %*% mu2[l]`,
  `mean_w = alpha %*% mu[l]` (proper law of total variance)
- inverse-DWT projects `var_w` to position-space variance
  via `W_inv^2` (`R/utils_wavelet.R:544-566`)
- band uses `qnorm((1+level)/2)` (proper `1-α` CI)

**Verdict.** Section 5 of the OpenSpec change can shrink:
- mfsusieR's HMM band formula is fine; no derivation needed.
- The "suppress the band" gating that mfsusieR currently
  applies should be **removed** (the band is correct;
  hiding it is what's wrong).
- Update the HMM band roxygen with a one-paragraph note
  explaining that the band uses proper mixture variance,
  unlike fsusieR's DWT-band path.

This drops Section 5 from "needs derivation note + code
change" to "delete the gating + add doc paragraph" — much
cheaper.

## C. Heavy-fixture re-measurement

Heavy fixture: `n=84, p=3500, M=6, T=c(128,128,128,128,128,128),
L=10, max_iter=30`. Three causal SNPs at indices 100, 1500,
2800 with smooth bell-shape effects; sigma = 0.5; standard
seed.

**Result: timed out at 10 minutes (exit code 143)** on
current HEAD with the cache + subsetting + cpp11 perf work
already landed.

**Diagnosis** (from code-shape reasoning, not profile yet):
the prior perf work optimized the M-step (mixsqp input
shrunk via `mixsqp_alpha_eps` from p=3500 down to ~tens of
SNPs). The **SER step path was untouched**:
`compute_residuals.mf_individual` and
`compute_ser_statistics.mf_individual` build
`bhat = X^T R / xtx_diag` across all `p` SNPs every
(effect, outcome, IBSS iter):

  per call: ~n × p × T_basis = 84 × 3500 × 128 ≈ 38M FP ops
  calls per IBSS iter: L × M = 10 × 6 = 60
  per IBSS iter: ~2.3B FP ops
  pure-R matrix code dominates wall-clock for this shape

**Action** (Section 6 of the OpenSpec change): cpp11 port
of `mf_per_outcome_bhat_shat()` core loop, with pure-R
reference at `tol = 1e-12`. Expected speedup 5-10× on the
SER step. Profile first to confirm; then port; then
re-measure.

## D. Implementation order — refined per audit findings

### Tier A (no upstream coordination, immediate)

1. §3 Tier-A deletions (~9 mfsusieR overrides, ~50 LOC
   removed). Each deletion gets a fresh-context auditor
   review per the new humanize-approach rule.
2. §1.1 retire `V` field + `summary.mfsusie` `pi_V` line.
3. §1.3 document `sigma2` shape in `@return`.
4. §7 mechanical comment polish + helper promotion.

### Tier B (ride on upstream patches, bundled)

5. Coordinate the **`format_sigma2_summary` + `format_extra_diag`** susieR patch (user-pre-approved). Once landed:
   - delete `check_convergence.mf_individual`
   - implement §4a, §4b, §4c (verbose, convergence args, agreement test)
   - implement §4g extra-diag columns
6. Coordinate the **`get_objective` `na.rm` patch + `cleanup_extra_fields` hook** (one PR if scope allows). Once landed:
   - delete `get_objective.mf_individual`, `cleanup_model.mf_individual`

### Tier C (in-mfsusieR work)

7. §4d/§4e renames: `mixture_weight_method` →
   `estimate_prior_variance`, `lbf_min` → `greed_lbf_cutoff`.
8. §5 HMM band — drop the gating, add a doc paragraph (much
   smaller than originally scoped per audit B3).
9. §2 affected_reg port (~30 LOC).
10. §6 perf — re-measurement first (§6.1 in progress); fix
    only if measurement justifies.

### Tier D (low priority)

11. Remaining patch-susieR-and-delete items (P4–P6) if
    bandwidth allows.

Each Tier-A deletion and each Tier-C / Tier-B
implementation step is followed by a 3-auditor
fresh-context review per the user's new
humanize-approach instruction:

- Auditor 1: numerical correctness (read the deleted
  override + susieR default; confirm byte-equivalence).
- Auditor 2: API surface (read public-facing roxygen +
  vignettes; confirm no silent semantic shift).
- Auditor 3: test-suite integrity (read modified tests;
  confirm tolerances and assertions still mean what they
  claim).

Each auditor returns under 200 words. Findings are
addressed before the next section starts.
