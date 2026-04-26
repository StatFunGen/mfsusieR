# PLAN — please review before further implementation

## What's already in working tree (reviewable now)

- `mfsusie_plot()` (single unified function, base R, Okabe-Ito palette,
  M=1 -> 2-panel, M>1 -> grid). Working.
- `mf_post_smooth(fit, method = c("scalewise", "TI", "HMM"))` —
  all three methods run end-to-end. Operate on `fit` only.
- Always-save: `fit$residuals[[m]]` and `fit$lead_X[[l]]` populated
  for every fit (the `save_residuals` argument is removed).
- `vignettes/post_processing.Rmd` (scalewise demo only, before/after).

## Open decisions I want your sign-off on before more code

### D1. cpp11 acceleration — which loops, how much speedup, when?

**Hot loops in the smoothers**:

- TI inner loop: per-row stationary wavelet transform via
  `wavethresh::wd(..., type = "station")` (~`n` calls per
  outcome × per effect). Bottleneck for large `n`.
- TI ash loop: `ashr::ash` per scale (~`log2(T_basis)` calls). Already
  in compiled code via `ashr` C++.
- HMM forward-backward: `T_basis × (2 halfK + 1)` matrix.
  Currently pure R.

**My estimate**: at the typical sizes we run (`n = 200`,
`T_basis = 64`, `L = 5`, `M = 5`), TI takes ~1–2 s end-to-end and HMM
takes <1 s. Both are dominated by `wavethresh::wd`'s C internals
already. cpp11 ports of the *outer* loops would save < 50% in the
worst case.

**Recommendation**: defer cpp11 to a separate change
`add-cpp11-smoothers`. **Decide if the smoothers run too slow on
your real data**, then prioritize. We have evidence-based
prioritization that way (profile first).

If you want it now: TI's outer "for each effect / for each row of
`Y_pos`" loop is the candidate; HMM's emission matrix is the
other.

### D2. Variance helpers `wd.var` / `AvBasis.var` / `convert.var`

The TI credible band currently uses an approximation:
apply `av.basis` to the per-coefficient posterior SD and treat
the result as the pointwise SD. This is **conservative but not
exact** — variances are not linear under `av.basis`.

The proper port is fsusieR's three helpers (~300 lines in
`R/wavelet_utils.R`). Porting gives bit-identity to fsusieR's TI
band.

**Recommendation**:
- (a) port now into mfsusieR's `R/utils_wavelet.R` so
  `mf_post_smooth(method = "TI")` is fully faithful, OR
- (b) keep the approximation for v1 and label TI's band as
  "approximate" in roxygen, port helpers in a follow-up.

(a) takes ~1–2h of careful work + tests. (b) ships now. Your call.

### D3. Test fidelity vs fsusieR

For bit-identity testing of TI / HMM against fsusieR (when
installed), what tolerance?

**Options**:
- `<= 1e-12` strict — only achievable when we port the variance
  helpers (D2 = a) and use identical seeds + ash settings.
- `<= 1e-8` loose — works regardless of whether helpers are
  ported, as long as the algorithms match in the limit.

**Recommendation**: `<= 1e-8` to start; tighten when D2 = (a).

### D4. Data fixtures (CR1/CR2 + CASS4)

**Source files** (12 MB and 25 MB):
- `~/GIT/fsusie-experiments/plot/CR1_CR2/CR1_CR2_obj.RData`
- (CASS4 source TBD; Start_CASS4.R loads from external paths)

**Trim to package-friendly size (~250 KB each)**:
- Subsample to `n = 200` individuals, `p = 100` variables, `T_m =
  64` positions per region. Keep one true causal locus.
- Strip rownames, sample-IDs, all PHI fields. Replace with
  anonymous `S001..S200`.
- Save as `data/fsusie_methyl_example.rda` (CR1/CR2) and
  `data/mfsusie_joint_example.rda` (CASS4 — multi-cell-type).

**Build scripts under `data-raw/`**, `.Rbuildignore`d.

**Recommendation**: build now. The trimming + de-identification is
straightforward.

### D5. Vignette refresh — which to refresh now?

Eight vignettes. Refreshing all means substantial diff.

**Recommendation**: prioritize three:
- `getting_started.Rmd` — show `mfsusie_plot()` with the packaged
  data.
- `post_processing.Rmd` — extend to demonstrate all three smoother
  methods (scalewise / TI / HMM) with before/after plots and an
  LFSR overlay panel for HMM.
- `fsusie_dnam_case_study.Rmd` — use `data(fsusie_methyl_example)`
  and `mfsusie_plot()` to reproduce the CR1/CR2 figure.

The other four (`fsusie_intro`, `fsusie_covariates_and_coloc`,
`fsusie_why_functional`, `mfsusie_intro`) can be refreshed in a
follow-up if you want to keep this change tight.

### D6. Unit tests

**Recommendation**: add `tests/testthat/test_post_smooth.R` with:
- Smoke for `scalewise` / `TI` / `HMM` (each runs without error;
  output shapes correct).
- Bit-identity tests skipped if `fsusieR` not installed; assert
  TI's point estimate vs `fsusieR::univariate_TI_regression` at
  D3's tolerance; HMM's point estimate vs `univariate_HMM_regression`.

## Concrete proposed scope for THIS change

Given the decisions above, I'd pick:

- D1: **defer cpp11** -> `add-cpp11-smoothers` (separate change).
- D2: **option (b)** for now — approximate band, label clearly.
- D3: `<= 1e-8` for any bit-identity test.
- D4: **build CR1/CR2 + CASS4 fixtures now**.
- D5: refresh **getting_started + post_processing + dnam_case_study**;
  defer the four motivation/intro vignettes.
- D6: smoke + bit-identity tests for all three smoothers.

That's a tight, shippable change. Tell me yes/no on each decision
and I'll execute against the agreed scope. No further code changes
until then.
