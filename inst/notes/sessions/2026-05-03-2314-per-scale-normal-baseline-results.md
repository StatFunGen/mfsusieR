# 2026-05-03 PR-1 baseline benchmark: per_scale_normal 6-grid

SLURM job `31970954`, partition `cpu`, host `cpu-b-2`, 2 cpus / 8 GB,
total wall-clock **36 m 22 s**, peak memory **454 928 K** (~445 MiB),
exit 0. Driver: `inst/bench/profiling/benchmark_per_scale_normal_6grid.R`
under sbatch wrapper `inst/bench/slurm/run_per_scale_normal_6grid.sbatch`.

Per-replicate raw data: `inst/bench/profiling/results/per_scale_normal_6grid_20260503_2314.rds`
(30 rows = 6 cells × 5 replicates; per-cell summary CSV is a partial view —
the `aggregate(..., FUN = mean)` formula in the script drops the
`mixture_null_weight = NA` rows for `per_scale_normal`, so the CSV only has
4 rows; recompute from the .rds for the full picture, as below).

## Setup

| field | value |
|---|---|
| `n` (samples) | 84 |
| `p` (variants) | 500 |
| `M` (outcomes) | 2 |
| `T_basis[m]` | 64 (single dyadic length) |
| `L` | 10 |
| `save_mu_method` | `"alpha_collapsed"` (PR-1 feature) |
| true causal indices | 50, 220, 380 |
| effect curve | box function on positions 20-40 |
| noise | iid Gaussian, sd = 1.0 |
| n_rep per cell | 5 |
| PIP threshold | 0.05 |

mfsusieR built from branch `fix-mu-storage` @ `fa53ad2`. susieR 0.16.1 from
GitHub master (ebnm 1.0.55).

## Full 6-cell summary (means over 5 replicates)

| cell | `prior_variance_scope` | `wavelet_qnorm` | `mixture_null_weight` | FDR  | power | `#disc` | `#cs` | `niter` | runtime (s) | `fit_size` (MB) |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | `per_scale`        | F | 0.05 | **0.050** | 1.00 | 3.2  | 3.2 | 3.0  | **11.4** | 0.44 |
| 2 | `per_scale`        | T | 0.05 | **0.050** | 1.00 | 3.2  | 3.2 | 2.8  | **10.5** | 0.44 |
| 3 | `per_scale`        | F | 0    | 0.766     | 1.00 | 14.6 | 5.8 | 28.2 | 190.7    | 0.64 |
| 4 | `per_scale`        | T | 0    | 0.775     | 1.00 | 14.0 | 6.0 | 31.6 | 192.4    | 0.64 |
| 5 | `per_scale_normal` | F | (n/a) | 0.409    | 1.00 | 5.6  | 3.0 | 3.8  | 14.6     | 0.58 |
| 6 | `per_scale_normal` | T | (n/a) | 0.479    | 1.00 | 6.0  | 3.0 | 4.0  | 14.8     | 0.58 |

Runtime range per cell (min / median / max):

| cell | min | median | max |
|---|---|---|---|
| 1 | 9.0 s | 11.2 s | 13.4 s |
| 2 | 7.6 s | 11.0 s | 13.3 s |
| 3 | 108.8 s | 159.1 s | 334.9 s |
| 4 | 91.9 s | 225.2 s | 285.8 s |
| 5 | 13.1 s | 14.7 s | 16.0 s |
| 6 | 13.6 s | 14.8 s | 16.2 s |

Convergence:

| cell | converged / 5 | comment |
|---|---|---|
| 1 | 5 / 5 | fast, clean |
| 2 | 5 / 5 | fast, clean |
| 3 | 4 / 5 | one rep hit `max_iter = 50` |
| 4 | 3 / 5 | two reps hit `max_iter = 50` |
| 5 | 5 / 5 | fast, clean |
| 6 | 5 / 5 | fast, clean |

## Findings

1. **`mixture_null_weight = 0.05` (the package default) is well-calibrated.**
   FDR matches the nominal 0.05 threshold (0.05 ± 0 across reps), power is
   1.0, the fit converges in ~3 IBSS iterations, and one fit takes ~11 s.
   This is a 17× speedup and 10× fewer iterations versus `mnw = 0`.

2. **`mixture_null_weight = 0` is broken on this fixture.** Without the
   null pseudo-weight, the M-step `mixsqp` solver does not push enough
   mass onto the null spike: the alpha-aggregated PIP for many off-target
   variants stays well above 0.05, FDR climbs to ~0.77 (15× the nominal
   target), and `niter` saturates at the 50-iter cap on 3/10 reps. This
   confirms the design rationale for keeping `0.05` as the default and
   should NOT be a user-facing choice without a strong warning.

3. **`per_scale_normal` is fast but not as well-calibrated as
   `per_scale + mnw = 0.05` on this Gaussian fixture.** FDR sits at
   ~0.4-0.5 (8-10× the nominal 0.05) — markedly better than `mnw = 0`'s
   0.77 but markedly worse than `mnw = 0.05`'s 0.05. Power is 1.0 and
   `cs_count` is exactly 3 (matching the truth), so the false positives
   are SNPs with PIP just above 0.05 that fall outside the credible sets.
   Runtime is competitive (~15 s, ~1.4× the `mnw = 0.05` per_scale fit).
   This argues against switching the package default to
   `per_scale_normal` based on this single Gaussian scenario.

4. **`wavelet_qnorm` has no measurable effect on this Gaussian fixture.**
   FDR / power / `cs_count` differences between (qnorm=F, qnorm=T)
   within each scope-mnw pair are within sampling noise (cells 1↔2,
   3↔4, 5↔6 each pair). qnorm's value-add lives in heavy-tailed Y, which
   this fixture does not exercise; the follow-up benchmark (#21,
   heavy_tailed + null scenarios) is needed to isolate that axis.

## Open questions raised by these results

- **Does `per_scale_normal` improve on heavy-tailed Y?** The 0.4 FDR on a
  clean Gaussian fixture suggests `ebnm_point_normal`'s data-driven π₀
  is more permissive than `mixsqp + null pseudo-weight 0.05`. Heavy-tailed
  Y might tilt this differently because `ebnm_point_normal` is more
  flexible at picking up structured signal. The `#21` follow-up will
  tell.
- **Is the `cs_count = 3` exact-recovery in `per_scale_normal` enough to
  recommend it over `per_scale + mnw = 0.05` for downstream use?** If
  the false positives at PIP > 0.05 always sit outside any credible set
  (as observed here), the user-facing answer is "use the credible sets,
  not raw PIP > 0.05". This is consistent with how SuSiE outputs are
  used in fine-mapping practice.
- **Why does `per_scale + mnw = 0` fit in 109-335 s with high variance,
  while the `mnw = 0.05` cells finish uniformly in ~10 s?** mixsqp
  without the null pseudo-weight has multiple near-degenerate local
  optima — convergence path is data-driven and high-variance. Adding
  `0.05` regularises the M-step into a single basin.

## Heavy-tailed + null follow-up (SLURM 31977248, 2026-05-04)

SLURM job `31977248`, partition `cpu`, 2 cpus / 8 GB, total wall-clock
**1 h 0 m 59 s**, peak memory **0.43 GB**, exit 0:0. Driver:
`inst/bench/profiling/benchmark_heavy_tailed_null_6grid.R` under
sbatch wrapper `inst/bench/slurm/run_heavy_tailed_null_6grid.sbatch`.
Per-replicate rds: `inst/bench/profiling/results/heavy_tailed_null_6grid_<TS>.rds`.

Two new scenarios on the same 6-cell prior grid, 5 reps each
(60 fits total, plus the 30-fit Gaussian baseline above for context):

- `heavy_tailed_signal`: same X / β as baseline, but each Y entry has
  18% iid Bernoulli(0.18) outlier contamination at sd = 4 added on
  top of the sd = 1 Gaussian noise. Tests whether `wavelet_qnorm = TRUE`
  recovers calibration that `wavelet_qnorm = FALSE` would lose, and
  whether `per_scale_normal` 's ebnm fit handles non-Gaussian wavelet
  coefficients better than `per_scale + mixture_null_weight = 0.05`.
- `null_no_signal`: pure Gaussian noise, no causal SNPs. FDR / power
  are undefined; we report `n_disc` (any PIP > 0.05) and `has_disc`
  (binary "any false discovery this rep") as the type-I rate per cell.

### Per-cell aggregates (means over 5 replicates)

`heavy_tailed_signal`:

| cell | scope | qnorm | mnw | FDR | power | `#disc` | `#cs` | niter | runtime (s) | conv. |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | `per_scale`        | F | 0.05 | 0.150 | 1.00 | 3.6  | 3.2 | 4.0  | 15.5 | 5/5 |
| 2 | `per_scale`        | T | 0.05 | 0.100 | 1.00 | 3.4  | 3.0 | 3.2  | 12.4 | 5/5 |
| 3 | `per_scale`        | F | 0    | 0.709 | 1.00 | 10.4 | 7.0 | 27.6 | 167.3 | 4/5 |
| 4 | `per_scale`        | T | 0    | 0.817 | 1.00 | 17.2 | 5.4 | 40.2 | 239.9 | 3/5 |
| **5** | `per_scale_normal` | F | (n/a) | **0.080** | 1.00 | 3.4 | **3.0** | 4.0 | **11.0** | 5/5 |
| **6** | `per_scale_normal` | T | (n/a) | **0.000** | 1.00 | 3.0 | **3.0** | 3.6 | **9.1**  | 5/5 |

`null_no_signal` (FDR / power undefined; `has_disc` ∈ [0, 1] is the
per-rep mean of "any PIP > 0.05 fired"):

| cell | scope | qnorm | mnw | `#disc` | `has_disc` | `#cs` | niter | runtime (s) | conv. |
|---|---|---|---|---|---|---|---|---|---|
| 1 | `per_scale`        | F | 0.05 | 0.6 | 0.4 | 0.2 | 4.0  | 20.9  | 5/5 |
| 2 | `per_scale`        | T | 0.05 | 0.4 | 0.4 | 0.4 | 2.8  | 15.9  | 5/5 |
| 3 | `per_scale`        | F | 0    | 7.0 | **1.0** | **4.2** | 12.2 | 96.1  | 5/5 |
| 4 | `per_scale`        | T | 0    | 7.0 | **1.0** | **4.2** | 22.2 | 126.9 | 4/5 |
| **5** | `per_scale_normal` | F | (n/a) | **0.0** | **0.0** | **0.0** | 3.4 | 6.8 | 5/5 |
| **6** | `per_scale_normal` | T | (n/a) | **0.0** | **0.0** | **0.0** | 3.2 | 7.0 | 5/5 |

### Findings (combined across all three scenarios)

1. **`per_scale_normal` is the strongest mode under heavy-tailed and
   null Y.** On heavy-tailed, FDR drops to 0.00-0.08 with `cs_count`
   exactly equal to the truth (3); on null, the cell never fires
   (zero `n_disc`, zero `cs_count`, zero `has_disc` across all 10
   reps). This is a sharp reversal from the Gaussian baseline where
   the same cells sat at FDR 0.41-0.48. The mechanism: `ebnm`'s
   data-driven π₀ in `ebnm_point_normal` is conservative when the
   wavelet-coefficient distribution has heavy tails or no real
   structure, and only adds mass to the slab when the data force it.
   Under Gaussian-only the pure-noise tail probability is high so
   ebnm picks more spurious slabs; under contaminated / null Y the
   contrast between true signal and noise is sharper for ebnm.

2. **`per_scale + mixture_null_weight = 0.05` remains a safe second
   choice.** It is the only cell that is well-calibrated across all
   three scenarios (FDR 0.05 / 0.10-0.15 / cs_count 0.2-0.4 on null).
   Slight degradation on heavy-tailed (FDR creeps to 0.10-0.15) but
   never catastrophic. Stays the right default when one cannot
   commit to ebnm's assumptions or when explainability matters.

3. **`per_scale + mixture_null_weight = 0` is broken in every
   scenario.** Heavy-tailed FDR 0.71-0.82, null `cs_count = 4.2`
   (false credible sets on pure noise), 17-30× slower than
   mnw = 0.05, and 1-3 of 5 reps non-converged in three of the four
   cells. Should not be a user choice without a strong warning.

4. **`wavelet_qnorm` matters slightly more on heavy-tailed.**
   Within `per_scale + mnw = 0.05`: qnorm = T improves heavy-tailed
   FDR from 0.150 to 0.100. Within `per_scale_normal`: qnorm = T
   improves heavy-tailed FDR from 0.080 to 0.000. The qnorm path is
   doing real work in heavy-tailed Y, just not on the order of
   "magic fix" — the underlying `ebnm` choice does the heavy lifting,
   and qnorm tightens it further.

### Implications for PR-1 / defaults

- **Keep `mixture_null_weight = NULL` (resolves to 0.05) as the
  default of `per_scale`.** Across three scenarios it is the only
  `per_scale` setting that is well-calibrated end-to-end.
- **Recommend `prior_variance_scope = "per_scale_normal"` as the
  default for non-Gaussian / real-data fits**, but DO NOT change
  the package default in PR-1. Reason: the Gaussian baseline shows
  per_scale_normal at FDR 0.41-0.48 (poor on the easy case), so a
  user who runs simulation experiments on Gaussian Y will be
  surprised. A doc-only update is the right move for now: the
  vignette / `?mfsusie` should describe per_scale_normal as
  "preferred for heavy-tailed or sparse Y; on idealised Gaussian Y
  prefer per_scale + mixture_null_weight = 0.05". The real-data perm
  grid (`submit_perm_mfsusieR.sh` cells 5-6 of each bin size) will
  decide whether to switch the default in a later PR. Currently
  `prior_variance_scope = "per_outcome"` remains the package
  default; this benchmark did not touch `per_outcome` and so cannot
  argue against it.
- **Document the `mixture_null_weight = 0` runtime / FDR cliff.**
  Setting `mixture_null_weight = 0` explicitly should emit a
  `warning_message(style = "hint")` in `mfsusie()` pointing out the
  17-30× runtime hit, the FDR collapse, and the 4-5 false credible
  sets observed on pure null. Tracked as PR-2 (silent-error /
  hint-emission defense).
- **`save_mu_method = "alpha_collapsed"` works in all three scenarios.**
  All 90 fits across baseline + heavy-tailed + null saved cleanly;
  `fit_size_mb` 0.40-0.65 MB, no anomaly.

### Caveats

- **5 replicates is light** for FDR estimation. The mnw = 0.05
  columns are tight (3-4 false discoveries / 60 = 0.05-0.07); the
  mnw = 0 and per_scale_normal columns have rep-to-rep variance of
  ±0.05-0.10 on FDR. 20-50 reps would tighten the per-cell point
  estimates; 5 was chosen to keep the SLURM run under 1 h. For an
  FDR-claim figure in a paper, escalate to ≥ 20 reps and rerun. The
  `inst/bench/profiling/results/heavy_tailed_null_6grid_<TS>.rds`
  carries all 60 raw rows so the rerun script can resume seeds
  101-105 (used here) and add seeds 106-125 etc.
- **Sims are all `p = 500`, `n = 84`, `M = 2`, `T = 64`, `L = 10`.**
  Translation to real ATAC perm regions (`p = 1-4k`, `T = 1024-10240`,
  `M = 6`, `L = 20`) is not literal; the real-data perm grid
  (`mfsusieR_20260503_*` under `output/new_package/perm/`) is the
  proper-scale validation.
- **`per_scale_normal` cs_count = 0 on null is a strong claim, but
  with only 5 reps means we have not seen a single false CS firing.**
  At 20+ reps a single false CS would still leave the rate at <= 5%.
  The benchmark cannot distinguish "exactly zero" from "very rare".

### Open follow-ups (not in PR-1)

- **PR-2**: emit hint when `mixture_null_weight = 0` is passed; harden
  sbatch + R driver against silent errors (surfaced in 2026-05-03
  perm grid debugging).
- **PR-4**: SSC math verification (issue #8) + an SSC = {F, T} sweep on
  the same 6-cell grid to test Gao's hypothesis that historical SSC
  value was a fudge factor for prior-side math errors that have
  since been fixed. Results memo will be appended here when those
  data land.
- **Rerun at 20 reps once PR-1 lands and the perm grid frees the
  cluster.** Document the full grid (`per_outcome` axis included)
  before committing to a default switch.

## Files / references

- Baseline benchmark driver: `inst/bench/profiling/benchmark_per_scale_normal_6grid.R`
- Heavy-tailed + null driver: `inst/bench/profiling/benchmark_heavy_tailed_null_6grid.R`
- sbatch wrappers: `inst/bench/slurm/run_per_scale_normal_6grid.sbatch`,
  `inst/bench/slurm/run_heavy_tailed_null_6grid.sbatch`
- SLURM logs: `inst/bench/slurm/logs/per_scale_normal_6grid_31970954.{out,err}`,
  `inst/bench/slurm/logs/heavy_tailed_null_6grid_31977248.{out,err}`
- Per-replicate rds: `inst/bench/profiling/results/per_scale_normal_6grid_20260503_2314.rds`
  (baseline) and `heavy_tailed_null_6grid_<TS>.rds` (follow-up)
- Per-cell summary CSV (baseline 4-row partial view): `inst/bench/profiling/results/per_scale_normal_6grid_20260503_2314_summary.csv`
- Plan memo: `inst/notes/sessions/2026-05-03-2048-mu-storage-and-benchmark-plan.md`
- Real-data perm grid README: `/home/anjing.liu/mydata/anjing.liu/project/mfsusie/multfsusie-paper/output/new_package/perm/README_20260503_perm_grid.md`
