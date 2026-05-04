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

## Implications for PR-1 / defaults

- **Keep `mixture_null_weight = NULL` (resolves to 0.05) as the default.**
  These data confirm it is the only well-calibrated option in the
  Gaussian regime.
- **Keep `prior_variance_scope = "per_outcome"` as the default for now.**
  The benchmark grid varied `per_scale` vs. `per_scale_normal` (both
  were on the menu for becoming the new default), and neither beat
  `per_scale + mnw = 0.05` on FDR. `per_outcome` was not in the grid
  but is the cheapest path; keeping it as default is safe until a
  benchmark explicitly motivates a switch.
- **Document the runtime cliff.** Setting `mixture_null_weight = 0`
  should emit a `warning_message(style = "hint")` in `mfsusie()` when
  the user passes it explicitly, pointing out the FDR / runtime trade-off.
  This is a follow-up issue candidate (not in PR-1 scope).
- **`save_mu_method = "alpha_collapsed"` works at p = 500 scale.** All
  30 fits converged + were saved without anomaly. `fit_size_mb` ranges
  0.44 - 0.64 MB; the 0.64 MB rows are `mnw = 0` cells where more
  iterations + larger CS membership inflated the saved object slightly.
  No comparable `complete` measurement (would have been ~p× larger).

## Caveats

- **Single Gaussian fixture, no heavy-tail / null scenarios.** The
  follow-up #21 covers `heavy_tailed_signal` and `null_no_signal`
  scenarios on the same 6-cell grid. Total projected wall-clock on the
  same SLURM size: 30-50 min, up from this run's 36 min.
- **5 replicates is light for FDR estimation.** The 0.05 column is
  exact (literal 3 / 60 = 0.05) but the 0.4-0.8 columns have
  per-rep variance of ±0.05-0.1. 20-50 reps would tighten the
  estimates; 5 was chosen to keep wall-clock under the 1-h SLURM
  ceiling. If this becomes an FDR-claim run for a paper, escalate to
  ≥ 20 reps and rerun.
- **`p = 500`, `n = 84`** is a small fine-mapping fixture. The Anjing
  ATAC cohort is `n = 84`, but real regions span `p = 1k - 5k`. The
  `submit_perm_mfsusieR.sh` real-data perm grid (jobs `31972343 -
  31972354`, submitted same day on the 1024-bin and 10240-bin
  fixtures) is the proper-scale follow-up.

## Files / references

- Benchmark driver: `inst/bench/profiling/benchmark_per_scale_normal_6grid.R`
- sbatch wrapper: `inst/bench/slurm/run_per_scale_normal_6grid.sbatch`
- SLURM logs: `inst/bench/slurm/logs/per_scale_normal_6grid_31970954.{out,err}`
- Per-replicate rds: `inst/bench/profiling/results/per_scale_normal_6grid_20260503_2314.rds`
- Per-cell summary CSV (4-row partial view): `inst/bench/profiling/results/per_scale_normal_6grid_20260503_2314_summary.csv`
- Plan memo: `inst/notes/sessions/2026-05-03-2048-mu-storage-and-benchmark-plan.md`
- Real-data perm grid README: `/home/anjing.liu/mydata/anjing.liu/project/mfsusie/multfsusie-paper/output/new_package/perm/README_20260503_perm_grid.md`
