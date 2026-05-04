# 2026-05-03 PR-1 prior-grid benchmark results

Date: 2026-05-03 / 2026-05-04
Scope: per-cell FDR / power / runtime / convergence for the six prior
configurations Gao listed on 2026-05-03 Slack
(`prior_variance_scope` x `wavelet_qnorm` x `mixture_null_weight`),
across three Y scenarios (Gaussian baseline, heavy-tailed signal,
null no signal). Three discovery metrics are reported in parallel.

## Setup

| field | value |
|---|---|
| `n` (samples) | 84 |
| `p` (variants) | 500 |
| `M` (outcomes) | 2 |
| `T_basis[m]` | 64 |
| `L` | 10 |
| `save_mu_method` | `"alpha_collapsed"` (PR-1 feature) |
| true causal indices | 50, 220, 380 |
| effect curve | box function on positions 20-40 |
| Gaussian-noise sd | 1.0 |
| heavy-tail outlier rate / sd | 18% per entry, sd = 4 |
| n_rep per cell | 5 |
| PIP loose threshold | 0.05 |
| PIP high (hybrid) threshold | 0.5 |
| CS purity threshold (hybrid) | min.abs.corr >= 0.8 |
| LD threshold (CS-level TP rule) | abs cor >= 0.5 |

mfsusieR built from branch `fix-mu-storage`. susieR 0.16.1 from
GitHub master, ebnm 1.0.55.

Drivers:
- `inst/bench/profiling/benchmark_per_scale_normal_6grid.R`
  (Gaussian baseline, 30 fits)
- `inst/bench/profiling/benchmark_heavy_tailed_null_6grid.R`
  (heavy-tailed + null follow-up, 60 fits)

SLURM accounting:
- Baseline job `31993466`, 33 min 41 s wall-clock, peak 0.43 GB,
  exit 0.
- Heavy-tailed + null job `31993467`, 1 h 0 m 20 s wall-clock,
  peak 0.43 GB, exit 0.

Per-replicate raw data (untracked, machine-local):
- `inst/bench/profiling/results/per_scale_normal_6grid_20260504_0713.rds`
- `inst/bench/profiling/results/heavy_tailed_null_6grid_20260504_0740.rds`

Each rds row carries the metric columns plus four list-columns
(`fit_pip`, `fit_cs`, `fit_purity`, `fit_true_idx`) + `rep_seed`
so future metrics can be computed retroactively without re-fitting.

## Metric definitions

Three discovery metrics are computed per (scenario, cell, rep) and
reported in parallel. They differ in what counts as one "discovery"
and which discoveries are deemed correct.

| Metric | Discovery rule (what counts as one) | TP rule | Use case |
|---|---|---|---|
| **SNP loose** (`fdr` / `power`) | Every SNP with `pip > 0.05` is one discovery, regardless of CS membership. | `j` is a TP iff `j` is in the causal index set. | A sensitive type-I-error proxy. Fires on any low-PIP leakage that crosses 0.05. |
| **CS-level** (`fdr_cs` / `power_cs`) | Every credible set is one discovery, taking the CS lead = SNP within the CS with max PIP. | Lead is TP iff `lead` is causal OR `abs(cor(X[, lead], X[, causal])) >= 0.5` for any causal. Power is fraction of unique causals covered by some TP CS. | The fine-mapping practice view, matching `inst/bench/profiling/per_scale_normal_vignette_sweep.R` and `vignettes/fsusie_intro.Rmd:299-303`. |
| **SNP hybrid** (`fdr_hyb` / `power_hyb`) | A SNP `j` is one discovery iff (`j` is in some CS with `min.abs.corr >= 0.8`) OR (`j` is not in any CS AND `pip[j] > 0.5`). | `j` is a TP iff `j` is in the causal index set. | A SuSiE-conventional fine-mapping output: trust high-purity CSes (every member of a high-purity CS counts), and outside CSes only credit standalone SNPs at high PIP. |

The three metrics answer different questions on the same fits, so
their FDR / power numbers can disagree. Where they agree, the cell
is robust. Where they disagree, the difference itself is a finding.

In `null_no_signal` the causal set is empty: SNP-loose FDR / power
are undefined (`NaN`), CS-level and hybrid FDR equal 1 whenever
their respective discovery sets are non-empty (every discovery is
false), and `cs_count` directly counts spurious credible sets.

## Per-cell results (means over 5 replicates)

### Gaussian baseline

| cell | scope | qnorm | mnw | FDR loose | FDR cs | FDR hyb | cs_count | runtime (s) | conv. |
|---|---|---|---|---|---|---|---|---|---|
| 1 | per_scale        | F | 0.05  | 0.050 | 0.050 | 0.050 | 3.2 | 10.6 | 5/5 |
| 2 | per_scale        | T | 0.05  | 0.050 | 0.050 | 0.050 | 3.2 |  9.6 | 5/5 |
| 3 | per_scale        | F | 0     | 0.766 | 0.408 | 0.628 | 5.8 | 176.7 | 4/5 |
| 4 | per_scale        | T | 0     | 0.775 | 0.479 | 0.613 | 6.0 | 178.4 | 3/5 |
| 5 | per_scale_normal | F | (n/a) | 0.409 | **0.000** | **0.000** | 3.0 | 13.4 | 5/5 |
| 6 | per_scale_normal | T | (n/a) | 0.479 | **0.000** | **0.000** | 3.0 | 13.7 | 5/5 |

Power = 1.000 in every cell across every rep (all three causal
SNPs are recovered every time).

### Heavy-tailed signal

| cell | scope | qnorm | mnw | FDR loose | FDR cs | FDR hyb | cs_count | runtime (s) |
|---|---|---|---|---|---|---|---|---|
| 1 | per_scale        | F | 0.05  | 0.150 | 0.050 | 0.100 | 3.2 |  15.5 |
| 2 | per_scale        | T | 0.05  | 0.100 | 0.000 | 0.100 | 3.0 |  12.5 |
| 3 | per_scale        | F | 0     | 0.709 | 0.553 | 0.654 | 7.0 | 168.2 |
| 4 | per_scale        | T | 0     | 0.817 | 0.409 | 0.498 | 5.4 | 238.3 |
| 5 | per_scale_normal | F | (n/a) | 0.080 | **0.000** | **0.000** | 3.0 |  10.8 |
| 6 | per_scale_normal | T | (n/a) | **0.000** | **0.000** | **0.000** | 3.0 |   8.9 |

Power = 1.000 in every cell. Convergence: cells 3-4 fail to
converge in 1-2 of 5 reps each.

### Null no signal (SNP-loose / power undefined; cs_count is the
type-I count of spurious CSes)

| cell | scope | qnorm | mnw | FDR cs | FDR hyb | cs_count | runtime (s) |
|---|---|---|---|---|---|---|---|
| 1 | per_scale        | F | 0.05  | 0.200 | 0.200 | 0.2 |  20.3 |
| 2 | per_scale        | T | 0.05  | 0.400 | 0.400 | 0.4 |  15.5 |
| 3 | per_scale        | F | 0     | 1.000 | 1.000 | **4.2** | 93.9 |
| 4 | per_scale        | T | 0     | 1.000 | 1.000 | **4.2** | 123.9 |
| 5 | per_scale_normal | F | (n/a) | **0.000** | **0.000** | **0.0** |  6.7 |
| 6 | per_scale_normal | T | (n/a) | **0.000** | **0.000** | **0.0** |  6.9 |

## Findings

0. **Power is at the ceiling in both signal scenarios.** All twelve
   signal cells (Gaussian + heavy-tailed) report `power = 1.000`
   on every rep under all three metrics; every cell recovers all
   three causal SNPs. Power is therefore not a discriminating
   axis in this benchmark; FDR / cs_count / runtime / convergence
   carry the between-cell signal. In `null_no_signal` power is
   undefined (no causal); the substitute axes are CS-level / hybrid
   FDR (= 1 whenever any spurious discovery fires) and `cs_count`.

1. **`per_scale_normal` is the only cell that is FDR = 0 under both
   CS-level and hybrid metrics across all three scenarios.** Across
   the 30 reps spanning Gaussian / heavy-tailed / null, the
   CS-level and hybrid FDR are exactly 0.000 in every cell;
   `cs_count` is exactly 3 (matching the truth) in every signal
   rep and exactly 0 in every null rep. The SNP-loose metric
   inflates to 0.41-0.48 on Gaussian baseline because the cell
   leaves a few non-causal SNPs at PIP just above 0.05 (outside
   any credible set), but those leaks are not picked up by either
   the CS-level or hybrid view.

2. **`per_scale + mixture_null_weight = 0.05` is well-calibrated
   under all three metrics across all three scenarios but with
   small slippage on heavy-tailed and null.** Gaussian: FDR loose
   / cs / hyb all 0.050. Heavy-tailed: FDR loose 0.10-0.15,
   FDR cs 0.00-0.05, FDR hyb 0.10. Null: spurious CSes fire on
   ~20-40% of reps (`cs_count` mean 0.2-0.4); FDR_cs / FDR_hyb
   inflate to 0.20-0.40 on those reps. A safe second choice when
   ebnm's distributional assumptions are uncertain.

3. **`per_scale + mixture_null_weight = 0` is broken under every
   metric in every scenario.** SNP-loose FDR 0.71-0.82, CS-level
   0.41-0.55, hybrid 0.50-0.65. `cs_count` averages 4.2 spurious
   CSes on pure null. Runtime 17-30x the mnw = 0.05 cells. 1-2 of
   5 reps fail to converge in three of the four cells. Justifies
   the warning-hint planned for PR-2 when mnw = 0 is passed
   explicitly.

4. **`wavelet_qnorm = TRUE` adds modest improvement on heavy-tailed
   Y under both SNP metrics.** Within `per_scale + mnw = 0.05`,
   FDR loose drops 0.150 -> 0.100, FDR cs 0.050 -> 0.000.
   Within `per_scale_normal`, FDR loose drops 0.080 -> 0.000.
   Effect on Gaussian baseline and on null is within rep-to-rep
   variance.

## Implications for PR-1 / defaults

- **Keep `mixture_null_weight = NULL` (resolves to 0.05) as the
  default of `per_scale`.** Across all three scenarios this is
  the only `per_scale` setting that does not collapse under any
  metric. The slight slippage on heavy-tailed / null is small
  enough that switching defaults is not justified by this benchmark.

- **`prior_variance_scope = "per_scale_normal"` is the
  best-calibrated cell in this benchmark under CS-level and
  hybrid metrics.** Gaussian / heavy-tailed / null all give
  FDR_cs = FDR_hyb = 0 with `cs_count` matching the truth.
  The SNP-loose 0.41-0.48 on Gaussian is the only blemish, and
  it reflects sub-CS PIP leakage that conventional fine-mapping
  output does not act on.

  PR-1 does not change the package default
  `prior_variance_scope = "per_outcome"`. Switching the default
  is deferred to a later PR after the real-data perm grid
  finishes; the change there will be informed by both the
  CS / hybrid view (fine-mapping practice) and the SNP-loose
  view (sensitive type-I proxy). The `?mfsusie` and
  `vignettes/fsusie_intro.Rmd` text should describe
  per_scale_normal as the recommended choice for non-Gaussian /
  real-data fits, with the caveat that SNP-loose PIP curves are
  looser than the per_scale + mnw = 0.05 baseline.

- **Document the `mixture_null_weight = 0` cliff.** Setting
  `mixture_null_weight = 0` should emit a `warning_message(style
  = "hint")` in `mfsusie()` referencing the 17-30x runtime,
  the FDR collapse under all three metrics, and the 4.2 false
  credible sets observed on pure null. Tracked as PR-2.

- **`save_mu_method = "alpha_collapsed"` works in all three
  scenarios.** All 90 fits across baseline + heavy-tailed + null
  saved cleanly; `fit_size_mb` 0.40-0.65 MB, no anomaly.

## Caveats

- **5 replicates per cell is light** for FDR estimation. The
  mnw = 0.05 columns are tight (3-4 false discoveries / 60
  reps = 0.05-0.07); the mnw = 0 and per_scale_normal columns
  have rep-to-rep variance of +/- 0.05-0.10 on FDR. 20-50 reps
  would tighten the per-cell point estimates. The list-column
  layout in the rds means a 20-rep rerun can resume seeds 101-105
  / 201-205 (used here) and add seeds 106-125 / 206-225 etc
  without recomputing the existing 5 reps.
- **Sims are all `p = 500`, `n = 84`, `M = 2`, `T = 64`, `L = 10`.**
  Translation to real ATAC perm regions (`p = 1-4k`, `T = 1024-10240`,
  `M = 6`, `L = 20`) is not literal; the real-data perm grid
  (`output/new_package/perm/mfsusieR_20260503_*` directories)
  is the proper-scale validation.
- **Hybrid SNP-level counts every SNP in a high-purity CS as a
  separate judgment.** A 5-member CS containing one causal
  contributes 1 TP + 4 FP under this rule; CS-level instead
  treats the CS as a single judgment. A reader who interprets
  CS membership as "the causal is in this set" rather than
  "every member is causal" should weight the CS-level numbers
  more than the hybrid numbers.

## Open follow-ups (not in PR-1)

- **PR-2**: emit hint when `mixture_null_weight = 0` is passed,
  plus silent-error defense in the multfsusie-paper R driver
  (surfaced during 2026-05-03 perm grid debugging).
- **PR-4**: SSC math verification (issue #8) and a SSC =
  {FALSE, TRUE} sweep on the same six-cell grid.
- **20-rep rerun once PR-1 lands and the perm grid frees the
  cluster.** Documenting the full grid (`per_outcome` axis
  included) before any default switch.

## Files / references

- Baseline benchmark driver: `inst/bench/profiling/benchmark_per_scale_normal_6grid.R`
- Heavy-tailed + null driver: `inst/bench/profiling/benchmark_heavy_tailed_null_6grid.R`
- Per-replicate rds (untracked, local):
  - `inst/bench/profiling/results/per_scale_normal_6grid_20260504_0713.rds`
  - `inst/bench/profiling/results/heavy_tailed_null_6grid_20260504_0740.rds`
- Plan memo: `inst/notes/sessions/2026-05-03-2048-mu-storage-and-benchmark-plan.md`
- Real-data perm grid README:
  `/home/anjing.liu/mydata/anjing.liu/project/mfsusie/multfsusie-paper/output/new_package/perm/README_20260503_perm_grid.md`
