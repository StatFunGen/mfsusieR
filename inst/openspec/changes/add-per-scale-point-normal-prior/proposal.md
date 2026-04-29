# Add per-scale point-normal and point-laplace priors

## Why

The current `prior_variance_scope = "per_scale"` flavor fits a
flexible scale-mixture-of-normals at every (outcome, scale) cell.
With K+1 ≈ 30 mixture weights and per-scale row counts in the
wavelet pyramid of `idx_size ∈ {1, 1, 2, 4, 8, 16, 32, 64}` (for
`T_basis = 128`), the M-step is severely under-determined on
coarse scales: 30 unknowns from 1-4 effective rows. Practical use
falls back to `per_outcome` (pools all scales, ~128 rows for the
same K+1 weights), but that forces a single mixture across
genuinely scale-heterogeneous wavelet coefficients.

A two-parameter spike-and-slab is the right tool inside a single
wavelet scale `s`. Within scale `s`, the coefficients have a
single characteristic effect-size scale; the truth is

```
beta_s ~ pi_0_s * delta_0  +  (1 - pi_0_s) * G_s(. ; sigma_s)
```

with `G_s` either Normal (`per_scale_normal`) or Laplace
(`per_scale_laplace`). Two parameters per (outcome, scale) cell
versus the K+1 of `mixture_normal_per_scale`. The data goes ~15x
further (~15x more rows per fitted parameter); the parametric
form is intrinsically regularized so the null pseudo-count knob
(`mixture_null_weight`) is unnecessary on these paths.

Both Normal and Laplace ship as parallel options through a single
`.opv_<class>` helper per slab choice sharing the same loop body.
Laplace gives the heavier tail when the wavelet-domain signal is
sparse within a scale (1-2 active positions per scale), the
common case for sharp transitions or localized features.

## What changes

Extend `prior_variance_scope` with two new public values:

```r
prior_variance_scope = c("per_outcome",         # current default
                         "per_scale",           # current alt
                         "per_scale_normal",    # NEW
                         "per_scale_laplace")   # NEW
```

Internal `prior_class` mapping:

| public scope         | internal `prior_class`            | M-step solver           |
|----------------------|-----------------------------------|-------------------------|
| `per_outcome`        | `mixture_normal`                  | `mixsqp` (current)      |
| `per_scale`          | `mixture_normal_per_scale`        | `mixsqp` (current)      |
| `per_scale_normal`   | `mixture_point_normal_per_scale`  | `ebnm::ebnm_point_normal`  |
| `per_scale_laplace`  | `mixture_point_laplace_per_scale` | `ebnm::ebnm_point_laplace` |

Both new prior classes carry G_prior in the same shape the
existing kernel reads: the slot is a `fitted_g` record returned
by ebnm. For Normal, `fitted_g` is a 2-component normalmix
(`pi = c(pi_0, 1 - pi_0)`, `sd = c(0, sigma)`, `mean = c(0, 0)`).
For Laplace, `fitted_g` is a `laplacemix` (`pi`, `scale`,
`mean`). The downstream loglik / posterior-moment kernels
iterate over `fitted_g` and remain unchanged for the Normal path;
the Laplace path adds a small posterior-moment helper for the
Laplace slab.

### ebnm mirrors the mixsqp call shape

The new M-step is a per-(outcome, scale) loop that calls ebnm,
identical in shape to the existing `.opv_mixsqp` per-(outcome,
scale) loop that calls `mixsqp`. The dispatch follows the
existing `.opv_<class>` helper-per-class pattern; no S3 generic
refactor, no parent-class hierarchy. One helper per ebnm-backed
class plus a shared body:

```r
.opv_ebnm_point        <- function(data, params, model, ser_stats,
                                    keep_idx, zeta_keep, ebnm_fn) { ... }
.opv_ebnm_point_normal  <- function(...) .opv_ebnm_point(..., ebnm_fn = ebnm::ebnm_point_normal)
.opv_ebnm_point_laplace <- function(...) .opv_ebnm_point(..., ebnm_fn = ebnm::ebnm_point_laplace)
```

The data shape ebnm sees is the same alpha-thinned
`(|keep_idx| × |idx_s|)` rectangle of (Bhat, Shat) that mixsqp
consumes, just unweighted: ebnm cannot accept observation
weights so the per-row alpha weighting drops. The unweighted-
many-variables design preserves mixsqp's structural property
that the per-(m, s) parametric MLE sees the bulk-noise
distribution, so `pi_0` settles at a sensible conservative
value rather than collapsing toward 0 (the lead-only
alternative would fit "this single variable is non-null" and
produce spurious CSes from the resulting too-liberal slab).

ebnm owns: optim transport, init heuristics, the small-n
freeze, the `g_init` warm-start path, convergence reporting,
and the `fitted_g` record format we store back into
`G_prior[[m]][[s]]$fitted_g`.

The `g_init = <previous fit>` warm-start mirrors mixsqp's
`pi_warm_start = pi_prev` pattern (always passed; no cold/warm
branch). The init helper writes `fitted_g` into
`G_prior[[m]][[s]]` at iter 0 by picking a marginal-data lead
per scale (`lead_s = which.max(rowMeans(bs$Bhat[, idx_s]^2))`)
and fitting ebnm on that lead's per-scale slice; the IBSS-loop
M-step then refits on the multi-variable rectangle.

`alpha_thin_eps` (renamed from `mixsqp_alpha_eps`; see
`mfsusie()` formal) is the unified per-effect alpha threshold
used by both M-step solvers: drop variables with `alpha[l, j]
< alpha_thin_eps` from the M-step input. Mixsqp uses it to cap
the L-matrix size; ebnm uses it to scope the multi-variable
ebnm call to non-negligible variables. Same threshold, same
truncation-error semantics.

Why ebnm:

- `ebnm_point_normal()` and `ebnm_point_laplace()` share the same
  `(x, s, g_init, fix_g)` interface, so one helper covers both
  slabs.
- `g_init = <previous fit>` warm-starts across IBSS iters without
  a bespoke cache, mirroring the mixsqp wrapper's
  `pi_warm_start`.
- `ashr` and `mixsqp` are already in `susieR`'s closure, so the
  marginal weight is `ebnm` itself.
- ebnm owns numerical edge cases (`s == 0`, degenerate `x`,
  `n <= 2`, `optim` non-convergence).

## Why this layered out the way the model does

Functional Y goes through a discrete wavelet transform that
decorrelates positions. Wavelet columns at different scales span
orders of magnitude in natural variance; the mixture-of-normals
adapts to that cross-scale heterogeneity (the `per_outcome`
section above). Within a single scale the variance is
approximately constant, and a two-component spike-and-slab is the
right Bayesian model for "most of these positions are noise; the
rest share a variance." Normal versus Laplace controls the
slab tail: Normal for diffuse activity within a scale, Laplace
for sparse activity (1-2 positions out of `idx_size`).

## Acceptance criteria

- `mfsusie(prior_variance_scope = "per_scale_normal", ...)` and
  `mfsusie(prior_variance_scope = "per_scale_laplace", ...)` both
  run end-to-end on standard fixtures and produce a fit with the
  documented shape (`G_prior[[m]][[s]]$fitted_g` is a `fitted_g`
  record; `lbf_variable`, `pip`, `sets$cs` populated).
- Three machine-precision degeneracy tests (T=1, fixed prior,
  no EB) bit-equivalent to `susieR::susie()` at `tol = 1e-12`:
  1. `per_outcome` + length-1 grid (already exists).
  2. `per_scale` + length-1 grid.
  3. `per_scale_normal` with `pi_0` clamped at 0 and a fixed
     `sigma` (no ebnm M-step on the IBSS loop). Primary
     correctness lock for the new path.
- `per_scale_normal` and `per_scale_laplace` each recover the
  planted causal on the synthetic sparse-wavelet fixture
  (`set.seed(2L); n = 200; p = 30; T_m = 64; signal at
  variables c(7, 18) and times c(20, 44)`) with PIP > 0.9 on the
  causal variables and `nCS = 2`. `per_outcome` is the
  side-by-side reference at the same fixture.
- The full unit-test suite passes.
- **Post-implementation audit, three independent agents from
  three perspectives (math / code / test), each running both
  checks below:**
  1. **Over-engineering / hallucination check.** Did the author
     add complexity, edge-case logic, optimizations, or tests
     that the proposal did not call for? Are there fictional or
     hallucinated APIs, mis-cited file:line references,
     fabricated convergence guarantees, or "AI fluff" present
     for thoroughness theatre rather than catching a real bug?
  2. **Missing-piece check.** What did the author not implement
     that the proposal called for? Are any acceptance criteria
     only partially covered? Are documented edge cases (cache
     gate on the new prior class, lead picker at init,
     warm-start through `g_init`) actually present in the code?
     Did the dispatch site forget to plumb a parameter?

  Three perspectives:
  - **Math agent**: marginal-likelihood correctness, ebnm
    output unpacking into our `fitted_g` slot, lead-variable
    coherence (init from marginal data, M-step from `alpha`),
    Laplace posterior-moment math, degeneracy-claim audit.
  - **Code agent**: `.opv_<class>` helper-per-class dispatch
    pattern matches `.opv_mixsqp`, cache gating one-liner is
    correct, covariate-adjust path untouched, no hidden
    coupling, ebnm `g_init` warm-start passed every iter.
  - **Test agent**: tolerance protocol holds, degeneracy suite
    covers all flavors at machine precision, sparse-fixture
    recovery test asserts both `per_scale_normal` and
    `per_scale_laplace`, vignette sweep CS-agreement metric is
    well-defined, benchmark fixtures exist.

  Implementation does not merge until all three report "no
  significant issues."
- **Vignette sweep**: every shipped vignette except
  `_practical_data_applications.Rmd` is re-run three times: at
  the default (`per_outcome`), at `per_scale_normal`, and at
  `per_scale_laplace`. For each pair against the default:
  (a) `|nCS_new - nCS_po| <= 1` on every workload; (b) for every
  matched CS, the new fit's lead variable equals the
  `per_outcome` CS's lead variable or appears inside it (LD-tag
  shift OK); (c) per-CS Jaccard against the matched
  `per_outcome` CS is `>= 0.5`; (d) median wallclock speedup is
  `>= 2x` on at least 4 of the 9 vignette workloads. The sweep
  script and TSV land in `inst/bench/profiling/`.
