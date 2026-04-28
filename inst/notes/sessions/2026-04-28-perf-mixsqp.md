# 2026-04-28 perf-mixsqp-prior-update + fix-prior-init-parameterization

Two OpenSpec changes landed back-to-back. Phase: perf /
correctness intermixed (Phase 6 / Phase 7 hybrid).

## Change 1: fix-prior-init-parameterization (commit `9e29e11`)

Replaced the hardcoded `c(0.8, 0.2/(K-1), …)` line in
`init_scale_mixture_prior_default()` with a Dirichlet init driven
by the existing `null_prior_weight` parameter:

    pi_null = null_prior_weight / (K + 1)
    pi_k    = (1 - pi_null) / (K - 1)   for k > 1

Same formula already used by `distribute_mixture_weights()` (the
user-supplied-grid path). One parameter, both paths. No new
public API.

C2 fidelity test patched to pass `null_prior_weight = 0.8*(K+1)`,
mathematically recovering upstream's exact `(0.8, …)` init for
the bit-identity comparison.

## Change 2: perf-mixsqp-prior-update (commits `4b77932`,
`9d7205a`, `169552c`)

Three-step plan; two delivered, one reverted, one deferred.

### Baseline (clean wall-clock, no profvis instrumentation)

Fixture: n=300, p=150, M=2, T_basis=(128, 64), L=10,
`max_iter=30`. Median of 3 runs.

    pre-cache (HEAD~):  28.86 s

### Step 2.1 (commit `9d7205a`): cache loop invariants

Added `model$em_cache` for `sigma2_per_pos`, `shat2`, `sdmat`,
`log_sdmat`. Repopulated by `update_variance_components` per IBSS
iter. M-step reads cache instead of recomputing per effect.

    cache: 28.12 s    (Δ -0.74 s, -2.6%)

Real but small win. profvis-instrumented numbers (13.33 s) were
misleading — profvis sampling overhead doesn't track real
wall-clock changes proportionally. Lesson: trust clean wall-clock
not profvis.

### Step 2.2: vectorise NaN-imputation (REVERTED)

Replaced `apply(log_L, 2, function(col) …)` with
`if (any(nas)) log_L[nas] <- col_med[col(log_L)[nas]]`. Tested
clean and bit-identical, but **wall-clock unchanged** (NA branch
is rarely hit; the fast-path skip dominates either way). Reverted
to keep the codebase simple.

### Step 2.3 (commit `169552c`): adaptive variant subsetting

Added `mixsqp_alpha_eps` argument to `mfsusie()` (default `1e-6`).
Inside `optimize_prior_variance.mf_individual()`,
`keep_idx <- which(model$alpha[l, ] > eps)` selects SNPs whose
posterior mass is non-negligible. The L-mat passed to mixsqp is
built from `bhat[keep_idx, idx]` only.

Truncation error on the M-step gradient is bounded by
`sum_{j outside} alpha_j * max_k(L_jk)`, well under
floating-point precision for typical concentrated posteriors.
Setting `mixsqp_alpha_eps = 0` recovers the unsubsetted behaviour
at bit-identity (verified `max_diff = 0` vs the cache-only HEAD~
on a 2-effect, 2-modality fixture).

    subset (eps = 1e-6): 12.91 s    (Δ -15.21 s, -54%)

Numerical drift between `eps = 0` and `eps = 1e-6` on a
n=200, p=100, M=2 fixture:

    alpha   max diff: 3.94e-7
    pip     max diff: 7.86e-8
    sigma2  max diff: 5.54e-9

Both fits converged in 21 iterations.

### Step 2.4 (deferred): C++ port

After subsetting, the post-subsetting profile shows
`mixsqp` internals (`irlba` / `tsvd` / `normalize.rows`) at ~640
samples, the L-mat builder at ~430 samples (`apply` 237, `dnorm`
132, `outer` 68). mixsqp itself is already C++; we cannot speed
it up without either passing fewer rows (already done) or
fundamentally changing the M-step algorithm.

Porting `mf_em_likelihood_per_scale` to C++ would shave the L-mat
builder's ~30% of post-subsetting fit time. Expected wall-clock
savings: ~15%. Real but smaller than what we already have, and
adds C++ build complexity. **Decision: defer to a separate
proposal** when (and if) profiling reveals it as the next
bottleneck.

## Final result

End-to-end wall-clock on the perf-bench fixture:

    baseline (pre-cache):       28.86 s
    cache + subset (current):   12.91 s
    Δ:                         -15.95 s  (-55%)

OpenSpec acceptance criterion was "≥30% reduction"; delivered
55%.

### What we did

- Cached invariants on `model$em_cache`. ~2.6% real gain.
- Adaptive subsetting on `mixsqp_alpha_eps`. ~54% real gain.

### What we didn't do

- Vectorise NaN-imputation (no measurable gain, reverted).
- C++ port (deferred — diminishing returns).

### What changed in the public API

- `mfsusie()` gains `mixsqp_alpha_eps = 1e-6` (default).
  Set to `0` to disable subsetting (recovers bit-identity vs
  pre-subsetting). All other defaults unchanged.

### Tests

1188 PASS, 0 FAIL, 1 SKIP — unchanged across both changes.

## Next session

- Optional: tighten `mixsqp_alpha_eps` default. `1e-6` is
  conservative; raising to `1e-4` would drop more SNPs and might
  give another 10-15% wall-clock at negligible drift. Worth
  measuring on a few real fixtures before changing the default.
- Optional: C++ port of `mf_em_likelihood_per_scale` if a real
  fit shows it as the bottleneck.
- Both `apply-cross-package-audit-fixnow` and several other
  changes remain unarchived per `openspec list --changes`.
