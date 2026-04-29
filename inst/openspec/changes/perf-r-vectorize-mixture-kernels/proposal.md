# Replace 3 cpp mixture kernels with vectorized-R + cached log_dens

## Why

The three cpp kernels in `src/posterior_mixture.cpp` each
independently evaluate the same per-(j, t, k) mixture-component
log density:
`log dnorm(bhat[j,t]; 0, sqrt(shat[j,t]^2 + sd_k^2))`.

```
mixture_log_bf_per_scale_cpp(bhat, shat, sd_grid, pi, V)        # length-J
mixture_posterior_per_scale_cpp(bhat, shat, sd_grid, pi, V)     # (pmean, pmean2)
mf_em_log_likelihood_per_scale_cpp(bvec, sdmat, log_sdmat)      # log_L (N x K)
```

profvis at p=500 (mfsusieR per_outcome, post `perf-skip-elbo-on-pip`):

| function | self time |
|---|---|
| `mf_em_log_likelihood_per_scale_cpp` | 13.4 % |
| `mixture_log_bf_per_scale_cpp` | 8.9 % |
| `mixture_posterior_per_scale_cpp` | 1.8 % |
| **3-kernel total** | **24 %** |

Three independent passes through the same `(p × T × K)` log /
arithmetic. Same operands, identical math.

mvf and fsusieR achieve competitive speed without per-cell cpp:
both implement the equivalent computation in vectorized R + BLAS
(`fsusieR::cal_Bhat_Shat`, `fsusieR::log_BF.multfsusie_prior`).
The math is straightforward dnorm + softmax + closed-form
posterior — vectorizable in R using `outer`, broadcasting, and
column-wise `colSums` / `pmax` / `exp`. R + OpenBLAS matches or
beats hand-rolled cpp on this kind of dense element-wise
arithmetic, while saving:

* the per-call cpp ABI marshalling overhead (small but real, with
  `M*L*S` calls per outer iter on per_scale)
* triplicate per-cell density evaluation (3 cpp passes → 1 R pass)
* the cpp source + `cpp11.cpp` regen + `src/Makevars` plumbing
* a maintenance surface that has no demonstrated cpp benefit

The em_cache already maintains `sdmat[m, s]` and `log_sdmat[m, s]`
per (m, scope group); we extend the cache (or stash on
`ser_stats`) to include `log_dens[r, k]` and `log_dens_null[r]`
computed once per (l, m, scope group). All three downstream
consumers (mixsqp L, BF aggregation, posterior moments) read
from this cache.

## What changes

### New: `mf_em_log_dens_per_scale(bhat_slice, shat_slice, sdmat, log_sdmat)` (R)

Returns `list(log_dens, log_dens_null)` computed via vectorized
R using cached `sdmat` / `log_sdmat`:

```r
log_dens     <- -0.5 * log(2 * pi) - log_sdmat - 0.5 * (bvec / sdmat)^2
log_dens_null <- -0.5 * log(2 * pi) - log(shat_vec) - 0.5 * (bvec / shat_vec)^2
```

`bvec` is `as.numeric(bhat_slice)` (column-major flatten);
broadcasting handles the `(N x K)` result against `bvec` of
length N.

### Modified: `compute_ser_statistics.mf_individual`

After computing `betahat` / `shat2` per outcome, also computes
the per-(m, scope group) `log_dens` cache and stashes on
`ser_stats$log_dens_per_m`. Reads `sdmat` / `log_sdmat` from
`model$em_cache` (already populated).

### Modified: `mf_em_likelihood_per_scale` (em_helpers.R)

Becomes pure R. When `log_dens_cache` is supplied:

```r
L <- exp(log_dens_cache)        # (N x K), the mixsqp likelihood
# NA imputation in vectorised R (no apply()):
nas <- is.na(L)
if (any(nas)) {
  col_med <- matrixStats::colMedians(L, na.rm = TRUE)
  L[nas]  <- col_med[col(L)[nas]]
}
if (is_ebmvfr) L else rbind(c(100, rep(0, K - 1)), L)
```

The cpp call to `mf_em_log_likelihood_per_scale_cpp` is
eliminated.

### Modified: `loglik.mf_individual`

Aggregates BF in vectorized R from cached `log_dens`:

```r
# log_dens: (N=p*Ti) x K, log_dens_null: length-N
log_w  <- log_dens + matrix(log_pi, N, K, byrow = TRUE) -
          log_dens_null                                   # (N x K)
m_max  <- matrixStats::rowMaxs(log_w)
lbf_jt <- m_max + log(rowSums(exp(log_w - m_max)))         # length-N
# lbf_j = sum over t of lbf_jt
lbf_m_per_outcome <- rowSums(matrix(lbf_jt, J, Ti))       # length-J
```

Replaces `mixture_log_bf_per_scale_cpp`.

### Modified: `calculate_posterior_moments.mf_individual`

Aggregates posterior moments via the closed form using cached
`log_dens` and `sd_grid`:

```r
# Component-wise weights (same softmax as BF):
log_w  <- log_dens + log_pi_row
m_max  <- rowMaxs(log_w)
w      <- exp(log_w - m_max)
w      <- w / rowSums(w)                                   # (N x K)
# Per-component shrinkage (closed form):
shrink_k <- (sd_k^2 * V) / (sd_k^2 * V + shat^2)           # (N x K)
post_mean_k <- shrink_k * b                                # (N x K)
post_var_k  <- shrink_k * shat^2                           # (N x K)
pmean  <- rowSums(w * post_mean_k)
pmean2 <- rowSums(w * (post_var_k + post_mean_k^2))
```

Replaces `mixture_posterior_per_scale_cpp`.

### Deleted

* `src/posterior_mixture.cpp` — the three kernels
  (`mixture_log_bf_per_scale_cpp`,
  `mixture_posterior_per_scale_cpp`,
  `mf_em_log_likelihood_per_scale_cpp`).
* Their generated entries in `src/cpp11.cpp` and `R/cpp11.R`
  (regenerated automatically by `cpp11::cpp_register()`).

The Johnson-t correction kernel `mixture_log_bf_per_scale_johnson`
remains (it's R-only already, but uses a pre-computed `log_dens`
from cpp internally — adapt that to consume the new R-cached
`log_dens` directly).

## Acceptance criteria

* All 1326 existing tests pass.
* Numerical parity at `tol = 1e-12` against the pre-change
  cpp implementation: PIPs, alpha, mu, mu2, sigma2, ELBO (when
  `convergence_method = "elbo"`).
* On the n=84, p=3500, M=6, T=128, M=6 0/1/2-genotype fixture,
  default config, per-iter wall time decreases by an additional
  **≥ 10%** vs the post-`perf-skip-elbo-on-pip` baseline.
* On the apples-to-apples mvf comparison (n=100, p=2000, M=2),
  the cumulative speedup from
  `perf-skip-elbo-on-pip` + `perf-r-vectorize-mixture-kernels`
  puts mfsusieR strictly faster per iter than mvf.
* `src/` shrinks: `posterior_mixture.cpp` deleted entirely (no
  remaining cpp kernels for the per-(j, t, k) mixture work).

## Impact

* `src/posterior_mixture.cpp` deleted (~210 LOC removed).
* `src/cpp11.cpp` and `R/cpp11.R` regenerated (auto).
* `R/posterior_mixture.R` becomes the only implementation; the
  internal R reference oracles
  (`mixture_log_bf_per_scale_R`,
  `mixture_posterior_per_scale_R`) are renamed without the `_R`
  suffix and become the production code path.
* `R/em_helpers.R` modified (mf_em_likelihood_per_scale loses
  cpp call).
* `R/individual_data_methods.R` modified
  (`compute_ser_statistics.mf_individual`,
  `loglik.mf_individual`,
  `calculate_posterior_moments.mf_individual`).

## Out of scope

* Removing other cpp kernels in the repo (smoothing wavelet ops
  in `R/utils_wavelet.R`; those have demonstrated cpp benefit per
  the wavelet-cpp11 audit).
* Any algorithmic change to the mixture-EB prior structure.
* Restructuring the per-effect SER orchestrator beyond plumbing
  the new cache through `ser_stats`.

## Risk + mitigation

* R may be slower than cpp on platforms with un-optimized BLAS.
  Mitigation: the bench fixture uses pixi's openblas (the conda
  default); we verify `≤ cpp` speed on this baseline. If a user
  has reference BLAS only, the slowdown is bounded by ~1.5×
  (measured by independent profvis comparing R-only vs cpp
  paths on a simplified fixture).
* Edge cases (Bhat = 0, Shat → 0, NaN from numerical underflow)
  may behave differently in R vs cpp. Mitigation: 1e-12 parity
  test on the existing fixtures, plus dedicated tests for the
  degenerate-input cases that the cpp kernels handle.
* Memory: the `log_dens` cache is `(N × K)` per (l, m, scope
  group), cached on ser_stats per (l, m). For per_outcome at
  p=3500, T=128, K=10: 4.48M doubles per (l, m) = 36 MB. Across
  L=10, M=6 = 2.16 GB peak. **This is significantly more memory
  than current.** Mitigation: discard the cache as soon as the
  (l, m) SER step completes; do not retain across effects.
  Quantified peak addition: only the current (l, m)'s cache, so
  ~36 MB peak addition, not 2.16 GB.
