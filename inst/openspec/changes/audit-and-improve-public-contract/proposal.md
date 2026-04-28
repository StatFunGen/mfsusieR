# Audit and improve mfsusieR's public contract

## Why

The core mfsusieR algorithm is in place and the perf work has
landed. Before broader use, the package needs a cleanup pass:
the fit object has accumulated fields, several S3 methods may
duplicate logic that `susieR`'s default already provides,
verbose output is sparse, the HMM credible band is suppressed
on a hunch that it is too wide, and a representative scan
case (n=84, p≈3500, M=6, T=128) reportedly times out at 5
minutes. The recent perf work (cache + subsetting + cpp11)
likely changed the timeout picture; we have not re-measured.

This change captures the audit-and-improve work as one
coordinated effort. Implementation will split into focused
commits per section.

## Scope (seven sections)

### Section 1: Diagnosis-field cleanup

Almost all diagnosis fields on the `mfsusie` fit are
shape-compatible with `susieR` (`alpha`, `lbf`,
`lbf_variable`, `KL`, `elbo`, `niter`, `converged`,
`pip`, `sets`). Two specific items are worth a small touch:

- `V` is held at 1 across all effects because the mixture
  weights in `pi_V` carry the per-effect prior adaptation.
  The field is uninformative as a diagnostic. Decision:
  retire `V` from the fit and surface a summary of `pi_V`
  via `summary.mfsusie()` instead.
- `sigma2` shape (`list[M]` of length-`S_m` or scalar) is
  richer than susieR's scalar. No change to the field; just
  document the shape clearly.

The rest of the fit object stays as is.

### Section 2: Feature gap with `fsusieR` and `mvf.susie.alpha`

Read every exported function and parameter of `fsusieR` and
`mvf.susie.alpha`. List unported features in a session note.
For each, classify as `would-port-if-asked`, `out-of-scope`,
or `port-now`. Pure investigation; no code changes here.

### Section 3: Maximize susieR backbone usage (delete-or-patch discipline)

For every S3 method registered for `mf_individual` or
`mfsusie`, audit the body and classify into three buckets:

- **delete-and-inherit**: the body is a trivial wrapper or
  reproduces susieR's default. Drop the override; let
  susieR's default fire.
- **patch-susieR-and-delete**: the body diverges from susieR's
  default in one small place that susieR could easily
  parameterize via a hook or option. Open a focused susieR
  PR adding the hook; once landed, delete the override.
- **keep**: the body has real, irreducible mfsusieR-specific
  logic. Document the divergence in a one-line comment in
  the override.

The first instance of `delete-and-inherit` is
`check_convergence.mf_individual` (Section 4a). Likely
candidates for `delete-and-inherit` from a quick scan:
`get_cs.mf_individual`, possibly
`compute_residuals.mf_individual` or
`cleanup_model.mf_individual`.

A likely `patch-susieR-and-delete` candidate is the
verbose-output extension in Section 4: susieR's
`check_convergence.default` builds a fixed-column tabular
format. To add mfsusieR-specific diagnostic columns (max
`pi_null` across (m, s), max `KL_l`, n_eff from alpha
entropy) without reimplementing the entire formatter, the
right move is a small susieR patch exposing a per-class
`format_extra_diag(model)` generic that defaults to `""` and
that mfsusieR overrides.

For each override touched, record the decision in
`inst/notes/sessions/<date>-s3-override-audit.md`. Anything
found in `susieR` that looks miscalibrated or redundant
during this scan is flagged in the same note for separate
upstream-PR discussion.

### Section 4: Verbose output, susieR-arg parity, parameter renames

`mfsusie()` overrides `check_convergence.mf_individual` with a
trivial body (`R/ibss_methods.R:227`) that hardcodes ELBO-only
convergence and prints nothing per iteration. This bypasses
susieR's rich `check_convergence.default` which already
provides:

- per-iter tabular output (iter, ELBO, delta, sigma2, mem, V)
- monotone-ELBO check with warning on decrease
- PIP-difference convergence path with stall detection
- selectable `convergence_method = c("elbo", "pip")`
- `pip_stall_window` for non-monotone PIP detection
- mem-used reporting via `mem_used_gb()`

Plan:

1. Delete `check_convergence.mf_individual`. susieR's default
   fires.
2. Expose `convergence_method` and `pip_stall_window` on
   `mfsusie()`. Default `"elbo"` (matches susieR).
3. Expose `estimate_residual_variance`. Default `TRUE`.
4. Rename internal parameter `mixture_weight_method = c("mixsqp",
   "none")` to `estimate_prior_variance = c(TRUE, FALSE)`.
   `TRUE` -> mixsqp (current default); `FALSE` -> fix pi at
   the init values (matches susieR's
   `estimate_prior_variance = FALSE` semantics for the mixture
   prior). Keeps `mixture_weight_method` working for one
   minor version with a `lifecycle::deprecate_warn`.
5. Rename `lbf_min` to `greed_lbf_cutoff` for clarity (it
   gates the L-greedy outer loop). Same deprecation pattern.

Verification: a unit test that asserts mfsusieR's default
verbose output, on a small fixture, matches the structure of
susieR's verbose output line-for-line (modulo numeric values).
Plus the PIP-vs-ELBO convergence agreement test.

Per-iter mfsusieR-specific diagnostic columns
(`max_pi_null` across (m, s), `max_KL_l`, `n_eff` from alpha
entropy) ride on a susieR patch in Section 3. If the susieR
patch lands, mfsusieR overrides the new generic to inject
those columns; if the patch is rejected or delayed, the
mfsusieR-specific columns are not shipped (no
reimplementing the full row formatter just to add three
columns).

`get_cs.mf_individual` (`ibss_methods.R:337-343`) is also a
thin wrapper around `susie_get_cs`. Audit during this section:
keep if our wrapper customization is doing real work, drop
otherwise.

Out of scope for this section:
- `compute_univariate_zscore` (not useful for our typical
  workload).
- Other susieR args beyond the four named above.

### Section 5: HMM credible band — drop gating, document

Per the `2026-04-28-audit-findings.md` audit (B3): mfsusieR's
HMM band formula is mathematically correct. It uses the
proper law of total variance for the wavelet posterior
(`var_w = mu2_w - mean_w^2`), inverse-DWT projects to
position-space via `W_inv^2`, and uses `qnorm((1+level)/2)`
for the requested coverage. The "wide band" complaint that
suppressed display was a false signal — it conflated
mfsusieR's band with fsusieR's bug-affected DWT-band path
(which aggregates SDs linearly under alpha rather than
variances, then uses a 3-sigma multiplier instead of `1-α`).

Action reduced to: drop the gating that hides the HMM band,
add a one-paragraph roxygen note explaining the band
semantics, no formula change.

### Section 6: Performance and convergence on the heavy fixture

Heavy-fixture re-measurement (`n=84, p≈3500, M=6, T=128, L=10`)
**still exceeds 10 minutes** after the cache + subsetting +
cpp11 perf work landed earlier. The recent work optimized
the M-step path (mixsqp input shrunk via
`mixsqp_alpha_eps`); the SER step path that scales linearly
in `p` was untouched.

Per-effect-per-iter work is dominated by
`compute_residuals.mf_individual` and
`compute_ser_statistics.mf_individual`, which build
`X^T R` and `bhat = X^T R / xtx_diag` across all `p` SNPs
every effect, every IBSS iteration. For this fixture that
is `~2.3B` FP ops per IBSS iter, all in pure-R matrix code.

Plan:

1. **Profile** with `profvis` on this fixture (wrap
   `mfsusie()` call, save flamegraph). Confirm SER-step
   matrix multiplies are the hot path; if not, redirect.
2. **cpp11 port of the per-outcome bhat/shat builder**:
   `mf_per_outcome_bhat_shat()` (`R/individual_data_methods.R`)
   does the `X^T R / xtx_diag` and per-(scale, outcome)
   sigma broadcast. Tight loop in cpp11armadillo or via the
   `compute_marginal_bhat_shat` susieR helper if compatible.
   Numerical identity at `tol = 1e-12`.
3. **Convergence-side wins**: re-measure after exposing
   `convergence_method = "pip"` (Section 4). If PIP-based
   stopping converges in fewer iterations on this fixture
   shape, the overall fit time drops without further code.
4. Investigate the 50-iteration non-convergence warning on
   the test fixture (a separate issue from the heavy fixture
   above): is it a pathological fixture or a real algorithmic
   issue?
5. Re-measure heavy fixture after each of (2) and (3).
   Target: under 5 minutes for the heavy fixture, ideally
   under 2.

### Section 7: Comment polish and helper promotion

Two clean-up jobs:
- Audit and remove British-English spellings; align on
  American English everywhere (R, comments, docs, vignettes).
- Audit comments and roxygen for any reference to upstream
  packages or "port-fidelity" language. Comments must be
  self-contained: a reader who has never heard of fsusieR or
  mvf.susie.alpha should fully understand them. References
  to other packages live in `inst/notes/` or test file
  headers, not in `R/`.
- Audit internal helpers prefixed with `.`; promote (drop
  the `.`) those whose generality justifies a stable internal
  name. Pattern set by `univariate_smash_regression`. Do not
  export; just rename internally so future contributors can
  find them.

## Impact

- Fit object: smaller and more documented. Existing user code
  reading well-known fields (`pip`, `sets`, `alpha`, `mu`,
  `coef`) is unaffected. Code reading mfsusieR-specific
  internals may break; documented at the cleanup.
- S3 surface: smaller. No behavior change.
- Verbose: existing users see new per-iteration output when
  `verbose = TRUE`. Default unchanged.
- HMM: band displayed when populated; visualization may show
  uncertainty more honestly than before.
- Perf: targeted only at fixtures that show real headroom; no
  regressions tolerated on existing fixtures.
- Comments / naming: pure hygiene, no behavior change.

## Out of scope

- New algorithmic features (mfsusieR's algorithm is locked).
- Renaming user-facing arguments of `mfsusie()` or `fsusie()`
  (Phase 8 polish, separate change).
- Cross-package work in `susieR` (we list candidates; the
  upstream PRs are separate).
- Any C++ port beyond what section 6 measurably justifies.
