# 2026-04-30 cross-package audit follow-ups (track-for-later bundle)

## Why

The 2026-04-30 cross-package audit (Agents A, B, C; see
`inst/openspec/changes/audit-cross-package-post-hooks/`)
surfaced five small polish items that are not blocking but
should be addressed before they decay into real bugs. Bundling
them under one OpenSpec change keeps the issue tracker simple;
each has its own task in section 2 and can land independently.

## Items

### A-2. Warm-start validator

mfsusieR's `ibss_initialize.mf_individual` does not call
`validate_init` (susieR's NA/Inf checker on warm-start fits)
because mfsusieR's `mu`/`mu2` are list-of-list and would crash
the `dim(mu) == dim(alpha)` check in `validate_init`. Net effect:
warm-start fits with NA/Inf in alpha or mismatched mu list
shapes fail later, deeper in the IBSS loop, with less helpful
errors. Add an mfsusieR-flavored validator alongside
`expand_model_init_to_L` covering the list-of-list shapes.

### A-3. Post-hook S3 dispatch cleanup

`post_loglik_prior_hook.mf_individual` calls
`optimize_prior_variance.mf_individual(...)` directly by name
(line 711) and passes a `moments = get_post(model, l)` argument
that the callee ignores. Also `loglik`,
`calculate_posterior_moments`, `compute_kl`,
`get_alpha_l`, `get_posterior_moments_l` are fetched via
`getFromNamespace("susieR")` (lines 703-707). Both work but bypass
the S3 dispatch table, making any future subclass override
silently inert. Replace direct calls with dispatched
`optimize_prior_variance(...)`. Drop the unused `moments` arg.

### B-5. CS purity / coverage parity smoke test

mfsusieR's CS construction routes through `susie_get_cs`
(susieR backbone); fsusieR uses its own
`update_cal_cs.susiF` + post-hoc purity filter. Both use the
same coverage and min-abs-corr threshold. The static read shows
the algebra agrees; observable parity depends on whether the
purity filter fires before or after CS aggregation in each path.
Add a C2 fidelity smoke test on a small fixture comparing
`fit$sets$cs` between mfsusieR (`fsusie(Y, X, pos)`) and
`fsusieR::susiF(Y, X, pos)` to lock in any per-CS rounding and
catch future drift.

### C-7. Combiner default method

`combine_outcome_lbfs` is a generic with one
`.mf_prior_cross_outcome_independent` method. Any future
combiner that subclasses `mf_prior_cross_outcome` without
inheriting from `_independent` will fall through to UseMethod
with no method and error opaquely. Add a
`combine_outcome_lbfs.default` method that errors with a helpful
message naming the registered combiner classes, OR document the
extension protocol in the file header so the seam's stated
purpose is defended.

### C-8. L_greedy expansion regression test

`expand_model_init_to_L` appends new effects with `alpha = 1/p`,
`mu = mu2 = 0`, and `pi_V`/`fitted_g_per_effect` seeded from
effect 1's prior state. Defensible warm-start convention but not
covered by a regression test. Add a test on L_greedy growth from
5 to 10 asserting the expanded fit converges within `1e-8` of a
cold-start fit at L=10.

## What changes

Each item gets its own task. PRs are independent; nothing in
this change blocks anything else.

## Impact

- **Severity**: low. Each item is a polish or defensive-
  programming improvement, none affect numerical correctness on
  the supported scenarios.
- **Source-code changes**: small per-item; total surface ~50
  lines of R + ~3 new test files.
- **Documentation**: optional NEWS entry per item.

## Out of scope

- Performance work.
- Anything outside the five items listed above.

## References

- Audit Agent A report (Findings A-2, A-3):
  `inst/notes/sessions/2026-04-30-cross-package-audit-a-susier-posthooks.md`
- Audit Agent B report (Finding B-5):
  `inst/notes/sessions/2026-04-30-cross-package-audit-b-fsusier-posthooks.md`
- Audit Agent C report (Findings C-7, C-8):
  `inst/notes/sessions/2026-04-30-cross-package-audit-c-mvf-posthooks.md`
- Audit summary:
  `inst/notes/cross-package-audit-summary-posthooks.md`
- Source change:
  `inst/openspec/changes/audit-cross-package-post-hooks/`
