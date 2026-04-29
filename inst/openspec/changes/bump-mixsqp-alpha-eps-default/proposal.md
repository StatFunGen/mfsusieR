# Bump `mixsqp_alpha_eps` default `1e-6 → 5e-5`

## Why

`mixsqp_alpha_eps` is the alpha-thinning threshold in the
mixsqp M-step (`R/individual_data_methods.R:579`). Variants with
`model$alpha[l, j] < mixsqp_alpha_eps` are dropped from the
likelihood-matrix input. The current default `1e-6` was chosen
to keep the truncation error well below floating-point precision;
empirically the bound is

```
truncation error / total likelihood ≈ eps × (max L_jk / mean L_jk)
                                    ≈ eps × O(10)
```

so `eps = 1e-6` ⇒ ~1e-5 relative error in the EB-fitted `pi`
(negligible).

For typical SuSiE posteriors, alpha is concentrated on a handful
of variants. With `eps = 1e-6` we keep ~50 variants out of
`p = 1000`; with `eps = 5e-5` we keep ~15. The smaller eps
spends compute on variants whose alpha is between `0.0005%` and
`0.005%` of total — by any practical definition of an effect,
indistinguishable from null. The relative numerical error climbs
from ~1e-5 to ~5e-4, still 200× below the IBSS `tol = 1e-4`
default, so the fit is unchanged at any tolerance a user cares
about.

## What changes

In `R/mfsusie.R`, change one default:

```diff
-                    mixsqp_alpha_eps          = 1e-6,
+                    mixsqp_alpha_eps          = 5e-5,
```

Update the roxygen `@param mixsqp_alpha_eps` to reflect the new
default and the truncation-error bound.

## Acceptance criteria

- Bumping the default to `5e-5` does not visibly change the fit
  on the canonical fixtures (`why_functional`,
  `susie_post_outcome_configuration` end-to-end). Specifically:
  - `max|alpha_new - alpha_old| <= 1e-3` and
    `max|pip_new - pip_old| <= 1e-3` on at least three
    fixtures (small p, medium p, large p).
  - All existing unit tests pass at their existing tolerances.
- A micro-benchmark on a `p = 5000` fixture shows the M-step
  runtime drops by at least 30% with the new default vs the old.
- The roxygen for `mixsqp_alpha_eps` documents the new default
  and its truncation-error rationale.
- The pre-existing escape hatch `mixsqp_alpha_eps = 0` (use
  every variant regardless of alpha) keeps working.
