# Refactor exceptions ledger

This file records every line of the two port sources that mfsusieR
intentionally does NOT carry over. Doctrine lives in
`inst/openspec/changes/add-mfsusier-s3-architecture/design.md`
D11e. Reviewer pass on every Phase 3 numerical PR confirms the
ledger is up to date; PRs that omit lines without a corresponding
entry SHALL be blocked.

## Port sources in scope

- `mvf.susie.alpha/R/multfsusie.R` and supporting
  `operation_on_multfsusie_*.R`, `EM.R`, `ELBO_mutlfsusie.R`,
  `computational_routine.R`, `utils_wavelet_transform.R`, and the
  parts of `utils.R` and `utils_formatting.R` that the IBSS path
  touches.
- `fsusieR/R/susiF.R` and `susiF_workhorse.R` and supporting
  `operation_on_susiF_obj.R`, `wavelet_utils.R`, the parts of
  `operation_on_prior.R` that init the scale-mixture prior, the
  per-method smoothers in `operation_on_susiF_obj.R`
  (`TI_regression`, `smash_regression`, `HMM_regression`), and
  the parts of `utils.R` needed by the susiF path.

## Out of scope (file-level entries; no per-line walk required)

- fsusieR/R/EBmvFR.R
  Behavior: EB multivariate functional regression, a different
    model with no SuSiE structure.
  Decision: out-of-scope-EBmvFR
  Reason: design.md D5 and CLAUDE.md hard rule #2 record EBmvFR as
    an entirely separate algorithm. Not in v1; possibly addressed
    in a follow-up OpenSpec change after the SuSiE-track ports
    stabilize.
- fsusieR/R/EBmvFR_workhorse.R
  Behavior: Inner loop for EBmvFR.
  Decision: out-of-scope-EBmvFR
  Reason: As above.
- fsusieR/R/operation_on_EBmvFR_obj.R
  Behavior: Operations on the EBmvFR fit object.
  Decision: out-of-scope-EBmvFR
  Reason: As above.

## Entry schema for in-scope omissions

```
- <port_source>/R/<file>.R:<lo>-<hi>
  Behavior: <one-line description of what the original lines do>
  Decision: omit | replaced-by-<upstream> | deferred-to-<phase> |
            out-of-scope-EBmvFR
  Reason: <one-paragraph justification, citing OpenSpec change or
          manuscript section>
```

## In-scope omissions

### PR group 2 (R/utils_wavelet.R, 2026-04-25)

- fsusieR/R/utils.R:113-119
  Behavior: `for (col in zero_var_cols) { csd[col] <- 1 }`. Loop
    over zero-variance columns to replace `csd[col]` with 1.
  Decision: replaced-by-vectorized (col_scale uses
    `csd[csd == 0] <- 1`).
  Reason: Trivially equivalent to the vectorized form; no
    behavioural change. Cosmetic style cleanup logged here per
    D11e.

- fsusieR/R/utils.R:135-137
  Behavior: `d = n * cm^2 + (n - 1) * csd^2; d = (d - n * cm^2) /
    csd^2`. Two-step computation that algebraically simplifies to
    a constant `(n - 1)` per column since `csd[csd == 0] <- 1`
    was applied upstream.
  Decision: replaced-by-simplified-math (port computes
    `attr(x, "d") <- rep(nrow(x) - 1, ncol(x))` directly).
  Reason: The two-step form is wasteful compute and introduces
    ULP-level rounding noise via the `(a + b) - a` pattern. The
    simplification is correct math (algebraically equivalent),
    and the C2 fidelity test on the `d` attribute is relaxed to
    `tolerance = 1e-12` (well within contract C2's `<= 1e-8`
    floor) to accept the resulting machine-epsilon difference.
    The test file documents the fsusieR-bug nature inline. Net
    effect: cleaner mfsusieR code, identical observable behaviour
    on the `d` attribute up to machine epsilon.

- fsusieR/R/wavelet_utils.R:295-297
  Behavior: `verbose` message text uses raw newlines:
    `"Response matrix dimensions not equal to nx 2^J \n or
    unevenly spaced data \n interpolation procedure used"`.
  Decision: replaced-by-cleaner-message-text.
  Reason: The message text is user-facing diagnostic output, not
    part of any numerical contract. mfsusieR replaces with a
    semicolon-separated single line for readability. No effect on
    C2 / C3 contracts.

- fsusieR/R/wavelet_utils.R:148
  Behavior: `for (s in 1:(lev_res-1))` in `gen_wavelet_indx`.
    For `lev_res = 1` this expands to `1:0 = c(1, 0)`, executing
    the loop twice with `s = 1` and `s = 0`, producing a
    nonsensical result.
  Decision: replaced-by-seq_len (port uses `seq_len(lev_res - 1)`
    which yields `integer(0)` for `lev_res = 1`).
  Reason: `lev_res = 1` corresponds to a 2-sample input which is
    too small for a meaningful wavelet decomposition and is not
    part of the C2 contract floor (`T_1 in c(64, 128)` per
    design.md D11b). The fsusieR pattern is a real bug at this
    edge case; the port adopts the safe `seq_len` idiom. The C2
    test range (lev_res 6-7) is bit-identical between fsusieR and
    mfsusieR; the bug-fixed edge case is below the contract.

- fsusieR/R/wavelet_utils.R:18-37 + R/wavelet_utils.R:284-297
  Behavior: `interpol_mat` returns `grid` rescaled to a `[0,
    length(grid)]` integer-like scale, and `remap_data` then
    rescales it back to position units via
    `outing_grid <- start_pos + (end_pos - start_pos) /
    length(grid) * grid`.
  Decision: replaced-by-single-step (port folds the rescaling
    into `interpol_mat`, returning `outing_grid` in position
    units directly).
  Reason: The two-step composition is mathematically identical to
    the one-step linear interpolation, but the order of floating-
    point operations differs, producing a ULP-level difference
    (~3e-15) in `outing_grid`. Y output is bit-identical (max
    diff = 0). The reordering is well within contract C2's `<=
    1e-8` floor. The simplification eliminates redundant
    `min(pos)`, `max(pos)`, `min(raw_grid)`, `max(raw_grid)`
    passes and is locked by a regression test in
    `tests/testthat/test_utils_wavelet.R`.

- fsusieR/R/wavelet_utils.R:102-105 (DWT2's `min.scale = 10`)
  Behavior: `wavethresh::wd(..., min.scale = 10)` hardcoded inside
    the DWT helper.
  Decision: replaced-by-parameter-unification (port adds
    `max_scale = 10` argument to `dwt_matrix`, forwarded as
    wavethresh's `min.scale`).
  Reason: The `10` is the same conceptual quantity as `max_scale`
    elsewhere in `R/utils_wavelet.R`. Hardcoding it inside the
    DWT helper is a magic-number style issue the audit flagged.
    The port lifts it to a function-level parameter named
    consistently with the rest of the file. No behavioural change
    when `max_scale = 10` (the default); behaviour preserved.

### PR group 2 (R/dwt.R + R/data_class.R, 2026-04-25)

- mvf.susie.alpha/R/utils_wavelet_transform.R:24-50 (DWT2)
  Behavior: Row-wise wavelet transform of a matrix of curves.
    Identical in shape to `fsusieR::DWT2` but without the
    `min.scale = 10` argument; the wavethresh default applies
    (decompose to scale 0).
  Decision: replaced-by-`R/utils_wavelet.R::dwt_matrix` (already
    landed in PR group 2 task 2.2).
  Reason: A single canonical DWT helper is enough; mvf.susie.alpha
    and fsusieR each shipped a near-duplicate. The C3 fidelity
    test in `tests/testthat/test_dwt.R` ("matches mvf.susie.alpha
    DWT2 row-wise") verifies bit-identity on the same input
    (after the same `col_scale` step).

- mvf.susie.alpha/R/utils_wavelet_transform.R:64-68 (pack_dwt)
  Behavior: One-line helper computing `cbind(W$D, W$C)` where W
    is `DWT2(Y)`.
  Decision: replaced-by-inline (the `cbind(D, C)` pack step is
    inlined into `R/dwt.R::mf_dwt`).
  Reason: The packing convention is exposed in `mf_dwt` directly;
    a separate one-line helper is unnecessary.

- mvf.susie.alpha/R/multfsusie.R:255-318 (per-modality DWT pipeline
  inside `multfsusie()`)
  Behavior: Inline loop in the `multfsusie()` body that calls
    `fsusieR::remap_data` -> `fsusieR::colScale` -> `DWT2` ->
    `cbind(D, C)` -> `fsusieR::gen_wavelet_indx` per modality.
  Decision: replaced-by-`R/dwt.R::mf_dwt` (per-modality wrapper)
    invoked from `R/data_class.R::create_mf_individual` (per-fit
    constructor).
  Reason: The original buries the DWT pipeline inside the public
    fitting function; mfsusieR factors it into a reusable
    per-modality helper (`mf_dwt`) and a separate constructor
    (`create_mf_individual`). The pipeline order and arguments
    are bit-identical to mvf.susie.alpha's; the C3 DWT fidelity
    test verifies this on the matrix output.

- fsusieR/R/operation_on_susiF_obj.R:1785-1796 (inverse-DWT
  reconstruction skeleton)
  Behavior: Constructs a zero `wavethresh::wd` object, sets `$D`
    and the last entry of `$C`, applies `wavethresh::wr` to
    reconstruct the position-space curve. Used by `update_cal_fit_func`
    and downstream coef / predict.
  Decision: replaced-by-`R/dwt.R::mf_invert_dwt` (PR group 2 task
    2.4).
  Reason: Same skeleton, factored into a reusable inverse-DWT
    helper. The 2e-9 wd/wr roundtrip noise inherent to
    `filter.number = 10` is documented in
    `tests/testthat/test_dwt.R` ("mf_dwt + mf_invert_dwt is
    identity to within wd/wr precision"), with tolerance set to
    `1e-8` (the contract C2 floor).
