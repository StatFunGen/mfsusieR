## ADDED Requirements

### Requirement: PIP is a function of the post-filter alpha state

`fit$pip` SHALL be computed from the alpha matrix only after credible
sets have been constructed and any user-requested CS filtering has been
applied. The PIP formula SHALL be `pip_j = 1 - prod_l (1 - alpha[l, j])`
where `alpha` is the post-filter matrix, matching
`methods/online_method.tex` PIP definition (line 144).

This requirement directly addresses the strongest FDR-bug candidate
identified in `inst/notes/paradigms/mvf-original.md` section 8: the
original code computes PIP at
`R/operation_on_multfsusie_obj.R:1447` and applies CS filtering via
`discard_cs` at `:1452-1454` *afterwards*, leaving stale PIPs in the
returned object.

#### Scenario: filtered effects do not contribute to PIP

- **WHEN** `mfsusie()` runs with `filter_credible_sets = TRUE` on a
  scenario where one of the L effects is dropped by purity filtering
- **THEN** `fit$pip[j]` SHALL equal `1 - prod_{l in retained} (1 -
  alpha[l, j])` and SHALL NOT include any contribution from the
  dropped effect

#### Scenario: no filtering preserves the marginal definition

- **WHEN** `mfsusie()` runs with `filter_credible_sets = FALSE`
- **THEN** `fit$pip[j]` SHALL equal `1 - prod_l (1 - alpha[l, j])` over
  all L effects

### Requirement: CS construction order

`ibss_finalize.mf_individual` SHALL execute the following ordered steps
in this exact sequence, each step strictly following the previous:

1. Finalize alpha and posterior moments after the last IBSS sweep.
2. Compute credible sets at coverage `coverage` from the current alpha.
3. Apply the purity filter (`min_abs_corr`) and any user-requested
   filtering, dropping rows from `alpha` for filtered effects.
4. Compute `pip` from the post-filter `alpha`.
5. Attach `cs` and `pip` to the fit object.

#### Scenario: ordered execution

- **WHEN** the finalize routine is invoked
- **THEN** the sequence of internal calls SHALL match steps 1-5 above,
  and `pip` SHALL be computed at most once and only at step 4

#### Scenario: removed effects are reflected in alpha

- **WHEN** purity filtering removes one effect at step 3
- **THEN** `nrow(fit$alpha)` SHALL equal `L_initial - 1` and `length(fit$cs)`
  SHALL equal `L_initial - 1`

### Requirement: credible-set object structure

`fit$cs` SHALL be a list of credible-set entries, each carrying:

- `cs`: integer vector of variant indices into `1:J`.
- `coverage`: numeric scalar, the requested coverage achieved.
- `purity`: numeric scalar, the minimum absolute pairwise correlation
  in the CS, or `NA_real_` if not computed.
- `cs_index`: integer scalar, the original effect index `1..L_initial`
  before any filtering.
- `lbf`: numeric scalar, the log Bayes factor of the corresponding
  effect.

#### Scenario: CS object has documented fields

- **WHEN** `fit$cs[[1]]` is inspected after a successful run
- **THEN** the entry SHALL contain at minimum the fields `cs`,
  `coverage`, `purity`, `cs_index`, `lbf` with the types above

### Requirement: agreement with susieR get accessors

`susieR::susie_get_pip(fit)` and `susieR::susie_get_cs(fit, X = X)` SHALL
return values consistent with `fit$pip` and `fit$cs` respectively, given
that `class(fit) = c("mfsusie", "susie")`.

#### Scenario: get_pip equals fit$pip

- **WHEN** `susieR::susie_get_pip(fit)` is called
- **THEN** the returned vector SHALL equal `fit$pip` element-wise

#### Scenario: get_cs returns equivalent CS structure

- **WHEN** `susieR::susie_get_cs(fit, X = X)` is called
- **THEN** the returned list SHALL contain credible-set indices that
  match `fit$cs`, possibly with additional purity diagnostics
  computed from X

### Requirement: CS / PIP invariants hold for the single-modality functional case

mfsusieR SHALL apply the same CS construction order, PIP definition,
and CS object structure to the single-modality case (`M = 1`,
contract C2) as to the multi-modality case (`M >= 1`, contract C3).
The PIP-after-CS-filter ordering rule SHALL fix a known bug in
`fsusieR::susiF` where PIP is computed before CS filtering; the
ported mfsusieR path SHALL apply the filter first.

#### Scenario: M = 1 fits produce same-shape CS structure as M > 1

- **WHEN** `mfsusie()` is called with `M = 1, T_1 > 1` (or
  equivalently via `fsusie()`) and again with `M >= 2`, both with
  `filter_credible_sets = TRUE`
- **THEN** `fit$cs[[i]]` SHALL have the same fields (`cs`,
  `coverage`, `purity`, `cs_index`, `lbf`) and the same dtypes in
  both cases; `fit$pip` SHALL be a length-J numeric vector in both
  cases

#### Scenario: M = 1 PIP-after-CS-filter fix vs fsusieR::susiF

- **WHEN** `mfsusie(X, list(Y_matrix), list(pos), filter_credible_sets
  = TRUE)` is run on a fixture where one of the L effects is
  dropped by purity filtering, and `fsusieR::susiF(Y_matrix, X,
  pos, ...)` is run on the same fixture and seed
- **THEN** the two `pip` vectors SHALL differ at indices that were
  contributed by the dropped effect; the contract C2 test in
  `mf-ibss/spec.md` SHALL assert this deviation explicitly with a
  comment citing this OpenSpec change (the mfsusieR side is the
  fixed side; the test does not "expect" fsusieR's pre-filter PIP)
