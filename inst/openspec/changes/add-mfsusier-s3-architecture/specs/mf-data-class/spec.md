## ADDED Requirements

### Requirement: mf_individual data class

The package SHALL provide an S3 class `mf_individual` that inherits from
`individual` and represents the full input to a multi-modality functional
SuSiE fit. The constructor `create_mf_individual(X, Y, pos, ...)` SHALL
accept `Y` as a list of M numeric matrices (one per modality), `pos` as a
list of M position vectors, and `X` as an `n x J` numeric matrix.

#### Scenario: ragged modality lengths

- **WHEN** Y has M modalities with different T_m values
  (e.g., M = 3, T = c(64, 128, 1024))
- **THEN** the constructor SHALL succeed and the resulting object SHALL
  store per-modality DWT caches, per-modality scale-index lists, and
  per-modality padded lengths without padding to a common T

#### Scenario: univariate trait collapses to T_m = 1

- **WHEN** any modality m has `T_m = 1` (a univariate trait)
- **THEN** the DWT pipeline for that modality SHALL short-circuit and
  store the data directly with no wavelet transform, the scale index
  list for that modality SHALL be `1L`, and the inverse-DWT helper SHALL
  pass the value through unchanged

#### Scenario: single modality collapses to M = 1

- **WHEN** Y is a list of length 1
- **THEN** the constructor SHALL succeed, all per-modality fields SHALL
  be length-1 lists, and downstream IBSS SHALL produce a fit equivalent
  to a single-modality functional regression

### Requirement: DWT cache stored at construction

The constructor SHALL perform forward discrete wavelet transform once at
construction and cache the resulting wavelet-coefficient matrices for
reuse across IBSS iterations. The IBSS loop SHALL NOT perform wavelet
transforms.

#### Scenario: cache fields exist after construction

- **WHEN** an `mf_individual` object is constructed from Y, X, pos
- **THEN** the object SHALL contain fields `D` (list of M matrices of
  wavelet coefficients), `scale_index` (list of M integer vectors),
  `T_padded` (integer vector of M padded lengths), `wavelet_meta` (list
  carrying `filter_number`, `family`, and per-modality column-scale
  factors)

### Requirement: shape uniformity for IBSS dispatch

`mf_individual` SHALL satisfy the same shape contract regardless of M
and T_m, so that the same S3 method implementations work for univariate,
single-modality functional, and multi-modality functional inputs.

#### Scenario: methods work uniformly across degenerate cases

- **WHEN** the same `compute_residuals.mf_individual` method is invoked
  on objects with `(M = 1, T_1 = 1)`, `(M = 1, T_1 = 64)`, and
  `(M = 6, T_m varies)`
- **THEN** the method SHALL return correctly-shaped residual outputs
  consistent with the per-modality storage shape, with no
  dimension-specific code branches in the public method body
