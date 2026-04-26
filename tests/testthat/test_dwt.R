# Apple-to-apple comparison against:
#   fsusieR/R/susiF_workhorse.R DWT pipeline (col_scale + DWT2 +
#     gen_wavelet_indx, packed as cbind(D, C))
#   mvf.susie.alpha/R/multfsusie.R:308-316 (per-modality DWT
#     pipeline, same shape)
#
# Contract C2 / C3 (design.md D11b/c) at the per-modality wavelet
# representation. The packed `D` matrix mfsusieR caches on
# `mf_individual` SHALL match the row-by-row output of fsusieR's
# DWT2 (after the same col_scale step) at machine epsilon.

test_that("mf_dwt forward shape on a power-of-two functional input", {
  set.seed(mfsusier_test_seed())
  Y_m <- matrix(rnorm(20 * 64), nrow = 20)
  pos <- seq_len(64)

  out <- mfsusieR:::mf_dwt(Y_m, pos, max_padded_log2 = 10)

  expect_identical(out$T_basis, 64L)
  expect_identical(dim(out$D), c(20L, 64L))
  expect_identical(length(out$pos), 64L)
  expect_identical(length(out$column_center), 64L)
  expect_identical(length(out$column_scale), 64L)
  expect_identical(out$family, "DaubLeAsymm")
  expect_identical(out$filter_number, 10)
})

test_that("mf_dwt univariate short-circuit (T_m = 1)", {
  set.seed(mfsusier_test_seed())
  y <- matrix(rnorm(20), ncol = 1)

  out <- mfsusieR:::mf_dwt(y, pos_m = 1)

  expect_identical(out$T_basis, 1L)
  expect_identical(out$scale_index, list(1L))
  expect_identical(dim(out$D), c(20L, 1L))
  # The "wavelet coefficient" is just the centered+scaled column.
  expect_equal(mean(out$D), 0, tolerance = 1e-12)
  expect_equal(stats::sd(as.numeric(out$D)), 1, tolerance = 1e-12)
})

test_that("mf_dwt + mf_invert_dwt is identity to within wd/wr precision", {
  set.seed(mfsusier_test_seed())
  Y_m <- matrix(rnorm(20 * 64), nrow = 20)
  pos <- seq_len(64)

  out <- mfsusieR:::mf_dwt(Y_m, pos)
  recon <- mfsusieR:::mf_invert_dwt(
    out$D, out$column_center, out$column_scale,
    filter_number = out$filter_number,
    family = out$family
  )

  # wavethresh wd/wr roundtrip is ~2e-9 for filter.number = 10
  # (Daubechies-10); the additional col_scale inversion stays well
  # within the contract C2 floor of 1e-8. The inherent noise comes
  # from float accumulation through the long filter, NOT our port;
  # the next test uses filter_number = 1 (Haar wavelet) to confirm
  # bit-identity in the absence of long-filter accumulation.
  expect_equal(recon, Y_m, tolerance = 1e-8)
})

test_that("mf_dwt + mf_invert_dwt is bit-identical with filter_number = 1 (Haar wavelet)", {
  # Haar wavelet (filter.number = 1, family = "DaubExPhase" or
  # "DaubLeAsymm") has filter length 2; wd/wr does not accumulate
  # roundoff through the pyramid the way Daubechies-10 does. Running
  # the same forward+inverse pipeline with the Haar wavelet
  # demonstrates that the per-row reconstruction skeleton in
  # mf_invert_dwt is exact at machine precision and that the ~2e-9
  # noise observed with filter_number = 10 is inherent to wavethresh's
  # long-filter wd/wr, not a bug in our port. Production code uses
  # filter_number = 10 by default; this test pins the edge case.
  set.seed(mfsusier_test_seed())
  Y_m <- matrix(rnorm(20 * 64), nrow = 20)
  pos <- seq_len(64)

  out <- mfsusieR:::mf_dwt(Y_m, pos,
                           filter_number = 1,
                           family        = "DaubExPhase")
  recon <- mfsusieR:::mf_invert_dwt(
    out$D, out$column_center, out$column_scale,
    filter_number = out$filter_number,
    family        = out$family
  )

  expect_equal(recon, Y_m, tolerance = 1e-12)
})

test_that("mf_dwt + mf_invert_dwt is bit-identical for univariate (T_m = 1)", {
  set.seed(mfsusier_test_seed())
  y <- matrix(rnorm(20), ncol = 1)

  out <- mfsusieR:::mf_dwt(y, pos_m = 1)
  recon <- mfsusieR:::mf_invert_dwt(out$D, out$column_center, out$column_scale)

  # Univariate path skips the wavelet machinery entirely; the
  # inverse is a single multiply-add and SHALL be bit-identical
  # up to machine epsilon.
  expect_equal(recon, y, tolerance = 1e-12)
})

test_that("mf_dwt forward output matches fsusieR's susiF DWT pipeline (C2 fidelity)", {
  skip_if_not_installed("fsusieR")
  set.seed(mfsusier_test_seed())
  Y_m <- matrix(rnorm(20 * 64), nrow = 20)
  pos <- seq_len(64)

  # mfsusieR's path.
  out <- mfsusieR:::mf_dwt(Y_m, pos, max_padded_log2 = 10)

  # fsusieR's susiF pipeline does colScale (centering+scaling) then
  # DWT2 then cbind(D, C). Reproduce that here from fsusieR's
  # exported helpers.
  Y_scaled <- fsusieR::colScale(Y_m)
  W <- asNamespace("fsusieR")$DWT2(Y_scaled,
                                    filter.number = 10,
                                    family        = "DaubLeAsymm")
  D_ref <- cbind(W$D, W$C)

  # mfsusieR's `dwt_matrix` uses `min.scale = max_scale = 10`
  # internally; fsusieR's `DWT2` uses `min.scale = 10`. Same value,
  # bit-identical wavelet output.
  expect_equal(out$D, D_ref, tolerance = 0)
  expect_identical(out$scale_index,
                   fsusieR::gen_wavelet_indx(log2(ncol(out$D))))

  # Per-modality metadata that downstream IBSS / inverse-DWT code
  # depends on. These must also be bit-identical (or ULP-level under
  # Pattern B for the per-position centering / scaling computed via
  # `col_scale`, which the d-attribute test in test_utils_wavelet
  # already pins).
  expect_identical(out$T_basis, ncol(Y_m))
  expect_equal(out$pos, pos, tolerance = 0)  # power-of-two regular
  expect_equal(out$column_center,
               attr(Y_scaled, "scaled:center"), tolerance = 0)
  expect_equal(out$column_scale,
               attr(Y_scaled, "scaled:scale"), tolerance = 0)
  expect_identical(out$family, "DaubLeAsymm")
  expect_identical(out$filter_number, 10)
})

test_that("mf_dwt forward output matches mvf.susie.alpha::DWT2 row-wise (C3 fidelity)", {
  skip_if_no_mvf_alpha()
  set.seed(mfsusier_test_seed())
  Y_m <- matrix(rnorm(20 * 64), nrow = 20)
  pos <- seq_len(64)

  # mvf.susie.alpha does the same pipeline (colScale + DWT2 +
  # cbind(D, C)) per modality.
  out <- mfsusieR:::mf_dwt(Y_m, pos, max_padded_log2 = 10)

  Y_scaled <- mfsusieR:::col_scale(Y_m)  # bit-identical to fsusieR::colScale
  W_ref <- asNamespace("mvf.susie.alpha")$DWT2(Y_scaled,
                                                filter.number = 10,
                                                family        = "DaubLeAsymm")
  D_ref <- cbind(W_ref$D, W_ref$C)

  expect_equal(out$D, D_ref, tolerance = 0)

  # Per-modality metadata mirrors the C2 contract checks above.
  expect_identical(out$T_basis, ncol(Y_m))
  expect_equal(out$pos, pos, tolerance = 0)
  expect_equal(out$column_center,
               attr(Y_scaled, "scaled:center"), tolerance = 0)
  expect_equal(out$column_scale,
               attr(Y_scaled, "scaled:scale"), tolerance = 0)
})

test_that("mf_dwt forwards verbose to remap_data only when interpolation is needed", {
  set.seed(mfsusier_test_seed())
  Y_m <- matrix(rnorm(20 * 30), nrow = 20)  # not a power of two
  pos <- seq_len(30)

  expect_message(
    mfsusieR:::mf_dwt(Y_m, pos, max_padded_log2 = 10, verbose = TRUE),
    "interpolated to a regular dyadic grid"
  )

  Y_p2 <- matrix(rnorm(20 * 64), nrow = 20)
  pos_p2 <- seq_len(64)
  expect_silent(mfsusieR:::mf_dwt(Y_p2, pos_p2, verbose = TRUE))
})
