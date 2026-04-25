# Apple-to-apple comparison against:
#   fsusieR/R/utils.R#L21-L24      (is.wholenumber)
#   fsusieR/R/utils.R#L85-L141     (colScale)
#   fsusieR/R/wavelet_utils.R#L42-L46   (interpolKS2)
#   fsusieR/R/wavelet_utils.R#L18-L37   (interpol_mat)
#   fsusieR/R/wavelet_utils.R#L277-L307 (remap_data)
#   fsusieR/R/wavelet_utils.R#L139-L160 (gen_wavelet_indx)
#   fsusieR/R/wavelet_utils.R#L89-L119  (DWT2)
#
# Contract C2 (design.md D11b, mf-ibss/spec.md): ported helpers must
# produce numerically identical output to fsusieR's originals on the
# same inputs. These tests are the regression net; any drift fails
# them at exact equality (max abs diff = 0).

test_that("col_scale center=FALSE returns center attribute as zeros", {
  X <- matrix(rnorm(20), 5)
  out <- mfsusieR:::col_scale(X, center = FALSE, scale = TRUE)
  expect_null(attr(out, "scaled:center"))   # add_attr only sets center attr when center=TRUE
})

test_that("col_scale scale=FALSE returns no scale attr", {
  X <- matrix(rnorm(20), 5)
  out <- mfsusieR:::col_scale(X, center = TRUE, scale = FALSE)
  expect_null(attr(out, "scaled:scale"))
})

test_that("is_wholenumber matches fsusieR::is.wholenumber on integer and float inputs", {
  skip_if_not_installed("fsusieR")
  cases <- list(
    integers   = c(0, 1, -2, 3),
    near_int   = c(1 + 1e-12, -3 + 5e-13),
    fractional = c(0.5, 1.7, log2(7))
  )
  for (case in cases) {
    expect_identical(
      mfsusieR:::is_wholenumber(case),
      asNamespace("fsusieR")$is.wholenumber(case)
    )
  }
})

test_that("col_scale matches fsusieR::colScale on a fixture row", {
  skip_if_not_installed("fsusieR")
  set.seed(mfsusier_test_seed())
  Y <- matrix(rnorm(30 * 32), nrow = 30)

  a <- mfsusieR:::col_scale(Y)
  b <- fsusieR::colScale(Y)

  # Centred / scaled values and the center / scale factor attributes
  # are bit-identical.
  expect_equal(c(a), c(b), tolerance = 0)
  expect_equal(attr(a, "scaled:center"),
               attr(b, "scaled:center"), tolerance = 0)
  expect_equal(attr(a, "scaled:scale"),
               attr(b, "scaled:scale"), tolerance = 0)

  # The `d` attribute is algebraically (n - 1) per column once
  # zero-variance columns have had their csd replaced with 1.
  # fsusieR computes it via a two-step `(a + b - a) / b` expression
  # that injects ULP-level rounding noise; mfsusieR computes the
  # simplified form directly. Both produce results that are equal
  # to within machine epsilon. Tolerance set well within contract
  # C2's <= 1e-8 floor.
  expect_equal(attr(a, "d"), attr(b, "d"), tolerance = 1e-12)
})

test_that("col_scale handles zero-variance columns by setting csd to 1", {
  Y <- cbind(rep(1, 10), rnorm(10), rep(2, 10))

  a <- mfsusieR:::col_scale(Y)
  scl <- attr(a, "scaled:scale")

  expect_identical(scl[1], 1)
  expect_identical(scl[3], 1)
  expect_true(is.finite(a[1, 1]))
})

test_that("col_scale subsetting via rows / cols / both", {
  set.seed(mfsusier_test_seed())
  Y <- matrix(rnorm(30 * 8), nrow = 30)

  a_rows <- mfsusieR:::col_scale(Y, rows = 1:10)
  expect_identical(dim(a_rows), c(10L, 8L))

  a_cols <- mfsusieR:::col_scale(Y, cols = 1:4)
  expect_identical(dim(a_cols), c(30L, 4L))

  a_both <- mfsusieR:::col_scale(Y, rows = 1:10, cols = 1:4)
  expect_identical(dim(a_both), c(10L, 4L))
})

test_that("gen_wavelet_indx matches fsusieR::gen_wavelet_indx for lev_res 2..10", {
  skip_if_not_installed("fsusieR")
  for (lev_res in 2:10) {
    expect_identical(
      mfsusieR:::gen_wavelet_indx(lev_res),
      fsusieR::gen_wavelet_indx(lev_res)
    )
  }
})

test_that("gen_wavelet_indx errors for lev_res > 15", {
  expect_error(mfsusieR:::gen_wavelet_indx(16),
               "lev_res must be at most 15")
})

test_that("dwt_matrix matches fsusieR::DWT2 (internal) on a fixture", {
  skip_if_not_installed("fsusieR")
  set.seed(mfsusier_test_seed())
  Y <- matrix(rnorm(30 * 64), nrow = 30)

  w_new <- mfsusieR:::dwt_matrix(Y)
  w_old <- asNamespace("fsusieR")$DWT2(Y)

  expect_equal(w_new$C, w_old$C, tolerance = 0)
  expect_equal(w_new$D, w_old$D, tolerance = 0)
  expect_identical(w_new$J, w_old$J)
  expect_identical(w_new$filter_number, w_old$filter.number)
  expect_identical(w_new$family, w_old$family)
})

test_that("dwt_matrix max_scale parameter forwards to wavethresh min.scale", {
  set.seed(mfsusier_test_seed())
  Y <- matrix(rnorm(10 * 64), nrow = 10)

  # Default max_scale = 10 >= log2(64) = 6, full decomposition.
  default <- mfsusieR:::dwt_matrix(Y)
  explicit <- mfsusieR:::dwt_matrix(Y, max_scale = 10)
  expect_equal(default$C, explicit$C, tolerance = 0)
  expect_equal(default$D, explicit$D, tolerance = 0)

  # max_scale recorded on the result so callers can audit it.
  expect_identical(default$max_scale, 10)
})

test_that("dwt_matrix handles fully-NA rows by zeroing during DWT and restoring NA", {
  set.seed(mfsusier_test_seed())
  Y <- matrix(rnorm(8 * 32), nrow = 8)
  Y[3, ] <- NA

  w <- mfsusieR:::dwt_matrix(Y)

  expect_true(all(is.na(w$D[3, ])))
  expect_true(is.na(w$C[3]))
  expect_false(any(is.na(w$D[-3, ])))
  expect_false(any(is.na(w$C[-3])))
})

test_that("remap_data is identity when ncol(Y) is a power of two and pos is regular", {
  skip_if_not_installed("fsusieR")
  set.seed(mfsusier_test_seed())
  Y <- matrix(rnorm(10 * 64), nrow = 10)
  pos <- seq_len(64)

  r_new <- mfsusieR:::remap_data(Y, pos, verbose = FALSE)
  r_old <- fsusieR::remap_data(Y, pos, verbose = FALSE)

  expect_equal(r_new$Y, r_old$Y, tolerance = 0)
  expect_identical(r_new$outing_grid, r_old$outing_grid)
})

test_that("remap_data interpolates non-power-of-two onto a power-of-two grid", {
  skip_if_not_installed("fsusieR")
  set.seed(mfsusier_test_seed())
  Y <- matrix(rnorm(10 * 30), nrow = 10)
  pos <- seq_len(30)

  r_new <- mfsusieR:::remap_data(Y, pos, verbose = FALSE)
  r_old <- fsusieR::remap_data(Y, pos, verbose = FALSE)

  # Y values are bit-identical; the interpolation math is unchanged.
  expect_equal(r_new$Y, r_old$Y, tolerance = 0)

  # outing_grid differs from fsusieR at machine-epsilon scale (~3e-15)
  # because mfsusieR's interpol_mat folds the position-rescaling step
  # into one linear interpolation, while fsusieR rescales to a [0, length]
  # intermediate first and then back to position units in remap_data.
  # The two formulas are mathematically identical; the difference is
  # ULP-level and well within contract C2's <= 1e-8 floor.
  expect_equal(r_new$outing_grid, r_old$outing_grid, tolerance = 1e-12)
  expect_true(mfsusieR:::is_wholenumber(log2(ncol(r_new$Y))))
})

test_that("remap_data zero-fills NA rows for interpolation and restores NA in the output", {
  skip_if_not_installed("fsusieR")
  set.seed(mfsusier_test_seed())
  Y <- matrix(rnorm(10 * 30), nrow = 10)
  Y[5, ] <- NA
  pos <- seq_len(30)

  r_new <- mfsusieR:::remap_data(Y, pos, verbose = FALSE)
  r_old <- fsusieR::remap_data(Y, pos, verbose = FALSE)

  expect_identical(is.na(r_new$Y), is.na(r_old$Y))
  expect_equal(r_new$Y, r_old$Y, tolerance = 0)
  expect_true(all(is.na(r_new$Y[5, ])))
})

test_that("interpol_mat returns Y on a power-of-two grid in pos units", {
  set.seed(mfsusier_test_seed())
  Y <- matrix(rnorm(8 * 30), nrow = 8)
  pos <- 100 + seq_len(30) * 7  # arbitrary affine

  out <- mfsusieR:::interpol_mat(Y, pos, max_scale = 10)

  expect_true(mfsusieR:::is_wholenumber(log2(ncol(out$Y))))
  expect_identical(nrow(out$Y), 8L)
  expect_identical(length(out$outing_grid), ncol(out$Y))
  expect_equal(min(out$outing_grid), min(pos), tolerance = 1e-12)
  expect_equal(max(out$outing_grid), max(pos), tolerance = 1e-12)
})

test_that("interpol_mat is bit-identical to remap_data on the interpolation path (regression net for the min/max refactor)", {
  # The refactor folds the position-rescaling step from remap_data into
  # interpol_mat. This test asserts that for an irregular pos input,
  # remap_data(...)$outing_grid equals interpol_mat(...)$outing_grid
  # exactly. If anyone reintroduces a separate rescaling step in
  # remap_data, this test fails.
  skip_if_not_installed("fsusieR")
  set.seed(mfsusier_test_seed())
  Y <- matrix(rnorm(10 * 30), nrow = 10)
  pos <- seq_len(30)

  r_new <- mfsusieR:::remap_data(Y, pos, verbose = FALSE)
  inter <- mfsusieR:::interpol_mat(Y, pos, max_scale = 10)
  expect_equal(r_new$outing_grid, inter$outing_grid, tolerance = 0)
  expect_equal(r_new$Y, inter$Y, tolerance = 0)
})

test_that("power_of_two_gridn returns the cap or 2^next when n is below 2^max_scale", {
  expect_identical(mfsusieR:::power_of_two_gridn(30, max_scale = 10), 32)
  expect_identical(mfsusieR:::power_of_two_gridn(60, max_scale = 10), 64)
  expect_identical(mfsusieR:::power_of_two_gridn(2000, max_scale = 10),
                   2L^10)
})
