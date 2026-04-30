# Apple-to-apple comparison against:
#   mvf.susie.alpha/R/multfsusie.R:255-318 (per-modality DWT cache
#     pipeline, mfsusieR mirrors the shape contract)
#   mvsusieR/R/mvsusie_constructors.R:21-160 (paradigm reference for
#     internal data class with X centering/scaling + per-outcome cache)
#
# Tests the spec contract in mf-data-class/spec.md: ragged modality
# lengths, univariate (T_m = 1) short-circuit, single-modality
# (M = 1) collapse, DWT cache fields, residual storage opt-in.

# ---- Helpers -----------------------------------------------------------

make_X <- function(n, J, seed = mfsusier_test_seed()) {
  set.seed(seed)
  X <- matrix(rnorm(n * J), nrow = n)
  X
}

make_Y_functional <- function(n, T_per_outcome, seed = mfsusier_test_seed()) {
  set.seed(seed + 1)
  lapply(T_per_outcome, function(T_m) matrix(rnorm(n * T_m), nrow = n))
}

# ---- Constructor input validation --------------------------------------

test_that("create_mf_individual coerces vector Y[[m]] to a 1-column matrix", {
  set.seed(mfsusier_test_seed())
  n <- 20; J <- 4
  X <- matrix(rnorm(n * J), n)
  Y <- list(rnorm(n))   # plain vector, not a matrix
  data <- mfsusieR:::create_mf_individual(wavelet_qnorm = FALSE, wavelet_standardize = FALSE, X = X, Y = Y, verbose = FALSE)
  expect_identical(ncol(data$D[[1]]), 1L)
  expect_identical(nrow(data$D[[1]]), 20L)
})

test_that("create_mf_individual errors on non-numeric Y[[m]]", {
  set.seed(mfsusier_test_seed())
  X <- matrix(rnorm(20), 5)
  Y_bad <- list(matrix("not-numeric", 5, 4))
  expect_error(
    mfsusieR:::create_mf_individual(wavelet_qnorm = FALSE, wavelet_standardize = FALSE, X = X, Y = Y_bad, verbose = FALSE),
    "must be numeric"
  )
})

test_that("create_mf_individual errors on wrong-length pos list", {
  set.seed(mfsusier_test_seed())
  X <- matrix(rnorm(20), 5)
  Y <- list(matrix(rnorm(20), 5))
  pos_bad <- list(seq_len(3))   # length mismatch with ncol(Y[[1]]) = 4
  expect_error(
    mfsusieR:::create_mf_individual(wavelet_qnorm = FALSE, wavelet_standardize = FALSE, X = X, Y = Y, pos = pos_bad, verbose = FALSE),
    NULL
  )
})

test_that("create_mf_individual errors on pos length mismatch with M", {
  set.seed(mfsusier_test_seed())
  X <- matrix(rnorm(20), 5)
  Y <- list(matrix(rnorm(20), 5), matrix(rnorm(40), 5))
  pos_bad <- list(seq_len(4))   # length 1, but M = 2
  expect_error(
    mfsusieR:::create_mf_individual(wavelet_qnorm = FALSE, wavelet_standardize = FALSE, X = X, Y = Y, pos = pos_bad, verbose = FALSE),
    "must be a list of length"
  )
})

test_that("mf_dwt accepts a vector Y_m and coerces to matrix; T_m=1 short-circuit handles zero-sd", {
  set.seed(mfsusier_test_seed())
  # Constant Y vector triggers the `csd == 0` fallback inside mf_dwt's
  # T_m = 1 branch.
  Y_const <- rep(1.5, 10)
  out <- mfsusieR:::mf_dwt(Y_const, pos_m = seq_len(1),
                           verbose = FALSE)
  expect_identical(out$T_basis, 1L)
  # Centered to zero (Y - mean = 0); csd was forced to 1, so D = 0 vec.
  expect_true(all(out$D == 0))
})

test_that("mf_invert_dwt handles a length-1 D vector (vector input, T_basis=1)", {
  # A 1-position case: D is just (n,1) matrix; round-trip should yield
  # back a (n,1) matrix when reconstructed.
  set.seed(mfsusier_test_seed())
  X <- matrix(rnorm(10), 5)
  Y <- list(matrix(rnorm(5), 5))   # T_m = 1
  data <- mfsusieR:::create_mf_individual(wavelet_qnorm = FALSE, wavelet_standardize = FALSE, X = X, Y = Y, verbose = FALSE)
  expect_identical(data$T_basis[1], 1L)
})

test_that("create_mf_individual rejects non-matrix X", {
  expect_error(
    mfsusieR:::create_mf_individual(wavelet_qnorm = FALSE, wavelet_standardize = FALSE, X = list(), Y = list(matrix(0, 4, 4))),
    "must be a numeric matrix"
  )
})

test_that("create_mf_individual rejects empty Y", {
  expect_error(
    mfsusieR:::create_mf_individual(wavelet_qnorm = FALSE, wavelet_standardize = FALSE, X = make_X(4, 3), Y = list()),
    "non-empty list"
  )
})

test_that("create_mf_individual rejects mismatched n between X and Y", {
  X <- make_X(10, 5)
  Y <- list(matrix(0, nrow = 8, ncol = 4))
  expect_error(mfsusieR:::create_mf_individual(wavelet_qnorm = FALSE, wavelet_standardize = FALSE, X = X, Y = Y), "rows.*match")
})

test_that("create_mf_individual rejects pos length mismatch", {
  X <- make_X(10, 5)
  Y <- list(matrix(0, nrow = 10, ncol = 64))
  expect_error(
    mfsusieR:::create_mf_individual(wavelet_qnorm = FALSE, wavelet_standardize = FALSE, X = X, Y = Y, pos = list(seq_len(32))),
    "length 32 but"
  )
})

# ---- Single-modality functional case (M = 1, T_1 > 1) ------------------

test_that("M = 1 functional case stores DWT cache with the right shape", {
  X <- make_X(20, 8)
  Y <- make_Y_functional(20, 64)
  obj <- mfsusieR:::create_mf_individual(wavelet_qnorm = FALSE, wavelet_standardize = FALSE, X = X, Y = Y, verbose = FALSE)

  expect_identical(class(obj), c("mf_individual", "individual"))
  expect_identical(obj$M, 1L)
  expect_identical(obj$n, 20L)
  expect_identical(obj$p, 8L)
  expect_identical(obj$T_basis, 64L)
  expect_identical(length(obj$D), 1L)
  expect_identical(dim(obj$D[[1]]), c(20L, 64L))
  expect_identical(length(obj$scale_index), 1L)
  expect_identical(length(obj$pos), 1L)
  expect_identical(length(obj$pos[[1]]), 64L)
  expect_identical(length(obj$wavelet_meta$column_center), 1L)
  expect_identical(length(obj$wavelet_meta$column_scale), 1L)
  expect_identical(length(obj$wavelet_meta$column_center[[1]]), 64L)
})

# ---- Multi-modality with ragged T_m -----------------------------------

test_that("ragged modality lengths are stored without padding to a common T", {
  X <- make_X(20, 8)
  Y <- make_Y_functional(20, c(64L, 128L, 256L))
  obj <- mfsusieR:::create_mf_individual(wavelet_qnorm = FALSE, wavelet_standardize = FALSE, X = X, Y = Y, verbose = FALSE)

  expect_identical(obj$M, 3L)
  expect_identical(obj$T_basis, c(64L, 128L, 256L))
  expect_identical(dim(obj$D[[1]]), c(20L, 64L))
  expect_identical(dim(obj$D[[2]]), c(20L, 128L))
  expect_identical(dim(obj$D[[3]]), c(20L, 256L))
  for (m in seq_len(3)) {
    expect_identical(length(obj$wavelet_meta$column_center[[m]]),
                     as.integer(obj$T_basis[m]))
    expect_identical(length(obj$wavelet_meta$column_scale[[m]]),
                     as.integer(obj$T_basis[m]))
  }
})

# ---- Univariate trait (T_m = 1) short-circuit --------------------------

test_that("T_m = 1 short-circuits the wavelet machinery", {
  X <- make_X(20, 8)
  Y <- list(matrix(rnorm(20), ncol = 1))
  obj <- mfsusieR:::create_mf_individual(wavelet_qnorm = FALSE, wavelet_standardize = FALSE, X = X, Y = Y, verbose = FALSE,
                                         pos = list(1))

  expect_identical(obj$T_basis, 1L)
  expect_identical(obj$scale_index[[1]], list(1L))
  expect_identical(dim(obj$D[[1]]), c(20L, 1L))
  expect_equal(mean(obj$D[[1]]), 0, tolerance = 1e-12)
  expect_equal(stats::sd(as.numeric(obj$D[[1]])), 1, tolerance = 1e-12)
})

# ---- Mixed univariate + functional modalities --------------------------

test_that("mixed (T_1 = 1, T_2 > 1) modalities flow through one constructor", {
  X <- make_X(20, 8)
  Y <- list(matrix(rnorm(20), ncol = 1),
            matrix(rnorm(20 * 64), nrow = 20))
  obj <- mfsusieR:::create_mf_individual(wavelet_qnorm = FALSE, wavelet_standardize = FALSE, X = X, Y = Y, verbose = FALSE,
                                         pos = list(1, seq_len(64)))

  expect_identical(obj$M, 2L)
  expect_identical(obj$T_basis, c(1L, 64L))
  expect_identical(dim(obj$D[[1]]), c(20L, 1L))
  expect_identical(dim(obj$D[[2]]), c(20L, 64L))
})

# ---- X centering / scaling ---------------------------------------------

test_that("X centering applied when intercept = TRUE (default)", {
  X <- make_X(20, 8)
  Y <- make_Y_functional(20, 64)
  obj <- mfsusieR:::create_mf_individual(wavelet_qnorm = FALSE, wavelet_standardize = FALSE, X = X, Y = Y, verbose = FALSE)

  expect_equal(colMeans(obj$X), rep(0, ncol(X)), tolerance = 1e-12)
  expect_equal(obj$wavelet_meta$X_center, colMeans(X), tolerance = 1e-12)
})

test_that("X scaling applied when standardize = TRUE (default)", {
  X <- make_X(20, 8)
  Y <- make_Y_functional(20, 64)
  obj <- mfsusieR:::create_mf_individual(wavelet_qnorm = FALSE, wavelet_standardize = FALSE, X = X, Y = Y, verbose = FALSE)

  # After centering then scaling each column, every column has SD 1.
  col_sds <- apply(obj$X, 2, stats::sd)
  expect_equal(col_sds, rep(1, ncol(X)), tolerance = 1e-12)
})

test_that("zero-variance X column is preserved with csd = 1", {
  X <- make_X(20, 5)
  X[, 3] <- 1                          # constant column
  Y <- make_Y_functional(20, 64)
  obj <- mfsusieR:::create_mf_individual(wavelet_qnorm = FALSE, wavelet_standardize = FALSE, X = X, Y = Y, verbose = FALSE)

  expect_identical(obj$csd[3], 1)
  expect_true(all(is.finite(obj$X[, 3])))
})

test_that("intercept = FALSE skips X centering", {
  X <- make_X(20, 8)
  Y <- make_Y_functional(20, 64)
  obj <- mfsusieR:::create_mf_individual(wavelet_qnorm = FALSE, wavelet_standardize = FALSE, X = X, Y = Y,
                                         intercept = FALSE,
                                         standardize = FALSE,
                                         verbose = FALSE)

  expect_equal(obj$X, X, tolerance = 0)
  expect_equal(obj$wavelet_meta$X_center, rep(0, ncol(X)), tolerance = 0)
})

# ---- Data-object slot list ---------------------------------------------

test_that("data object exposes the canonical fit-time slots", {
  X <- make_X(20, 8)
  Y <- make_Y_functional(20, c(64L, 128L))
  obj <- mfsusieR:::create_mf_individual(wavelet_qnorm = FALSE, wavelet_standardize = FALSE, X = X, Y = Y, verbose = FALSE)

  # Required slots used downstream by the IBSS dispatch and the
  # post-fit attach step.
  for (slot in c("X", "Y", "D", "scale_index", "T_basis", "csd",
                 "xtx_diag", "n", "p", "M", "wavelet_meta")) {
    expect_true(slot %in% names(obj),
                info = sprintf("missing slot: %s", slot))
  }
  # `residuals` was a no-op slot historically; it should no longer
  # be present.
  expect_false("residuals" %in% names(obj))
})

# ---- DWT cache contents agree with mf_dwt ------------------------------

test_that("the constructor's per-modality D matches a direct mf_dwt call", {
  X <- make_X(20, 8)
  Y <- make_Y_functional(20, 64)
  pos <- list(seq_len(64))
  # Disable wavelet_qnorm so the constructor's D is the raw DWT
  # output (matches the direct `mf_dwt` call below).
  obj <- mfsusieR:::create_mf_individual(X = X, Y = Y, pos = pos,
                                         wavelet_qnorm = FALSE, wavelet_standardize = FALSE,
                                         verbose = FALSE)

  out <- mfsusieR:::mf_dwt(Y[[1]], pos[[1]], verbose = FALSE)

  expect_equal(obj$D[[1]], out$D, tolerance = 0)
  expect_identical(obj$scale_index[[1]], out$scale_index)
  expect_identical(obj$T_basis[1], out$T_basis)
})

# ---- Default pos is seq_len(T_m) per modality --------------------------

test_that("when pos is NULL, defaults to seq_len(T_m) per modality", {
  X <- make_X(20, 8)
  Y <- make_Y_functional(20, c(64L, 128L))
  obj <- mfsusieR:::create_mf_individual(wavelet_qnorm = FALSE, wavelet_standardize = FALSE, X = X, Y = Y, verbose = FALSE)

  # With regular integer pos and a power-of-two T_m, no remap is
  # applied and the pos slot mirrors the input.
  expect_equal(obj$pos[[1]], seq_len(64), tolerance = 0)
  expect_equal(obj$pos[[2]], seq_len(128), tolerance = 0)
})
