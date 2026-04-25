# Smoke tests for the package skeleton (PR group 1).
#
# These tests verify that the test infrastructure itself works:
#   - the package loads,
#   - setup.R helpers are callable,
#   - the canonical minimal fixture loads with the documented shape.
#
# Numerical fidelity tests live in later PR groups; nothing here compares
# against `mvf.susie.alpha`.

test_that("package loads and namespace is wired", {
  expect_true("mfsusieR" %in% loadedNamespaces())
})

test_that("setup helpers are available and behave", {
  expect_identical(mfsusier_test_seed(), 1L)
  expect_type(skip_if_no_mvf_alpha, "closure")
})

test_that("minimal fixture loads with the documented shape", {
  fx <- mfsusier_load_fixture("scenario_minimal")

  expect_identical(fx$name, "scenario_minimal")
  expect_identical(fx$seed, 20260425L)
  expect_identical(fx$N, 200L)
  expect_identical(fx$J, 50L)
  expect_identical(fx$M, 2L)
  expect_identical(fx$T_per_modality, c(64L, 128L))
  expect_identical(fx$L_true, 3L)

  expect_true(is.matrix(fx$X))
  expect_identical(dim(fx$X), c(fx$N, fx$J))

  expect_type(fx$Y, "list")
  expect_named(fx$Y, c("Y_f", "Y_u"), ignore.order = TRUE)
  expect_length(fx$Y$Y_f, fx$M)
  for (m in seq_len(fx$M)) {
    expect_true(is.matrix(fx$Y$Y_f[[m]]))
    expect_identical(dim(fx$Y$Y_f[[m]]), c(fx$N, fx$T_per_modality[m]))
    expect_false(anyNA(fx$Y$Y_f[[m]]))
  }
  expect_null(fx$Y$Y_u)

  expect_length(fx$pos, fx$M)
  for (m in seq_len(fx$M)) {
    expect_identical(length(fx$pos[[m]]), as.integer(fx$T_per_modality[m]))
  }

  expect_length(fx$true_idx, fx$L_true)
  expect_true(all(fx$true_idx >= 1L & fx$true_idx <= fx$J))
})
