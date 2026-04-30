# Audit follow-up A-2: warm-start validator on `model_init` fields.

test_that("validate_model_init catches NA in alpha", {
  mi <- list(alpha = matrix(c(0.5, NA, 0.3, 0.2), 2, 2))
  expect_error(mfsusieR:::validate_model_init(mi), "alpha")
})

test_that("validate_model_init catches Inf in V", {
  mi <- list(alpha = matrix(0.25, 2, 4),
             V = c(1, Inf))
  expect_error(mfsusieR:::validate_model_init(mi), "V")
})

test_that("validate_model_init catches NA in mu list-of-list", {
  mi <- list(
    alpha = matrix(0.25, 2, 4),
    mu    = list(list(matrix(c(NA, 0, 0, 0), 2, 2)),
                 list(matrix(0, 2, 2))))
  expect_error(mfsusieR:::validate_model_init(mi), "mu")
})

test_that("validate_model_init catches mu length mismatch with alpha rows", {
  mi <- list(
    alpha = matrix(0.25, 2, 4),
    mu    = list(list(matrix(0, 2, 2))))
  expect_error(mfsusieR:::validate_model_init(mi), "alpha")
})

test_that("validate_model_init catches non-list mu", {
  mi <- list(
    alpha = matrix(0.25, 2, 4),
    mu    = matrix(0, 2, 4))
  expect_error(mfsusieR:::validate_model_init(mi), "list-of-list")
})

test_that("validate_model_init returns invisibly on a clean fixture", {
  mi <- list(
    alpha  = matrix(0.25, 2, 4),
    V      = c(1, 1),
    KL     = c(0, 0),
    sigma2 = c(0.5),
    mu     = list(list(matrix(0, 4, 8)),
                  list(matrix(0, 4, 8))),
    mu2    = list(list(matrix(0, 4, 8)),
                  list(matrix(0, 4, 8))))
  expect_silent(mfsusieR:::validate_model_init(mi))
})

test_that("expand_model_init_to_L errors on a corrupted model_init", {
  bad <- list(
    alpha = matrix(c(0.5, NA, 0.5, 0), 2, 2),
    mu    = list(list(matrix(0, 2, 2)), list(matrix(0, 2, 2))),
    mu2   = list(list(matrix(0, 2, 2)), list(matrix(0, 2, 2))))
  expect_error(
    mfsusieR:::expand_model_init_to_L(bad, L_new = 4L, p = 2L, M = 1L,
                                       T_basis = 2L),
    "alpha")
})
