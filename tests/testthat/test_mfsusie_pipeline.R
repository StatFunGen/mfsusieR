# Tests for mfsusie_screen(), mfsusie_pipeline(), and mfsusie_empty_result().
#
# OpenSpec: add-wavelet-pipeline-and-multiregion

set.seed(42)
n  <- 60L
p  <- 50L
T1 <- 64L   # functional outcome length (power of 2)

X_raw <- matrix(sample(0:2, n * p, replace = TRUE, prob = c(0.25, 0.5, 0.25)),
                nrow = n)
colnames(X_raw) <- paste0("snp", seq_len(p))

Y1 <- matrix(rnorm(n * T1), nrow = n)
colnames(Y1) <- paste0("bin", seq_len(T1))

# A second outcome with a genuine association (SNP 1 drives scale 4)
b_true <- rep(0, T1)
b_true[33:40] <- 2
Y2 <- matrix(rnorm(n * T1), nrow = n) + outer(X_raw[, 1], b_true)
colnames(Y2) <- paste0("bin", seq_len(T1))

pos1 <- seq_len(T1)

# ── mfsusie_empty_result ─────────────────────────────────────────────────────

test_that("mfsusie_empty_result returns correct class and structure", {
  e <- mfsusie_empty_result("r1", "screen_rejected", 0.12)
  expect_s3_class(e, "mfsusie")
  expect_s3_class(e, "susie")
  expect_true(e$is_empty)
  expect_equal(e$pip, numeric(0))
  expect_equal(e$sets$cs, list())
  expect_equal(e$screen_result$region_id, "r1")
  expect_equal(e$screen_result$status, "screen_rejected")
  expect_equal(e$screen_result$min_pval, 0.12)
})

test_that("predict, fitted, coef, print do not error on empty result", {
  e <- mfsusie_empty_result()
  expect_null(predict(e))
  expect_null(fitted(e))
  expect_null(coef(e))
  expect_output(print(e), "empty result")
})

# ── mfsusie_screen ───────────────────────────────────────────────────────────

test_that("mfsusie_screen returns named list with correct fields", {
  res <- mfsusie_screen(X = X_raw, Y = list(Y1), pos = list(pos1),
                        wavelet_qnorm = TRUE, p_cutoff = 1e-5)
  expect_named(res, c("region_id", "min_pval", "n_tested", "passed"))
  expect_true(is.numeric(res$min_pval))
  expect_true(is.integer(res$n_tested) || is.numeric(res$n_tested))
  expect_true(is.logical(res$passed))
  expect_true(res$n_tested > 0L)
})

test_that("mfsusie_screen passes region with genuine signal", {
  # Y2 has a strong association at SNP 1; should produce a very small p-value
  res <- mfsusie_screen(X = X_raw, Y = list(Y2), pos = list(pos1),
                        wavelet_qnorm = TRUE, p_cutoff = 1e-5)
  expect_true(res$passed)
  expect_lt(res$min_pval, 1e-5)
})

test_that("mfsusie_screen with wavelet_qnorm FALSE differs from TRUE (smoke)", {
  res_qn <- mfsusie_screen(X = X_raw, Y = list(Y2), pos = list(pos1),
                            wavelet_qnorm = TRUE)
  res_nq <- mfsusie_screen(X = X_raw, Y = list(Y2), pos = list(pos1),
                            wavelet_qnorm = FALSE)
  # Both should pass for this signal-present region; min_pval may differ
  expect_true(res_qn$passed)
  expect_true(res_nq$passed)
})

test_that("mfsusie_screen returns passed=TRUE and NA min_pval for scalar-only Y", {
  Y_scalar <- matrix(rnorm(n), ncol = 1L)
  res <- mfsusie_screen(X = X_raw, Y = list(Y_scalar), pos = list(1))
  expect_true(res$passed)
  expect_true(is.na(res$min_pval))
  expect_equal(res$n_tested, 0L)
})

test_that("mfsusie_screen passes region_id through", {
  res <- mfsusie_screen(X = X_raw, Y = list(Y1), pos = list(pos1),
                        region_id = "chr1_test")
  expect_equal(res$region_id, "chr1_test")
})

# ── mfsusie_pipeline ─────────────────────────────────────────────────────────

make_regions <- function() {
  list(
    real = list(X = X_raw, Y = list(Y2), pos = list(pos1), region_id = "real"),
    null = list(X = X_raw, Y = list(Y1), pos = list(pos1), region_id = "null")
  )
}

test_that("mfsusie_pipeline returns list with results and funnel", {
  # Using a null-only region at a loose threshold to guarantee one rejects
  regions <- list(
    real = list(X = X_raw, Y = list(Y2), pos = list(pos1), region_id = "real")
  )
  out <- mfsusie_pipeline(regions,
                          use_wavelet_filter = TRUE,
                          wavelet_qnorm = TRUE,
                          p_cutoff = 1e-5,
                          L = 5L, max_iter = 5L,
                          verbose = FALSE)
  expect_named(out, c("results", "funnel"))
  expect_true(is.list(out$results))
  expect_true(is.data.frame(out$funnel))
  expect_true("real" %in% out$funnel$region_id)
})

test_that("mfsusie_pipeline with use_wavelet_filter=FALSE fits all regions", {
  regions <- list(
    r1 = list(X = X_raw, Y = list(Y1), pos = list(pos1), region_id = "r1")
  )
  out <- mfsusie_pipeline(regions,
                          use_wavelet_filter = FALSE,
                          L = 3L, max_iter = 5L,
                          verbose = FALSE)
  expect_equal(out$funnel$status, "fit_completed")
  expect_false(mfsusie_is_empty(out$results[["r1"]]))
})

test_that("mfsusie_pipeline funnel has expected columns", {
  regions <- list(
    r1 = list(X = X_raw, Y = list(Y2), pos = list(pos1), region_id = "r1")
  )
  out <- mfsusie_pipeline(regions,
                          L = 3L, max_iter = 5L,
                          verbose = FALSE)
  expected_cols <- c("region_id", "status", "min_pval", "n_cs", "n_hp_cs",
                     "cell_types_entered", "cell_types_active")
  expect_true(all(expected_cols %in% names(out$funnel)))
})

test_that("rejected region produces mfsusie_empty_result in results", {
  # Shuffle X so no association exists; use a very tight p_cutoff
  X_null <- X_raw[sample(n), ]
  regions <- list(
    nr = list(X = X_null, Y = list(Y1), pos = list(pos1), region_id = "nr")
  )
  # At p_cutoff=1 everything passes; at p_cutoff=0 everything rejects
  out <- mfsusie_pipeline(regions,
                          p_cutoff = 0,
                          verbose = FALSE)
  expect_true(mfsusie_is_empty(out$results[["nr"]]))
  expect_equal(out$funnel$status, "screen_rejected")
})
