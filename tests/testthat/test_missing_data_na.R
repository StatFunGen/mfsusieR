# Tests for per-trait complete-case NA handling (commit be0ce13).
#
# Before this fix, Y[[m]] rows of NA propagated through the DWT into
# D[[m]], making compute_marginal_bhat_shat produce all-NA Bhat/Shat
# and causing ash to throw "all input values are missing". The fix
# pre-computes na_idx[[m]] / xtx_diag_list[[m]] at data-class
# construction and threads them through every row-wise computation.
#
# Two paths are tested for each scenario:
#   qnorm = FALSE  -- previously crashed; now the primary regression target
#   qnorm = TRUE   -- was accidentally working (rank() treats NA as a value)
#                     but used wrong n and xtx; now correct for both paths

# ---- helpers -----------------------------------------------------------

make_na_XY <- function(n = 50, p = 20, T_m = 8, M = 3,
                       causal_j = 2L, effect = 2,
                       miss_frac = 0.2, seed = mfsusier_test_seed()) {
  set.seed(seed)
  X <- matrix(rbinom(n * p, 2, 0.3), n, p)
  X <- scale(X)
  signal <- as.numeric(X[, causal_j]) * effect
  Y <- lapply(seq_len(M), function(m)
    matrix(signal, n, T_m) + matrix(rnorm(n * T_m, sd = 0.5), n, T_m))
  # Each trait gets a different set of missing rows.
  set.seed(seed + 10)
  for (m in seq_len(M)) {
    na_rows <- sort(sample(n, round(n * miss_frac)))
    Y[[m]][na_rows, ] <- NA
  }
  list(X = X, Y = Y, causal_j = causal_j, n = n, p = p)
}

mfsusie_args_noqnorm <- list(L = 5, max_iter = 40,
                              prior_variance_scope    = "per_outcome",
                              residual_variance_scope = "per_outcome",
                              quantile_norm = FALSE,
                              verbose = FALSE)

mfsusie_args_qnorm   <- list(L = 5, max_iter = 40,
                              prior_variance_scope    = "per_outcome",
                              residual_variance_scope = "per_outcome",
                              quantile_norm = TRUE,
                              verbose = FALSE)

# ---- create_mf_individual: na_idx and xtx_diag_list contract ----------

test_that("na_idx = seq_len(n) and xtx_diag_list == xtx_diag when Y has no NA", {
  set.seed(mfsusier_test_seed())
  n <- 30; p <- 6
  X <- matrix(rnorm(n * p), n, p)
  Y <- list(matrix(rnorm(n * 8), n, 8),
            matrix(rnorm(n * 4), n, 4))
  data <- mfsusieR:::create_mf_individual(X, Y, verbose = FALSE)

  expect_identical(data$na_idx[[1]], seq_len(n))
  expect_identical(data$na_idx[[2]], seq_len(n))
  expect_equal(data$xtx_diag_list[[1]], data$xtx_diag, tolerance = 1e-12)
  expect_equal(data$xtx_diag_list[[2]], data$xtx_diag, tolerance = 1e-12)
})

test_that("na_idx excludes NA rows and xtx_diag_list < xtx_diag for trait with NAs", {
  set.seed(mfsusier_test_seed())
  n <- 40; p <- 5
  X <- matrix(rnorm(n * p), n, p)
  Y <- list(matrix(rnorm(n * 8), n, 8),
            matrix(rnorm(n * 8), n, 8))
  na_rows <- c(3L, 7L, 15L, 20L)
  Y[[2]][na_rows, ] <- NA

  data <- mfsusieR:::create_mf_individual(X, Y, verbose = FALSE)

  # Trait 1: complete, no change.
  expect_identical(data$na_idx[[1]], seq_len(n))
  expect_equal(data$xtx_diag_list[[1]], data$xtx_diag, tolerance = 1e-12)

  # Trait 2: na_idx excludes the 4 NA rows.
  expect_identical(sort(data$na_idx[[2]]), sort(setdiff(seq_len(n), na_rows)))
  expect_equal(data$xtx_diag_list[[2]],
               colSums(data$X[data$na_idx[[2]], , drop = FALSE]^2),
               tolerance = 1e-12)
  # xtx_diag_list is strictly smaller than global xtx_diag.
  expect_true(all(data$xtx_diag_list[[2]] < data$xtx_diag))
})

test_that("each trait gets an independent na_idx when missing rows differ", {
  set.seed(mfsusier_test_seed())
  n <- 30; p <- 4
  X <- matrix(rnorm(n * p), n, p)
  Y <- list(matrix(rnorm(n * 8), n, 8),
            matrix(rnorm(n * 8), n, 8),
            matrix(rnorm(n * 8), n, 8))
  Y[[1]][1:5,  ] <- NA
  Y[[2]][6:12, ] <- NA  # different rows
  # Y[[3]] complete

  data <- mfsusieR:::create_mf_individual(X, Y, verbose = FALSE)

  expect_false(identical(data$na_idx[[1]], data$na_idx[[2]]))
  expect_identical(data$na_idx[[3]], seq_len(n))
  expect_equal(length(data$na_idx[[1]]), n - 5L)
  expect_equal(length(data$na_idx[[2]]), n - 7L)
})

# ---- mfsusie end-to-end: no crash with NA rows -------------------------

test_that("mfsusie quantile_norm=FALSE does not crash with 20% NA rows per trait", {
  d <- make_na_XY()
  expect_no_error(
    do.call(mfsusie, c(list(X = d$X, Y = d$Y), mfsusie_args_noqnorm))
  )
})

test_that("mfsusie quantile_norm=TRUE does not crash with 20% NA rows per trait", {
  d <- make_na_XY()
  expect_no_error(
    do.call(mfsusie, c(list(X = d$X, Y = d$Y), mfsusie_args_qnorm))
  )
})

test_that("mfsusie quantile_norm=FALSE returns valid pip with 20% NA", {
  d <- make_na_XY()
  fit <- do.call(mfsusie, c(list(X = d$X, Y = d$Y), mfsusie_args_noqnorm))
  expect_true(is.numeric(fit$pip))
  expect_equal(length(fit$pip), d$p)
  expect_true(all(fit$pip >= 0 & fit$pip <= 1))
})

test_that("mfsusie quantile_norm=TRUE returns valid pip with 20% NA", {
  d <- make_na_XY()
  fit <- do.call(mfsusie, c(list(X = d$X, Y = d$Y), mfsusie_args_qnorm))
  expect_true(is.numeric(fit$pip))
  expect_equal(length(fit$pip), d$p)
  expect_true(all(fit$pip >= 0 & fit$pip <= 1))
})

# ---- causal signal recovery under NA -----------------------------------

test_that("causal pip > 0.5 under quantile_norm=FALSE with 20% NA (strong signal)", {
  d <- make_na_XY(effect = 2)
  fit <- do.call(mfsusie, c(list(X = d$X, Y = d$Y), mfsusie_args_noqnorm))
  expect_gt(fit$pip[d$causal_j], 0.5)
})

test_that("causal pip > 0.5 under quantile_norm=TRUE with 20% NA (strong signal)", {
  d <- make_na_XY(effect = 2)
  fit <- do.call(mfsusie, c(list(X = d$X, Y = d$Y), mfsusie_args_qnorm))
  expect_gt(fit$pip[d$causal_j], 0.5)
})

test_that("causal pip ranks #1 under both qnorm paths with 20% NA", {
  d <- make_na_XY(effect = 2)
  fit_noq <- do.call(mfsusie, c(list(X = d$X, Y = d$Y), mfsusie_args_noqnorm))
  fit_q   <- do.call(mfsusie, c(list(X = d$X, Y = d$Y), mfsusie_args_qnorm))
  expect_identical(which.max(fit_noq$pip), d$causal_j)
  expect_identical(which.max(fit_q$pip),   d$causal_j)
})

# ---- edge cases --------------------------------------------------------

test_that("mfsusie quantile_norm=FALSE does not crash: univariate T_m=1 with NA", {
  set.seed(mfsusier_test_seed())
  n <- 50; p <- 15
  X <- scale(matrix(rbinom(n * p, 2, 0.3), n, p))
  Y <- list(matrix(2 * X[, 3] + rnorm(n), n, 1),
            matrix(2 * X[, 3] + rnorm(n), n, 1))
  Y[[1]][sample(n, 10), ] <- NA
  expect_no_error(
    mfsusie(X, Y, L = 3, max_iter = 30, quantile_norm = FALSE, verbose = FALSE)
  )
})

test_that("mfsusie quantile_norm=TRUE does not crash: univariate T_m=1 with NA", {
  set.seed(mfsusier_test_seed())
  n <- 50; p <- 15
  X <- scale(matrix(rbinom(n * p, 2, 0.3), n, p))
  Y <- list(matrix(2 * X[, 3] + rnorm(n), n, 1),
            matrix(2 * X[, 3] + rnorm(n), n, 1))
  Y[[1]][sample(n, 10), ] <- NA
  expect_no_error(
    mfsusie(X, Y, L = 3, max_iter = 30, quantile_norm = TRUE, verbose = FALSE)
  )
})

test_that("mfsusie quantile_norm=FALSE: null simulation with NA stays well-behaved", {
  set.seed(mfsusier_test_seed() + 5L)
  n <- 50; p <- 20
  X <- scale(matrix(rbinom(n * p, 2, 0.3), n, p))
  Y <- lapply(1:3, function(m) {
    y <- matrix(rnorm(n * 8), n, 8)
    y[sample(n, 12), ] <- NA
    y
  })
  fit <- mfsusie(X, Y, L = 5, max_iter = 30, quantile_norm = FALSE, verbose = FALSE)
  expect_true(all(fit$pip >= 0 & fit$pip <= 1))
  expect_true(max(fit$pip) < 0.99)
})

test_that("mfsusie quantile_norm=TRUE: null simulation with NA stays well-behaved", {
  set.seed(mfsusier_test_seed() + 5L)
  n <- 50; p <- 20
  X <- scale(matrix(rbinom(n * p, 2, 0.3), n, p))
  Y <- lapply(1:3, function(m) {
    y <- matrix(rnorm(n * 8), n, 8)
    y[sample(n, 12), ] <- NA
    y
  })
  fit <- mfsusie(X, Y, L = 5, max_iter = 30, quantile_norm = TRUE, verbose = FALSE)
  expect_true(all(fit$pip >= 0 & fit$pip <= 1))
  expect_true(max(fit$pip) < 0.99)
})

test_that("mfsusie handles 40% NA per trait without crash (qnorm=FALSE)", {
  d <- make_na_XY(miss_frac = 0.4)
  expect_no_error(
    do.call(mfsusie, c(list(X = d$X, Y = d$Y), mfsusie_args_noqnorm))
  )
})

test_that("mfsusie handles 40% NA per trait without crash (qnorm=TRUE)", {
  d <- make_na_XY(miss_frac = 0.4)
  expect_no_error(
    do.call(mfsusie, c(list(X = d$X, Y = d$Y), mfsusie_args_qnorm))
  )
})

# ---- per_scale scope also works with NA --------------------------------

test_that("mfsusie per_scale scope does not crash with NA (qnorm=FALSE)", {
  d <- make_na_XY(T_m = 16)
  expect_no_error(
    mfsusie(d$X, d$Y, L = 3, max_iter = 30,
            prior_variance_scope    = "per_scale",
            residual_variance_scope = "per_scale",
            quantile_norm = FALSE, verbose = FALSE)
  )
})

test_that("mfsusie per_scale scope does not crash with NA (qnorm=TRUE)", {
  d <- make_na_XY(T_m = 16)
  expect_no_error(
    mfsusie(d$X, d$Y, L = 3, max_iter = 30,
            prior_variance_scope    = "per_scale",
            residual_variance_scope = "per_scale",
            quantile_norm = TRUE, verbose = FALSE)
  )
})
