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
  set.seed(seed + 10)
  for (m in seq_len(M)) {
    na_rows <- sort(sample(n, round(n * miss_frac)))
    Y[[m]][na_rows, ] <- NA
  }
  list(X = X, Y = Y, causal_j = causal_j, n = n, p = p)
}

# Run `expr`, record wall time, message "[RUNTIME] label: X.XXs", return
# list(result = <expr value>, elapsed = <seconds>). The generous 120 s
# budget catches only catastrophic regressions on the small fixtures used
# here (n=50, p=20, T=8); real-data runs are covered by submit_runtime_test.sh.
with_runtime <- function(label, expr, budget_sec = 120) {
  t <- system.time(result <- expr)
  elapsed <- unname(t["elapsed"])
  message(sprintf("[RUNTIME] %s: %.2f s", label, elapsed))
  list(result = result, elapsed = elapsed, budget_sec = budget_sec)
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
  rt <- with_runtime("create_mf_individual (no NA)", {
    mfsusieR:::create_mf_individual(X, Y, verbose = FALSE)
  })
  data <- rt$result
  expect_lt(rt$elapsed, rt$budget_sec)

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

  rt <- with_runtime("create_mf_individual (NA trait)", {
    mfsusieR:::create_mf_individual(X, Y, verbose = FALSE)
  })
  data <- rt$result
  expect_lt(rt$elapsed, rt$budget_sec)

  expect_identical(data$na_idx[[1]], seq_len(n))
  expect_equal(data$xtx_diag_list[[1]], data$xtx_diag, tolerance = 1e-12)
  expect_identical(sort(data$na_idx[[2]]), sort(setdiff(seq_len(n), na_rows)))
  expect_equal(data$xtx_diag_list[[2]],
               colSums(data$X[data$na_idx[[2]], , drop = FALSE]^2),
               tolerance = 1e-12)
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
  Y[[2]][6:12, ] <- NA
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
  rt <- with_runtime("mfsusie qnorm=FALSE 20% NA", {
    do.call(mfsusie, c(list(X = d$X, Y = d$Y), mfsusie_args_noqnorm))
  })
  expect_lt(rt$elapsed, rt$budget_sec)
  expect_s3_class(rt$result, "susie")
})

test_that("mfsusie quantile_norm=TRUE does not crash with 20% NA rows per trait", {
  d <- make_na_XY()
  rt <- with_runtime("mfsusie qnorm=TRUE 20% NA", {
    do.call(mfsusie, c(list(X = d$X, Y = d$Y), mfsusie_args_qnorm))
  })
  expect_lt(rt$elapsed, rt$budget_sec)
  expect_s3_class(rt$result, "susie")
})

test_that("mfsusie quantile_norm=FALSE returns valid pip with 20% NA", {
  d <- make_na_XY()
  rt <- with_runtime("mfsusie qnorm=FALSE pip check", {
    do.call(mfsusie, c(list(X = d$X, Y = d$Y), mfsusie_args_noqnorm))
  })
  fit <- rt$result
  expect_lt(rt$elapsed, rt$budget_sec)
  expect_true(is.numeric(fit$pip))
  expect_equal(length(fit$pip), d$p)
  expect_true(all(fit$pip >= 0 & fit$pip <= 1))
})

test_that("mfsusie quantile_norm=TRUE returns valid pip with 20% NA", {
  d <- make_na_XY()
  rt <- with_runtime("mfsusie qnorm=TRUE pip check", {
    do.call(mfsusie, c(list(X = d$X, Y = d$Y), mfsusie_args_qnorm))
  })
  fit <- rt$result
  expect_lt(rt$elapsed, rt$budget_sec)
  expect_true(is.numeric(fit$pip))
  expect_equal(length(fit$pip), d$p)
  expect_true(all(fit$pip >= 0 & fit$pip <= 1))
})

# ---- causal signal recovery under NA -----------------------------------

test_that("causal pip > 0.5 under quantile_norm=FALSE with 20% NA (strong signal)", {
  d <- make_na_XY(effect = 2)
  rt <- with_runtime("mfsusie qnorm=FALSE signal recovery", {
    do.call(mfsusie, c(list(X = d$X, Y = d$Y), mfsusie_args_noqnorm))
  })
  expect_lt(rt$elapsed, rt$budget_sec)
  expect_gt(rt$result$pip[d$causal_j], 0.5)
})

test_that("causal pip > 0.5 under quantile_norm=TRUE with 20% NA (strong signal)", {
  d <- make_na_XY(effect = 2)
  rt <- with_runtime("mfsusie qnorm=TRUE signal recovery", {
    do.call(mfsusie, c(list(X = d$X, Y = d$Y), mfsusie_args_qnorm))
  })
  expect_lt(rt$elapsed, rt$budget_sec)
  expect_gt(rt$result$pip[d$causal_j], 0.5)
})

test_that("causal pip ranks #1 under both qnorm paths with 20% NA", {
  d <- make_na_XY(effect = 2)
  rt_noq <- with_runtime("mfsusie qnorm=FALSE rank-1 causal", {
    do.call(mfsusie, c(list(X = d$X, Y = d$Y), mfsusie_args_noqnorm))
  })
  rt_q <- with_runtime("mfsusie qnorm=TRUE rank-1 causal", {
    do.call(mfsusie, c(list(X = d$X, Y = d$Y), mfsusie_args_qnorm))
  })
  message(sprintf("[RUNTIME] qnorm=FALSE vs qnorm=TRUE: %.2f s vs %.2f s",
                  rt_noq$elapsed, rt_q$elapsed))
  expect_identical(which.max(rt_noq$result$pip), d$causal_j)
  expect_identical(which.max(rt_q$result$pip),   d$causal_j)
})

# ---- runtime comparison: qnorm=FALSE vs qnorm=TRUE across missing rates --

test_that("runtime comparison: qnorm=FALSE vs qnorm=TRUE at 5%, 20%, 40% NA", {
  for (miss_frac in c(0.05, 0.20, 0.40)) {
    d <- make_na_XY(miss_frac = miss_frac)
    rt_noq <- with_runtime(
      sprintf("qnorm=FALSE miss=%.0f%%", miss_frac * 100),
      do.call(mfsusie, c(list(X = d$X, Y = d$Y), mfsusie_args_noqnorm))
    )
    rt_q <- with_runtime(
      sprintf("qnorm=TRUE  miss=%.0f%%", miss_frac * 100),
      do.call(mfsusie, c(list(X = d$X, Y = d$Y), mfsusie_args_qnorm))
    )
    expect_lt(rt_noq$elapsed, rt_noq$budget_sec)
    expect_lt(rt_q$elapsed,   rt_q$budget_sec)
    # Both paths must produce valid pips.
    expect_true(all(rt_noq$result$pip >= 0 & rt_noq$result$pip <= 1))
    expect_true(all(rt_q$result$pip   >= 0 & rt_q$result$pip   <= 1))
  }
})

# ---- edge cases --------------------------------------------------------

test_that("mfsusie quantile_norm=FALSE does not crash: univariate T_m=1 with NA", {
  set.seed(mfsusier_test_seed())
  n <- 50; p <- 15
  X <- scale(matrix(rbinom(n * p, 2, 0.3), n, p))
  Y <- list(matrix(2 * X[, 3] + rnorm(n), n, 1),
            matrix(2 * X[, 3] + rnorm(n), n, 1))
  Y[[1]][sample(n, 10), ] <- NA
  rt <- with_runtime("mfsusie qnorm=FALSE T_m=1 NA", {
    mfsusie(X, Y, L = 3, max_iter = 30, quantile_norm = FALSE, verbose = FALSE)
  })
  expect_lt(rt$elapsed, rt$budget_sec)
  expect_s3_class(rt$result, "susie")
})

test_that("mfsusie quantile_norm=TRUE does not crash: univariate T_m=1 with NA", {
  set.seed(mfsusier_test_seed())
  n <- 50; p <- 15
  X <- scale(matrix(rbinom(n * p, 2, 0.3), n, p))
  Y <- list(matrix(2 * X[, 3] + rnorm(n), n, 1),
            matrix(2 * X[, 3] + rnorm(n), n, 1))
  Y[[1]][sample(n, 10), ] <- NA
  rt <- with_runtime("mfsusie qnorm=TRUE T_m=1 NA", {
    mfsusie(X, Y, L = 3, max_iter = 30, quantile_norm = TRUE, verbose = FALSE)
  })
  expect_lt(rt$elapsed, rt$budget_sec)
  expect_s3_class(rt$result, "susie")
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
  rt <- with_runtime("mfsusie qnorm=FALSE null+NA", {
    mfsusie(X, Y, L = 5, max_iter = 30, quantile_norm = FALSE, verbose = FALSE)
  })
  expect_lt(rt$elapsed, rt$budget_sec)
  expect_true(all(rt$result$pip >= 0 & rt$result$pip <= 1))
  expect_true(max(rt$result$pip) < 0.99)
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
  rt <- with_runtime("mfsusie qnorm=TRUE null+NA", {
    mfsusie(X, Y, L = 5, max_iter = 30, quantile_norm = TRUE, verbose = FALSE)
  })
  expect_lt(rt$elapsed, rt$budget_sec)
  expect_true(all(rt$result$pip >= 0 & rt$result$pip <= 1))
  expect_true(max(rt$result$pip) < 0.99)
})

test_that("mfsusie handles 40% NA per trait without crash (qnorm=FALSE)", {
  d <- make_na_XY(miss_frac = 0.4)
  rt <- with_runtime("mfsusie qnorm=FALSE 40% NA", {
    do.call(mfsusie, c(list(X = d$X, Y = d$Y), mfsusie_args_noqnorm))
  })
  expect_lt(rt$elapsed, rt$budget_sec)
  expect_s3_class(rt$result, "susie")
})

test_that("mfsusie handles 40% NA per trait without crash (qnorm=TRUE)", {
  d <- make_na_XY(miss_frac = 0.4)
  rt <- with_runtime("mfsusie qnorm=TRUE 40% NA", {
    do.call(mfsusie, c(list(X = d$X, Y = d$Y), mfsusie_args_qnorm))
  })
  expect_lt(rt$elapsed, rt$budget_sec)
  expect_s3_class(rt$result, "susie")
})

# ---- per_scale scope also works with NA --------------------------------

test_that("mfsusie per_scale scope does not crash with NA (qnorm=FALSE)", {
  d <- make_na_XY(T_m = 16)
  rt <- with_runtime("mfsusie per_scale qnorm=FALSE NA", {
    mfsusie(d$X, d$Y, L = 3, max_iter = 30,
            prior_variance_scope    = "per_scale",
            residual_variance_scope = "per_scale",
            quantile_norm = FALSE, verbose = FALSE)
  })
  expect_lt(rt$elapsed, rt$budget_sec)
  expect_s3_class(rt$result, "susie")
})

test_that("mfsusie per_scale scope does not crash with NA (qnorm=TRUE)", {
  d <- make_na_XY(T_m = 16)
  rt <- with_runtime("mfsusie per_scale qnorm=TRUE NA", {
    mfsusie(d$X, d$Y, L = 3, max_iter = 30,
            prior_variance_scope    = "per_scale",
            residual_variance_scope = "per_scale",
            quantile_norm = TRUE, verbose = FALSE)
  })
  expect_lt(rt$elapsed, rt$budget_sec)
  expect_s3_class(rt$result, "susie")
})
