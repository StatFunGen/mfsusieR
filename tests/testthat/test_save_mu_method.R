# Locks the contract for save_mu_method = c("complete", "alpha_collapsed", "lead")
# and the mf_thin() post-fit helper.
#
# Invariants enforced:
#   1. complete is the default and produces fits with attr "complete".
#   2. alpha_collapsed reduces fit$mu / fit$mu2 to 1 x T per (l, m) and adds
#      fit$coef_wavelet (1 x T) for the raw-X coef path.
#   3. lead reduces fit$mu / fit$mu2 to 1 x T per (l, m) and stores
#      fit$top_index = which.max(alpha[l, ]).
#   4. coef.mfsusie complete == alpha_collapsed numerically (1e-12).
#   5. coef.mfsusie complete != lead in general (cheap-coef, biased).
#   6. mf_post_smooth scalewise complete == alpha_collapsed numerically (1e-12).
#   7. mf_thin(complete, "alpha_collapsed") == mfsusie(..., "alpha_collapsed").
#   8. predict.mfsusie(newx) errors on 1D modes; fitted() works on 1D modes.
#   9. mfsusie(model_init = thinned_fit) errors with a save_mu_method message.
#  10. mf_thin() rejects re-thinning a 1D fit.

# Each test_that block in this file is wrapped to print a [RUNTIME]
# line on completion, matching the convention used in test_missing_data_na.R
# and test_qnorm_*. The macro is intentionally local to this file.
.timed <- function(label, expr) {
  t0 <- Sys.time()
  force(expr)
  cat(sprintf("[RUNTIME] %s: %.2f s\n", label,
              as.numeric(difftime(Sys.time(), t0, units = "secs"))))
}

build_toy_fit <- function(seed = 1L, save_mu_method = "complete",
                          L = 5L, max_iter = 6L) {
  set.seed(seed)
  n <- 50L; p <- 30L; Tlen <- 32L
  X <- matrix(rnorm(n * p), n, p)
  beta <- numeric(p); beta[c(3, 17)] <- c(2, -1.5)
  Y <- list(matrix(rnorm(n * Tlen), n, Tlen) + as.numeric(X %*% beta))
  fit <- mfsusie(X, Y, L = L, max_iter = max_iter, verbose = FALSE,
                 save_mu_method = save_mu_method)
  list(fit = fit, X = X, Y = Y)
}

test_that("save_mu_method default is complete and recorded as attr", {
  .timed("save_mu_method default", {
    obj <- build_toy_fit()
    expect_equal(attr(obj$fit, "save_mu_method"), "complete")
    expect_true(nrow(obj$fit$mu[[1]][[1]]) == ncol(obj$X))
    expect_null(obj$fit$top_index)
    expect_null(obj$fit$coef_wavelet)
  })
})

test_that("alpha_collapsed shrinks mu/mu2 to 1 x T and adds coef_wavelet", {
  .timed("alpha_collapsed shape", {
    obj <- build_toy_fit(save_mu_method = "alpha_collapsed")
    expect_equal(attr(obj$fit, "save_mu_method"), "alpha_collapsed")
    expect_equal(nrow(obj$fit$mu[[1]][[1]]), 1L)
    expect_equal(nrow(obj$fit$mu2[[1]][[1]]), 1L)
    expect_false(is.null(obj$fit$coef_wavelet))
    expect_equal(nrow(obj$fit$coef_wavelet[[1]][[1]]), 1L)
  })
})

test_that("lead shrinks mu/mu2 to 1 x T and stores top_index", {
  .timed("lead shape", {
    obj <- build_toy_fit(save_mu_method = "lead")
    expect_equal(attr(obj$fit, "save_mu_method"), "lead")
    expect_equal(nrow(obj$fit$mu[[1]][[1]]), 1L)
    expect_equal(nrow(obj$fit$mu2[[1]][[1]]), 1L)
    expect_equal(length(obj$fit$top_index), nrow(obj$fit$alpha))
    for (l in seq_along(obj$fit$top_index)) {
      expect_equal(obj$fit$top_index[l], which.max(obj$fit$alpha[l, ]))
    }
  })
})

test_that("coef.mfsusie alpha_collapsed is numerically equivalent to complete", {
  .timed("coef alpha_collapsed equiv", {
    obj_full <- build_toy_fit(save_mu_method = "complete")
    obj_ac   <- build_toy_fit(save_mu_method = "alpha_collapsed")
    expect_equal(coef(obj_full$fit), coef(obj_ac$fit), tolerance = 1e-12)
  })
})

test_that("coef.mfsusie lead diverges from complete (cheap-coef, biased)", {
  .timed("coef lead divergence", {
    obj_full <- build_toy_fit(save_mu_method = "complete")
    obj_ld   <- build_toy_fit(save_mu_method = "lead")
    diff_max <- max(abs(coef(obj_full$fit)[[1]] - coef(obj_ld$fit)[[1]]))
    expect_gt(diff_max, 1e-3)
  })
})

test_that("mf_post_smooth scalewise alpha_collapsed equals complete (1e-12)", {
  .timed("post_smooth alpha_collapsed equiv", {
    obj_full <- build_toy_fit(save_mu_method = "complete")
    obj_ac   <- build_toy_fit(save_mu_method = "alpha_collapsed")
    sm_full <- mf_post_smooth(obj_full$fit, method = "scalewise")
    sm_ac   <- mf_post_smooth(obj_ac$fit,   method = "scalewise")
    expect_equal(sm_full$smoothed$scalewise$effect_curves,
                 sm_ac$smoothed$scalewise$effect_curves,
                 tolerance = 1e-12)
  })
})

test_that("mf_thin(complete, 'alpha_collapsed') matches a fresh alpha_collapsed fit", {
  .timed("mf_thin alpha_collapsed equiv", {
    obj_full <- build_toy_fit(save_mu_method = "complete")
    obj_ac   <- build_toy_fit(save_mu_method = "alpha_collapsed")
    thinned  <- mf_thin(obj_full$fit, method = "alpha_collapsed")
    expect_equal(attr(thinned, "save_mu_method"), "alpha_collapsed")
    expect_equal(coef(thinned), coef(obj_ac$fit), tolerance = 1e-12)
  })
})

test_that("mf_thin(complete, 'lead') matches a fresh lead fit", {
  .timed("mf_thin lead equiv", {
    obj_full <- build_toy_fit(save_mu_method = "complete")
    obj_ld   <- build_toy_fit(save_mu_method = "lead")
    thinned  <- mf_thin(obj_full$fit, method = "lead")
    expect_equal(attr(thinned, "save_mu_method"), "lead")
    expect_equal(thinned$top_index, obj_ld$fit$top_index)
    expect_equal(coef(thinned), coef(obj_ld$fit), tolerance = 1e-12)
  })
})

test_that("predict.mfsusie(newx) errors on 1D modes; fitted() works", {
  .timed("predict guard 1D", {
    obj_ac <- build_toy_fit(save_mu_method = "alpha_collapsed")
    obj_ld <- build_toy_fit(save_mu_method = "lead")
    expect_error(predict(obj_ac$fit, newx = obj_ac$X), "save_mu_method")
    expect_error(predict(obj_ld$fit, newx = obj_ld$X), "save_mu_method")
    # fitted() reads fit$fitted directly, no per-variant mu needed.
    expect_silent(fitted(obj_ac$fit))
    expect_silent(fitted(obj_ld$fit))
  })
})

test_that("mfsusie(model_init = thinned_fit) errors with save_mu_method message", {
  .timed("model_init guard 1D", {
    obj_ac <- build_toy_fit(save_mu_method = "alpha_collapsed")
    expect_error(
      mfsusie(obj_ac$X, obj_ac$Y, L = 5L, max_iter = 1L, verbose = FALSE,
              model_init = obj_ac$fit),
      "save_mu_method")
  })
})

test_that("mf_thin() rejects re-thinning a fit that is already 1D", {
  .timed("mf_thin rejects 1D", {
    obj_ac <- build_toy_fit(save_mu_method = "alpha_collapsed")
    expect_error(mf_thin(obj_ac$fit, method = "lead"),
                 "save_mu_method")
  })
})
