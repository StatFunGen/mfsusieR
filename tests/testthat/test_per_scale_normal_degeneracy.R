# Tests for `prior_variance_scope = "per_scale_normal"` and
# `"per_scale_laplace"` (ebnm-backed point-* priors).
#
# Section layout mirrors
# `inst/openspec/changes/add-per-scale-point-normal-prior/tasks.md`
# sections 4 (init helper), 5 (M-step dispatch), 6 (cache gate),
# 7 (degeneracy locks), 7e (ebnm wrapper forwarding), 7f
# (numerical-property tests), 7g (end-to-end shape), and 8
# (end-to-end power tests).

# ---- Shared fixtures --------------------------------------------

# C1 fixture: scalar Y, the susie-degenerate locks reuse this
# shape. Y is mean-centered + sd-scaled so mfsusie's univariate
# scaling is a no-op vs susieR's.
make_c1_fixture <- function() {
  set.seed(1L); n <- 200L; p <- 50L
  X <- matrix(rnorm(n * p), n)
  beta <- numeric(p); beta[c(7L, 18L)] <- c(1.5, -1.0)
  y_raw <- as.numeric(X %*% beta + rnorm(n, sd = 1))
  y <- (y_raw - mean(y_raw)) / stats::sd(y_raw)
  list(X = X, y = y, signal_idx = c(7L, 18L))
}

# Sparse-wavelet fixture from the spec acceptance criteria.
# Signal at variables c(7, 18), positions c(20, 44).
make_sparse_fixture <- function() {
  set.seed(2L); n <- 200L; p <- 30L; T_m <- 64L
  X <- matrix(rnorm(n * p), n)
  shape <- numeric(T_m); shape[c(20L, 44L)] <- 1
  beta  <- numeric(p);    beta[c(7L, 18L)]  <- c(1.0, -1.0)
  Y <- X %*% (matrix(beta, p, 1) %*% matrix(shape, 1, T_m)) +
       matrix(rnorm(n * T_m, sd = 0.5), n)
  list(X = X, Y = Y, T_m = T_m, signal_idx = c(7L, 18L),
       signal_pos = c(20L, 44L))
}

# Run mfsusie under the susie-degenerate parameter point: scalar
# Y, fixed slab variance via length-1 prior_variance_grid,
# `null_prior_init = 0`, EB off, ELBO convergence to match
# susieR's stopping rule.
fit_susie_degen <- function(scope, X, y, sigma2 = 0.2, L = 5L) {
  mfsusieR::mfsusie(
    wavelet_qnorm           = FALSE,
    X, list(matrix(y, ncol = 1L)),
    L                       = L,
    prior_variance_scope    = scope,
    prior_variance_grid     = sigma2,
    null_prior_init         = 0,
    residual_variance_scope = "per_outcome",
    estimate_prior_variance = FALSE,
    convergence_method      = "elbo",
    L_greedy                = NULL,
    max_iter = 100L, tol = 1e-8,
    verbose = FALSE)
}

# ============================================================
# Section 4. Init helper (`init_ebnm_prior_per_scale`)
# ============================================================

test_that("4: init helper returns the documented G_prior shape (Normal)", {
  fx <- make_sparse_fixture()
  data <- mfsusieR:::create_mf_individual(
    X = fx$X, Y = list(fx$Y), pos = list(seq_len(fx$T_m)),
    standardize = TRUE, intercept = TRUE,
    max_padded_log2 = 10, wavelet_basis_order = 10,
    wavelet_family = "DaubLeAsymm",
    wavelet_magnitude_cutoff = 0,
    wavelet_qnorm = TRUE, verbose = FALSE)
  out <- mfsusieR:::init_ebnm_prior_per_scale(
    Y_m         = data$D[[1L]],
    X           = data$X,
    prior_class = "mixture_point_normal_per_scale",
    groups      = data$scale_index[[1L]],
    na_idx      = data$na_idx[[1L]])
  expect_equal(class(out$G_prior), "mixture_point_normal_per_scale")
  expect_equal(length(out$G_prior), length(data$scale_index[[1L]]))
  expect_setequal(names(out$G_prior[[1L]]),
                   c("fitted_g", "idx", "lead_init_s"))
  expect_equal(class(out$G_prior[[1L]]$fitted_g), "normalmix")
  expect_length(out$G_prior[[1L]]$fitted_g$pi, 2L)
  expect_equal(sum(out$G_prior[[1L]]$fitted_g$pi), 1, tolerance = 1e-12)
})

test_that("4: init helper returns the documented G_prior shape (Laplace)", {
  fx <- make_sparse_fixture()
  data <- mfsusieR:::create_mf_individual(
    X = fx$X, Y = list(fx$Y), pos = list(seq_len(fx$T_m)),
    standardize = TRUE, intercept = TRUE,
    max_padded_log2 = 10, wavelet_basis_order = 10,
    wavelet_family = "DaubLeAsymm",
    wavelet_magnitude_cutoff = 0,
    wavelet_qnorm = TRUE, verbose = FALSE)
  out <- mfsusieR:::init_ebnm_prior_per_scale(
    Y_m         = data$D[[1L]],
    X           = data$X,
    prior_class = "mixture_point_laplace_per_scale",
    groups      = data$scale_index[[1L]],
    na_idx      = data$na_idx[[1L]])
  expect_equal(class(out$G_prior), "mixture_point_laplace_per_scale")
  expect_equal(class(out$G_prior[[1L]]$fitted_g), "laplacemix")
  expect_length(out$G_prior[[1L]]$fitted_g$scale, 2L)
  expect_equal(out$G_prior[[1L]]$fitted_g$scale[1L], 0,
                tolerance = 1e-12)
})

test_that("4: marginal-data lead picker selects a signal-bearing variable on sparse fixture", {
  fx   <- make_sparse_fixture()
  data <- mfsusieR:::create_mf_individual(
    X = fx$X, Y = list(fx$Y), pos = list(seq_len(fx$T_m)),
    standardize = TRUE, intercept = TRUE,
    max_padded_log2 = 10, wavelet_basis_order = 10,
    wavelet_family = "DaubLeAsymm",
    wavelet_magnitude_cutoff = 0, wavelet_qnorm = TRUE,
    verbose = FALSE)
  out <- mfsusieR:::init_ebnm_prior_per_scale(
    Y_m         = data$D[[1L]],
    X           = data$X,
    prior_class = "mixture_point_normal_per_scale",
    groups      = data$scale_index[[1L]],
    na_idx      = data$na_idx[[1L]])
  leads   <- vapply(out$G_prior, function(g) g$lead_init_s, integer(1L))
  groups  <- data$scale_index[[1L]]
  signal_pos <- fx$signal_pos
  overlap <- vapply(groups,
                     function(idx) any(idx %in% signal_pos),
                     logical(1L))
  # Every scale whose `idx_s` overlaps the signal positions must
  # pick a signal-bearing variable.
  expect_true(all(leads[overlap] %in% fx$signal_idx))
})

# ============================================================
# Section 5. M-step dispatch (`.opv_ebnm_point_*`)
# ============================================================

test_that("5: end-to-end mfsusie() builds a fit on per_scale_normal", {
  fx  <- make_sparse_fixture()
  set.seed(1L)
  fit <- mfsusieR::mfsusie(
    fx$X, list(fx$Y), pos = list(seq_len(fx$T_m)),
    L = 5L,
    prior_variance_scope = "per_scale_normal",
    estimate_prior_variance = TRUE,
    L_greedy = NULL, max_iter = 50L, tol = 1e-3,
    verbose = FALSE)
  expect_s3_class(fit, "mfsusie")
  expect_s3_class(fit, "susie")
  expect_equal(class(fit$G_prior[[1L]]),
                "mixture_point_normal_per_scale")
  S_m <- length(fit$G_prior[[1L]])
  expect_equal(dim(fit$pi_V[[1L]]), c(S_m, 2L))
  expect_true(is.matrix(fit$alpha))
  expect_true(!is.null(fit$pip))
  expect_true(!is.null(fit$sets$cs))
})

test_that("5: end-to-end mfsusie() builds a fit on per_scale_laplace", {
  fx <- make_sparse_fixture()
  set.seed(1L)
  fit <- suppressWarnings(mfsusieR::mfsusie(
    fx$X, list(fx$Y), pos = list(seq_len(fx$T_m)),
    L = 5L,
    prior_variance_scope = "per_scale_laplace",
    estimate_prior_variance = TRUE,
    L_greedy = NULL, max_iter = 50L, tol = 1e-3,
    verbose = FALSE))
  expect_s3_class(fit, "mfsusie")
  expect_equal(class(fit$G_prior[[1L]]),
                "mixture_point_laplace_per_scale")
  expect_equal(class(fit$G_prior[[1L]][[1L]]$fitted_g), "laplacemix")
  S_m <- length(fit$G_prior[[1L]])
  expect_equal(dim(fit$pi_V[[1L]]), c(S_m, 2L))
})

test_that("5: ebnm M-step writes both fitted_g and pi_V", {
  # Dispatch wires `fit$fitted_g` -> `G_prior$fitted_g` and
  # `fit$fitted_g$pi` -> `pi_V[[m]][s, ]`. Pin the equality.
  fx  <- make_sparse_fixture()
  set.seed(1L)
  fit <- mfsusieR::mfsusie(
    fx$X, list(fx$Y), pos = list(seq_len(fx$T_m)),
    L = 5L,
    prior_variance_scope = "per_scale_normal",
    estimate_prior_variance = TRUE,
    L_greedy = NULL, max_iter = 50L, tol = 1e-3,
    verbose = FALSE)
  for (s in seq_along(fit$G_prior[[1L]])) {
    expect_equal(fit$pi_V[[1L]][s, ],
                 fit$G_prior[[1L]][[s]]$fitted_g$pi,
                 tolerance = 0)
  }
})

# ============================================================
# Section 5e. ebnm wrapper forwarding (mocked)
# ============================================================

test_that("5e: .opv_ebnm_point_normal forwards (x, s, g_init, fix_g) to ebnm::ebnm_point_normal", {
  fx <- make_sparse_fixture()
  data <- mfsusieR:::create_mf_individual(
    X = fx$X, Y = list(fx$Y), pos = list(seq_len(fx$T_m)),
    standardize = TRUE, intercept = TRUE,
    max_padded_log2 = 10, wavelet_basis_order = 10,
    wavelet_family = "DaubLeAsymm",
    wavelet_magnitude_cutoff = 0, wavelet_qnorm = TRUE,
    verbose = FALSE)
  prior <- mfsusieR:::mf_prior_scale_mixture(data,
              prior_variance_scope = "per_scale_normal",
              null_prior_init = 0)
  ser_stats <- list(
    betahat = list(matrix(rnorm(data$p * data$T_basis[1L]),
                           nrow = data$p)),
    shat2   = list(matrix(0.04, nrow = data$p,
                           ncol = data$T_basis[1L])))
  model <- list(
    M = data$M,
    G_prior = prior$G_prior,
    pi_V    = prior$pi,
    alpha   = matrix(1 / data$p, nrow = 1L, ncol = data$p))

  call_log <- list()
  fake_ebnm <- function(x, s, g_init, fix_g, ...) {
    call_log[[length(call_log) + 1L]] <<- list(
      x = x, s = s, g_init = g_init, fix_g = fix_g)
    g_init  # echo back so the M-step write succeeds
    list(fitted_g = g_init)
  }
  testthat::with_mocked_bindings(
    {
      mfsusieR:::.opv_ebnm_point_normal(
        data, list(estimate_prior_variance = TRUE,
                   alpha_thin_eps = 1e-6),
        model, ser_stats,
        keep_idx = seq_len(data$p),
        zeta_keep = rep(1 / data$p, data$p))
    },
    ebnm_point_normal = fake_ebnm,
    .package = "ebnm"
  )
  expect_equal(length(call_log),
                length(prior$G_prior[[1L]]))   # one call per (m=1, s)
  first <- call_log[[1L]]
  idx1  <- prior$G_prior[[1L]][[1L]]$idx
  expect_length(first$x, length(idx1) * data$p)   # all kept variables
  expect_length(first$s, length(idx1) * data$p)
  expect_false(first$fix_g)                        # estimate_prior_variance=TRUE
  expect_equal(first$g_init, prior$G_prior[[1L]][[1L]]$fitted_g)
})

test_that("5e: fix_g = TRUE when estimate_prior_variance = FALSE", {
  fx <- make_sparse_fixture()
  data <- mfsusieR:::create_mf_individual(
    X = fx$X, Y = list(fx$Y), pos = list(seq_len(fx$T_m)),
    standardize = TRUE, intercept = TRUE,
    max_padded_log2 = 10, wavelet_basis_order = 10,
    wavelet_family = "DaubLeAsymm",
    wavelet_magnitude_cutoff = 0, wavelet_qnorm = TRUE,
    verbose = FALSE)
  prior <- mfsusieR:::mf_prior_scale_mixture(data,
              prior_variance_scope = "per_scale_normal",
              null_prior_init = 0)
  ser_stats <- list(
    betahat = list(matrix(rnorm(data$p * data$T_basis[1L]),
                           nrow = data$p)),
    shat2   = list(matrix(0.04, nrow = data$p,
                           ncol = data$T_basis[1L])))
  model <- list(
    M = data$M, G_prior = prior$G_prior, pi_V = prior$pi,
    alpha = matrix(1 / data$p, nrow = 1L, ncol = data$p))

  fix_g_seen <- NA
  fake_ebnm <- function(x, s, g_init, fix_g, ...) {
    fix_g_seen <<- fix_g
    list(fitted_g = g_init)
  }
  testthat::with_mocked_bindings(
    mfsusieR:::.opv_ebnm_point_normal(
      data, list(estimate_prior_variance = FALSE,
                 alpha_thin_eps = 1e-6),
      model, ser_stats,
      keep_idx = seq_len(data$p),
      zeta_keep = rep(1 / data$p, data$p)),
    ebnm_point_normal = fake_ebnm,
    .package = "ebnm")
  expect_true(isTRUE(fix_g_seen))
})

# ============================================================
# Section 6. Cache gate (`refresh_iter_cache`)
# ============================================================

test_that("6: iter_cache skips sdmat / log_sdmat on the ebnm path", {
  fx  <- make_sparse_fixture()
  set.seed(1L)
  fit <- suppressWarnings(mfsusieR::mfsusie(
    fx$X, list(fx$Y), pos = list(seq_len(fx$T_m)),
    L = 5L, prior_variance_scope = "per_scale_normal",
    estimate_prior_variance = TRUE, L_greedy = NULL,
    max_iter = 50L, tol = 1e-3, verbose = FALSE))
  expect_true(!is.null(fit$iter_cache$shat2[[1L]]))
  # shat2 is `p x T_basis` per outcome.
  expect_equal(nrow(fit$iter_cache$shat2[[1L]]), ncol(fx$X))
  expect_null(fit$iter_cache$sdmat)
  expect_null(fit$iter_cache$log_sdmat)
})

test_that("6: iter_cache keeps sdmat / log_sdmat on the mixsqp path", {
  fx  <- make_sparse_fixture()
  set.seed(1L)
  fit <- suppressWarnings(mfsusieR::mfsusie(
    fx$X, list(fx$Y), pos = list(seq_len(fx$T_m)),
    L = 5L, prior_variance_scope = "per_outcome",
    estimate_prior_variance = TRUE, L_greedy = NULL,
    max_iter = 50L, tol = 1e-3, verbose = FALSE))
  expect_false(is.null(fit$iter_cache$shat2[[1L]]))
  expect_false(is.null(fit$iter_cache$sdmat[[1L]][[1L]]))
  expect_false(is.null(fit$iter_cache$log_sdmat[[1L]][[1L]]))
})

# ============================================================
# Section 7a. Susie-degenerate locks (T = 1, EB off)
# ============================================================

test_that("7a.3: per_scale_normal bit-matches susieR::susie at machine precision", {
  fx <- make_c1_fixture()
  fit_s <- susieR::susie(fx$X, fx$y, L = 5L,
                          scaled_prior_variance      = 0.2,
                          estimate_prior_variance    = FALSE,
                          estimate_residual_variance = TRUE,
                          max_iter = 100L, tol = 1e-8)
  fit_psn <- fit_susie_degen("per_scale_normal", fx$X, fx$y)
  expect_equal(fit_psn$alpha, fit_s$alpha,        tolerance = 1e-12)
  expect_equal(unname(fit_psn$pip),
                unname(fit_s$pip),                  tolerance = 1e-12)
  expect_equal(fit_psn$sigma2[[1L]], fit_s$sigma2, tolerance = 1e-12)
  expect_identical(fit_psn$niter, fit_s$niter)
  expect_equal(lapply(fit_psn$sets$cs, sort),
                lapply(fit_s$sets$cs, sort))
})

test_that("7a.2: per_scale + length-1 grid bit-matches susieR::susie", {
  fx <- make_c1_fixture()
  fit_s   <- susieR::susie(fx$X, fx$y, L = 5L,
                            scaled_prior_variance      = 0.2,
                            estimate_prior_variance    = FALSE,
                            estimate_residual_variance = TRUE,
                            max_iter = 100L, tol = 1e-8)
  fit_pms <- fit_susie_degen("per_scale", fx$X, fx$y)
  expect_equal(fit_pms$alpha, fit_s$alpha,    tolerance = 1e-12)
  expect_equal(unname(fit_pms$pip),
                unname(fit_s$pip),              tolerance = 1e-12)
})

# ============================================================
# Section 7b. Cross-flavor consistency at the degenerate point
# ============================================================

test_that("7b.2: per_outcome and per_scale_normal bit-match at degenerate point", {
  fx <- make_c1_fixture()
  fit_po  <- fit_susie_degen("per_outcome",      fx$X, fx$y)
  fit_psn <- fit_susie_degen("per_scale_normal", fx$X, fx$y)
  expect_equal(fit_psn$alpha, fit_po$alpha, tolerance = 1e-12)
  expect_equal(fit_psn$pip,   fit_po$pip,   tolerance = 1e-12)
})

test_that("7b.3: per_scale and per_scale_normal bit-match at degenerate point", {
  fx <- make_c1_fixture()
  fit_pms <- fit_susie_degen("per_scale",        fx$X, fx$y)
  fit_psn <- fit_susie_degen("per_scale_normal", fx$X, fx$y)
  expect_equal(fit_psn$alpha, fit_pms$alpha, tolerance = 1e-12)
  expect_equal(fit_psn$pip,   fit_pms$pip,   tolerance = 1e-12)
})

# ============================================================
# Section 7d. L invariance under the degenerate setup
# ============================================================

test_that("7d: per_scale_normal bit-matches susie across L = 1, 5, 10", {
  fx <- make_c1_fixture()
  for (L_value in c(1L, 5L, 10L)) {
    fit_s <- susieR::susie(fx$X, fx$y, L = L_value,
                            scaled_prior_variance      = 0.2,
                            estimate_prior_variance    = FALSE,
                            estimate_residual_variance = TRUE,
                            max_iter = 100L, tol = 1e-8)
    fit_psn <- fit_susie_degen("per_scale_normal", fx$X, fx$y,
                                L = L_value)
    expect_equal(fit_psn$alpha, fit_s$alpha, tolerance = 1e-12,
                 info = sprintf("L = %d", L_value))
    expect_equal(unname(fit_psn$pip),
                 unname(fit_s$pip),            tolerance = 1e-12,
                 info = sprintf("L = %d", L_value))
  }
})

# ============================================================
# Section 7f. Numerical-property tests
# ============================================================

test_that("7f.1: M-step idempotence at machine precision (Normal + Laplace)", {
  set.seed(1L)
  n <- 500L
  x <- rnorm(n) * stats::rbinom(n, 1L, 0.4)
  s <- rep(0.3, n)
  for (ebnm_fn in list(ebnm::ebnm_point_normal,
                        ebnm::ebnm_point_laplace)) {
    fit1 <- ebnm_fn(x = x, s = s)
    fit2 <- ebnm_fn(x = x, s = s, g_init = fit1$fitted_g, fix_g = FALSE)
    expect_equal(fit2$fitted_g, fit1$fitted_g, tolerance = 1e-12)
  }
})

test_that("7f.2 / 7f.3: probability constraints + slab non-negativity", {
  set.seed(1L)
  n <- 500L
  x <- rnorm(n) * stats::rbinom(n, 1L, 0.4); s <- rep(0.3, n)
  fit_n <- ebnm::ebnm_point_normal(x = x, s = s)
  fit_l <- ebnm::ebnm_point_laplace(x = x, s = s)
  expect_length(fit_n$fitted_g$pi, 2L)
  expect_equal(sum(fit_n$fitted_g$pi), 1, tolerance = 1e-12)
  expect_true(all(fit_n$fitted_g$pi >= 0))
  expect_equal(fit_n$fitted_g$sd[1L], 0, tolerance = 1e-12)
  expect_gte(fit_n$fitted_g$sd[2L], 0)
  expect_equal(fit_l$fitted_g$scale[1L], 0, tolerance = 1e-12)
  expect_gte(fit_l$fitted_g$scale[2L], 0)
})

test_that("7f.4 / 7f.5: signal recovery + noise recovery sanity", {
  set.seed(1L)
  n <- 1000L
  z <- stats::rbinom(n, 1L, 0.5)
  x <- ifelse(z == 1L, rnorm(n, sd = 1), 0) + rnorm(n, sd = 0.3)
  s <- rep(0.3, n)
  fit <- ebnm::ebnm_point_normal(x = x, s = s)
  expect_gte(fit$fitted_g$pi[1L], 0.4)
  expect_lte(fit$fitted_g$pi[1L], 0.6)
  expect_gte(fit$fitted_g$sd[2L], 0.85)
  expect_lte(fit$fitted_g$sd[2L], 1.15)

  # All-null recovery: pi_0 should be substantially above 0.5.
  # ebnm cannot push pi_0 above ~0.85 here because, with
  # `s = 0.3` matching the data's noise sd, the spike (sd=0)
  # and slab (sd ~ 0) marginals coincide and ebnm settles on
  # an interior optimum. Stricter recovery requires a fixture
  # where the noise sd is well below s.
  set.seed(1L)
  x_null <- rnorm(n, sd = 0.3); s_null <- rep(0.3, n)
  fit_null <- ebnm::ebnm_point_normal(x = x_null, s = s_null)
  expect_gte(fit_null$fitted_g$pi[1L], 0.8)
})

# ============================================================
# Section 7g. End-to-end shape contract (downstream methods)
# ============================================================

test_that("7g.4: predict / coef / fitted / summary / print work on per_scale_normal", {
  fx  <- make_sparse_fixture()
  set.seed(1L)
  fit <- mfsusieR::mfsusie(
    fx$X, list(fx$Y), pos = list(seq_len(fx$T_m)),
    L = 3L,
    prior_variance_scope = "per_scale_normal",
    estimate_prior_variance = TRUE,
    L_greedy = NULL, max_iter = 30L, tol = 1e-3,
    verbose = FALSE)
  expect_no_error(predict(fit))
  expect_no_error(coef(fit))
  expect_no_error(fitted(fit))
  expect_no_error(summary(fit))
  expect_no_error(print(fit))
})

test_that("7g.4: predict / coef / fitted / summary / print work on per_scale_laplace", {
  fx <- make_sparse_fixture()
  set.seed(1L)
  fit <- suppressWarnings(mfsusieR::mfsusie(
    fx$X, list(fx$Y), pos = list(seq_len(fx$T_m)),
    L = 3L,
    prior_variance_scope = "per_scale_laplace",
    estimate_prior_variance = TRUE,
    L_greedy = NULL, max_iter = 30L, tol = 1e-3,
    verbose = FALSE))
  expect_no_error(predict(fit))
  expect_no_error(coef(fit))
  expect_no_error(fitted(fit))
  expect_no_error(summary(fit))
  expect_no_error(print(fit))
})

# ============================================================
# Section 8. End-to-end power tests
# ============================================================

test_that("8.1: per_scale_normal recovers planted causals on sparse fixture", {
  fx  <- make_sparse_fixture()
  set.seed(1L)
  fit <- suppressWarnings(mfsusieR::mfsusie(
    fx$X, list(fx$Y), pos = list(seq_len(fx$T_m)),
    L = 5L, prior_variance_scope = "per_scale_normal",
    estimate_prior_variance = TRUE,
    L_greedy = NULL, max_iter = 50L, tol = 1e-3,
    verbose = FALSE))
  expect_gte(fit$pip[fx$signal_idx[1L]], 0.9)
  expect_gte(fit$pip[fx$signal_idx[2L]], 0.9)
  expect_equal(length(fit$sets$cs), 2L)
})

test_that("8.1: per_scale_laplace recovers planted causals on sparse fixture", {
  fx  <- make_sparse_fixture()
  set.seed(1L)
  fit <- suppressWarnings(mfsusieR::mfsusie(
    fx$X, list(fx$Y), pos = list(seq_len(fx$T_m)),
    L = 5L, prior_variance_scope = "per_scale_laplace",
    estimate_prior_variance = TRUE,
    L_greedy = NULL, max_iter = 50L, tol = 1e-3,
    verbose = FALSE))
  expect_gte(fit$pip[fx$signal_idx[1L]], 0.9)
  expect_gte(fit$pip[fx$signal_idx[2L]], 0.9)
  expect_equal(length(fit$sets$cs), 2L)
})

test_that("8.4: null locus returns no spurious CSes (per_scale_normal)", {
  set.seed(3L)
  n <- 200L; p <- 30L; T_m <- 32L
  X <- matrix(rnorm(n * p), n)
  Y <- matrix(rnorm(n * T_m, sd = 0.5), n)
  set.seed(1L)
  fit <- suppressWarnings(mfsusieR::mfsusie(
    X, list(Y), pos = list(seq_len(T_m)),
    L = 5L, prior_variance_scope = "per_scale_normal",
    estimate_prior_variance = TRUE,
    L_greedy = NULL, max_iter = 50L, tol = 1e-3,
    verbose = FALSE))
  expect_equal(length(fit$sets$cs), 0L)
})

# ============================================================
# Section 9. Laplace lbf kernel agreement with ebnm internals
# ============================================================

test_that("9: mixture_log_bf_laplace_per_scale matches ebnm internal vloglik", {
  set.seed(1L); p <- 30L; T_idx <- 8L
  bhat <- matrix(rnorm(p * T_idx), p, T_idx)
  shat <- matrix(0.3, p, T_idx)
  pi_0 <- 0.7; lambda <- 0.5
  fitted_g <- structure(
    list(pi = c(pi_0, 1 - pi_0), mean = c(0, 0),
         scale = c(0, lambda)),
    class = "laplacemix")

  # Ours: per-j sum_t log mixture density - log null density.
  ours <- mfsusieR:::mixture_log_bf_laplace_per_scale(
            bhat, shat, fitted_g)

  # Reference: vectorize, call ebnm internal, reshape, sum.
  log_dens_mix  <- ebnm:::vloglik_point_laplace(
    x  = as.vector(bhat),
    s  = as.vector(shat),
    w  = 1 - pi_0, a = 1 / lambda, mu = 0)
  log_dens_null <- dnorm(as.vector(bhat), 0,
                         as.vector(shat), log = TRUE)
  ref <- rowSums(matrix(log_dens_mix - log_dens_null,
                         nrow = p, ncol = T_idx))
  expect_equal(ours, ref, tolerance = 1e-12)
})

test_that("9: mixture_log_bf_laplace_per_scale handles pi_0 = 1 / lambda = 0 boundaries", {
  set.seed(1L); p <- 10L; T_idx <- 4L
  bhat <- matrix(rnorm(p * T_idx), p, T_idx)
  shat <- matrix(0.3, p, T_idx)
  # All-null prior: BF = log(pi_0) per cell, summed over T_idx.
  fg_null <- structure(list(pi = c(1, 0), mean = c(0, 0),
                             scale = c(0, 0.5)), class = "laplacemix")
  out <- mfsusieR:::mixture_log_bf_laplace_per_scale(
            bhat, shat, fg_null)
  expect_equal(out, rep(T_idx * log(1), p), tolerance = 1e-12)
  # lambda = 0 same effect.
  fg_lam0 <- structure(list(pi = c(0.5, 0.5), mean = c(0, 0),
                             scale = c(0, 0)), class = "laplacemix")
  out2 <- mfsusieR:::mixture_log_bf_laplace_per_scale(
            bhat, shat, fg_lam0)
  expect_equal(out2, rep(T_idx * log(0.5), p), tolerance = 1e-12)
})

# ============================================================
# Section 10. mixture_null_weight is a no-op on the ebnm paths
# ============================================================

test_that("10: mixture_null_weight is silently ignored on per_scale_normal", {
  fx <- make_sparse_fixture()
  set.seed(1L)
  fit_a <- suppressWarnings(mfsusieR::mfsusie(
    fx$X, list(fx$Y), pos = list(seq_len(fx$T_m)),
    L = 5L, prior_variance_scope = "per_scale_normal",
    estimate_prior_variance = TRUE,
    L_greedy = NULL, max_iter = 50L, tol = 1e-3,
    verbose = FALSE))
  set.seed(1L)
  fit_b <- suppressWarnings(mfsusieR::mfsusie(
    fx$X, list(fx$Y), pos = list(seq_len(fx$T_m)),
    L = 5L, prior_variance_scope = "per_scale_normal",
    mixture_null_weight = 0.5,                # very different value
    estimate_prior_variance = TRUE,
    L_greedy = NULL, max_iter = 50L, tol = 1e-3,
    verbose = FALSE))
  expect_equal(fit_a$alpha, fit_b$alpha, tolerance = 1e-12)
  expect_equal(fit_a$pip,   fit_b$pip,   tolerance = 1e-12)
})
