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
    wavelet_qnorm = FALSE, wavelet_standardize = FALSE,
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

# Build the data class on the sparse fixture (used by every test
# that needs the wavelet-decomposed `mf_individual` view).
make_sparse_data <- function() {
  fx   <- make_sparse_fixture()
  data <- mfsusieR:::create_mf_individual(
    X = fx$X, Y = list(fx$Y), pos = list(seq_len(fx$T_m)),
    standardize = TRUE, intercept = TRUE,
    max_padded_log2 = 10, wavelet_basis_order = 10,
    wavelet_family = "DaubLeAsymm",
    wavelet_magnitude_cutoff = 0, wavelet_qnorm = TRUE,
    verbose = FALSE)
  list(fx = fx, data = data)
}

# Synthetic ser_stats matching the data fixture's (p, T_basis[1])
# shape; used by the dispatch / wrapper-level tests.
make_ser_stats <- function(data) list(
  betahat = list(matrix(rnorm(data$p * data$T_basis[1L]),
                         nrow = data$p)),
  shat2   = list(matrix(0.04, nrow = data$p,
                         ncol = data$T_basis[1L])))

# Minimal model state for the dispatch tests: one outcome,
# uniform alpha, prior just constructed. pi_V and
# fitted_g_per_effect are per-effect (list[L]) sidecars; the
# dispatch tests use L = 1.
make_seed_model <- function(data, prior, L = 1L) list(
  M       = data$M,
  G_prior = prior$G_prior,
  pi_V    = lapply(seq_len(L), function(.) prior$pi),
  fitted_g_per_effect = lapply(seq_len(L), function(.)
    lapply(prior$G_prior, function(G_m)
      lapply(G_m, function(g_s) g_s$fitted_g))),
  alpha   = matrix(1 / data$p, nrow = L, ncol = data$p))

# Run mfsusie on the sparse fixture under the named scope. Used
# by the end-to-end shape and power tests.
fit_sparse <- function(scope, fx = make_sparse_fixture(), L = 5L,
                       max_iter = 50L) {
  set.seed(1L)
  suppressWarnings(mfsusieR::mfsusie(
    fx$X, list(fx$Y), pos = list(seq_len(fx$T_m)),
    L = L, prior_variance_scope = scope,
    estimate_prior_variance = TRUE,
    L_greedy = NULL, max_iter = max_iter, tol = 1e-3,
    verbose = FALSE))
}

# ============================================================
# Section 4. Init helper (`init_ebnm_prior_per_scale`)
# ============================================================

test_that("4: init helper returns the documented G_prior shape (Normal)", {
  data <- make_sparse_data()$data
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
  data <- make_sparse_data()$data
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
  sp   <- make_sparse_data(); fx <- sp$fx; data <- sp$data
  out  <- mfsusieR:::init_ebnm_prior_per_scale(
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
  fit <- fit_sparse("per_scale_normal")
  expect_s3_class(fit, "mfsusie")
  expect_s3_class(fit, "susie")
  expect_equal(class(fit$G_prior[[1L]]),
                "mixture_point_normal_per_scale")
  S_m <- length(fit$G_prior[[1L]])
  expect_equal(dim(fit$pi_V[[1L]][[1L]]), c(S_m, 2L))
  expect_true(is.matrix(fit$alpha))
  expect_true(!is.null(fit$pip))
  expect_true(!is.null(fit$sets$cs))
})

test_that("5: end-to-end mfsusie() builds a fit on per_scale_laplace", {
  fit <- fit_sparse("per_scale_laplace")
  expect_s3_class(fit, "mfsusie")
  expect_equal(class(fit$G_prior[[1L]]),
                "mixture_point_laplace_per_scale")
  expect_equal(class(fit$G_prior[[1L]][[1L]]$fitted_g), "laplacemix")
  S_m <- length(fit$G_prior[[1L]])
  expect_equal(dim(fit$pi_V[[1L]][[1L]]), c(S_m, 2L))
})

test_that("5: ebnm M-step writes both fitted_g and pi_V", {
  # Dispatch writes to per-effect sidecars in lockstep:
  #   fit$pi_V[[l]][[m]][s, ] == fit$fitted_g_per_effect[[l]][[m]][[s]]$pi.
  fit <- fit_sparse("per_scale_normal")
  for (l in seq_along(fit$pi_V)) {
    for (s in seq_along(fit$G_prior[[1L]])) {
      expect_equal(fit$pi_V[[l]][[1L]][s, ],
                   fit$fitted_g_per_effect[[l]][[1L]][[s]]$pi,
                   tolerance = 0)
    }
  }
})

# ============================================================
# Section 5e. ebnm wrapper forwarding (mocked)
# ============================================================

test_that("5e: .opv_ebnm_point forwards (x, s, g_init, fix_g) to ebnm correctly", {
  data      <- make_sparse_data()$data
  prior     <- mfsusieR:::mf_prior_scale_mixture(data,
                  prior_variance_scope = "per_scale_normal",
                  null_prior_init = 0)
  ser_stats <- make_ser_stats(data)
  model     <- make_seed_model(data, prior)

  call_log <- list()
  fake_ebnm <- function(x, s, mode, g_init, fix_g, ...) {
    call_log[[length(call_log) + 1L]] <<- list(
      x = x, s = s, g_init = g_init, fix_g = fix_g)
    list(fitted_g = g_init)
  }
  testthat::with_mocked_bindings(
    mfsusieR:::.opv_ebnm_point(
      data, list(estimate_prior_variance = TRUE,
                 alpha_thin_eps = 1e-6),
      model, ser_stats,
      keep_idx = seq_len(data$p),
      ebnm_fn  = ebnm::ebnm_point_normal),
    ebnm_point_normal = fake_ebnm,
    .package = "ebnm")
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
  data      <- make_sparse_data()$data
  prior     <- mfsusieR:::mf_prior_scale_mixture(data,
                  prior_variance_scope = "per_scale_normal",
                  null_prior_init = 0)
  ser_stats <- make_ser_stats(data)
  model     <- make_seed_model(data, prior)

  fix_g_seen <- NA
  fake_ebnm <- function(x, s, mode, g_init, fix_g, ...) {
    fix_g_seen <<- fix_g
    list(fitted_g = g_init)
  }
  testthat::with_mocked_bindings(
    mfsusieR:::.opv_ebnm_point(
      data, list(estimate_prior_variance = FALSE,
                 alpha_thin_eps = 1e-6),
      model, ser_stats,
      keep_idx = seq_len(data$p),
      ebnm_fn  = ebnm::ebnm_point_normal),
    ebnm_point_normal = fake_ebnm,
    .package = "ebnm")
  expect_true(isTRUE(fix_g_seen))
})

test_that("5e: g_init warm-starts across consecutive M-step calls", {
  data      <- make_sparse_data()$data
  prior     <- mfsusieR:::mf_prior_scale_mixture(data,
                  prior_variance_scope = "per_scale_normal",
                  null_prior_init = 0)
  ser_stats <- make_ser_stats(data)
  model     <- make_seed_model(data, prior)
  S_m       <- length(prior$G_prior[[1L]])

  # Each call returns a sentinel `fitted_g` so the next call's
  # `g_init` is detectable as the previous call's return value.
  call_log <- list()
  fake_ebnm <- function(x, s, mode, g_init, fix_g, ...) {
    new_g <- g_init
    new_g$pi <- c(0.42, 0.58)
    call_log[[length(call_log) + 1L]] <<- list(g_init  = g_init,
                                                returned = new_g)
    list(fitted_g = new_g)
  }
  testthat::with_mocked_bindings(
    {
      m1 <- mfsusieR:::.opv_ebnm_point(
        data, list(estimate_prior_variance = TRUE,
                   alpha_thin_eps = 1e-6),
        model, ser_stats,
        keep_idx = seq_len(data$p),
        ebnm_fn  = ebnm::ebnm_point_normal)
      mfsusieR:::.opv_ebnm_point(
        data, list(estimate_prior_variance = TRUE,
                   alpha_thin_eps = 1e-6),
        m1, ser_stats,
        keep_idx = seq_len(data$p),
        ebnm_fn  = ebnm::ebnm_point_normal)
    },
    ebnm_point_normal = fake_ebnm,
    .package = "ebnm")
  for (s in seq_len(S_m)) {
    expect_equal(call_log[[S_m + s]]$g_init,
                 call_log[[s]]$returned)
  }
})

# ============================================================
# Section 6. Cache gate (`refresh_iter_cache`)
# ============================================================

# The iter_cache slot is IBSS hot-path scratch only and is
# stripped by `cleanup_extra_fields.mf_individual` before the fit
# is returned, so these tests exercise
# `refresh_iter_cache.mf_individual` directly on a freshly-built
# model object rather than the cleaned-up fit.

test_that("6: iter_cache skips sdmat / log_sdmat on the ebnm path", {
  ds <- make_sparse_data()
  fit <- fit_sparse("per_scale_normal", fx = ds$fx)
  cached <- mfsusieR:::refresh_iter_cache.mf_individual(ds$data, fit)
  expect_true(!is.null(cached$iter_cache$shat2[[1L]]))
  expect_equal(nrow(cached$iter_cache$shat2[[1L]]), ncol(ds$fx$X))
  expect_null(cached$iter_cache$sdmat)
  expect_null(cached$iter_cache$log_sdmat)
})

test_that("6: iter_cache keeps sdmat / log_sdmat on the mixsqp path", {
  ds  <- make_sparse_data()
  fit <- fit_sparse("per_outcome", fx = ds$fx)
  cached <- mfsusieR:::refresh_iter_cache.mf_individual(ds$data, fit)
  expect_false(is.null(cached$iter_cache$shat2[[1L]]))
  expect_false(is.null(cached$iter_cache$sdmat[[1L]][[1L]]))
  expect_false(is.null(cached$iter_cache$log_sdmat[[1L]][[1L]]))
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

expect_downstream_methods <- function(fit) {
  pr <- predict(fit)
  expect_true(is.numeric(pr) || is.list(pr))
  expect_true(all(is.finite(unlist(pr))))
  cf <- coef(fit)
  expect_true(is.list(cf) || is.numeric(cf))
  expect_no_error(fitted(fit))
  expect_s3_class(summary(fit), "summary.mfsusie")
  expect_no_error(print(fit))
}

test_that("7g.4: predict / coef / fitted / summary / print return non-trivial output on per_scale_normal", {
  expect_downstream_methods(fit_sparse("per_scale_normal",
                                        L = 3L, max_iter = 30L))
})

test_that("7g.4: predict / coef / fitted / summary / print return non-trivial output on per_scale_laplace", {
  expect_downstream_methods(fit_sparse("per_scale_laplace",
                                        L = 3L, max_iter = 30L))
})

# ============================================================
# Section 8.4. Multi-outcome (M >= 2) end-to-end (per spec scenarios)
# ============================================================

test_that("8.4: per_scale_normal recovers shared causal across M=2 outcomes", {
  set.seed(3L)
  n <- 200L; p <- 30L; T_per <- c(32L, 32L)
  X <- matrix(rnorm(n * p), n)
  beta <- numeric(p); beta[c(7L, 18L)] <- c(1.0, -1.0)
  Y_list <- lapply(T_per, function(T_m)
    X %*% matrix(rep(beta, T_m), nrow = p) +
      matrix(rnorm(n * T_m, sd = 0.4), nrow = n))
  set.seed(1L)
  fit <- suppressWarnings(mfsusieR::mfsusie(
    X, Y_list, pos = lapply(T_per, seq_len),
    L = 5L, prior_variance_scope = "per_scale_normal",
    estimate_prior_variance = TRUE,
    L_greedy = NULL, max_iter = 50L, tol = 1e-3,
    verbose = FALSE))
  expect_gte(fit$pip[7L],  0.9)
  expect_gte(fit$pip[18L], 0.9)
  expect_equal(length(fit$sets$cs), 2L)
})

# ============================================================
# Section 7e.4. Wrapper-level fix_g short-circuit at machine precision
# ============================================================

test_that("7e.4: estimate_prior_variance = FALSE leaves G_prior$fitted_g unchanged", {
  data      <- make_sparse_data()$data
  prior     <- mfsusieR:::mf_prior_scale_mixture(data,
                  prior_variance_scope = "per_scale_normal",
                  null_prior_init = 0)
  ser_stats <- make_ser_stats(data)
  model     <- make_seed_model(data, prior)
  before <- lapply(model$G_prior[[1L]],
                    function(g) g$fitted_g)
  m_after <- mfsusieR:::.opv_ebnm_point(
    data, list(estimate_prior_variance = FALSE,
               alpha_thin_eps = 1e-6),
    model, ser_stats,
    keep_idx = seq_len(data$p),
    ebnm_fn  = ebnm::ebnm_point_normal)
  after <- lapply(m_after$G_prior[[1L]],
                   function(g) g$fitted_g)
  for (s in seq_along(before)) {
    expect_equal(after[[s]], before[[s]], tolerance = 1e-12)
  }
})

# ============================================================
# Section 7e.5. Wrapper-level idempotence on identical inputs
# ============================================================

test_that("7e.5: M-step is bit-idempotent on identical inputs (Normal + Laplace)", {
  data <- make_sparse_data()$data
  for (scope in c("per_scale_normal", "per_scale_laplace")) {
    prior     <- mfsusieR:::mf_prior_scale_mixture(data,
                    prior_variance_scope = scope,
                    null_prior_init = 0)
    ser_stats <- make_ser_stats(data)
    model     <- make_seed_model(data, prior)
    ebnm_fn <- if (scope == "per_scale_normal") ebnm::ebnm_point_normal
               else                              ebnm::ebnm_point_laplace
    m1 <- mfsusieR:::.opv_ebnm_point(
      data, list(estimate_prior_variance = TRUE,
                 alpha_thin_eps = 1e-6),
      model, ser_stats,
      keep_idx = seq_len(data$p), ebnm_fn = ebnm_fn)
    m2 <- mfsusieR:::.opv_ebnm_point(
      data, list(estimate_prior_variance = TRUE,
                 alpha_thin_eps = 1e-6),
      model, ser_stats,   # same input
      keep_idx = seq_len(data$p), ebnm_fn = ebnm_fn)
    for (s in seq_along(m1$G_prior[[1L]])) {
      expect_equal(m1$G_prior[[1L]][[s]]$fitted_g,
                   m2$G_prior[[1L]][[s]]$fitted_g,
                   tolerance = 1e-12,
                   info = sprintf("scope=%s, s=%d", scope, s))
    }
  }
})

# ============================================================
# Section 8. End-to-end power tests
# ============================================================

expect_recovers_sparse_signals <- function(scope) {
  fx  <- make_sparse_fixture()
  fit <- fit_sparse(scope, fx = fx)
  expect_gte(fit$pip[fx$signal_idx[1L]], 0.9)
  expect_gte(fit$pip[fx$signal_idx[2L]], 0.9)
  expect_equal(length(fit$sets$cs), 2L)
}

test_that("8.1: per_scale_normal recovers planted causals on sparse fixture", {
  expect_recovers_sparse_signals("per_scale_normal")
})

test_that("8.1: per_scale_laplace recovers planted causals on sparse fixture", {
  expect_recovers_sparse_signals("per_scale_laplace")
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
  # pi_0 = 1: mixture density = dnull, BF per cell = log(pi_0) = 0,
  # sum to 0 over T_idx.
  fg_null <- structure(list(pi = c(1, 0), mean = c(0, 0),
                             scale = c(0, 0.5)), class = "laplacemix")
  expect_equal(
    mfsusieR:::mixture_log_bf_laplace_per_scale(bhat, shat, fg_null),
    rep(0, p), tolerance = 1e-12)
  # lambda = 0 with pi_0 < 1: slab degenerates to delta at zero,
  # so mixture density = (pi_0 + pi_1) * dnull = dnull, BF = 0.
  fg_lam0 <- structure(list(pi = c(0.5, 0.5), mean = c(0, 0),
                             scale = c(0, 0)), class = "laplacemix")
  expect_equal(
    mfsusieR:::mixture_log_bf_laplace_per_scale(bhat, shat, fg_lam0),
    rep(0, p), tolerance = 1e-12)
})

test_that("9b: mixture_log_bf_laplace_per_scale handles single-component laplacemix collapse", {
  # ebnm's MLE collapses to a 1-element laplacemix when pi_0 -> 0
  # (pure slab) or pi_0 -> 1 (pure null). Defend against the
  # `scale[2L] = NA` indexing failure.
  set.seed(1L); p <- 10L; T_idx <- 4L
  bhat <- matrix(rnorm(p * T_idx), p, T_idx)
  shat <- matrix(0.3, p, T_idx)
  # Pure-null collapse: pi = 1, scale = 0.
  fg_null1 <- structure(list(pi = 1, mean = 0, scale = 0),
                         class = "laplacemix")
  expect_equal(
    mfsusieR:::mixture_log_bf_laplace_per_scale(bhat, shat, fg_null1),
    rep(0, p), tolerance = 1e-12)
  # Pure-slab collapse: pi = 1, scale = lambda > 0.
  fg_slab1 <- structure(list(pi = 1, mean = 0, scale = 0.5),
                         class = "laplacemix")
  out <- mfsusieR:::mixture_log_bf_laplace_per_scale(bhat, shat, fg_slab1)
  expect_length(out, p)
  expect_true(all(is.finite(out)))
})

# ============================================================
# Section 10. mixture_null_weight is a no-op on the ebnm paths
# ============================================================

test_that("10: mixture_null_weight has no numerical effect on per_scale_normal (and warns)", {
  fx <- make_sparse_fixture()
  fit_a <- fit_sparse("per_scale_normal", fx = fx)
  fit_b <- suppressWarnings(mfsusieR::mfsusie(
    fx$X, list(fx$Y), pos = list(seq_len(fx$T_m)),
    L = 5L, prior_variance_scope = "per_scale_normal",
    mixture_null_weight = 0.5,
    estimate_prior_variance = TRUE,
    L_greedy = NULL, max_iter = 50L, tol = 1e-3,
    verbose = FALSE))
  expect_equal(fit_a$alpha, fit_b$alpha, tolerance = 1e-12)
  expect_equal(fit_a$pip,   fit_b$pip,   tolerance = 1e-12)

  # The user-facing warning fires at the front-door call.
  expect_message(
    suppressWarnings(mfsusieR::mfsusie(
      fx$X, list(fx$Y), pos = list(seq_len(fx$T_m)),
      L = 2L, prior_variance_scope = "per_scale_normal",
      mixture_null_weight = 0.5,
      estimate_prior_variance = TRUE,
      L_greedy = NULL, max_iter = 5L, tol = 1e-3,
      verbose = FALSE)),
    regexp = "mixture_null_weight.*ignored")
})
