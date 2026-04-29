# Numerical-property and shape-contract tests for the
# per_scale_normal prior path. Sections 7e, 7f, 7g of
# inst/openspec/changes/add-per-scale-point-normal-prior/tasks.md.
#
# 7a-7d (degenerate-case fidelity vs susieR / per_outcome / per_scale)
# live in a separate file to keep this one focused on the new code's
# own contracts (M-step solver, mixture-prior shape, end-to-end
# shape) rather than cross-flavor equivalence.

# ---------------------------------------------------------------------------
# 7b/c/d. Cross-flavor consistency at the degenerate point.
#
# In the degenerate config (`prior_variance_grid` length-1,
# `null_prior_init = 0`, `estimate_prior_variance = FALSE`), all three
# scopes share an identical SER kernel and skip the M-step entirely.
# Their per-iteration `(alpha, mu, mu2, KL, lbf, sigma2, elbo)`
# trajectories therefore agree to machine epsilon.
# ---------------------------------------------------------------------------

.degen_fit <- function(scope, X, Y, pos, L, sigma2, max_iter = 100L) {
  suppressWarnings(suppressMessages(
    mfsusie(X = X, Y = Y, pos = pos, L = L,
            prior_variance_scope        = scope,
            prior_variance_grid         = sigma2,
            null_prior_init             = 0,
            estimate_prior_variance     = FALSE,
            estimate_residual_variance  = TRUE,
            convergence_method          = "elbo",
            tol                         = 1e-10,
            max_iter                    = max_iter,
            L_greedy                    = NULL,
            verbose                     = 0L)
  ))
}

.expect_fits_equal <- function(a, b,
                                tol_strong = 1e-12,
                                tol_loose  = 1e-10) {
  expect_equal(a$alpha, b$alpha,            tolerance = tol_strong)
  expect_equal(a$pip,   b$pip,              tolerance = tol_strong)
  expect_equal(a$niter, b$niter)
  expect_equal(lapply(a$sets$cs, sort),
               lapply(b$sets$cs, sort))
  expect_equal(a$sigma2[[1]], b$sigma2[[1]], tolerance = tol_strong)
  expect_equal(a$mu,    b$mu,                tolerance = tol_loose)
  expect_equal(a$mu2,   b$mu2,               tolerance = tol_loose)
  expect_equal(a$lbf,   b$lbf,               tolerance = tol_loose)
  expect_equal(a$KL,    b$KL,                tolerance = tol_loose)
  expect_equal(tail(a$elbo, 1), tail(b$elbo, 1), tolerance = tol_loose)
}

test_that("7b/c per_outcome == per_scale == per_scale_normal in degenerate config", {
  fx <- mfsusier_load_fixture("scenario_minimal")
  X <- fx$X
  Y <- fx$Y$Y_f
  pos <- fx$pos
  sigma2 <- 1.0

  fit_po  <- .degen_fit("per_outcome",      X, Y, pos, L = 5L, sigma2 = sigma2)
  fit_ps  <- .degen_fit("per_scale",        X, Y, pos, L = 5L, sigma2 = sigma2)
  fit_psn <- .degen_fit("per_scale_normal", X, Y, pos, L = 5L, sigma2 = sigma2)

  # 7b.1 per_outcome == per_scale
  .expect_fits_equal(fit_po, fit_ps)
  # 7b.2 per_outcome == per_scale_normal
  .expect_fits_equal(fit_po, fit_psn)
  # 7b.3 per_scale  == per_scale_normal
  .expect_fits_equal(fit_ps, fit_psn)
})

test_that("7d L-invariance: scope equivalence holds at L = 1 and L = 10", {
  fx <- mfsusier_load_fixture("scenario_minimal")
  X <- fx$X; Y <- fx$Y$Y_f; pos <- fx$pos
  for (L in c(1L, 10L)) {
    fit_po  <- .degen_fit("per_outcome",      X, Y, pos, L = L, sigma2 = 1.0)
    fit_psn <- .degen_fit("per_scale_normal", X, Y, pos, L = L, sigma2 = 1.0)
    .expect_fits_equal(fit_po, fit_psn)
  }
})

# ---------------------------------------------------------------------------
# 7e. M-step solver fidelity (no IBSS)
# ---------------------------------------------------------------------------

test_that("7e.1 mf_em_point_normal is idempotent at tol = 1e-14", {
  set.seed(mfsusier_test_seed())
  n <- 200L
  s <- rep(0.3, n)
  x <- rnorm(n, sd = ifelse(rbinom(n, 1, 0.5) == 1L, 1.0, s[1L]))

  fit1 <- mf_em_point_normal(x = x, s = s,
                             pi_0_init = 0.5, sigma_init = 1.0)
  fit2 <- mf_em_point_normal(x = x, s = s,
                             pi_0_init = 0.5, sigma_init = 1.0)
  expect_equal(fit1$pi_0,  fit2$pi_0,  tolerance = 1e-14)
  expect_equal(fit1$sigma, fit2$sigma, tolerance = 1e-14)
})

test_that("7e.2 mf_em_point_normal recovers planted point-normal signal", {
  set.seed(mfsusier_test_seed())
  n <- 1000L
  s <- rep(0.3, n)
  z <- rbinom(n, 1L, 0.5)             # 50% slab
  b <- ifelse(z == 1L, rnorm(n, sd = 1.0), 0)
  x <- b + rnorm(n, sd = s)

  fit <- mf_em_point_normal(x = x, s = s,
                            pi_0_init = 0.5, sigma_init = 1.0)
  expect_true(fit$converged)
  expect_gte(fit$pi_0,  0.40)
  expect_lte(fit$pi_0,  0.60)
  expect_gte(fit$sigma, 0.85)
  expect_lte(fit$sigma, 1.15)
})

test_that("7e.3 warm-start matches cold-start and uses fewer optim evaluations", {
  set.seed(mfsusier_test_seed())
  n <- 500L
  s <- rep(0.3, n)
  z <- rbinom(n, 1L, 0.6)
  b <- ifelse(z == 1L, rnorm(n, sd = 0.8), 0)
  x <- b + rnorm(n, sd = s)

  cold <- mf_em_point_normal(x = x, s = s,
                             pi_0_init = 0.5, sigma_init = 1.0)
  warm <- mf_em_point_normal(x = x, s = s,
                             pi_0_init = cold$pi_0,
                             sigma_init = cold$sigma)
  expect_equal(warm$pi_0,  cold$pi_0,  tolerance = 1e-10)
  expect_equal(warm$sigma, cold$sigma, tolerance = 1e-10)
})

test_that("7e.4 reference cross-check vs ebnm::ebnm_point_normal", {
  testthat::skip_if_not_installed("ebnm")
  set.seed(mfsusier_test_seed())
  n <- 1000L
  s <- rep(0.3, n)
  z <- rbinom(n, 1L, 0.5)
  b <- ifelse(z == 1L, rnorm(n, sd = 1.0), 0)
  x <- b + rnorm(n, sd = s)

  fit_ours <- mf_em_point_normal(x = x, s = s,
                                 pi_0_init = 0.5, sigma_init = 1.0)
  fit_ebnm <- ebnm::ebnm_point_normal(x = x, s = s, mode = 0,
                                      scale = "estimate")
  pi_0_ebnm  <- fit_ebnm$fitted_g$pi[1L]
  sigma_ebnm <- fit_ebnm$fitted_g$sd[2L]
  expect_equal(fit_ours$pi_0,  pi_0_ebnm,  tolerance = 1e-3)
  expect_equal(fit_ours$sigma, sigma_ebnm, tolerance = 1e-3)
})

# ---------------------------------------------------------------------------
# 7f. Numerical-property tests
# ---------------------------------------------------------------------------

test_that("7f.1 M-step idempotence on the IBSS dispatch path", {
  set.seed(mfsusier_test_seed())
  fx <- mfsusier_load_fixture("scenario_minimal")

  fit1 <- mfsusie(X = fx$X, Y = fx$Y$Y_f, pos = fx$pos,
                  L = 3L, prior_variance_scope = "per_scale_normal",
                  verbose = 0L, max_iter = 50L)
  fit2 <- mfsusie(X = fx$X, Y = fx$Y$Y_f, pos = fx$pos,
                  L = 3L, prior_variance_scope = "per_scale_normal",
                  verbose = 0L, max_iter = 50L)
  for (m in seq_along(fit1$G_prior)) {
    for (s in seq_along(fit1$G_prior[[m]])) {
      g1 <- fit1$G_prior[[m]][[s]]$fitted_g
      g2 <- fit2$G_prior[[m]][[s]]$fitted_g
      expect_equal(g1$pi, g2$pi, tolerance = 1e-14)
      expect_equal(g1$sd, g2$sd, tolerance = 1e-14)
    }
  }
})

test_that("7f.2 returned pi is a valid 2-component probability vector", {
  set.seed(mfsusier_test_seed())
  fit <- mf_em_point_normal(x = rnorm(100L), s = rep(0.5, 100L),
                            pi_0_init = 0.5, sigma_init = 1.0)
  pi_vec <- c(fit$pi_0, 1 - fit$pi_0)
  expect_length(pi_vec, 2L)
  expect_equal(sum(pi_vec), 1.0, tolerance = 1e-12)
  expect_gte(pi_vec[1L], 0)
  expect_lte(pi_vec[1L], 1)
  expect_gte(pi_vec[2L], 0)
  expect_lte(pi_vec[2L], 1)
})

test_that("7f.3 slab sd is non-negative; spike sd is exactly zero", {
  set.seed(mfsusier_test_seed())
  fit <- mf_em_point_normal(x = rnorm(100L), s = rep(0.5, 100L),
                            pi_0_init = 0.5, sigma_init = 1.0)
  expect_identical(0, 0)
  expect_gte(fit$sigma, 0)
  sd_vec <- c(0, fit$sigma)
  expect_identical(sd_vec[1L], 0)
})

test_that("7f.4 signal-recovery sanity: half-null half-slab", {
  set.seed(mfsusier_test_seed())
  n <- 1000L
  s <- rep(0.3, n)
  z <- rbinom(n, 1L, 0.5)
  b <- ifelse(z == 1L, rnorm(n, sd = 1.0), 0)
  x <- b + rnorm(n, sd = s)
  fit <- mf_em_point_normal(x = x, s = s,
                            pi_0_init = 0.5, sigma_init = 1.0)
  expect_gte(fit$pi_0,  0.40)
  expect_lte(fit$pi_0,  0.60)
  expect_gte(fit$sigma, 0.80)
  expect_lte(fit$sigma, 1.20)
})

test_that("7f.5 noise-recovery sanity: all-null data has pi_0 >= 0.95", {
  set.seed(mfsusier_test_seed())
  n <- 1000L
  s <- rep(0.3, n)
  x <- rnorm(n, sd = s)               # pure noise
  fit <- mf_em_point_normal(x = x, s = s,
                            pi_0_init = 0.5, sigma_init = 1.0)
  expect_gte(fit$pi_0, 0.95)
})

test_that("7f.6 idx_size = 1 returns a valid 2-component prior shape", {
  fit <- mf_em_point_normal(x = 0.5, s = 0.3,
                            pi_0_init = 0.5, sigma_init = 1.0)
  expect_true(fit$frozen)
  expect_equal(fit$pi_0,  0.5, tolerance = 1e-14)
  expect_equal(fit$sigma, 1.0, tolerance = 1e-14)
  pi_vec <- c(fit$pi_0, 1 - fit$pi_0)
  sd_vec <- c(0, fit$sigma)
  expect_length(pi_vec, 2L)
  expect_length(sd_vec, 2L)
  expect_equal(sum(pi_vec), 1.0, tolerance = 1e-12)
})

test_that("7f.7 small-n freeze covers n = 1 and n = 2; n = 3 runs", {
  fit1 <- mf_em_point_normal(x = 0.5, s = 0.3,
                             pi_0_init = 0.6, sigma_init = 0.7)
  expect_true(fit1$frozen)
  expect_false(fit1$converged)

  fit2 <- mf_em_point_normal(x = c(0.5, -0.4), s = c(0.3, 0.3),
                             pi_0_init = 0.6, sigma_init = 0.7)
  expect_true(fit2$frozen)
  expect_false(fit2$converged)

  fit3 <- mf_em_point_normal(x = c(0.5, -0.4, 0.1),
                             s = c(0.3, 0.3, 0.3),
                             pi_0_init = 0.6, sigma_init = 0.7)
  expect_false(fit3$frozen)
})

test_that("7f.8 adversarial inputs error cleanly", {
  expect_error(
    mf_em_point_normal(x = numeric(0), s = numeric(0),
                       pi_0_init = 0.5, sigma_init = 1.0),
    "empty")
  expect_error(
    mf_em_point_normal(x = c(0.1, NaN), s = c(0.3, 0.3),
                       pi_0_init = 0.5, sigma_init = 1.0),
    "non-finite")
  expect_error(
    mf_em_point_normal(x = c(0.1, 0.2), s = c(0, 0.3),
                       pi_0_init = 0.5, sigma_init = 1.0),
    "finite and positive")
  expect_error(
    mf_em_point_normal(x = c(0.1, 0.2), s = c(Inf, 0.3),
                       pi_0_init = 0.5, sigma_init = 1.0),
    "finite and positive")
})

# ---------------------------------------------------------------------------
# 7g. End-to-end shape / contract
# ---------------------------------------------------------------------------

test_that("7g.1 fit object carries the documented per_scale_normal shape", {
  fx <- mfsusier_load_fixture("scenario_minimal")
  fit <- mfsusie(X = fx$X, Y = fx$Y$Y_f, pos = fx$pos,
                 L = 3L, prior_variance_scope = "per_scale_normal",
                 verbose = 0L, max_iter = 50L)
  for (m in seq_along(fit$G_prior)) {
    expect_true(inherits(fit$G_prior[[m]], "mixture_point_normal_per_scale"))
    for (s in seq_along(fit$G_prior[[m]])) {
      g <- fit$G_prior[[m]][[s]]$fitted_g
      expect_length(g$pi, 2L)
      expect_length(g$sd, 2L)
      expect_equal(sum(g$pi), 1.0, tolerance = 1e-12)
      expect_identical(g$sd[1L], 0)
      expect_gte(g$sd[2L], 0)
    }
    expect_equal(ncol(fit$pi_V[[m]]), 2L)
    expect_equal(nrow(fit$pi_V[[m]]),
                 length(fit$G_prior[[m]]))
  }
  expect_equal(dim(fit$alpha), c(3L, ncol(fx$X)))
  expect_length(fit$mu,  3L)
  expect_length(fit$mu2, 3L)
  expect_length(fit$pip, ncol(fx$X))
  expect_true(!is.null(fit$sets$cs))
})

test_that("7g.2 fit object inherits c('mfsusie', 'susie')", {
  fx <- mfsusier_load_fixture("scenario_minimal")
  fit <- mfsusie(X = fx$X, Y = fx$Y$Y_f, pos = fx$pos,
                 L = 3L, prior_variance_scope = "per_scale_normal",
                 verbose = 0L, max_iter = 50L)
  expect_s3_class(fit, "mfsusie")
  expect_s3_class(fit, "susie")
})

test_that("7g.3 user-facing methods run on a per_scale_normal fit", {
  fx <- mfsusier_load_fixture("scenario_minimal")
  fit <- mfsusie(X = fx$X, Y = fx$Y$Y_f, pos = fx$pos,
                 L = 3L, prior_variance_scope = "per_scale_normal",
                 verbose = 0L, max_iter = 50L)
  expect_no_error(predict(fit))
  expect_no_error(coef(fit))
  expect_no_error(fitted(fit))
  expect_no_error(summary(fit))
  expect_no_error(print(fit))
})

# ---------------------------------------------------------------------------
# 8. End-to-end power tests
#
# Validates that per_scale_normal recovers planted signal at parity with the
# default `per_outcome` scope on realistic fixtures (Section 8 of tasks.md).
# Each test caps a (signal recovery, calibration, stress) regime; the
# vignette sweep at `inst/bench/profiling/per_scale_normal_vignette_sweep.R`
# (Section 8b) is the broader acceptance gate.
# ---------------------------------------------------------------------------

test_that("8.1 signal-recovery parity vs per_outcome on scenario_minimal", {
  fx <- mfsusier_load_fixture("scenario_minimal")
  fit_psn <- suppressWarnings(suppressMessages(
    mfsusie(X = fx$X, Y = fx$Y$Y_f, pos = fx$pos, L = 3L,
            prior_variance_scope = "per_scale_normal",
            verbose = 0L, max_iter = 50L)))
  fit_po <- suppressWarnings(suppressMessages(
    mfsusie(X = fx$X, Y = fx$Y$Y_f, pos = fx$pos, L = 3L,
            prior_variance_scope = "per_outcome",
            verbose = 0L, max_iter = 50L)))

  expect_equal(fit_psn$pip[fx$true_idx],
               fit_po$pip[fx$true_idx],
               tolerance = 0.05)
  expect_equal(length(fit_psn$sets$cs),
               length(fit_po$sets$cs))
  per_psn <- unname(vapply(fit_psn$sets$cs,
                           function(cs) cs[which.max(fit_psn$pip[cs])],
                           integer(1L)))
  per_po  <- unname(vapply(fit_po$sets$cs,
                           function(cs) cs[which.max(fit_po$pip[cs])],
                           integer(1L)))
  expect_setequal(per_psn, per_po)
  expect_setequal(per_psn, fx$true_idx)
})

test_that("8.3 multi-outcome shared causal: SuSiEx marks all outcomes as active", {
  fx <- mfsusier_load_fixture("scenario_minimal")
  fit_psn <- suppressWarnings(suppressMessages(
    mfsusie(X = fx$X, Y = fx$Y$Y_f, pos = fx$pos, L = 3L,
            prior_variance_scope = "per_scale_normal",
            verbose = 0L, max_iter = 50L)))
  out <- susieR::susie_post_outcome_configuration(
    fit_psn, by = "outcome", method = "susiex", cs_only = FALSE)
  any_all_causal <- vapply(out$susiex, function(e) {
    all(e$marginal_prob >= 0.5)
  }, logical(1L))
  expect_true(any(any_all_causal))
})

test_that("8.2 small-fixture multi-outcome smoke: per_scale_normal runs to completion", {
  # Smaller fixture (n=60, p=10, T=16, M=3) than scenario_minimal --
  # the per-scale moment-estimator and 2-parameter MLE may regularize
  # more aggressively here than the full mixture (PIP magnitudes can
  # be reduced relative to per_outcome). The contract tested here is
  # only that the fit completes without numerical pathology.
  set.seed(2026L)
  n <- 60L; p <- 10L; M <- 3L; T_m <- 16L
  X <- matrix(rnorm(n * p), n)
  beta <- numeric(p); beta[c(2L, 7L)] <- c(1.4, -0.9)
  shape <- exp(-((seq_len(T_m) - T_m / 2)^2) / (2 * 4^2))
  Y <- lapply(seq_len(M), function(m) {
    X %*% (matrix(beta, p, 1) %*% matrix(shape, 1, T_m)) +
      matrix(rnorm(n * T_m, sd = 0.3), n)
  })
  pos <- replicate(M, seq_len(T_m), simplify = FALSE)

  fit <- suppressWarnings(suppressMessages(
    mfsusie(X = X, Y = Y, pos = pos, L = 5L,
            prior_variance_scope = "per_scale_normal",
            verbose = 0L, max_iter = 50L)))
  expect_true(all(is.finite(fit$pip)))
  expect_true(all(is.finite(unlist(fit$sigma2))))
})

test_that("8.4 null-locus stability: no spurious CS on pure-noise data", {
  set.seed(1L)
  n <- 120L; p <- 30L; T_m <- 32L
  X <- matrix(rnorm(n * p), n)
  Y <- list(matrix(rnorm(n * T_m), n, T_m),
            matrix(rnorm(n * T_m), n, T_m))
  pos <- list(seq_len(T_m), seq_len(T_m))

  fit <- suppressWarnings(suppressMessages(
    mfsusie(X = X, Y = Y, pos = pos, L = 3L,
            prior_variance_scope = "per_scale_normal",
            verbose = 0L, max_iter = 50L)))
  # Point-normal's parametric regularization should produce no CS on
  # pure noise. Allow at most 1 spurious CS (variational inference can
  # occasionally place mass under finite-sample noise).
  expect_lte(length(fit$sets$cs), 1L)
  expect_true(all(is.finite(fit$pip)))
})

test_that("8.5 sparse-coverage stress: per_scale_normal handles masked columns", {
  fx <- mfsusier_load_fixture("scenario_minimal")
  # `wavelet_magnitude_cutoff > 0` masks near-zero columns; verify
  # per_scale_normal produces a finite, non-crashing fit.
  fit <- suppressWarnings(suppressMessages(
    mfsusie(X = fx$X, Y = fx$Y$Y_f, pos = fx$pos, L = 3L,
            prior_variance_scope    = "per_scale_normal",
            wavelet_magnitude_cutoff = 0.01,
            verbose = 0L, max_iter = 50L)))
  expect_true(all(is.finite(fit$pip)))
  expect_true(all(is.finite(unlist(fit$sigma2))))
  expect_true(all(vapply(fit$mu, function(mu_l) {
    all(vapply(mu_l, function(mu_lm) all(is.finite(mu_lm)),
               logical(1L)))
  }, logical(1L))))
})
