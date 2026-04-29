# PR group 6: variance components, mixture-weight update, ELBO.
#
# Apple-to-apple comparisons:
#   - per_outcome variance update vs the per-position formula
#     in mvf.susie.alpha::estimate_residual_variance (R/computational_routine.R#L395-L431).
#   - mixsqp output vs fsusieR::scale_m_step on a canonical fixture
#     (machine precision; the cpp11 kernels and the ported helpers
#     evaluate the same closed form mixsqp).
#
# Manuscript references:
#   methods/derivation.tex eq:residual_variance_per_scale
#   methods/derivation.tex line 216 (mixture-weight update)
#   methods/online_method.tex (Eloglik, ELBO aggregate)

# ---- ER2 helper ----------------------------------------------------

make_data_for_ev <- function(n = 30, J = 8, T_per_outcome = c(64L, 128L)) {
  set.seed(mfsusier_test_seed())
  X <- matrix(rnorm(n * J), nrow = n)
  Y <- lapply(T_per_outcome, function(T_m) matrix(rnorm(n * T_m), nrow = n))
  mfsusieR:::create_mf_individual(wavelet_qnorm = FALSE, X = X, Y = Y, verbose = FALSE)
}

make_model_with_post <- function(data, L = 2, sigma2_scalar = 1) {
  prior <- mfsusieR:::mf_prior_scale_mixture(
    data, prior_variance_grid = c(0.1, 0.5, 2), null_prior_init = 0.5
  )
  params <- list(L = L, prior_weights = NULL, prior = prior,
                 cross_outcome_prior = NULL,
                 residual_variance = lapply(seq_len(data$M),
                                            function(m) sigma2_scalar))
  var_y <- lapply(seq_len(data$M),
                  function(m) rep(sigma2_scalar, data$T_basis[m]))
  model <- mfsusieR:::initialize_susie_model.mf_individual(data, params, var_y)
  # Inject deterministic non-zero alpha + mu so ER2 / Eloglik are non-trivial.
  for (l in seq_len(L)) {
    model$alpha[l, ] <- runif(data$p)
    model$alpha[l, ] <- model$alpha[l, ] / sum(model$alpha[l, ])
    for (m in seq_len(data$M)) {
      model$mu[[l]][[m]]  <- matrix(rnorm(data$p * data$T_basis[m], sd = 0.1),
                                    nrow = data$p)
      model$mu2[[l]][[m]] <- model$mu[[l]][[m]]^2 +
                             matrix(runif(data$p * data$T_basis[m], 0.05, 0.2),
                                    nrow = data$p)
    }
  }
  # Sync the running fit cache with the injected posterior so
  # `mf_get_ER2_per_position` (which reads `model$fitted[[m]]`) is
  # consistent with `alpha * mu`.
  for (m in seq_len(data$M)) {
    postF <- matrix(0, nrow = data$p, ncol = data$T_basis[m])
    for (l in seq_len(L)) {
      postF <- postF + model$alpha[l, ] * model$mu[[l]][[m]]
    }
    model$fitted[[m]] <- data$X %*% postF
  }
  model
}

test_that("mf_get_ER2_per_position matches the susieR-style formula", {
  data  <- make_data_for_ev()
  model <- make_model_with_post(data)

  for (m in seq_len(data$M)) {
    er2 <- mfsusieR:::mf_get_ER2_per_position(data, model, m)

    # Hand-recompute.
    pw     <- data$xtx_diag
    postF  <- matrix(0, nrow = data$p, ncol = data$T_basis[m])
    for (l in seq_len(model$L)) {
      postF <- postF + model$alpha[l, ] * model$mu[[l]][[m]]
    }
    res    <- data$D[[m]] - data$X %*% postF
    rss_t  <- colSums(res^2)
    bias_t <- numeric(data$T_basis[m])
    for (l in seq_len(model$L)) {
      a   <- model$alpha[l, ]
      Bl  <- a * model$mu[[l]][[m]]
      XBl <- data$X %*% Bl
      bias_t <- bias_t + colSums((a * pw) * model$mu2[[l]][[m]]) -
                colSums(XBl^2)
    }
    expect_equal(er2, rss_t + bias_t, tolerance = 1e-12)
  }
})

# ---- update_variance_components --------------------------------

test_that("update_variance_components per_outcome returns a scalar per modality", {
  data  <- make_data_for_ev()
  model <- make_model_with_post(data)
  params <- list(residual_variance_scope = "per_outcome")

  out <- mfsusieR:::update_variance_components.mf_individual(data, params, model)
  for (m in seq_len(data$M)) {
    expect_length(out$sigma2[[m]], 1L)
    er2 <- mfsusieR:::mf_get_ER2_per_position(data, model, m)
    expect_equal(out$sigma2[[m]],
                 sum(er2) / (data$n * data$T_basis[m]),
                 tolerance = 1e-12)
  }
})

test_that("update_variance_components per_scale returns length-S_m vector per modality", {
  data  <- make_data_for_ev()
  model <- make_model_with_post(data)
  params <- list(residual_variance_scope = "per_scale")

  out <- mfsusieR:::update_variance_components.mf_individual(data, params, model)
  for (m in seq_len(data$M)) {
    indx <- data$scale_index[[m]]
    expect_length(out$sigma2[[m]], length(indx))
    er2 <- mfsusieR:::mf_get_ER2_per_position(data, model, m)
    for (s in seq_along(indx)) {
      expect_equal(out$sigma2[[m]][s],
                   sum(er2[indx[[s]]]) / (data$n * length(indx[[s]])),
                   tolerance = 1e-12)
    }
  }
})

test_that("update_variance_components defaults to per_scale when method is unset", {
  data  <- make_data_for_ev()
  model <- make_model_with_post(data)
  out_default <- mfsusieR:::update_variance_components.mf_individual(data, list(), model)
  for (m in seq_len(data$M)) {
    expect_length(out_default$sigma2[[m]], length(data$scale_index[[m]]))
  }
})

# ---- update_model_variance ------------------------------------

test_that("update_model_variance returns valid pi vectors per (m, s) for one effect", {
  data  <- make_data_for_ev()
  model <- make_model_with_post(data)
  model <- mfsusieR:::compute_residuals.mf_individual(data, NULL, model, l = 1)
  ser   <- mfsusieR:::compute_ser_statistics.mf_individual(data, NULL, model, l = 1)

  out <- mfsusieR:::update_model_variance.mf_individual(data, list(), model, ser, l = 1)

  for (m in seq_len(data$M)) {
    for (s in seq_along(data$scale_index[[m]])) {
      pi_ms <- out$G_prior[[m]][[s]]$fitted_g$pi
      expect_equal(sum(pi_ms), 1, tolerance = 1e-8)
      expect_true(all(pi_ms >= 0))
      expect_equal(out$pi_V[[m]][s, ], pi_ms, tolerance = 0)
    }
  }
})

# ---- Eloglik / get_objective ---------------------------------

test_that("Eloglik aggregates per-modality, per-position log-likelihood", {
  data  <- make_data_for_ev()
  model <- make_model_with_post(data)

  out <- mfsusieR:::Eloglik.mf_individual(data, model)

  expected <- 0
  for (m in seq_len(data$M)) {
    er2 <- mfsusieR:::mf_get_ER2_per_position(data, model, m)
    s2  <- mfsusieR:::mf_sigma2_per_position(data, model, m)
    expected <- expected + sum(-0.5 * data$n * log(2 * pi * s2) - 0.5 * er2 / s2)
  }
  expect_equal(out, expected, tolerance = 1e-12)
})

test_that("get_objective dispatches to susieR's default and returns Eloglik - sum(KL)", {
  # mfsusieR no longer overrides get_objective; the standard branch in
  # susieR's `get_objective.default` computes
  # `Eloglik(data, model) - sum(model$KL, na.rm = TRUE)` and dispatches
  # `Eloglik` per class, so `Eloglik.mf_individual` still fires.
  data  <- make_data_for_ev()
  model <- make_model_with_post(data)
  model$KL <- runif(model$L, 0.1, 1)
  class(data) <- c("mf_individual", "individual")

  out <- susieR:::get_objective(data, list(use_NIG = FALSE), model)

  el <- mfsusieR:::Eloglik.mf_individual(data, model)
  expect_equal(out, el - sum(model$KL), tolerance = 1e-12)
})

# ---- S3 registration on susieR namespace ---------------------

test_that("update_variance_components / update_model_variance / Eloglik registered on susieR generics", {
  # `get_objective` was deliberately removed from the per-class
  # registration list when the trivial mfsusieR override was deleted in
  # favour of dispatching to susieR's `get_objective.default`.
  for (g in c("update_variance_components", "update_model_variance",
              "Eloglik")) {
    method_fn <- getS3method(g, "mf_individual",
                             optional = TRUE,
                             envir = asNamespace("susieR"))
    expect_true(is.function(method_fn),
                info = paste("Missing registration:", g, ".mf_individual"))
  }
})

# ---- Comprehensive: per-scale sigma2 path through Eloglik / get_objective ----

test_that("Eloglik with per-scale sigma2 broadcasts via mf_sigma2_per_position", {
  data  <- make_data_for_ev()
  model <- make_model_with_post(data)
  for (m in seq_len(data$M)) {
    S_m <- length(data$scale_index[[m]])
    model$sigma2[[m]] <- seq_len(S_m) * 0.1
  }

  out <- mfsusieR:::Eloglik.mf_individual(data, model)

  expected <- 0
  for (m in seq_len(data$M)) {
    er2 <- mfsusieR:::mf_get_ER2_per_position(data, model, m)
    s2  <- mfsusieR:::mf_sigma2_per_position(data, model, m)
    expected <- expected + sum(-0.5 * data$n * log(2 * pi * s2) - 0.5 * er2 / s2)
  }
  expect_equal(out, expected, tolerance = 1e-12)
  expect_true(is.finite(out))
})

test_that("get_objective with per-scale sigma2 + non-uniform model$pi", {
  data  <- make_data_for_ev()
  model <- make_model_with_post(data)
  class(data) <- c("mf_individual", "individual")
  for (m in seq_len(data$M)) {
    S_m <- length(data$scale_index[[m]])
    model$sigma2[[m]] <- seq_len(S_m) * 0.1
  }
  # Non-uniform variable-selection prior.
  model$pi <- runif(data$p)
  model$pi <- model$pi / sum(model$pi)
  model$KL <- runif(model$L, 0.1, 1)

  out <- susieR:::get_objective(data, list(use_NIG = FALSE), model)
  expect_true(is.finite(out))

  el <- mfsusieR:::Eloglik.mf_individual(data, model)
  expect_equal(out, el - sum(model$KL), tolerance = 1e-12)
})

# ---- Edge cases ------------------------------------------------

test_that("M = 1 single modality: variance / Eloglik / objective all return finite scalars", {
  data  <- make_data_for_ev(T_per_outcome = 64L)
  expect_identical(data$M, 1L)
  model <- make_model_with_post(data, L = 2)

  v_shared <- mfsusieR:::update_variance_components.mf_individual(
    data, list(residual_variance_scope = "per_outcome"), model)
  expect_length(v_shared$sigma2[[1]], 1L)

  v_per <- mfsusieR:::update_variance_components.mf_individual(
    data, list(residual_variance_scope = "per_scale"), model)
  expect_length(v_per$sigma2[[1]], length(data$scale_index[[1]]))

  el <- mfsusieR:::Eloglik.mf_individual(data, model)
  expect_true(is.finite(el))

  class(data) <- c("mf_individual", "individual")
  obj <- susieR:::get_objective(data, list(use_NIG = FALSE), model)
  expect_true(is.finite(obj))
})

test_that("L = 1 single effect: ER2 reduces to RSS + single bias correction", {
  data  <- make_data_for_ev()
  model <- make_model_with_post(data, L = 1)
  for (m in seq_len(data$M)) {
    er2 <- mfsusieR:::mf_get_ER2_per_position(data, model, m)
    res <- data$D[[m]] - model$fitted[[m]]
    rss <- colSums(res * res)
    a   <- model$alpha[1, ]
    Bl  <- a * model$mu[[1]][[m]]
    XBl <- data$X %*% Bl
    bias <- colSums((a * data$xtx_diag) * model$mu2[[1]][[m]]) -
            colSums(XBl * XBl)
    expect_equal(er2, rss + bias, tolerance = 1e-12)
  }
})

test_that("All-zero alpha (no effects): ER2 = ||D||^2 per position", {
  # When alpha = 0 (or mu = 0), bias correction is 0 and ER2 = ||D||^2 per t.
  data  <- make_data_for_ev()
  model <- make_model_with_post(data)
  for (l in seq_len(model$L)) model$alpha[l, ] <- 0
  for (m in seq_len(data$M)) {
    model$fitted[[m]] <- matrix(0, nrow = data$n, ncol = data$T_basis[m])
  }

  for (m in seq_len(data$M)) {
    er2 <- mfsusieR:::mf_get_ER2_per_position(data, model, m)
    expect_equal(er2, colSums(data$D[[m]]^2), tolerance = 1e-12)
  }
})

# ---- optimize_prior_variance: non-trivial pi update, fidelity vs hand-roll ----

test_that("optimize_prior_variance.mf_individual changes pi vs the initial uniform prior", {
  data  <- make_data_for_ev()
  model <- make_model_with_post(data)
  model <- mfsusieR:::compute_residuals.mf_individual(data, NULL, model, l = 1)
  ser   <- mfsusieR:::compute_ser_statistics.mf_individual(data, NULL, model, l = 1)

  # Inject a more peaked alpha at SNP 1 so mixsqp has a clear signal.
  model$alpha[1, ]    <- rep(1e-6, data$p)
  model$alpha[1, 1]   <- 1 - 1e-6 * (data$p - 1)

  init_pi <- vapply(seq_along(model$G_prior[[1]]),
                    function(s) model$G_prior[[1]][[s]]$fitted_g$pi,
                    numeric(length(model$G_prior[[1]][[1]]$fitted_g$pi)))
  out_list <- mfsusieR:::optimize_prior_variance.mf_individual(
    data, list(), model, ser, l = 1)
  expect_named(out_list, c("V", "model"))
  expect_equal(out_list$V, 1)
  out <- out_list$model
  new_pi <- vapply(seq_along(out$G_prior[[1]]),
                   function(s) out$G_prior[[1]][[s]]$fitted_g$pi,
                   numeric(length(out$G_prior[[1]][[1]]$fitted_g$pi)))
  expect_false(isTRUE(all.equal(init_pi, new_pi)))
})

test_that("update_model_variance.mf_individual orchestrates sigma2 update + derived quantities", {
  data  <- make_data_for_ev()
  model <- make_model_with_post(data)
  params <- list(estimate_residual_variance = TRUE,
                 residual_variance_scope = "per_outcome")

  out <- mfsusieR:::update_model_variance.mf_individual(data, params, model)

  # sigma2 was refreshed to a (length-1) scalar per modality.
  for (m in seq_len(data$M)) {
    expect_length(out$sigma2[[m]], 1L)
    expect_true(out$sigma2[[m]] > 0)
  }
  # Running fit was synced to current alpha * mu.
  for (m in seq_len(data$M)) {
    postF <- matrix(0, nrow = data$p, ncol = data$T_basis[m])
    for (l in seq_len(model$L)) {
      postF <- postF + model$alpha[l, ] * model$mu[[l]][[m]]
    }
    expect_equal(out$fitted[[m]], data$X %*% postF, tolerance = 1e-12)
  }
})

test_that("update_model_variance is a no-op when estimate_residual_variance = FALSE", {
  data  <- make_data_for_ev()
  model <- make_model_with_post(data)
  out <- mfsusieR:::update_model_variance.mf_individual(
    data, list(estimate_residual_variance = FALSE), model)
  expect_identical(out, model)
})

test_that("mf_em_m_step_per_scale matches a hand-rolled mixsqp invocation at <= 1e-12", {
  set.seed(mfsusier_test_seed())
  J        <- 30
  idx_size <- 8
  K        <- 4
  bhat <- matrix(rnorm(J * idx_size, sd = 0.5), nrow = J)
  shat <- matrix(runif(J * idx_size, 0.5, 1.5), nrow = J)
  sd_grid <- c(0, 0.5, 1.0, 2.0)
  zeta    <- runif(J); zeta <- zeta / sum(zeta)

  L_mat <- mfsusieR:::mf_em_likelihood_per_scale(bhat, shat, sd_grid)
  ours  <- mfsusieR:::mf_em_m_step_per_scale(L_mat, zeta, idx_size,
                                             null_weight = 0.7)

  # Hand mixsqp (same arguments, same control).
  w <- c(0.7 * idx_size, rep(zeta, idx_size))
  hand <- mixsqp::mixsqp(L_mat, w,
                        x0 = c(0.5, rep(1e-6, K - 1)),
                        log = FALSE,
                        control = list(verbose = FALSE))$x
  if (hand[1] > 1 - 0.001) {
    hand    <- numeric(K)
    hand[1] <- 1
  }
  expect_equal(ours, hand, tolerance = 1e-12)
})

# ---- update_variance_components on per-scale sigma2 already initialized ----

test_that("update_variance_components per_scale preserves shape when sigma2 was already a vector", {
  data  <- make_data_for_ev()
  model <- make_model_with_post(data)
  for (m in seq_len(data$M)) {
    model$sigma2[[m]] <- rep(0.5, length(data$scale_index[[m]]))
  }
  out <- mfsusieR:::update_variance_components.mf_individual(
    data, list(residual_variance_scope = "per_scale"), model)
  for (m in seq_len(data$M)) {
    expect_length(out$sigma2[[m]], length(data$scale_index[[m]]))
    expect_true(all(out$sigma2[[m]] > 0))
  }
})

# ---- mf_em_likelihood_per_scale shape and properties ----------

test_that("mf_em_likelihood_per_scale returns positive matrix with penalty row", {
  bhat <- matrix(rnorm(20), 5, 4)
  shat <- matrix(runif(20, 0.5, 1.5), 5, 4)
  sd_grid <- c(0, 0.3, 1.0)
  L <- mfsusieR:::mf_em_likelihood_per_scale(bhat, shat, sd_grid)

  expect_identical(dim(L), c(5L * 4L + 1L, 3L))     # (J*|idx|+1) x K
  expect_equal(L[1, ], c(100, 0, 0))                 # null-component penalty row
  expect_true(all(L[-1, ] > 0))                      # all real rows positive (likelihoods)
})

test_that("mf_em_likelihood_per_scale handles edge case sd_grid = 0 (pure null)", {
  bhat <- matrix(rnorm(8), 4, 2)
  shat <- matrix(runif(8, 0.5, 1.5), 4, 2)
  L <- mfsusieR:::mf_em_likelihood_per_scale(bhat, shat, sd_grid = 0)
  expect_identical(dim(L), c(4L * 2L + 1L, 1L))
  expect_equal(L[1, 1], 100)
})

test_that("mf_em_likelihood_per_scale NA-imputation path (Shat = 0 produces Inf, then median fill)", {
  set.seed(mfsusier_test_seed())
  bhat <- matrix(rnorm(20), 5, 4)
  shat <- matrix(runif(20, 0.5, 1.5), 5, 4)
  shat[1, 1] <- 0   # Force NA in dnorm by zero Shat
  sd_grid <- c(0, 0.5, 1)
  L <- mfsusieR:::mf_em_likelihood_per_scale(bhat, shat, sd_grid)
  expect_true(all(is.finite(L)))
  expect_identical(dim(L), c(5L * 4L + 1L, 3L))
})

# ---- Input validation -------------------------------------------

test_that("update_variance_components errors on unrecognized method", {
  data  <- make_data_for_ev()
  model <- make_model_with_post(data)
  expect_error(
    mfsusieR:::update_variance_components.mf_individual(
      data, list(residual_variance_scope = "bogus_mode"), model),
    "must be one of"
  )
})

# ---- ELBO monotonicity at machine precision (C-4.4 watertight) --
#
# Derivation: inst/notes/cross-package-audit-derivations.md section 2.
# The susieR / mfsusieR ELBO `Eloglik - sum(KL)` is the standard
# variational decomposition. Any correctly-implemented variational EM
# must be monotone non-decreasing across iterations. We assert this
# at machine precision (with a small floating-point summation slack)
# on a multi-iteration fit. Post-iter-1 only: the initial sigma2
# (var(Y) guess) is replaced by its first closed-form update, which
# is a one-shot non-monotone step. The audit test floor is `1e-10`,
# four orders of magnitude tighter than the smoke test's `1e-6`.

test_that("mfsusie ELBO is non-decreasing post-iter-1 at machine precision (<= 1e-10)", {
  set.seed(mfsusier_test_seed())
  n <- 40; J <- 10; T_per_outcome <- c(64L, 32L)
  X <- matrix(rnorm(n * J), nrow = n)
  beta <- numeric(J); beta[1] <- 1; beta[3] <- -0.5
  Y <- lapply(T_per_outcome, function(T_m) {
    eta <- X %*% beta
    matrix(rep(eta, T_m), nrow = n) + matrix(rnorm(n * T_m, sd = 0.3), nrow = n)
  })

  # max_iter = 200 because this small-n, under-determined fixture
  # (n=40, J=10, L=5, only 2 true signals) is the slow-convergence
  # regime: each iter improves ELBO by O(0.03) so settling below
  # tol = 1e-4 takes ~120 iters. The slow geometric decay is inherent
  # to the mixsqp-driven mixture-prior M-step on small data; cf the
  # ELBO trace `0.0317, 0.0305, 0.0293, 0.0282, ...`. Convergence is
  # correct, just slow; we want max_iter > niter so the convergence
  # warning doesn't fire and so the monotonicity check has real data
  # past the iter-1 transient.
  fit <- mfsusie(wavelet_qnorm = FALSE, X, Y, L = 5, max_iter = 200, verbose = FALSE)

  # Need at least three iterations to make the post-iter-1 check
  # meaningful (drop the first diff, then test the remainder).
  expect_true(length(fit$elbo) >= 3L)
  diffs_post1 <- diff(fit$elbo)[-1]
  # Machine-precision tolerance with a small slack for FP summation
  # noise. A real ELBO bug would show up as a downward step orders
  # of magnitude larger than 1e-10.
  expect_true(all(diffs_post1 >= -1e-10),
              info = paste("ELBO non-monotone past iter 1:",
                           paste(format(diffs_post1, digits = 3),
                                 collapse = " ")))
})

# ---- loglik.mf_individual zero-pw branch ------------------------

test_that("init_scale_mixture_prior_default errors on NULL groups", {
  set.seed(mfsusier_test_seed())
  expect_error(
    mfsusieR:::init_scale_mixture_prior_default(
      Y_m = matrix(rnorm(20), 5),
      X = matrix(rnorm(20), 5),
      prior_class = "mixture_normal_per_scale",
      groups = NULL),
    "groups"
  )
})

test_that("initialize_susie_model.mf_individual reads cross_outcome_prior from params when supplied", {
  data <- make_data_for_ev()
  prior <- mfsusieR:::mf_prior_scale_mixture(
    data, prior_variance_grid = c(0.1, 0.5), null_prior_init = 0.5)
  custom_xmod <- mfsusieR:::cross_outcome_prior_independent()
  class(custom_xmod) <- c("custom_combiner", class(custom_xmod))

  params <- list(L = 2, prior_weights = NULL, prior = prior,
                 cross_outcome_prior = custom_xmod,
                 residual_variance = lapply(seq_len(data$M), function(m) 1))
  var_y <- lapply(seq_len(data$M),
                  function(m) rep(1, data$T_basis[m]))
  model <- mfsusieR:::initialize_susie_model.mf_individual(data, params, var_y)
  expect_identical(model$cross_outcome_combiner, custom_xmod)
})

test_that("mf_prior_scale_mixture errors when data is not an mf_individual object", {
  expect_error(
    mfsusieR:::mf_prior_scale_mixture(data = list(M = 1)),
    "must be an mf_individual"
  )
})

test_that("mf_prior_scale_mixture per_outcome scope on data-driven path returns single-row pi", {
  set.seed(mfsusier_test_seed())
  data <- mfsusieR:::create_mf_individual(wavelet_qnorm = FALSE, 
    X = matrix(rnorm(20 * 4), 20),
    Y = list(matrix(rnorm(20 * 16), 20)),
    verbose = FALSE
  )
  prior <- mfsusieR:::mf_prior_scale_mixture(
    data, prior_variance_grid = NULL,
    prior_variance_scope = "per_outcome",
    null_prior_init = 0.5
  )
  for (m in seq_len(data$M)) {
    expect_identical(nrow(prior$pi[[m]]), 1L)
  }
})

test_that("mf_invert_dwt accepts a vector D_packed (single-row coercion)", {
  # Vector input -> matrix(D_packed, nrow = 1) coercion path.
  D_packed <- rep(0.5, 8)   # length-8 vector, no dim attribute
  out <- mfsusieR:::mf_invert_dwt(D_packed,
                                  column_center = rep(0, 8),
                                  column_scale  = rep(1, 8))
  # Either a numeric vector or a 1xT matrix is acceptable; key
  # contract is that the vector input did not error out.
  expect_true(is.numeric(out))
})

test_that("init_scale_mixture_prior_default with one group returns single-element G_prior", {
  set.seed(mfsusier_test_seed())
  out <- mfsusieR:::init_scale_mixture_prior_default(
    Y_m = matrix(rnorm(50), 10),
    X = matrix(rnorm(50), 10),
    prior_class = "mixture_normal",
    groups      = list(seq_len(5)))
  expect_identical(length(out$G_prior), 1L)
  expect_identical(attr(out$G_prior, "class"), "mixture_normal")
})

test_that("loglik handles zero-predictor-weight SNPs by zeroing their lbf", {
  data <- make_data_for_ev()
  # Inject a constant column in X so its colSum-of-squares = 0 after centering.
  data$X[, 1] <- 0
  data$xtx_diag <- colSums(data$X^2)
  expect_equal(data$xtx_diag[1], 0)

  model <- make_model_with_post(data)
  model <- mfsusieR:::compute_residuals.mf_individual(data, NULL, model, 1)
  ser   <- mfsusieR:::compute_ser_statistics.mf_individual(data, NULL, model, 1)

  out <- mfsusieR:::loglik.mf_individual(data, NULL, model, V = 1, ser, l = 1)
  # SNP 1 has zero predictor weight; lbf_variable should be 0 there.
  expect_equal(out$lbf_variable[1, 1], 0, tolerance = 1e-12)
})
