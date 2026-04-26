# Tests for `R/model_class.R::initialize_susie_model.mf_individual`
# and `R/model_methods.R` per-effect accessors.
#
# Apple-to-apple comparison context:
#   mvsusieR/R/mvsusie_constructors.R:341-427 (paradigm reference for
#     the model-class constructor signature and field shape).
#   susieR/R/individual_data_methods.R:initialize_susie_model.individual
#     (paradigm reference for the dispatch contract).
#
# This PR group does NOT run the IBSS loop. Tests cover the initial
# state of the model and the per-effect getter / setter contracts.
# End-to-end fit fidelity (contracts C1, C2, C3) lands in PR groups
# 5-7.

# ---- Helpers -----------------------------------------------------------

make_data <- function(n = 20, J = 8, T_per_outcome = c(64L, 128L)) {
  set.seed(mfsusier_test_seed())
  X <- matrix(rnorm(n * J), nrow = n)
  Y <- lapply(T_per_outcome, function(T_m) matrix(rnorm(n * T_m), nrow = n))
  mfsusieR:::create_mf_individual(X = X, Y = Y, verbose = FALSE)
}

make_params <- function(L = 3, J = 8) {
  list(
    L                 = L,
    prior_weights     = rep(1 / J, J),
    prior             = NULL,             # populated in PR group 4
    residual_variance = NULL              # falls back to var_y
  )
}

make_var_y <- function(M, T_basis) {
  # Per-modality variance of y, used as sigma2 fallback.
  lapply(seq_len(M), function(m) rep(1, T_basis[m]))
}

# ---- Initialization shape contract -------------------------------------

test_that("initialize_susie_model.mf_individual produces correct class", {
  data <- make_data()
  params <- make_params(L = 3, J = data$p)
  var_y <- make_var_y(data$M, data$T_basis)

  model <- mfsusieR:::initialize_susie_model.mf_individual(data, params, var_y)

  expect_identical(class(model), c("mfsusie", "susie"))
})

test_that("alpha is initialized uniformly to 1/J for every effect", {
  data <- make_data(n = 20, J = 8)
  params <- make_params(L = 3, J = data$p)
  var_y <- make_var_y(data$M, data$T_basis)

  model <- mfsusieR:::initialize_susie_model.mf_individual(data, params, var_y)

  expect_identical(dim(model$alpha), c(3L, 8L))
  expect_equal(model$alpha, matrix(1 / 8, 3, 8), tolerance = 0)
})

test_that("mu and mu2 are nested lists L x M of zero matrices with the right shape", {
  data <- make_data(n = 20, J = 8, T_per_outcome = c(64L, 128L))
  params <- make_params(L = 3, J = data$p)
  var_y <- make_var_y(data$M, data$T_basis)

  model <- mfsusieR:::initialize_susie_model.mf_individual(data, params, var_y)

  expect_identical(length(model$mu), 3L)
  expect_identical(length(model$mu2), 3L)
  for (l in 1:3) {
    expect_identical(length(model$mu[[l]]), 2L)
    expect_identical(length(model$mu2[[l]]), 2L)
    expect_identical(dim(model$mu[[l]][[1]]), c(8L, 64L))
    expect_identical(dim(model$mu[[l]][[2]]), c(8L, 128L))
    expect_true(all(model$mu[[l]][[1]] == 0))
    expect_true(all(model$mu2[[l]][[2]] == 0))
  }
})

test_that("KL, lbf, lbf_variable initialized to NA at the right shape", {
  data <- make_data(n = 20, J = 8)
  params <- make_params(L = 3, J = data$p)
  var_y <- make_var_y(data$M, data$T_basis)

  model <- mfsusieR:::initialize_susie_model.mf_individual(data, params, var_y)

  expect_identical(length(model$KL), 3L)
  expect_true(all(is.na(model$KL)))
  expect_identical(length(model$lbf), 3L)
  expect_true(all(is.na(model$lbf)))
  expect_identical(dim(model$lbf_variable), c(3L, 8L))
  expect_true(all(is.na(model$lbf_variable)))
})

test_that("V (per-effect prior scaling) initialized to ones", {
  data <- make_data()
  params <- make_params(L = 5, J = data$p)
  var_y <- make_var_y(data$M, data$T_basis)

  model <- mfsusieR:::initialize_susie_model.mf_individual(data, params, var_y)

  expect_equal(model$V, rep(1, 5), tolerance = 0)
})

test_that("sigma2 falls back to var_y when params$residual_variance is NULL", {
  data <- make_data()
  params <- make_params(L = 3, J = data$p)
  var_y <- make_var_y(data$M, data$T_basis)

  model <- mfsusieR:::initialize_susie_model.mf_individual(data, params, var_y)

  expect_identical(model$sigma2, var_y)
})

test_that("sigma2 takes params$residual_variance when supplied", {
  data <- make_data()
  T_basis <- data$T_basis
  custom <- lapply(seq_along(T_basis),
                   function(m) rep(0.5, T_basis[m]))
  params <- make_params(L = 3, J = data$p)
  params$residual_variance <- custom
  var_y <- make_var_y(data$M, T_basis)

  model <- mfsusieR:::initialize_susie_model.mf_individual(data, params, var_y)

  expect_identical(model$sigma2, custom)
})

test_that("variable-selection prior pi defaults to uniform 1/J when prior_weights is NULL", {
  data <- make_data(J = 6)
  params <- make_params(L = 3, J = data$p)
  params$prior_weights <- NULL
  var_y <- make_var_y(data$M, data$T_basis)

  model <- mfsusieR:::initialize_susie_model.mf_individual(data, params, var_y)

  expect_equal(model$pi, rep(1 / 6, 6), tolerance = 0)
})

test_that("variable-selection prior pi takes params$prior_weights when supplied", {
  data <- make_data(J = 5)
  custom_prior <- c(0.1, 0.2, 0.3, 0.2, 0.2)
  params <- make_params(L = 3, J = data$p)
  params$prior_weights <- custom_prior
  var_y <- make_var_y(data$M, data$T_basis)

  model <- mfsusieR:::initialize_susie_model.mf_individual(data, params, var_y)

  expect_equal(model$pi, custom_prior, tolerance = 0)
})

test_that("fitted is a per-modality list of zero matrices at the right shape", {
  data <- make_data(n = 25, J = 6, T_per_outcome = c(64L, 128L))
  params <- make_params(L = 3, J = data$p)
  var_y <- make_var_y(data$M, data$T_basis)

  model <- mfsusieR:::initialize_susie_model.mf_individual(data, params, var_y)

  expect_identical(length(model$fitted), 2L)
  expect_identical(dim(model$fitted[[1]]), c(25L, 64L))
  expect_identical(dim(model$fitted[[2]]), c(25L, 128L))
  expect_true(all(model$fitted[[1]] == 0))
  expect_true(all(model$fitted[[2]] == 0))
})

test_that("intercept is a length-M zero vector", {
  data <- make_data(T_per_outcome = c(64L, 128L, 256L))
  params <- make_params(L = 3, J = data$p)
  var_y <- make_var_y(data$M, data$T_basis)

  model <- mfsusieR:::initialize_susie_model.mf_individual(data, params, var_y)

  expect_identical(model$intercept, rep(0, 3))
})

test_that("L, J, M are stored on the model for downstream dispatch", {
  data <- make_data(n = 17, J = 11, T_per_outcome = c(64L, 128L))
  params <- make_params(L = 4, J = data$p)
  var_y <- make_var_y(data$M, data$T_basis)

  model <- mfsusieR:::initialize_susie_model.mf_individual(data, params, var_y)

  expect_identical(model$L, 4L)
  expect_identical(model$p, 11L)
  expect_identical(model$M, 2L)
})

# ---- Degenerate cases (M = 1 functional, M = 1 univariate) -------------

test_that("M = 1 functional case produces a length-1 nested list with one modality", {
  data <- make_data(T_per_outcome = 64L)
  params <- make_params(L = 3, J = data$p)
  var_y <- make_var_y(data$M, data$T_basis)

  model <- mfsusieR:::initialize_susie_model.mf_individual(data, params, var_y)

  expect_identical(model$M, 1L)
  expect_identical(length(model$mu[[1]]), 1L)
  expect_identical(dim(model$mu[[1]][[1]]), c(8L, 64L))
})

test_that("T_m = 1 univariate case produces a J x 1 mu matrix", {
  set.seed(mfsusier_test_seed())
  X <- matrix(rnorm(20 * 8), nrow = 20)
  Y <- list(matrix(rnorm(20), ncol = 1))
  data <- mfsusieR:::create_mf_individual(X = X, Y = Y, pos = list(1),
                                          verbose = FALSE)
  params <- make_params(L = 3, J = 8)
  var_y <- list(rep(1, 1))

  model <- mfsusieR:::initialize_susie_model.mf_individual(data, params, var_y)

  expect_identical(model$M, 1L)
  expect_identical(data$T_basis, 1L)
  expect_identical(dim(model$mu[[1]][[1]]), c(8L, 1L))
  expect_identical(dim(model$fitted[[1]]), c(20L, 1L))
})

test_that("ragged multi-modality (M = 3) preserves per-modality shape", {
  data <- make_data(n = 20, J = 8, T_per_outcome = c(64L, 128L, 256L))
  params <- make_params(L = 2, J = data$p)
  var_y <- make_var_y(data$M, data$T_basis)

  model <- mfsusieR:::initialize_susie_model.mf_individual(data, params, var_y)

  expect_identical(model$M, 3L)
  for (l in 1:2) {
    expect_identical(dim(model$mu[[l]][[1]]), c(8L, 64L))
    expect_identical(dim(model$mu[[l]][[2]]), c(8L, 128L))
    expect_identical(dim(model$mu[[l]][[3]]), c(8L, 256L))
  }
})

# ---- Per-effect accessors (model_methods.R) ---------------------------

test_that("get_alpha_l returns the l-th row of alpha", {
  data <- make_data()
  params <- make_params(L = 3, J = data$p)
  var_y <- make_var_y(data$M, data$T_basis)
  model <- mfsusieR:::initialize_susie_model.mf_individual(data, params, var_y)

  # Mutate alpha to verify the accessor returns the actual row.
  model$alpha[2, ] <- seq_len(data$p)
  expect_equal(mfsusieR:::get_alpha_l.mfsusie(model, 2),
               seq_len(data$p), tolerance = 0)
})

test_that("get_posterior_moments_l returns mu and mu2 lists for effect l", {
  data <- make_data()
  params <- make_params(L = 3, J = data$p)
  var_y <- make_var_y(data$M, data$T_basis)
  model <- mfsusieR:::initialize_susie_model.mf_individual(data, params, var_y)

  out <- mfsusieR:::get_posterior_moments_l.mfsusie(model, 1)

  expect_named(out, c("mu", "mu2"))
  expect_identical(length(out$mu), 2L)
  expect_identical(length(out$mu2), 2L)
  expect_identical(dim(out$mu[[1]]), c(data$p, data$T_basis[1]))
})

test_that("get_posterior_mean_l multiplies mu element-wise by alpha[l, ]", {
  data <- make_data()
  params <- make_params(L = 3, J = data$p)
  var_y <- make_var_y(data$M, data$T_basis)
  model <- mfsusieR:::initialize_susie_model.mf_individual(data, params, var_y)

  # Inject deterministic mu and alpha so we can verify the product.
  set.seed(mfsusier_test_seed())
  for (m in seq_len(model$M)) {
    model$mu[[1]][[m]] <- matrix(seq_len(data$p * data$T_basis[m]),
                                 nrow = data$p)
  }
  model$alpha[1, ] <- seq_len(data$p) / sum(seq_len(data$p))

  out <- mfsusieR:::get_posterior_mean_l.mfsusie(model, 1)

  for (m in seq_len(model$M)) {
    expected <- sweep(model$mu[[1]][[m]], 1, model$alpha[1, ], "*")
    expect_equal(out[[m]], expected, tolerance = 0)
  }
})

test_that("get_posterior_mean_sum sums get_posterior_mean_l across L", {
  data <- make_data()
  params <- make_params(L = 3, J = data$p)
  var_y <- make_var_y(data$M, data$T_basis)
  model <- mfsusieR:::initialize_susie_model.mf_individual(data, params, var_y)

  # Inject deterministic mu / alpha across all L.
  for (l in seq_len(model$L)) {
    for (m in seq_len(model$M)) {
      model$mu[[l]][[m]] <-
        matrix(l + seq_len(data$p * data$T_basis[m]),
               nrow = data$p)
    }
    model$alpha[l, ] <- rep(1 / data$p, data$p)
  }

  out <- mfsusieR:::get_posterior_mean_sum.mfsusie(model)

  expected <- vector("list", model$M)
  for (m in seq_len(model$M)) {
    expected[[m]] <- matrix(0, nrow = data$p, ncol = data$T_basis[m])
    for (l in seq_len(model$L)) {
      expected[[m]] <- expected[[m]] +
        sweep(model$mu[[l]][[m]], 1, model$alpha[l, ], "*")
    }
  }
  for (m in seq_len(model$M)) {
    expect_equal(out[[m]], expected[[m]], tolerance = 0)
  }
})

test_that("get_prior_variance_l returns the l-th entry of V", {
  data <- make_data()
  params <- make_params(L = 3, J = data$p)
  var_y <- make_var_y(data$M, data$T_basis)
  model <- mfsusieR:::initialize_susie_model.mf_individual(data, params, var_y)

  expect_identical(mfsusieR:::get_prior_variance_l.mfsusie(model, 1), 1)
  model$V[2] <- 7
  expect_identical(mfsusieR:::get_prior_variance_l.mfsusie(model, 2), 7)
})

test_that("set_prior_variance_l updates V[l] and returns the model", {
  data <- make_data()
  params <- make_params(L = 3, J = data$p)
  var_y <- make_var_y(data$M, data$T_basis)
  model <- mfsusieR:::initialize_susie_model.mf_individual(data, params, var_y)

  m1 <- mfsusieR:::set_prior_variance_l.mfsusie(model, 2, 5)
  expect_identical(m1$V[2], 5)
  expect_identical(m1$V[c(1, 3)], rep(1, 2))     # other effects unchanged
  expect_identical(class(m1), c("mfsusie", "susie"))
})
