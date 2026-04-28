# Coverage-driven tests for mfsusie_plot.R.
# Goal: exercise every public branch and most internal branches
# of the plotting code so coverage is well above 95%.

# ---- shared fixture builders -------------------------------------------------

make_M1 <- function(L = 5, n = 80, p = 25, T_m = 64) {
  set.seed(7L)
  X <- matrix(rnorm(n * p), n)
  beta <- numeric(p); beta[3] <- 1; beta[10] <- -0.6
  Y <- list(X %*% matrix(rep(beta, T_m), nrow = p) +
              matrix(rnorm(n * T_m, sd = 0.5), n))
  fit <- mfsusie(X, Y, L = L, verbose = FALSE)
  list(fit = fit, X = X, Y = Y, beta = beta)
}

make_M2 <- function(L = 5, n = 80, p = 25, T_per = c(32L, 32L)) {
  set.seed(11L)
  X <- matrix(rnorm(n * p), n)
  beta <- numeric(p); beta[5] <- 1; beta[12] <- -0.5
  Y <- lapply(T_per, function(T_m)
    X %*% matrix(rep(beta, T_m), nrow = p) +
      matrix(rnorm(n * T_m, sd = 0.5), n))
  fit <- mfsusie(X, Y, L = L,
                 prior_variance_scope = "per_outcome",
                 verbose = FALSE)
  list(fit = fit, X = X, Y = Y, beta = beta)
}

# Render to a PNG file under tempdir; verifies the call ran without
# errors AND produced a non-trivial output. Centralised here so we
# don't have to wrap every test in `png()` + `dev.off()`.
render <- function(expr) {
  out <- tempfile(fileext = ".png")
  png(out, width = 900, height = 700)
  res <- tryCatch(force(expr), finally = dev.off())
  expect_gt(file.info(out)$size, 1000)
  invisible(res)
}

# ---- mfsusie_plot dispatch surface ------------------------------------------

test_that("mfsusie_plot rejects non-mfsusie input", {
  expect_error(mfsusie_plot(list()),
               "must be an `mfsusie`")
  expect_error(mfsusie_plot_lfsr(list()),
               "must be an `mfsusie`")
  expect_error(mfsusie_plot_dimensions(list()),
               "must be an `mfsusie`")
})

test_that("mfsusie_plot rejects out-of-range m", {
  fx <- make_M1()
  expect_error(mfsusie_plot(fx$fit, m = 0L), "must be in 1")
  expect_error(mfsusie_plot(fx$fit, m = 5L), "must be in 1")
  expect_error(mfsusie_plot_dimensions(fx$fit, m = 0L),
               "must be in 1")
})

test_that("mfsusie_plot M = 1 default + smoothed paths", {
  fx <- make_M1()
  render(mfsusie_plot(fx$fit))                          # raw, no smoother
  fit_s <- mf_post_smooth(fx$fit, method = "TI")
  render(mfsusie_plot(fit_s))                           # auto smoother
  render(mfsusie_plot(fit_s, smooth_method = "TI"))     # explicit
})

test_that("mfsusie_plot facet modes for M = 1", {
  fx    <- make_M1()
  fit_s <- mf_post_smooth(fx$fit, method = "TI")
  render(mfsusie_plot(fit_s, facet_cs = "auto"))
  render(mfsusie_plot(fit_s, facet_cs = "overlay"))
  render(mfsusie_plot(fit_s, facet_cs = "stack"))
  # m focus
  render(mfsusie_plot(fit_s, m = 1L, facet_cs = "stack"))
  render(mfsusie_plot(fit_s, m = 1L, facet_cs = "overlay"))
})

test_that("mfsusie_plot effect_style errorbar requires bands", {
  fx <- make_M1()
  expect_error(
    mfsusie_plot(fx$fit, effect_style = "errorbar"),
    "errorbar")
  fit_s <- mf_post_smooth(fx$fit, method = "TI")
  render(mfsusie_plot(fit_s, effect_style = "errorbar"))
  render(mfsusie_plot(fit_s, effect_style = "errorbar",
                      facet_cs = "stack"))
})

test_that("mfsusie_plot M > 1 layout + stack", {
  fx    <- make_M2()
  fit_s <- mf_post_smooth(fx$fit, method = "TI")
  render(mfsusie_plot(fit_s))                       # default tile
  render(mfsusie_plot(fit_s, facet_cs = "stack"))   # M*K stack
  render(mfsusie_plot(fit_s, m = 2L))               # focus
})

test_that("mfsusie_plot truth overlay for M = 1 and M > 1", {
  fx <- make_M1(T_m = 32L)
  fit_s <- mf_post_smooth(fx$fit, method = "TI")
  T_m <- ncol(fx$Y[[1]])
  truth_vec <- runif(T_m, -1, 1)
  render(mfsusie_plot(fit_s, truth = truth_vec))
  K <- length(fit_s$sets$cs)
  if (K >= 1L)
    render(mfsusie_plot(fit_s,
                         truth = replicate(K, truth_vec, simplify = FALSE)))

  fx2 <- make_M2()
  fit2_s <- mf_post_smooth(fx2$fit, method = "TI")
  truth_M <- list(runif(32L), runif(32L))
  render(mfsusie_plot(fit2_s, truth = truth_M))
})

test_that("mfsusie_plot truth normalize errors on bad shape", {
  fx <- make_M2()
  fit_s <- mf_post_smooth(fx$fit, method = "TI")
  expect_error(mfsusie_plot(fit_s, truth = "not-a-list"),
               "length-M list")
  expect_error(mfsusie_plot(fit_s, truth = list(runif(32))),
               "length-M list")
})

test_that("mfsusie_plot show_grid_dots, show_affected_region, lfsr", {
  fx    <- make_M1()
  fit_s <- mf_post_smooth(fx$fit, method = "TI")
  render(mfsusie_plot(fit_s, show_grid_dots = TRUE))
  render(mfsusie_plot(fit_s, show_affected_region = FALSE))

  # lfsr overlay only triggered when has_lfsr is TRUE; TI smoother
  # populates lfsr_curves so lfsr should overlay by default.
  render(mfsusie_plot(fit_s, show_lfsr_curve = TRUE))
  render(mfsusie_plot(fit_s, show_lfsr_curve = FALSE))
})

test_that("mfsusie_plot legend toggling", {
  fx    <- make_M1()
  fit_s <- mf_post_smooth(fx$fit, method = "TI")
  render(mfsusie_plot(fit_s, add_legend = FALSE))
})

test_that("mfsusie_plot scalar-outcome path (T_m = 1)", {
  set.seed(13L)
  n <- 60; p <- 15
  X <- matrix(rnorm(n * p), n)
  Y_scalar <- list(matrix(rnorm(n), ncol = 1L))
  fit <- mfsusie(X, Y_scalar, L = 3, verbose = FALSE)
  render(mfsusie_plot(fit))
})

test_that("mfsusie_plot save = file.* writes the requested format", {
  fx    <- make_M1()
  fit_s <- mf_post_smooth(fx$fit, method = "TI")

  for (ext in c("pdf", "png", "jpg", "svg")) {
    out <- tempfile(fileext = paste0(".", ext))
    res <- mfsusie_plot(fit_s, save = out)
    expect_true(file.exists(out))
    expect_gt(file.info(out)$size, 100)
  }

  expect_error(mfsusie_plot(fit_s, save = "/tmp/foo.bogus"),
               "must end in")
})

# ---- mfsusie_plot_dimensions surface ----------------------------------------

test_that("mfsusie_plot_dimensions returns sensible sizes", {
  fx <- make_M1()
  fit_s <- mf_post_smooth(fx$fit, method = "TI")

  d_default <- mfsusie_plot_dimensions(fit_s)
  expect_named(d_default, c("width", "height", "n_cells"))
  expect_gt(d_default$height, 0)
  expect_gt(d_default$width, 0)

  d_stack <- mfsusie_plot_dimensions(fit_s, facet_cs = "stack")
  expect_gte(d_stack$height, d_default$height)

  fx2 <- make_M2()
  fit2_s <- mf_post_smooth(fx2$fit, method = "TI")
  d_M2_default <- mfsusie_plot_dimensions(fit2_s)
  d_M2_stack   <- mfsusie_plot_dimensions(fit2_s, facet_cs = "stack")
  expect_gt(d_M2_stack$height, d_M2_default$height)

  # m focus
  d_m <- mfsusie_plot_dimensions(fit_s, m = 1L)
  expect_gt(d_m$height, 0)

  d_m_stack <- mfsusie_plot_dimensions(fit_s, m = 1L,
                                        facet_cs = "stack")
  expect_gte(d_m_stack$height, d_m$height)
})

# ---- mfsusie_plot_lfsr ------------------------------------------------------

test_that("mfsusie_plot_lfsr requires HMM", {
  fx <- make_M1()
  expect_error(mfsusie_plot_lfsr(fx$fit),
               "HMM-smoothed lfsr")
  fit_TI <- mf_post_smooth(fx$fit, method = "TI")
  expect_error(mfsusie_plot_lfsr(fit_TI, smooth_method = "scalewise"),
               "HMM-smoothed lfsr")
})

test_that("mfsusie_plot_lfsr renders for M = 1 and M > 1", {
  fx <- make_M1()
  fit_h <- mf_post_smooth(fx$fit, method = "HMM")
  render(mfsusie_plot_lfsr(fit_h))
  render(mfsusie_plot_lfsr(fit_h, add_legend = FALSE))

  fx2 <- make_M2()
  fit2_h <- mf_post_smooth(fx2$fit, method = "HMM")
  render(mfsusie_plot_lfsr(fit2_h))
})

test_that("mfsusie_plot_lfsr truth coloring", {
  fx <- make_M1(T_m = 32L)
  fit_h <- mf_post_smooth(fx$fit, method = "HMM")
  T_m <- 32L
  K <- length(fit_h$sets$cs)
  if (K >= 1L) {
    truth_per_cs <- replicate(K, sample(c(TRUE, FALSE), T_m, TRUE),
                                simplify = FALSE)
    render(mfsusie_plot_lfsr(fit_h, truth = truth_per_cs))
    # Single bool vec (replicated)
    render(mfsusie_plot_lfsr(fit_h, truth = sample(c(TRUE, FALSE),
                                                   T_m, TRUE)))
  }

  # Wrong-length truth vector errors at draw time.
  if (K >= 1L) {
    expect_error(
      render(mfsusie_plot_lfsr(fit_h,
                                truth = list(rep(TRUE, T_m + 5L)))),
      "does not match T_m")
  }
})

test_that("mfsusie_plot_lfsr truth M > 1 list shape", {
  fx2 <- make_M2()
  fit2_h <- mf_post_smooth(fx2$fit, method = "HMM")
  M <- length(fit2_h$dwt_meta$T_basis)
  T_m <- ncol(fx2$Y[[1]])
  truth_M <- replicate(M,
                       sample(c(TRUE, FALSE), T_m, TRUE),
                       simplify = FALSE)
  render(mfsusie_plot_lfsr(fit2_h, truth = truth_M))

  expect_error(mfsusie_plot_lfsr(fit2_h, truth = "bad"),
               "length-M list")
})

test_that("mfsusie_plot_lfsr save = file.*", {
  fx <- make_M1()
  fit_h <- mf_post_smooth(fx$fit, method = "HMM")
  for (ext in c("pdf", "png")) {
    out <- tempfile(fileext = paste0(".", ext))
    mfsusie_plot_lfsr(fit_h, save = out)
    expect_true(file.exists(out))
    expect_gt(file.info(out)$size, 100)
  }
})

# ---- plot.mfsusie S3 method -------------------------------------------------

test_that("plot.mfsusie dispatches to mfsusie_plot", {
  fx <- make_M1()
  fit_s <- mf_post_smooth(fx$fit, method = "TI")
  render(plot(fit_s))
})

# ---- edge cases for coverage gaps -------------------------------------------

test_that("internal helpers handle the corner cases", {
  # .normalize_truth: per-CS list, length-1 list, errors
  norm <- mfsusieR:::.normalize_truth
  v <- runif(64)
  expect_length(norm(v, M = 1L, K = 3L), 1L)              # vec replicated
  expect_length(norm(list(v, v, v), M = 1L, K = 3L), 1L)  # length-K list
  expect_length(norm(list(v),       M = 1L, K = 3L), 1L)  # length-1 list
  expect_error(norm("nope", M = 1L, K = 1L),
               "must be NULL")
  expect_error(norm(list(v, v, v), M = 2L, K = 3L),       # length != M
               "length-M list")
  expect_error(norm(list(v, list(v, v)), M = 2L, K = 3L), # entry !K
               "length-K list")

  # .resolve_facet returns "stack"/"overlay" exactly
  fx <- make_M1()
  fit_s <- mf_post_smooth(fx$fit, method = "TI")
  facet_resolved <- mfsusieR:::.resolve_facet(
    "auto", fit_s, 1L, fit_s$smoothed$TI)
  expect_true(facet_resolved %in% c("stack", "overlay"))
  # Forced stack / overlay short-circuit the resolver.
  expect_identical(
    mfsusieR:::.resolve_facet("stack", fit_s, 1L,
                               fit_s$smoothed$TI),
    "stack")
  expect_identical(
    mfsusieR:::.resolve_facet("overlay", fit_s, 1L,
                               fit_s$smoothed$TI),
    "overlay")

  # .affected_mask handles NULL
  expect_identical(mfsusieR:::.affected_mask(NULL), logical(0))
  # affected_runs handles all-zero flags
  band_no_flags <- cbind(rep(-1, 5L), rep(1, 5L))
  expect_identical(mfsusieR:::affected_runs(band_no_flags), list())
  expect_identical(mfsusieR:::affected_runs(NULL), list())
})

test_that("mfsusie_plot save with svg path opens svg device", {
  fx <- make_M1()
  fit_s <- mf_post_smooth(fx$fit, method = "TI")
  out <- tempfile(fileext = ".svg")
  mfsusie_plot(fit_s, save = out)
  expect_true(file.exists(out))
  expect_gt(file.info(out)$size, 100)
})

test_that("mfsusie_plot_lfsr renders >2 outcome grid with leftover cell", {
  # M = 3 means ceil(sqrt(3)) = 2 cols, 2 rows -> 1 leftover plot.new()
  set.seed(17L)
  n <- 80; p <- 20
  X <- matrix(rnorm(n * p), n)
  beta <- numeric(p); beta[3] <- 1
  Y <- lapply(c(32L, 32L, 32L), function(T_m)
    X %*% matrix(rep(beta, T_m), nrow = p) +
      matrix(rnorm(n * T_m, sd = 0.5), n))
  fit <- mfsusie(X, Y, L = 3,
                 prior_variance_scope = "per_outcome",
                 verbose = FALSE)
  fit_h <- mf_post_smooth(fit, method = "HMM")
  render(mfsusie_plot_lfsr(fit_h))
})

test_that("mfsusie_plot_dimensions for tile path (M > 1 overlay)", {
  fx <- make_M2()
  fit_s <- mf_post_smooth(fx$fit, method = "TI")
  d <- mfsusie_plot_dimensions(fit_s, facet_cs = "overlay")
  expect_named(d, c("width", "height", "n_cells"))
  expect_gt(d$height, 0)
})

test_that("low-level helpers cover the small-K and edge branches", {
  # .open_save_device with NULL is a no-op.
  expect_null(mfsusieR:::.open_save_device(NULL, list(width = 1, height = 1)))

  # .resolve_facet K < 2 path: build an empty CS list.
  fx <- make_M1()
  fit <- fx$fit
  fit_one <- fit
  fit_one$sets$cs <- fit$sets$cs[1L]   # exactly 1 CS
  expect_identical(
    mfsusieR:::.resolve_facet("auto", fit_one, 1L, NULL),
    "overlay")

  fit_three <- fit
  fit_three$sets$cs <- replicate(3L, 1L, simplify = FALSE)
  fit_three$sets$cs_index <- c(1L, 2L, 3L)
  fit_three$alpha <- rbind(fit$alpha[1L, , drop = FALSE],
                           fit$alpha[1L, , drop = FALSE],
                           fit$alpha[1L, , drop = FALSE])
  expect_identical(
    mfsusieR:::.resolve_facet("auto", fit_three, 1L, NULL),
    "stack")

  # .overlay_truth_dots NULL early-return.
  expect_null(mfsusieR:::.overlay_truth_dots(1:5, NULL))

  # .normalize_truth: list-of-NULL entries via full-list path.
  v <- runif(64)
  out <- mfsusieR:::.normalize_truth(list(NULL, v), M = 2L, K = 1L)
  expect_length(out, 2L)
  expect_length(out[[1L]], 1L)
  expect_null(out[[1L]][[1L]])

  # mfsusie_plot_dimensions: K < 2 short-circuits the resolver.
  d <- mfsusie_plot_dimensions(fit_one, m = 1L, facet_cs = "stack")
  expect_named(d, c("width", "height", "n_cells"))
})

test_that("mfsusie_plot_lfsr error message is precise", {
  fx <- make_M1()
  fit_TI <- mf_post_smooth(fx$fit, method = "TI")
  expect_error(
    mfsusie_plot_lfsr(fit_TI, smooth_method = "smash"),
    "HMM-smoothed lfsr")
})

test_that("mfsusie_plot_lfsr bubble draws the no-CS subpanel", {
  # Build a fit-shaped object whose sets$cs is empty for the
  # second outcome by mutating a real fit.
  fx2 <- make_M2()
  fit_h <- mf_post_smooth(fx2$fit, method = "HMM")
  # Force the lfsr_curves[[2]] entry to be a list of length 0
  # so the panel exits via plot.new() + title().
  bad <- fit_h
  bad$sets$cs <- list()
  expect_error(mfsusie_plot_lfsr(bad), NA)  # should run silently
})
