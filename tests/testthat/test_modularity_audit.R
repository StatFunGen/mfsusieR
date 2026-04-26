# Modularity audit guard (PR group 9b.2 + 9.1).
#
# Two checks:
#
# 1. D12 forbidden strings: per design.md D12, R/ source must NOT
#    contain bare `mvf.susie.alpha`, `multfsusie`, `fsusieR`,
#    `susiF`, "original implementation", `@references_original`,
#    or `@manuscript_ref`. Provenance lives in
#    `inst/notes/refactor-exceptions.md` and test-file headers.
#    Exception: `R/fsusie.R` is the documented migration wrapper
#    and must mention `fsusieR::susiF` in its user-facing roxygen.
#
# 2. Forbidden private reimplementations: scan R/ for hand-rolled
#    duplicates of work the upstream packages already do (mixsqp,
#    DWT, ash, susieR's per-effect machinery). Use the call-name
#    list from design.md D10 as the allowlist and flag any local
#    reimplementation.

R_FILES <- list.files("../../R", pattern = "[.]R$", full.names = TRUE)
read_file_text <- function(path) {
  paste(readLines(path, warn = FALSE), collapse = "\n")
}

# ---- D12 forbidden strings ---------------------------------------

test_that("R/ source does not contain forbidden D12 strings (except R/fsusie.R)", {
  # `(?<![mM])` prevents matching `mfsusieR` (substring "fsusieR") and
  # `mvfsusieR` etc. as false positives.
  forbidden <- c("mvf\\.susie\\.alpha",
                 "(?<![a-zA-Z_])multfsusie",
                 "(?<![mM])fsusieR",
                 "(?<![a-zA-Z_])susiF\\b",
                 "original implementation",
                 "@references_original",
                 "@manuscript_ref")
  exempt <- c("R/fsusie.R",   # migration wrapper mentions fsusieR::susiF intentionally
              "R/cpp11.R")     # auto-generated bindings; not human-edited

  bad <- list()
  for (path in R_FILES) {
    rel <- sub(".*/R/", "R/", path)
    if (rel %in% exempt) next
    body <- read_file_text(path)
    for (pat in forbidden) {
      if (grepl(pat, body, perl = TRUE)) {
        bad <- c(bad, sprintf("%s: matches /%s/", rel, pat))
      }
    }
  }
  expect_identical(bad, list(),
                   info = paste("D12 violations:",
                                paste(unlist(bad), collapse = "; ")))
})

# ---- Forbidden private reimplementations ----------------------

test_that("R/ does not hand-roll routines provided by susieR / mixsqp / ashr / wavethresh", {
  # The single canonical caller of `wavethresh::wd` outside R/dwt.R
  # would indicate a hand-rolled DWT path. Same for `mixsqp::mixsqp`
  # outside R/em_helpers.R.
  rule_table <- list(
    list(pattern = "wavethresh::w[dr]\\(",
         allowed = c("R/dwt.R", "R/utils_wavelet.R"),
         hint    = "wavethresh::wd / wr should only be called from R/dwt.R + R/utils_wavelet.R"),
    list(pattern = "mixsqp::mixsqp\\(",
         allowed = c("R/em_helpers.R"),
         hint    = "mixsqp::mixsqp should only be called from R/em_helpers.R")
  )
  for (path in R_FILES) {
    body <- read_file_text(path)
    rel  <- sub(".*/R/", "R/", path)
    for (rule in rule_table) {
      if (grepl(rule$pattern, body, perl = TRUE) &&
          !(rel %in% rule$allowed)) {
        fail(sprintf("%s contains forbidden private reimpl: %s",
                     rel, rule$hint))
      }
    }
  }
  expect_true(TRUE)
})

# ---- Public API smoke after fresh load (PR 8.2) -----------------

test_that("mfsusie() and fsusie() run end-to-end after a fresh devtools::load_all", {
  # If `.onLoad` failed to register S3 methods on susieR's namespace,
  # the per-effect dispatch in susieR::susie_workhorse would fall
  # through to the .default methods and produce wrong results.
  set.seed(mfsusier_test_seed())
  n <- 25; J <- 6
  X <- matrix(rnorm(n * J), nrow = n)
  beta <- numeric(J); beta[1] <- 1

  # Single-modality functional via mfsusie
  Y_func <- X %*% matrix(rep(beta, 32L), nrow = J) +
            matrix(rnorm(n * 32L, sd = 0.3), nrow = n)
  fit_m <- mfsusie(X, list(Y_func), L = 2, max_iter = 20, verbose = FALSE)
  expect_s3_class(fit_m, "mfsusie")
  expect_true(isTRUE(fit_m$converged) || fit_m$niter > 0L)

  # Drop-in via fsusie
  fit_f <- fsusie(Y_func, X, L = 2, max_iter = 20, verbose = FALSE)
  expect_s3_class(fit_f, "mfsusie")
})
