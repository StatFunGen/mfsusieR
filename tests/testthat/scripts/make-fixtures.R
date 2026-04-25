# Fixture generator for mfsusieR test suite.
#
# Run from the package root:
#   Rscript tests/testthat/scripts/make-fixtures.R
#
# Writes one .rds per scenario into tests/testthat/fixtures/. Fixtures
# are committed to git when small (< 1 MB each), per the Phase 3 / 4
# convention recorded in CLAUDE.md. The script is the single source of
# truth for fixture contents; the .rds files are generated artefacts.
#
# Each fixture contains a self-contained list with at minimum:
#   - X        : N x J SNP matrix (centered, not scaled),
#   - Y        : list with $Y_f (list of N x T_m functional matrices)
#                and optional $Y_u (N x q_u univariate matrix), matching
#                the input shape expected by `mvf.susie.alpha::multfsusie`,
#   - pos      : list of length M, sampling positions per modality,
#   - L_true   : number of true causal effects,
#   - true_idx : integer vector of true causal SNP indices,
#   - seed     : the recorded RNG seed used to construct this fixture.
#
# This file does NOT depend on `mvf.susie.alpha`. Apple-to-apple fidelity
# tests in later PR groups will load a fixture here and call
# `mvf.susie.alpha::multfsusie` themselves; if that package is absent, the
# test is skipped via `skip_if_no_mvf_alpha()` in tests/testthat/setup.R.

suppressPackageStartupMessages({
  # No package deps; we use only base R + stats.
})

fixture_dir <- file.path("tests", "testthat", "fixtures")
dir.create(fixture_dir, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------------------
# Scenario: minimal (M = 2, T = c(64, 128), N = 200, J = 50, L = 3)
# ---------------------------------------------------------------------------
# Two functional modalities of different lengths (ragged T_m), small N and
# J so the fixture stays under 1 MB. Three true causal SNPs with sparse
# wavelet-domain effects added to Gaussian noise.

make_scenario_minimal <- function() {
  seed <- 20260425L
  set.seed(seed)

  N <- 200L
  J <- 50L
  T_per_modality <- c(64L, 128L)
  M <- length(T_per_modality)
  L_true <- 3L

  # Genotype matrix: 0/1/2 dosages with mild minor-allele frequency.
  maf <- runif(J, 0.1, 0.4)
  X <- vapply(maf, function(p) {
    rbinom(N, size = 2L, prob = p)
  }, FUN.VALUE = numeric(N))
  storage.mode(X) <- "double"
  X <- scale(X, center = TRUE, scale = FALSE)
  attr(X, "scaled:center") <- NULL

  # True causal SNPs.
  true_idx <- sort(sample.int(J, L_true))

  # Smooth functional effect for each (effect l, modality m). The signal
  # is a sum of two Gaussian bumps with modality-specific centres, scaled
  # so that the per-modality signal-to-noise ratio is moderate.
  make_func_effect <- function(T_m, l, m) {
    grid <- seq_len(T_m) / T_m
    centre1 <- ((l - 1L) %% 3L) / 3 + 0.15
    centre2 <- centre1 + 0.4 + 0.05 * m
    width <- 0.05
    amp <- 0.6 * (1 + 0.2 * (m - 1))
    amp * (exp(-((grid - centre1) ^ 2) / (2 * width ^ 2)) -
             0.7 * exp(-((grid - centre2) ^ 2) / (2 * width ^ 2)))
  }

  Y_f <- vector("list", M)
  for (m in seq_len(M)) {
    T_m <- T_per_modality[m]
    Y_f[[m]] <- matrix(rnorm(N * T_m, sd = 1.0), nrow = N, ncol = T_m)
    for (l in seq_len(L_true)) {
      eff <- make_func_effect(T_m, l, m)
      x_l <- X[, true_idx[l]]
      Y_f[[m]] <- Y_f[[m]] + tcrossprod(x_l, eff)
    }
  }

  pos <- lapply(T_per_modality, function(T_m) seq_len(T_m))

  list(
    name           = "scenario_minimal",
    seed           = seed,
    N              = N,
    J              = J,
    M              = M,
    T_per_modality = T_per_modality,
    L_true         = L_true,
    true_idx       = true_idx,
    X              = X,
    Y              = list(Y_f = Y_f, Y_u = NULL),
    pos            = pos
  )
}

scenario_minimal <- make_scenario_minimal()
out_path <- file.path(fixture_dir, paste0(scenario_minimal$name, ".rds"))
saveRDS(scenario_minimal, file = out_path, compress = "xz")
cat(sprintf("wrote %s (%.1f KB)\n",
            out_path,
            file.info(out_path)$size / 1024))
