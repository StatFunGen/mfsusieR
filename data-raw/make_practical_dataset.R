# Build `data/practical_fits.rda` for the
# `practical_data_applications` vignette.
#
# Rationale: a real-world ATAC-seq dataset by Anjing Liu and
# William Denault motivates the three round-3 features
# (`low_count_filter`, `quantile_norm`, `small_sample_correction`).
# The raw genotype, cell-count, and read-count files are not
# redistributable, so the vignette uses a simulated dataset
# constructed to share the relevant characteristics: small `n`,
# multiple outcomes (cell types), heavy-tailed wavelet
# coefficients (heavy-tailed log1p of sparse counts), and a
# subset of wavelet columns with near-zero median (sparse
# coverage). The framing in the vignette maps this back to the
# motivating ATAC-seq use case.
#
# Run from the package root with:
#   Rscript data-raw/make_practical_dataset.R

suppressPackageStartupMessages({
  library(devtools)
  load_all(quiet = TRUE)
})

set.seed(2026L)

# --- Simulation parameters ----------------------------------------
# n = 40 lands in the regime where the small-sample correction
# matters (the Anjing cohort is n = 84; here we shrink further
# to make the Johnson-t correction visible). M = 2 picks two
# representative cell types. T_m = 32 keeps the fits fast.
# p = 25 with one causal SNP at index 7 and 24 nulls.
n     <- 40L
M     <- 2L
T_m   <- 32L
p     <- 25L
sig   <- 7L          # signal SNP
cell_types <- c("CellTypeA", "CellTypeB")

# Genotypes from a binomial(2, MAF) with MAF drawn per SNP.
maf  <- runif(p, 0.10, 0.40)
X    <- vapply(maf, function(f) rbinom(n, 2L, f), numeric(n))
storage.mode(X) <- "double"

# Effect surface for the causal SNP: one localized bump at
# positions 16-24 (mid-curve).
beta_curve <- numeric(T_m)
beta_curve[16:24] <- 1

# Build per-outcome Y in the count-like scale, then log1p.
# Two characteristics carry over from the real data:
#   (i)  wavelet coefficients are heavy-tailed (some cell types
#        have rare-event spikes), which motivates `quantile_norm`.
#   (ii) some bins have median count ~ 0 across individuals,
#        which motivates `low_count_filter`.
# --- Build per-outcome Y -------------------------------------------
# Two outcomes carry the same signal SNP at index `sig`. Each
# outcome carries a different non-Gaussianity that one of the
# three options addresses:
#   - Outcome 1: a block of bins with median(|.|) below the
#     `low_count_filter` threshold but with non-zero noise that
#     pollutes the SER step. `low_count_filter` masks them.
#   - Outcome 2: heavy-tailed contamination but no near-zero
#     bins. `quantile_norm` Gaussianizes the wavelet
#     coefficients without destroying signal.
# Both outcomes are at small `n = 40` so the
# `small_sample_correction` (Johnson-t) is the right fit.
Y_f <- vector("list", M)

# Outcome 1: clean Gaussian noise on the signal bins (16-24)
# but a block of low-coverage bins at 1-8 with small but
# non-zero noise that the Wakefield kernel mistakes for
# information.
{
  signal_1 <- X[, sig, drop = FALSE] %*%
              matrix(0.55 * beta_curve, 1, T_m)
  noise_1  <- matrix(rnorm(n * T_m, mean = 0, sd = 0.8),
                     nrow = n, ncol = T_m)
  Y_f[[1]] <- signal_1 + noise_1
  # Low-coverage block: bins 1-8 are at very low SNR. Their
  # absolute-value median is small enough that
  # `low_count_filter = 0.02` masks them but `low_count_filter
  # = 0` keeps them in the fit, where they contribute noise.
  Y_f[[1]][, 1:8] <- rnorm(n * 8L, mean = 0, sd = 0.03)
}

# Outcome 2: heavy-tailed contamination only (no near-zero
# bins). 18% of cells receive a large-variance perturbation.
# The wavelet coefficients inherit heavy tails;
# `quantile_norm` Gaussianizes them and the SuSiE Bayes-factor
# assumption is restored.
{
  signal_2 <- X[, sig, drop = FALSE] %*%
              matrix(0.45 * beta_curve, 1, T_m)
  noise_2  <- matrix(rnorm(n * T_m, mean = 0, sd = 0.8),
                     nrow = n, ncol = T_m)
  outl_2   <- matrix(rbinom(n * T_m, 1L, 0.18),
                     nrow = n, ncol = T_m) *
              rnorm(n * T_m, mean = 0, sd = 4.0)
  Y_f[[2]] <- signal_2 + noise_2 + outl_2
}

pos <- lapply(seq_len(M), function(k) seq_len(T_m))

# --- Run four fits ------------------------------------------------
verbose <- FALSE
L_eff   <- 8L

cat(sprintf("Fitting n = %d, p = %d, M = %d, T_m = %d, L = %d\n",
            n, p, M, T_m, L_eff))

t_default <- system.time(
  fit_default <- mfsusie(X, Y_f, pos = pos, L = L_eff,
                         verbose = verbose)
)
cat(sprintf("default (Wakefield):              %.1f s, niter = %d\n",
            t_default[3L], fit_default$niter))

t_johnson <- system.time(
  fit_johnson <- mfsusie(X, Y_f, pos = pos, L = L_eff,
                         small_sample_correction = TRUE,
                         verbose = verbose)
)
cat(sprintf("small_sample_correction = TRUE:   %.1f s, niter = %d\n",
            t_johnson[3L], fit_johnson$niter))

t_lowcount <- system.time(
  fit_lowcount <- mfsusie(X, Y_f, pos = pos, L = L_eff,
                          low_count_filter = 0.02,
                          verbose = verbose)
)
cat(sprintf("low_count_filter = 0.02:          %.1f s, niter = %d\n",
            t_lowcount[3L], fit_lowcount$niter))

t_qn <- system.time(
  fit_qn <- mfsusie(X, Y_f, pos = pos, L = L_eff,
                    quantile_norm = TRUE,
                    verbose = verbose)
)
cat(sprintf("quantile_norm = TRUE:             %.1f s, niter = %d\n",
            t_qn[3L], fit_qn$niter))

t_combined <- system.time(
  fit_combined <- mfsusie(X, Y_f, pos = pos, L = L_eff,
                          low_count_filter        = 0.02,
                          quantile_norm           = TRUE,
                          small_sample_correction = TRUE,
                          verbose = verbose)
)
cat(sprintf("all three combined:               %.1f s, niter = %d\n",
            t_combined[3L], fit_combined$niter))

# --- Strip raw inputs and large internal state --------------------
strip_fit <- function(fit) {
  fit$residuals <- NULL
  fit$XtY       <- NULL
  fit$tracking  <- NULL
  fit$X         <- NULL
  fit$Y         <- NULL
  fit
}

practical_fits <- list(
  n          = n,
  p          = p,
  M          = M,
  T_m        = T_m,
  signal_snp = sig,
  cell_types = cell_types,
  default    = strip_fit(fit_default),
  johnson    = strip_fit(fit_johnson),
  lowcount   = strip_fit(fit_lowcount),
  qn         = strip_fit(fit_qn),
  combined   = strip_fit(fit_combined)
)

# --- Save ---------------------------------------------------------
usethis::use_data(practical_fits, overwrite = TRUE, compress = "xz")
cat("Saved practical_fits to data/practical_fits.rda\n")
