# 6-grid benchmark for the per_scale_normal vs per_scale + mixture_null_weight
# combinations Gao listed on 2026-05-03 Slack. Runs FDR / power / runtime /
# memory / convergence per cell over a small replicate count, writes a tidy
# data.frame to inst/bench/profiling/results/.
#
# Estimated wall-clock: ~10-25 minutes total on a single core at the default
# sim size (n = 84, p = 500, T = 64, M = 2, n_rep = 5, L = 10). Scaling p, T,
# n_rep, or L past these multiplies the wall-clock approximately linearly.
# Per CLAUDE.md hard rule 3, ask the user before pushing past the 30-minute
# threshold.
#
# Usage from the package root:
#   Rscript inst/bench/profiling/benchmark_per_scale_normal_6grid.R
#
# Output:
#   inst/bench/profiling/results/per_scale_normal_6grid_<timestamp>.rds
#   - tidy data.frame with one row per (grid cell, replicate)
#   inst/bench/profiling/results/per_scale_normal_6grid_<timestamp>_summary.csv
#   - per-cell aggregates ready for the results memo

suppressPackageStartupMessages({
  library(devtools)
  load_all(quiet = TRUE)
})

# ---- Simulation parameters (binding for this benchmark) -------------

set.seed(2026L * 503L)
n     <- 84L           # mirrors Anjing's ATAC cohort size
p     <- 500L          # large enough to exercise FDR; not so large the run
                       # blows past 30 min.
M     <- 2L            # two outcomes
Tlen  <- 64L           # single dyadic length, both outcomes
n_rep <- 5L            # replicates per grid cell
L     <- 10L           # IBSS effect cap
true_indices <- c(50L, 220L, 380L)  # three causal SNPs
pip_thresh   <- 0.05               # PIP cutoff for FDR / power summaries

# ---- Grid (6 cells per Gao 2026-05-03) ------------------------------

bench_grid <- rbind(
  expand.grid(wavelet_qnorm        = c(FALSE, TRUE),
              prior_variance_scope = "per_scale",
              mixture_null_weight  = c(0.05, 0),
              stringsAsFactors     = FALSE),
  data.frame  (wavelet_qnorm        = c(FALSE, TRUE),
               prior_variance_scope = "per_scale_normal",
               mixture_null_weight  = NA_real_,
               stringsAsFactors     = FALSE)
)
stopifnot(nrow(bench_grid) == 6L)

# ---- Per-replicate data builder -------------------------------------

build_dataset <- function(rep_seed) {
  set.seed(rep_seed)
  X <- matrix(rnorm(n * p), n, p)
  effect_curve <- numeric(Tlen)
  effect_curve[20:40] <- 1
  beta_per_snp <- matrix(0, nrow = p, ncol = Tlen)
  for (j in true_indices) beta_per_snp[j, ] <- effect_curve
  Y <- vector("list", M)
  for (m in seq_len(M)) {
    sig <- X %*% beta_per_snp
    Y[[m]] <- sig + matrix(rnorm(n * Tlen, sd = 1.0), n, Tlen)
  }
  list(X = X, Y = Y)
}

# ---- Per-cell metrics ------------------------------------------------

eval_fit <- function(fit, true_idx, pip_thresh) {
  pip <- fit$pip
  selected <- which(pip > pip_thresh)
  n_disc   <- length(selected)
  n_true   <- length(true_idx)
  n_tp     <- sum(selected %in% true_idx)
  n_fp     <- n_disc - n_tp
  fdr      <- if (n_disc > 0L) n_fp / n_disc else 0
  power    <- n_tp / n_true
  cs_count <- if (!is.null(fit$sets$cs)) length(fit$sets$cs) else 0L
  cs_purity <- if (!is.null(fit$sets$purity))
    mean(fit$sets$purity[, "min.abs.corr"]) else NA_real_
  list(n_disc = n_disc, n_tp = n_tp, n_fp = n_fp,
       fdr = fdr, power = power, cs_count = cs_count,
       cs_purity = cs_purity, niter = fit$niter,
       converged = isTRUE(fit$converged))
}

# ---- Driver ---------------------------------------------------------

results <- list()
row_idx <- 0L

t0 <- Sys.time()
for (rep_i in seq_len(n_rep)) {
  d <- build_dataset(rep_seed = 100L + rep_i)
  for (g in seq_len(nrow(bench_grid))) {
    cell <- bench_grid[g, ]
    fit_args <- list(
      X                     = d$X,
      Y                     = d$Y,
      L                     = L,
      verbose               = FALSE,
      prior_variance_scope  = cell$prior_variance_scope,
      wavelet_qnorm         = cell$wavelet_qnorm,
      save_mu_method        = "alpha_collapsed"  # benchmark uses thinned fits
    )
    if (!is.na(cell$mixture_null_weight)) {
      fit_args$mixture_null_weight <- cell$mixture_null_weight
    }
    t_start  <- Sys.time()
    mem_before <- gc(reset = TRUE, verbose = FALSE)[2L, 6L]  # Mb used post-gc
    fit <- tryCatch(
      do.call(mfsusie, fit_args),
      error = function(e) {
        warning(sprintf("rep=%d cell=%d errored: %s", rep_i, g,
                        conditionMessage(e)))
        NULL
      }
    )
    t_elapsed <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
    mem_after <- gc(verbose = FALSE)[2L, 6L]
    if (is.null(fit)) next
    metrics <- eval_fit(fit, true_indices, pip_thresh)
    fit_size_mb <- as.numeric(object.size(fit)) / (1024 * 1024)
    row_idx <- row_idx + 1L
    results[[row_idx]] <- data.frame(
      rep                   = rep_i,
      cell                  = g,
      wavelet_qnorm         = cell$wavelet_qnorm,
      prior_variance_scope  = cell$prior_variance_scope,
      mixture_null_weight   = cell$mixture_null_weight,
      n_disc    = metrics$n_disc,
      n_tp      = metrics$n_tp,
      n_fp      = metrics$n_fp,
      fdr       = metrics$fdr,
      power     = metrics$power,
      cs_count  = metrics$cs_count,
      cs_purity = metrics$cs_purity,
      niter     = metrics$niter,
      converged = metrics$converged,
      runtime_s = t_elapsed,
      fit_size_mb = fit_size_mb,
      mem_used_mb = mem_after - mem_before,
      stringsAsFactors = FALSE
    )
    cat(sprintf("rep=%d cell=%d (%s, qnorm=%s, mnw=%s)  fdr=%.3f power=%.3f t=%.1fs\n",
                rep_i, g, cell$prior_variance_scope,
                cell$wavelet_qnorm,
                if (is.na(cell$mixture_null_weight)) "NA"
                else format(cell$mixture_null_weight),
                metrics$fdr, metrics$power, t_elapsed))
  }
}
t_total <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
cat(sprintf("\nTotal wall-clock: %.1f s\n", t_total))

results_df <- do.call(rbind, results)

# ---- Persist outputs ------------------------------------------------

out_dir <- "inst/bench/profiling/results"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
ts  <- format(Sys.time(), "%Y%m%d_%H%M")
out_rds <- file.path(out_dir, sprintf("per_scale_normal_6grid_%s.rds", ts))
out_csv <- file.path(out_dir,
                     sprintf("per_scale_normal_6grid_%s_summary.csv", ts))
saveRDS(results_df, out_rds)
cat(sprintf("Saved per-replicate results to %s\n", out_rds))

summary_df <- aggregate(
  cbind(fdr, power, n_disc, cs_count, cs_purity, niter, runtime_s, fit_size_mb)
    ~ wavelet_qnorm + prior_variance_scope + mixture_null_weight,
  data = results_df, FUN = mean, na.action = na.pass)
write.csv(summary_df, out_csv, row.names = FALSE)
cat(sprintf("Saved per-cell summary to %s\n", out_csv))
print(summary_df)
