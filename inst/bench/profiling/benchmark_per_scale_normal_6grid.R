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
#
# Three metrics, all reported per (cell, rep):
#
#   SNP-level loose (PIP > pip_thresh, default 0.05)
#     Every SNP with pip above threshold is one judgment regardless
#     of CS membership. Sensitive to spurious low-PIP leakage.
#       n_disc / n_tp / n_fp / fdr / power
#
#   SNP-level hybrid (CS purity > purity_thresh) OR (no CS AND PIP > pip_high_thresh)
#     A SNP counts as a discovery if it sits in any credible set whose
#     min.abs.corr >= purity_thresh (default 0.8), OR if it is outside
#     every CS but has PIP above pip_high_thresh (default 0.5). This
#     is the fine-mapping practice view: trust high-purity CSes,
#     and only accept stand-alone SNPs at high PIP.
#       hyb_n_disc / hyb_n_tp / hyb_n_fp / fdr_hyb / power_hyb
#
#   CS-level (vignette / per_scale_normal_vignette_sweep view)
#     For each fit$sets$cs, take the lead = SNP within the CS with
#     max pip. Classify TP if lead is a causal SNP OR
#     |cor(X[, lead], X[, causal])| >= ld_thresh for some causal.
#       cs_tp / cs_fp / fdr_cs / power_cs

eval_fit <- function(fit, X, true_idx, pip_thresh,
                     ld_thresh         = 0.5,
                     purity_thresh     = 0.8,
                     pip_high_thresh   = 0.5) {
  pip      <- fit$pip
  cs_count <- if (!is.null(fit$sets$cs)) length(fit$sets$cs) else 0L
  cs_purities <- if (!is.null(fit$sets$purity))
    fit$sets$purity[, "min.abs.corr"] else rep(NA_real_, cs_count)
  cs_purity <- if (cs_count > 0L) mean(cs_purities, na.rm = TRUE) else NA_real_

  # ---- SNP-level (loose, PIP threshold)
  selected <- which(pip > pip_thresh)
  n_disc   <- length(selected)
  n_true   <- length(true_idx)
  n_tp     <- sum(selected %in% true_idx)
  n_fp     <- n_disc - n_tp
  fdr      <- if (n_disc > 0L) n_fp / n_disc else 0
  power    <- if (n_true > 0L) n_tp / n_true else NA_real_

  # ---- CS-level
  cs_tp <- 0L; cs_fp <- 0L
  causal_covered <- integer(0L)
  high_purity_cs_members <- integer(0L)
  in_any_cs <- integer(0L)
  if (cs_count > 0L) {
    for (i in seq_along(fit$sets$cs)) {
      members <- fit$sets$cs[[i]]
      in_any_cs <- union(in_any_cs, members)
      if (!is.na(cs_purities[i]) && cs_purities[i] >= purity_thresh) {
        high_purity_cs_members <- union(high_purity_cs_members, members)
      }
      lead <- members[which.max(pip[members])]
      hit_causal <- lead %in% true_idx
      if (!hit_causal && length(true_idx) > 0L) {
        cors <- abs(suppressWarnings(stats::cor(X[, lead], X[, true_idx])))
        if (any(cors >= ld_thresh, na.rm = TRUE)) {
          hit_causal <- TRUE
          covered <- true_idx[which.max(cors)]
          causal_covered <- union(causal_covered, covered)
        }
      } else if (hit_causal) {
        causal_covered <- union(causal_covered, lead)
      }
      if (hit_causal) cs_tp <- cs_tp + 1L else cs_fp <- cs_fp + 1L
    }
  }
  fdr_cs   <- if (cs_count > 0L) cs_fp / cs_count else 0
  power_cs <- if (length(true_idx) > 0L)
    length(causal_covered) / length(true_idx) else NA_real_

  # ---- SNP-level hybrid: high-purity CS members OR not-in-CS high PIP
  not_in_cs_high_pip <- which(pip > pip_high_thresh &
                              !(seq_along(pip) %in% in_any_cs))
  hyb_disc   <- union(high_purity_cs_members, not_in_cs_high_pip)
  hyb_n_disc <- length(hyb_disc)
  hyb_n_tp   <- sum(hyb_disc %in% true_idx)
  hyb_n_fp   <- hyb_n_disc - hyb_n_tp
  fdr_hyb    <- if (hyb_n_disc > 0L) hyb_n_fp / hyb_n_disc else 0
  power_hyb  <- if (length(true_idx) > 0L)
    sum(true_idx %in% hyb_disc) / length(true_idx) else NA_real_

  list(n_disc = n_disc, n_tp = n_tp, n_fp = n_fp,
       fdr = fdr, power = power,
       cs_count = cs_count, cs_purity = cs_purity,
       cs_tp = cs_tp, cs_fp = cs_fp,
       fdr_cs = fdr_cs, power_cs = power_cs,
       hyb_n_disc = hyb_n_disc, hyb_n_tp = hyb_n_tp, hyb_n_fp = hyb_n_fp,
       fdr_hyb = fdr_hyb, power_hyb = power_hyb,
       niter = fit$niter, converged = isTRUE(fit$converged))
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
    metrics <- eval_fit(fit, d$X, true_indices, pip_thresh)
    fit_size_mb <- as.numeric(object.size(fit)) / (1024 * 1024)
    row_idx <- row_idx + 1L
    # Save the minimal fit shape needed to recompute any future metric
    # (pip + sets$cs + sets$purity + true_idx + rep_seed). Stored as
    # list-columns so the resulting data.frame still aggregates on the
    # numeric columns; aggregation formulas just skip these fields.
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
      cs_tp     = metrics$cs_tp,
      cs_fp     = metrics$cs_fp,
      fdr_cs    = metrics$fdr_cs,
      power_cs  = metrics$power_cs,
      hyb_n_disc = metrics$hyb_n_disc,
      hyb_n_tp   = metrics$hyb_n_tp,
      hyb_n_fp   = metrics$hyb_n_fp,
      fdr_hyb    = metrics$fdr_hyb,
      power_hyb  = metrics$power_hyb,
      niter     = metrics$niter,
      converged = metrics$converged,
      runtime_s = t_elapsed,
      fit_size_mb = fit_size_mb,
      mem_used_mb = mem_after - mem_before,
      rep_seed  = 100L + rep_i,
      fit_pip    = I(list(fit$pip)),
      fit_cs     = I(list(fit$sets$cs)),
      fit_purity = I(list(fit$sets$purity)),
      fit_true_idx = I(list(true_indices)),
      stringsAsFactors = FALSE
    )
    cat(sprintf("rep=%d cell=%d (%s, qnorm=%s, mnw=%s)  fdr=%.3f power=%.3f  fdr_cs=%.3f power_cs=%.3f  fdr_hyb=%.3f power_hyb=%.3f  cs_tp=%d cs_fp=%d  t=%.1fs\n",
                rep_i, g, cell$prior_variance_scope,
                cell$wavelet_qnorm,
                if (is.na(cell$mixture_null_weight)) "NA"
                else format(cell$mixture_null_weight),
                metrics$fdr, metrics$power,
                metrics$fdr_cs, metrics$power_cs,
                metrics$fdr_hyb, metrics$power_hyb,
                metrics$cs_tp, metrics$cs_fp,
                t_elapsed))
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

# Aggregate by cell directly (cell uniquely identifies a row of the
# bench grid, including the per_scale_normal cells whose mnw is NA).
# The list-columns (fit_pip / fit_cs / fit_purity / fit_true_idx) are
# dropped from the formula so aggregate stays on the numeric columns.
summary_df <- aggregate(
  cbind(fdr, power, n_disc, cs_count, cs_purity,
        cs_tp, cs_fp, fdr_cs, power_cs,
        hyb_n_disc, hyb_n_tp, hyb_n_fp, fdr_hyb, power_hyb,
        niter, runtime_s, fit_size_mb)
    ~ cell + prior_variance_scope + wavelet_qnorm,
  data = transform(results_df, mixture_null_weight = NULL,
                   fit_pip = NULL, fit_cs = NULL,
                   fit_purity = NULL, fit_true_idx = NULL),
  FUN = mean, na.action = na.pass)
summary_df <- summary_df[order(summary_df$cell), ]
write.csv(summary_df, out_csv, row.names = FALSE)
cat(sprintf("Saved per-cell summary to %s\n", out_csv))
print(summary_df, digits = 3)
