# Follow-up to benchmark_per_scale_normal_6grid.R: same 6-cell prior grid
# but two new scenarios meant to expose what the Gaussian baseline could not:
#
#   1. "heavy_tailed_signal": same causal structure as the baseline plus
#      18% per-cell outlier contamination on each Y (sd = 4 vs sd = 1
#      Gaussian noise floor). Tests whether wavelet_qnorm = TRUE recovers
#      power that wavelet_qnorm = FALSE loses on heavy-tailed Y.
#
#   2. "null_no_signal": no causal SNPs, pure Gaussian noise. Tests the
#      type-I rate of each cell at the 0.05 PIP threshold and whether any
#      cell spuriously declares a credible set. Power and FDR are not
#      meaningful in this scenario; we report n_disc, has_disc (any
#      discovery T/F), cs_count, and runtime.
#
# Per Gao 2026-05-03 Slack the prior grid is held fixed at the 6 cells of
# the baseline (per_scale x {mnw 0.05, 0} x {qnorm F, T} plus
# per_scale_normal x {qnorm F, T}). Only the scenario axis is new.
#
# Estimated wall-clock based on the 36 min baseline (one Gaussian
# scenario, 30 fits): ~60-90 min for 60 fits across two scenarios.
# Submitted with 3 h SLURM time-limit for headroom.
#
# Usage:
#   Rscript inst/bench/profiling/benchmark_heavy_tailed_null_6grid.R
#
# Output:
#   inst/bench/profiling/results/heavy_tailed_null_6grid_<timestamp>.rds
#   inst/bench/profiling/results/heavy_tailed_null_6grid_<timestamp>_summary.csv

suppressPackageStartupMessages({
  library(devtools)
  load_all(quiet = TRUE)
})

# ---- Parameters (mirror the baseline so cells are directly comparable) ----

set.seed(2026L * 503L + 21L)
n     <- 84L
p     <- 500L
M     <- 2L
Tlen  <- 64L
n_rep <- 5L
L     <- 10L
true_indices_signal <- c(50L, 220L, 380L)
pip_thresh          <- 0.05

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

scenarios <- c("heavy_tailed_signal", "null_no_signal")

# ---- Per-scenario data builder -------------------------------------

build_dataset <- function(scenario, rep_seed) {
  set.seed(rep_seed)
  X <- matrix(rnorm(n * p), n, p)

  if (scenario == "heavy_tailed_signal") {
    effect_curve <- numeric(Tlen); effect_curve[20:40] <- 1
    beta_per_snp <- matrix(0, p, Tlen)
    for (j in true_indices_signal) beta_per_snp[j, ] <- effect_curve
    Y <- vector("list", M)
    for (m in seq_len(M)) {
      sig    <- X %*% beta_per_snp
      noise  <- matrix(rnorm(n * Tlen, sd = 1.0), n, Tlen)
      # 18% outlier contamination, sd = 4 (mirrors data-raw outcome 2).
      mask   <- matrix(rbinom(n * Tlen, 1L, 0.18), n, Tlen)
      outl   <- mask * matrix(rnorm(n * Tlen, sd = 4.0), n, Tlen)
      Y[[m]] <- sig + noise + outl
    }
    return(list(X = X, Y = Y, true_idx = true_indices_signal,
                scenario = scenario))
  }

  if (scenario == "null_no_signal") {
    Y <- vector("list", M)
    for (m in seq_len(M)) {
      Y[[m]] <- matrix(rnorm(n * Tlen, sd = 1.0), n, Tlen)
    }
    return(list(X = X, Y = Y, true_idx = integer(0L),
                scenario = scenario))
  }

  stop(sprintf("Unknown scenario: %s", scenario))
}

# ---- Per-cell metrics ------------------------------------------------
#
# Two metrics, both reported per (scenario, cell, rep). See
# benchmark_per_scale_normal_6grid.R::eval_fit for the SNP-level vs
# CS-level definitions; this version adds the null-scenario branch
# where FDR / power are undefined and per-CS judgments are all FP
# by construction.

eval_fit <- function(fit, X, true_idx, pip_thresh, scenario,
                     ld_thresh       = 0.5,
                     purity_thresh   = 0.8,
                     pip_high_thresh = 0.5) {
  pip      <- fit$pip
  selected <- which(pip > pip_thresh)
  n_disc   <- length(selected)
  n_true   <- length(true_idx)
  cs_count <- if (!is.null(fit$sets$cs)) length(fit$sets$cs) else 0L
  cs_purities <- if (!is.null(fit$sets$purity))
    fit$sets$purity[, "min.abs.corr"] else rep(NA_real_, cs_count)
  cs_purity <- if (cs_count > 0L) mean(cs_purities, na.rm = TRUE) else NA_real_

  # ---- CS-level + collect high-purity CS members for hybrid SNP-level
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

  # ---- Hybrid SNP-level: high-purity CS members OR not-in-CS PIP > pip_high_thresh
  not_in_cs_high_pip <- which(pip > pip_high_thresh &
                              !(seq_along(pip) %in% in_any_cs))
  hyb_disc   <- union(high_purity_cs_members, not_in_cs_high_pip)
  hyb_n_disc <- length(hyb_disc)
  hyb_n_tp   <- sum(hyb_disc %in% true_idx)
  hyb_n_fp   <- hyb_n_disc - hyb_n_tp
  fdr_hyb    <- if (hyb_n_disc > 0L) hyb_n_fp / hyb_n_disc else 0
  power_hyb  <- if (length(true_idx) > 0L)
    sum(true_idx %in% hyb_disc) / length(true_idx) else NA_real_

  if (scenario == "null_no_signal") {
    # SNP-level FDR / power undefined under no signal; every CS and
    # every hybrid discovery is a false positive by construction.
    return(list(n_disc = n_disc, n_tp = NA_integer_, n_fp = n_disc,
                fdr = NA_real_, power = NA_real_,
                has_disc = (n_disc > 0L),
                cs_count = cs_count, cs_purity = cs_purity,
                cs_tp = 0L, cs_fp = cs_fp,
                fdr_cs = if (cs_count > 0L) 1 else 0,
                power_cs = NA_real_,
                hyb_n_disc = hyb_n_disc,
                hyb_n_tp = 0L, hyb_n_fp = hyb_n_disc,
                fdr_hyb = if (hyb_n_disc > 0L) 1 else 0,
                power_hyb = NA_real_,
                niter = fit$niter, converged = isTRUE(fit$converged)))
  }

  n_tp  <- sum(selected %in% true_idx)
  n_fp  <- n_disc - n_tp
  fdr   <- if (n_disc > 0L) n_fp / n_disc else 0
  power <- n_tp / n_true
  fdr_cs   <- if (cs_count > 0L) cs_fp / cs_count else 0
  power_cs <- length(causal_covered) / length(true_idx)
  list(n_disc = n_disc, n_tp = n_tp, n_fp = n_fp, fdr = fdr, power = power,
       has_disc = (n_disc > 0L),
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

for (sc in scenarios) {
  for (rep_i in seq_len(n_rep)) {
    d <- build_dataset(sc, rep_seed = 200L + rep_i)
    for (g in seq_len(nrow(bench_grid))) {
      cell <- bench_grid[g, ]
      fit_args <- list(
        X                     = d$X,
        Y                     = d$Y,
        L                     = L,
        verbose               = FALSE,
        prior_variance_scope  = cell$prior_variance_scope,
        wavelet_qnorm         = cell$wavelet_qnorm,
        save_mu_method        = "alpha_collapsed"
      )
      if (!is.na(cell$mixture_null_weight)) {
        fit_args$mixture_null_weight <- cell$mixture_null_weight
      }
      t_start    <- Sys.time()
      mem_before <- gc(reset = TRUE, verbose = FALSE)[2L, 6L]
      fit <- tryCatch(
        do.call(mfsusie, fit_args),
        error = function(e) {
          warning(sprintf("scenario=%s rep=%d cell=%d errored: %s",
                          sc, rep_i, g, conditionMessage(e)))
          NULL
        }
      )
      t_elapsed <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
      mem_after <- gc(verbose = FALSE)[2L, 6L]
      if (is.null(fit)) next
      m <- eval_fit(fit, d$X, d$true_idx, pip_thresh, sc)
      fit_size_mb <- as.numeric(object.size(fit)) / (1024 * 1024)
      row_idx <- row_idx + 1L
      results[[row_idx]] <- data.frame(
        scenario              = sc,
        rep                   = rep_i,
        cell                  = g,
        wavelet_qnorm         = cell$wavelet_qnorm,
        prior_variance_scope  = cell$prior_variance_scope,
        mixture_null_weight   = cell$mixture_null_weight,
        n_disc    = m$n_disc, n_tp = m$n_tp, n_fp = m$n_fp,
        fdr = m$fdr, power = m$power, has_disc = m$has_disc,
        cs_count = m$cs_count, cs_purity = m$cs_purity,
        cs_tp = m$cs_tp, cs_fp = m$cs_fp,
        fdr_cs = m$fdr_cs, power_cs = m$power_cs,
        hyb_n_disc = m$hyb_n_disc, hyb_n_tp = m$hyb_n_tp, hyb_n_fp = m$hyb_n_fp,
        fdr_hyb = m$fdr_hyb, power_hyb = m$power_hyb,
        niter = m$niter, converged = m$converged,
        runtime_s = t_elapsed, fit_size_mb = fit_size_mb,
        mem_used_mb = mem_after - mem_before,
        rep_seed   = 200L + rep_i,
        fit_pip    = I(list(fit$pip)),
        fit_cs     = I(list(fit$sets$cs)),
        fit_purity = I(list(fit$sets$purity)),
        fit_true_idx = I(list(d$true_idx)),
        stringsAsFactors = FALSE
      )
      cat(sprintf("scn=%s rep=%d cell=%d (%s, qnorm=%s, mnw=%s)  fdr=%s power=%s  fdr_cs=%s power_cs=%s  fdr_hyb=%s power_hyb=%s  cs_tp=%d cs_fp=%d  t=%.1fs\n",
                  sc, rep_i, g, cell$prior_variance_scope, cell$wavelet_qnorm,
                  if (is.na(cell$mixture_null_weight)) "NA"
                  else format(cell$mixture_null_weight),
                  if (is.na(m$fdr)) "NA" else sprintf("%.3f", m$fdr),
                  if (is.na(m$power)) "NA" else sprintf("%.3f", m$power),
                  if (is.na(m$fdr_cs)) "NA" else sprintf("%.3f", m$fdr_cs),
                  if (is.na(m$power_cs)) "NA" else sprintf("%.3f", m$power_cs),
                  if (is.na(m$fdr_hyb)) "NA" else sprintf("%.3f", m$fdr_hyb),
                  if (is.na(m$power_hyb)) "NA" else sprintf("%.3f", m$power_hyb),
                  m$cs_tp, m$cs_fp, t_elapsed))
    }
  }
}

t_total <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
cat(sprintf("\nTotal wall-clock: %.1f s\n", t_total))

results_df <- do.call(rbind, results)

# ---- Persist outputs ------------------------------------------------

out_dir <- "inst/bench/profiling/results"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
ts  <- format(Sys.time(), "%Y%m%d_%H%M")
out_rds <- file.path(out_dir, sprintf("heavy_tailed_null_6grid_%s.rds", ts))
out_csv <- file.path(out_dir,
                     sprintf("heavy_tailed_null_6grid_%s_summary.csv", ts))
saveRDS(results_df, out_rds)
cat(sprintf("Saved per-replicate results to %s\n", out_rds))

# Aggregate per (scenario, cell). Drop mixture_null_weight from the formula
# so per_scale_normal rows (NA mnw) are not silently dropped. Drop the
# fit_* list-columns since aggregate cannot reduce them.
agg_df <- aggregate(
  cbind(fdr, power, n_disc, has_disc, cs_count, cs_purity,
        cs_tp, cs_fp, fdr_cs, power_cs,
        hyb_n_disc, hyb_n_tp, hyb_n_fp, fdr_hyb, power_hyb,
        niter, runtime_s, fit_size_mb)
    ~ scenario + cell + prior_variance_scope + wavelet_qnorm,
  data = transform(results_df, fit_pip = NULL, fit_cs = NULL,
                   fit_purity = NULL, fit_true_idx = NULL),
  FUN = mean, na.action = na.pass)
agg_df <- agg_df[order(agg_df$scenario, agg_df$cell), ]
write.csv(agg_df, out_csv, row.names = FALSE)
cat(sprintf("Saved per-cell summary to %s\n", out_csv))
print(agg_df, digits = 3)
