#!/usr/bin/env Rscript
# Vignette sweep for `per_scale_normal` vs default `per_outcome`.
#
# For each user-facing vignette in the package, reproduce the
# headline `fsusie()` / `mfsusie()` call twice (once per scope)
# and report nCS, CS leads, TP/FP per known causal set, and
# wallclock time. Output goes to
# `inst/bench/profiling/per_scale_normal_vignette_sweep.tsv`.
#
# A CS is counted as TP if its lead variable is in the causal
# set OR is LD-correlated (|cor(X[, lead], X[, causal])| >= 0.5)
# with any causal. Otherwise it's an FP. A truly missed CS would
# show up as an under-count of TP relative to the true number
# of causals.

suppressPackageStartupMessages({
  devtools::load_all(quiet = TRUE)
  library(susieR)
})

# Helpers ------------------------------------------------------------------

# CS classifier: TP if lead matches causal or LD >= 0.5 with causal.
classify_cs <- function(fit, X, causal, ld_thresh = 0.5) {
  cs <- fit$sets$cs
  if (length(cs) == 0L) return(list(tp = 0L, fp = 0L, leads = integer(0L)))
  leads <- vapply(cs, function(s) s[which.max(fit$pip[s])], integer(1L))
  is_tp <- vapply(leads, function(l) {
    if (l %in% causal) return(TRUE)
    if (length(causal) == 0L) return(FALSE)
    cors <- abs(suppressWarnings(stats::cor(X[, l], X[, causal])))
    any(cors >= ld_thresh, na.rm = TRUE)
  }, logical(1L))
  list(tp = sum(is_tp), fp = sum(!is_tp), leads = leads, tp_mask = is_tp)
}

# One workload: run both scopes, classify, print one row.
run_workload <- function(name, X, fit_fn, causal) {
  out_rows <- list()
  for (scope in c("per_outcome", "per_scale_normal")) {
    set.seed(1)
    t0 <- Sys.time()
    fit <- suppressWarnings(fit_fn(scope))
    dt  <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    cls <- classify_cs(fit, X, causal)
    cat(sprintf("  %-22s nCS=%d  TP=%d  FP=%d  niter=%2d  time=%6.2fs   leads=[%s]\n",
                scope, length(fit$sets$cs), cls$tp, cls$fp,
                fit$niter, dt,
                paste(cls$leads, collapse = ",")))
    out_rows[[scope]] <- data.frame(
      workload = name, scope = scope,
      n_causal = length(causal),
      nCS = length(fit$sets$cs),
      TP = cls$tp, FP = cls$fp,
      niter = fit$niter, time_s = round(dt, 3),
      leads = paste(cls$leads, collapse = ","),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, out_rows)
}

# Workloads ----------------------------------------------------------------

all_rows <- list()
cat("\n--- 1. fsusie_intro (coverage_example) ---\n")
data(coverage_example)
all_rows[["fsusie_intro"]] <- run_workload(
  "fsusie_intro",
  coverage_example$X,
  function(scope) fsusie(coverage_example$Y, coverage_example$X,
                          prior_variance_scope = scope, verbose = FALSE),
  coverage_example$causal_snps)

cat("\n--- 2. fsusie_why_functional (synthetic) ---\n")
local({
  set.seed(1); data(N3finemapping)
  X_full <- N3finemapping$X[seq_len(100), ]; T_m <- 128L
  effect <- 1.2 * cos(seq_len(T_m) / T_m * 3 * pi)
  effect[effect < 0] <- 0; effect[seq_len(40)] <- 0
  y <- matrix(0, 100L, T_m)
  for (i in seq_len(100L)) y[i, ] <- X_full[i, 700L] * effect + rnorm(T_m)
  keep <- which(apply(X_full, 2, var) > 0)
  X <- X_full[, keep]
  causal <- which(keep == 700L)
  all_rows[["why_functional"]] <<- run_workload(
    "why_functional", X,
    function(scope) fsusie(y, X, prior_variance_scope = scope, verbose = FALSE),
    causal)
})

cat("\n--- 3. fsusie_colocalization (Y1, Y2, Y3) ---\n")
local({
  set.seed(1); data(N3finemapping)
  X <- N3finemapping$X[, seq_len(150)]
  n <- nrow(X); p <- ncol(X); T_m <- 64
  pos <- seq_len(T_m)
  shape <- exp(-((pos - T_m / 2)^2) / (2 * (T_m / 6)^2))
  beta1 <- matrix(0, p, T_m); beta1[42, ] <- 1.2 * shape
  beta2 <- matrix(0, p, T_m); beta2[42, ] <- 0.9 * shape
  beta3 <- matrix(0, p, T_m); beta3[88, ] <- 0.9 * shape
  Y1 <- X %*% beta1 + matrix(rnorm(n * T_m, sd = 0.4), n)
  Y2 <- X %*% beta2 + matrix(rnorm(n * T_m, sd = 0.4), n)
  Y3 <- X %*% beta3 + matrix(rnorm(n * T_m, sd = 0.4), n)
  workloads <- list(coloc_Y1 = list(Y = Y1, causal = 42L),
                    coloc_Y2 = list(Y = Y2, causal = 42L),
                    coloc_Y3 = list(Y = Y3, causal = 88L))
  for (tag in names(workloads)) {
    Y      <- workloads[[tag]]$Y
    causal <- workloads[[tag]]$causal
    all_rows[[tag]] <<- run_workload(tag, X,
      function(scope) fsusie(Y, X, pos = pos,
                             prior_variance_scope = scope, verbose = FALSE),
      causal)
  }
})

cat("\n--- 4. fsusie_covariates_adjustment (raw + adjusted) ---\n")
local({
  set.seed(2L)
  rsnr <- 1; pos1 <- 25L; pos2 <- 75L; lev_res <- 7L; T_m <- 2L^lev_res
  f1 <- mf_simu_ibss_per_level(lev_res)$sim_func
  f2 <- mf_simu_ibss_per_level(lev_res)$sim_func
  f1_cov <- mf_simu_ibss_per_level(lev_res)$sim_func
  f2_cov <- mf_simu_ibss_per_level(lev_res)$sim_func
  f3_cov <- mf_simu_ibss_per_level(lev_res)$sim_func
  data(N3finemapping)
  Geno <- N3finemapping$X[, seq_len(100)]
  n <- nrow(Geno)
  Cov <- matrix(rnorm(3L * n, sd = 2), ncol = 3L)
  genetic_signal <- matrix(0, n, T_m); covariate_signal <- matrix(0, n, T_m)
  for (i in seq_len(n)) {
    noise_i <- rnorm(T_m, sd = (1 / rsnr) * var(f1))
    genetic_signal[i, ]   <- Geno[i, pos1] * f1 + Geno[i, pos2] * f2 + noise_i
    covariate_signal[i, ] <- Cov[i, 1L] * f1_cov + Cov[i, 2L] * f2_cov + Cov[i, 3L] * f3_cov
  }
  Y <- covariate_signal + genetic_signal
  adj <- mf_adjust_for_covariates(Y, Cov)
  causal <- c(pos1, pos2)
  all_rows[["covar_adjusted"]] <<- run_workload("covar_adjusted", Geno,
    function(scope) fsusie(adj$Y_adjusted, Geno,
                            prior_variance_scope = scope, verbose = FALSE),
    causal)
  all_rows[["covar_unadjusted"]] <<- run_workload("covar_unadjusted", Geno,
    function(scope) fsusie(Y, Geno,
                            prior_variance_scope = scope, verbose = FALSE),
    causal)
})

cat("\n--- 5. fsusie_dnam_case_study (dnam_example) ---\n")
data(dnam_example)
all_rows[["dnam"]] <- run_workload("dnam",
  dnam_example$X,
  function(scope) fsusie(dnam_example$Y, dnam_example$X,
                          prior_variance_scope = scope, verbose = FALSE),
  dnam_example$causal_snps)

cat("\n--- 6. fsusie_gtex_case_study (gtex_example) ---\n")
data(gtex_example)
all_rows[["gtex"]] <- run_workload("gtex",
  gtex_example$X,
  function(scope) fsusie(gtex_example$Y, gtex_example$X,
                          pos = gtex_example$pos,
                          prior_variance_scope = scope, verbose = FALSE),
  gtex_example$causal_snps)

cat("\n--- 7. getting_started (fsusie + mfsusie) ---\n")
local({
  set.seed(1)
  n <- 200; p <- 50; T_m <- 32
  X <- matrix(rnorm(n * p), nrow = n)
  beta <- numeric(p); beta[c(5L, 17L)] <- c(1.2, -0.8)
  Y <- X %*% matrix(rep(beta, T_m), nrow = p) +
       matrix(rnorm(n * T_m, sd = 0.4), nrow = n)
  all_rows[["getting_started_f"]] <<- run_workload("getting_started_f",
    X,
    function(scope) fsusie(Y, X, prior_variance_scope = scope, verbose = FALSE),
    c(5L, 17L))
  T_per <- c(32L, 64L)
  Y_func <- lapply(T_per, function(T_m)
    X %*% matrix(rep(beta, T_m), nrow = p) +
      matrix(rnorm(n * T_m, sd = 0.4), nrow = n))
  Y_scalar <- as.numeric(X %*% beta + rnorm(n, sd = 0.4))
  Y_list <- c(Y_func, list(matrix(Y_scalar, ncol = 1)))
  all_rows[["getting_started_m"]] <<- run_workload("getting_started_m",
    X,
    function(scope) mfsusie(X, Y_list, prior_variance_scope = scope, verbose = FALSE),
    c(5L, 17L))
})

cat("\n--- 8. mfsusie_intro (multiomic_example) ---\n")
data(multiomic_example)
Y_list <- multiomic_example$Y_list[c("DNAm", "RNA-seq")]
pos_list <- multiomic_example$pos_list[c("DNAm", "RNA-seq")]
all_rows[["mfsusie_intro"]] <- run_workload("mfsusie_intro",
  multiomic_example$X,
  function(scope) mfsusie(multiomic_example$X, Y_list, pos = pos_list,
                          prior_variance_scope = scope, verbose = FALSE),
  multiomic_example$causal_snps)

cat("\n--- 9. mfsusie_long_running_fits (signal + null loci) ---\n")
local({
  set.seed(1); data(N3finemapping)
  X <- N3finemapping$X[, seq_len(120)]
  n <- nrow(X); p <- ncol(X); T_func <- c(32L, 32L)
  beta_signal <- numeric(p); beta_signal[c(37L, 88L)] <- c(1, -0.6)
  Y_signal <- lapply(T_func, function(T_m)
    X %*% matrix(rep(beta_signal, T_m), nrow = p) +
      matrix(rnorm(n * T_m, sd = 0.4), nrow = n))
  Y_null <- lapply(T_func, function(T_m)
    matrix(rnorm(n * T_m, sd = 0.4), nrow = n))
  all_rows[["long_signal"]] <<- run_workload("long_signal", X,
    function(scope) mfsusie(X, Y_signal,
                             prior_variance_scope = scope, verbose = FALSE),
    c(37L, 88L))
  all_rows[["long_null"]] <<- run_workload("long_null", X,
    function(scope) mfsusie(X, Y_null,
                             prior_variance_scope = scope, verbose = FALSE),
    integer(0L))
})

cat("\n--- 10. post_processing (synthetic) ---\n")
local({
  set.seed(1); data(N3finemapping)
  X <- N3finemapping$X[, seq_len(150)]
  n <- nrow(X); p <- ncol(X); T_m <- 64
  pos <- seq_len(T_m)
  shape <- exp(-((pos - T_m / 2)^2) / (2 * (T_m / 8)^2))
  beta <- matrix(0, p, T_m); beta[42, ] <- 1.5 * shape; beta[97, ] <- -1.0 * shape
  Y <- X %*% beta + matrix(rnorm(n * T_m, sd = 0.5), n)
  all_rows[["post_processing"]] <<- run_workload("post_processing",
    X,
    function(scope) fsusie(Y, X, pos = pos,
                            prior_variance_scope = scope, verbose = FALSE),
    c(42L, 97L))
})

# Output --------------------------------------------------------------------
out <- do.call(rbind, all_rows)
out_path <- file.path("inst", "bench", "profiling",
                       "per_scale_normal_vignette_sweep.tsv")
write.table(out, out_path, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("\nWrote %s (%d rows)\n", out_path, nrow(out)))

# Per-workload summary: nCS / TP / FP / time per scope
cat("\n=== Summary ===\n")
po  <- out[out$scope == "per_outcome",      ]
psn <- out[out$scope == "per_scale_normal", ]
m   <- merge(po, psn, by = "workload", suffixes = c("_po", "_psn"))
m <- m[, c("workload", "nCS_po", "nCS_psn",
            "TP_po", "TP_psn", "FP_po", "FP_psn",
            "time_s_po", "time_s_psn")]
m$speedup <- round(m$time_s_po / m$time_s_psn, 2)
print(m, row.names = FALSE)

cat(sprintf("\nTotal time per scope: per_outcome=%.1fs, per_scale_normal=%.1fs (overall speedup %.2fx)\n",
            sum(po$time_s), sum(psn$time_s),
            sum(po$time_s) / sum(psn$time_s)))
cat(sprintf("Total TP: per_outcome=%d, per_scale_normal=%d\n",
            sum(po$TP), sum(psn$TP)))
cat(sprintf("Total FP: per_outcome=%d, per_scale_normal=%d\n",
            sum(po$FP), sum(psn$FP)))
