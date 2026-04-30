# Sweep all vignette workloads under three settings:
#   1. per_outcome,      mixture_null_weight = 0.05 (current default)
#   2. per_outcome,      mixture_null_weight = 0
#   3. per_scale_normal  (no mixture_null_weight, ebnm-driven)
#
# Reports nCS, TP, FP, FP CS members, runtime per configuration.
#
# Usage: Rscript inst/bench/profiling/sweep_vignettes_nullweight.R

suppressPackageStartupMessages({
  library(mfsusieR)
  library(susieR)
})

score <- function(fit, causals) {
  cs <- if (!is.null(fit$sets$cs)) fit$sets$cs else fit$cs
  hits <- vapply(cs, function(idx) any(idx %in% causals), logical(1L))
  tp <- sum(hits); fp <- sum(!hits)
  fp_str <- if (fp == 0L) "none" else
    paste0(sapply(which(!hits), function(k) paste(cs[[k]], collapse=",")),
           collapse="; ")
  list(ncs = length(cs), tp = tp, fp = fp, fp_str = fp_str,
       pip = round(fit$pip[causals], 3L))
}

run_three <- function(label, X, Y, causals, pos_list = NULL,
                      multi = FALSE, ...) {
  fit_fn <- if (multi) mfsusie else fsusie
  base   <- if (multi) list(X = X, Y = Y, pos = pos_list, ...)
            else        list(Y = Y, X = X, pos = pos_list[[1L]], ...)
  configs <- list(
    list(tag = "po nw=0.05",      args = list(prior_variance_scope = "per_outcome",
                                              mixture_null_weight = 0.05)),
    list(tag = "po nw=0",         args = list(prior_variance_scope = "per_outcome",
                                              mixture_null_weight = 0)),
    list(tag = "per_scale_normal", args = list(prior_variance_scope = "per_scale_normal"))
  )
  for (cfg in configs) {
    t0 <- Sys.time()
    fit <- tryCatch(
      do.call(fit_fn, c(base, cfg$args,
                        list(verbose = T))),
      error = function(e) { cat(sprintf("  %-30s | %-18s | <ERROR: %s>\n",
                                          label, cfg$tag, conditionMessage(e))); NULL })
    t_el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    if (is.null(fit)) next
    s <- score(fit, causals)
    cat(sprintf("  %-30s | %-18s | %5.1fs niter=%2d | nCS=%d TP=%d FP=%d  FP_CS=%-15s | PIP=%s\n",
                label, cfg$tag, t_el, fit$niter,
                s$ncs, s$tp, s$fp, substr(s$fp_str, 1L, 15L),
                paste(s$pip, collapse=",")))
  }
}

# ----- vignette workloads (same as sweep_vignettes.R) -----

cat("\n--- 1. getting_started.Rmd (single, n=200, p=50, T=32) ---\n")
set.seed(1L)
n <- 200L; p <- 50L; T_m <- 32L
X <- matrix(rnorm(n * p), n)
beta <- numeric(p); beta[c(5L, 17L)] <- c(1.2, -0.8)
Y <- X %*% matrix(rep(beta, T_m), p) + matrix(rnorm(n * T_m, sd = 0.4), n)
run_three("getting_started/single", X, Y, c(5L, 17L), pos_list = list(seq_len(T_m)))

cat("\n--- 2. getting_started.Rmd (multi, M=3) ---\n")
set.seed(1L)
n <- 200L; p <- 50L
X <- matrix(rnorm(n * p), n)
beta <- numeric(p); beta[c(5L, 17L)] <- c(1.2, -0.8)
T_per <- c(32L, 64L)
Y_func <- lapply(T_per, function(T_m)
  X %*% matrix(rep(beta, T_m), p) + matrix(rnorm(n * T_m, sd = 0.4), n))
Y_scalar <- as.numeric(X %*% beta + rnorm(n, sd = 0.4))
Y_list <- c(Y_func, list(matrix(Y_scalar, ncol = 1L)))
run_three("getting_started/multi", X, Y_list, c(5L, 17L),
          pos_list = lapply(Y_list, function(y) seq_len(ncol(y))), multi = TRUE)

cat("\n--- 3. fsusie_intro.Rmd (coverage_example) ---\n")
data(coverage_example)
run_three("fsusie_intro/coverage", coverage_example$X, coverage_example$Y,
          coverage_example$causal_snps,
          pos_list = list(seq_len(ncol(coverage_example$Y))))

cat("\n--- 4. fsusie_why_functional.Rmd ---\n")
set.seed(1L)
data(N3finemapping)
X_full <- N3finemapping$X[seq_len(100L), ]
keep <- which(apply(X_full, 2L, var) > 0)
X_w <- X_full[, keep]
true_pos_filt <- which(keep == 700L)
T_m <- 128L
positions <- seq_len(T_m)
effect <- 1.2 * cos(2 * pi * (positions - 41L) / 80L); effect[positions < 41L] <- 0
beta_f <- matrix(0, ncol(X_w), T_m); beta_f[true_pos_filt, ] <- effect
Y_f <- X_w %*% beta_f + matrix(rnorm(nrow(X_w) * T_m, sd = 0.4), nrow(X_w))
run_three("fsusie_why_functional", X_w, Y_f, true_pos_filt,
          pos_list = list(positions))

cat("\n--- 5. fsusie_covariates_adjustment.Rmd ---\n")
set.seed(2L)
lev_res <- 7L; T_m <- 2L^lev_res
f1     <- mf_simu_ibss_per_level(lev_res)$sim_func
f2     <- mf_simu_ibss_per_level(lev_res)$sim_func
f1_cov <- mf_simu_ibss_per_level(lev_res)$sim_func
f2_cov <- mf_simu_ibss_per_level(lev_res)$sim_func
f3_cov <- mf_simu_ibss_per_level(lev_res)$sim_func
data(N3finemapping)
Geno <- N3finemapping$X[, seq_len(100L)]
n_c <- nrow(Geno)
Cov <- matrix(rnorm(3L * n_c, sd = 2), ncol = 3L)
gen <- matrix(0, n_c, T_m); covsig <- matrix(0, n_c, T_m)
for (i in seq_len(n_c)) {
  noise_i <- rnorm(T_m, sd = var(f1))
  gen[i, ] <- Geno[i, 25L] * f1 + Geno[i, 75L] * f2 + noise_i
  covsig[i, ] <- Cov[i, 1L] * f1_cov + Cov[i, 2L] * f2_cov + Cov[i, 3L] * f3_cov
}
adj <- mf_adjust_for_covariates(covsig + gen, Cov)
run_three("covariates_adj", Geno, adj$Y_adjusted, c(25L, 75L),
          pos_list = list(seq_len(T_m)))

cat("\n--- 6. fsusie_colocalization.Rmd ---\n")
set.seed(1L)
data(N3finemapping)
X_co <- N3finemapping$X[, seq_len(150L)]
n_co <- nrow(X_co); p_co <- ncol(X_co); T_m <- 64L
positions <- seq_len(T_m)
shape <- exp(-((positions - T_m / 2)^2) / (2 * (T_m / 6)^2))
b <- function(j, a) { m <- matrix(0, p_co, T_m); m[j, ] <- a * shape; m }
for (case in list(list(tag="Y1 (causal=42)", b = b(42L, 1.2),  c = 42L),
                  list(tag="Y2 (causal=42)", b = b(42L, 0.9),  c = 42L),
                  list(tag="Y3 (causal=88)", b = b(88L, 0.9),  c = 88L))) {
  Y <- X_co %*% case$b + matrix(rnorm(n_co * T_m, sd = 0.4), n_co)
  run_three(paste0("coloc/", case$tag), X_co, Y, case$c, pos_list = list(positions))
}

cat("\n--- 7. fsusie_dnam_case_study.Rmd ---\n")
data(dnam_example)
run_three("dnam_example", dnam_example$X, dnam_example$Y,
          dnam_example$causal_snps,
          pos_list = list(if (!is.null(dnam_example$pos)) dnam_example$pos
                          else seq_len(ncol(dnam_example$Y))))

cat("\n--- 8. fsusie_gtex_case_study.Rmd ---\n")
data(gtex_example)
run_three("gtex_example", gtex_example$X, gtex_example$Y,
          gtex_example$causal_snps,
          pos_list = list(if (!is.null(gtex_example$pos)) gtex_example$pos
                          else seq_len(ncol(gtex_example$Y))))

cat("\n--- 9. post_processing.Rmd (single) ---\n")
set.seed(1L)
data(N3finemapping)
X_p <- N3finemapping$X[, seq_len(150L)]
n_p <- nrow(X_p); p_p <- ncol(X_p); T_m <- 64L
positions <- seq_len(T_m)
shape <- exp(-((positions - T_m / 2)^2) / (2 * (T_m / 8)^2))
beta_p <- matrix(0, p_p, T_m)
beta_p[42L, ] <-  1.5 * shape; beta_p[97L, ] <- -1.0 * shape
Y_p <- X_p %*% beta_p + matrix(rnorm(n_p * T_m, sd = 0.5), n_p)
run_three("post_processing/single", X_p, Y_p, c(42L, 97L),
          pos_list = list(positions))

cat("\n--- 10. mfsusie_intro.Rmd (multiomic_example) ---\n")
data(multiomic_example)
Y_list_mi <- multiomic_example$Y_list[c("DNAm", "RNA-seq")]
pos_list_mi <- multiomic_example$pos_list[c("DNAm", "RNA-seq")]
run_three("mfsusie_intro/2outcome", multiomic_example$X, Y_list_mi,
          multiomic_example$causal_snps,
          pos_list = pos_list_mi, multi = TRUE)

cat("\n--- 11. mfsusie_long_running_fits.Rmd ---\n")
set.seed(1L)
data(N3finemapping)
X_lr <- N3finemapping$X[, seq_len(120L)]
n_lr <- nrow(X_lr); p_lr <- ncol(X_lr); T_per <- c(32L, 32L)
beta_a <- numeric(p_lr); beta_a[c(37L, 88L)] <- c(1.0, -0.6)
Y_list_lr <- lapply(T_per, function(T_m)
  X_lr %*% matrix(rep(beta_a, T_m), p_lr) + matrix(rnorm(n_lr * T_m, sd = 0.4), n_lr))
run_three("long_running/signal", X_lr, Y_list_lr, c(37L, 88L),
          pos_list = lapply(Y_list_lr, function(y) seq_len(ncol(y))),
          multi = TRUE)

cat("\n--- 12. failing case (n=200 p=200 T=128) ---\n")
set.seed(1L)
n <- 200L; p <- 200L; T_m <- 128L
X <- matrix(rnorm(n * p), n)
shape <- exp(-((seq_len(T_m) - T_m / 2)^2) / (2 * (T_m / 8)^2))
beta_M <- matrix(0, p, T_m)
beta_M[40L, ]  <-  1.5 * shape; beta_M[120L, ] <- -1.0 * shape
Y <- X %*% beta_M + matrix(rnorm(n * T_m, sd = 0.5), n)
run_three("failing_case", X, Y, c(40L, 120L), pos_list = list(seq_len(T_m)))
