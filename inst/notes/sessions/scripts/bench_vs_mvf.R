## Apples-to-apples per-iteration bench: mfsusieR vs mvf.susie.alpha
## (multfsusie). Both run prior_variance_scope = per_outcome on identical
## data. Reports s/iter and a Rprof-attributed breakdown of where each
## spends its time.

suppressPackageStartupMessages({
  library(mfsusieR)
  library(mvf.susie.alpha)
})

set.seed(1)

bench_one <- function(n, p, M, T_, L = 10, max_iter = 50) {
  cat(sprintf("\n=========== n=%d p=%d M=%d T=%d L=%d ===========\n",
              n, p, M, T_, L))
  X <- matrix(rnorm(n * p), n, p)
  true_idx <- c(11, 200, min(p, 1500))
  beta <- matrix(0, p, T_)
  for (j in true_idx) beta[j, ] <- rnorm(T_, sd = 0.4)
  Y <- vector("list", M)
  for (m in seq_len(M)) {
    E <- matrix(rnorm(n * T_, sd = 0.6), n, T_)
    Y[[m]] <- X %*% beta + E
  }
  ## mvf wants Y as list(Y_f = list(M matrices)). Our mfsusieR takes
  ## Y as a flat list of M matrices.
  Y_mvf <- list(Y_f = Y)

  ## ---- mfsusieR (per_outcome, PIP convergence) ----
  message("> mfsusieR run")
  Rprof("/tmp/prof_mfsusie.out", interval = 0.01, line.profiling = FALSE)
  t_mf <- system.time({
    fit_mf <- mfsusieR::mfsusie(
      X = X, Y = Y, L = L,
      prior_variance_scope    = "per_outcome",
      residual_variance_scope = "per_outcome",
      max_iter                = max_iter,
      tol                     = 1e-4,
      verbose                 = FALSE
    )
  })
  Rprof(NULL)

  ## ---- mvf.susie.alpha (prior = mixture_normal = per_outcome) ----
  message("> mvf run")
  Rprof("/tmp/prof_mvf.out", interval = 0.01, line.profiling = FALSE)
  t_mvf <- system.time({
    fit_mvf <- mvf.susie.alpha::multfsusie(
      Y = Y_mvf, X = X, L = L, pos = NULL,
      prior            = "mixture_normal",
      maxit            = max_iter,
      tol              = 1e-4,
      max_SNP_EM       = p,
      cal_obj          = FALSE,    ## skip ELBO computation
      greedy           = FALSE,    ## match our L_greedy = NULL default
      backfit          = FALSE,
      post_processing  = "none",
      verbose          = FALSE
    )
  })
  Rprof(NULL)

  ## --- niter ---
  niter_mf  <- fit_mf$niter
  niter_mvf <- fit_mvf$iter %||% (length(fit_mvf$alpha_hist) - 1L)
  if (is.null(niter_mvf) || niter_mvf < 1L) niter_mvf <- 1L

  cat(sprintf("\n[mfsusieR]  wall %.3fs niter %d -> %.3fs/iter\n",
              t_mf[["elapsed"]],  niter_mf,  t_mf[["elapsed"]]  / niter_mf))
  cat(sprintf("[mvf]       wall %.3fs niter %d -> %.3fs/iter\n",
              t_mvf[["elapsed"]], niter_mvf, t_mvf[["elapsed"]] / niter_mvf))
  ratio <- (t_mf[["elapsed"]] / niter_mf) / (t_mvf[["elapsed"]] / niter_mvf)
  cat(sprintf("[gap per iter] mfsusieR / mvf = %.2fx %s\n",
              ratio, if (ratio > 1) "(slower)" else "(faster)"))

  list(p = p, n = n, M = M, T = T_, L = L,
       mf_wall  = t_mf[["elapsed"]],  mf_niter  = niter_mf,
       mvf_wall = t_mvf[["elapsed"]], mvf_niter = niter_mvf,
       mf_per_iter  = t_mf[["elapsed"]]  / niter_mf,
       mvf_per_iter = t_mvf[["elapsed"]] / niter_mvf,
       gap_ratio    = ratio)
}

results <- list()

## Single moderate-p case the user cares about; p=500 is small enough that
## both packages run in <2 min and large enough to surface the per-iter gap.
results[[1]] <- bench_one(n = 84, p = 500,  M = 6, T_ = 128, L = 10, max_iter = 30)
results[[2]] <- bench_one(n = 84, p = 1000, M = 6, T_ = 128, L = 10, max_iter = 30)

cat("\n========================= summary =========================\n")
cat(sprintf("%-6s %-6s %-10s %-10s %-12s %-12s %-8s\n",
            "p", "iters", "mf_per_it", "mvf_per_it", "mf_total", "mvf_total", "gap"))
for (r in results) {
  cat(sprintf("%-6d mf:%-3d mvf:%-3d %-10.3f %-10.3f %-12.3f %-12.3f %.2fx\n",
              r$p, r$mf_niter, r$mvf_niter,
              r$mf_per_iter, r$mvf_per_iter,
              r$mf_wall, r$mvf_wall, r$gap_ratio))
}

cat("\n=== mfsusieR Rprof summary (last run) ===\n")
print(head(summaryRprof("/tmp/prof_mfsusie.out")$by.self, 20))

cat("\n=== mvf Rprof summary (last run) ===\n")
print(head(summaryRprof("/tmp/prof_mvf.out")$by.self, 20))
