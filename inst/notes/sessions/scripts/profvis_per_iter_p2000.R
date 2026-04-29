## Profvis instrumentation for mfsusieR per-iter hot path at p=2000.
suppressPackageStartupMessages({
  library(mfsusieR)
  library(profvis)
})

set.seed(20250428)
n  <- 84
p  <- 2000
M  <- 6
T_ <- 128

X <- matrix(rnorm(n * p), n, p)
true_idx <- c(11, 200, 1500)
beta <- matrix(0, p, T_)
for (j in true_idx) beta[j, ] <- rnorm(T_, sd = 0.4)
Y <- vector("list", M)
for (m in seq_len(M)) {
  E <- matrix(rnorm(n * T_, sd = 0.6), n, T_)
  Y[[m]] <- X %*% beta + E
}

run_one <- function(scope, max_iter = 4) {
  message("--- profile run: prior_variance_scope = ", scope, " ---")
  pv <- profvis::profvis({
    fit <- mfsusieR::mfsusie(
      X = X, Y = Y, L = 10,
      prior_variance_scope     = scope,
      residual_variance_scope  = scope,
      max_iter                 = max_iter,
      tol                      = 1e-4,
      verbose                  = FALSE
    )
    invisible(fit)
  }, interval = 0.005)
  out_path <- sprintf("/tmp/profvis_%s.rds", scope)
  saveRDS(pv$x$message$prof, out_path)
  cat(sprintf("saved %s rows = %d\n", out_path, length(pv$x$message$prof$ref$ref)))

  ## Wall time per iter via system.time + max_iter limit
  t1 <- system.time({
    fit <- mfsusieR::mfsusie(
      X = X, Y = Y, L = 10,
      prior_variance_scope     = scope,
      residual_variance_scope  = scope,
      max_iter                 = max_iter,
      tol                      = 1e-4,
      verbose                  = FALSE
    )
  })
  cat(sprintf("[%s] wall %.3fs over %d iters -> %.3fs/iter (converged=%s, niter=%d)\n",
              scope, t1[["elapsed"]], fit$niter, t1[["elapsed"]] / fit$niter,
              fit$converged, fit$niter))
  invisible(NULL)
}

run_one("per_outcome")
run_one("per_scale")
