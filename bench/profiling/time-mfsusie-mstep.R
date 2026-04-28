# Wall-clock measurement of `mfsusie()` without profvis overhead.
# Same fixture as profile-mfsusie-mstep.R; runs N times and prints
# median + spread.

suppressPackageStartupMessages({ library(mfsusieR) })

args <- commandArgs(trailingOnly = TRUE)
n_runs <- if (length(args) >= 1L) as.integer(args[[1L]]) else 3L

set.seed(1)
n  <- 300L
p  <- 150L
X  <- matrix(rnorm(n * p), nrow = n)
T_per <- c(128L, 64L)
beta_idx <- c(40L, 100L)
Y  <- lapply(T_per, function(T_m) {
  shape <- exp(-((seq_len(T_m) - T_m / 2)^2) / (2 * (T_m / 8)^2))
  beta  <- matrix(0, p, T_m)
  beta[beta_idx[1L], ] <-  1.5 * shape
  beta[beta_idx[2L], ] <- -1.0 * shape
  X %*% beta + matrix(rnorm(n * T_m, sd = 0.5), n)
})

# Warm-up.
invisible(mfsusie(X, Y, L = 3, max_iter = 3, verbose = FALSE))

elapsed <- numeric(n_runs)
niter   <- integer(n_runs)
for (i in seq_len(n_runs)) {
  t0   <- proc.time()
  fit  <- mfsusie(X, Y, L = 10, max_iter = 30, verbose = FALSE)
  elapsed[i] <- (proc.time() - t0)[["elapsed"]]
  niter[i]   <- fit$niter
}

cat(sprintf("runs=%d niter=%s elapsed (s) min=%.2f median=%.2f max=%.2f\n",
            n_runs,
            paste(niter, collapse = ","),
            min(elapsed), median(elapsed), max(elapsed)))
