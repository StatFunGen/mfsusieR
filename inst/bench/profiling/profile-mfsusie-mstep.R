# Profile a representative `mfsusie()` fit to identify the M-step
# hot paths. Runs once with profvis; saves the flamegraph as HTML
# and prints a brief summary of the top hot paths.
#
# Estimated runtime: < 60 seconds on a single core.
#
# Usage:
#   Rscript bench/profiling/profile-mfsusie-mstep.R [TAG]
#
# `TAG` defaults to "baseline"; pass any short string to label the
# output (e.g., "cache", "vectorize", "subset", "cpp").

suppressPackageStartupMessages({
  library(profvis)
  library(mfsusieR)
  library(htmlwidgets)
})

args <- commandArgs(trailingOnly = TRUE)
tag  <- if (length(args) >= 1L) args[[1L]] else "baseline"

# A medium-size synthetic problem in the practical_dataset
# regime: M = 4 outcomes, n = 500 samples, p = 200 SNPs,
# T_basis = c(128, 128, 64, 1). Two scalar + two functional
# would test more code paths but four functional is closer to
# what the user said matters.
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

cat(sprintf(
  "fixture: n=%d p=%d M=%d T_basis=%s\n",
  n, p, length(T_per),
  paste(T_per, collapse = ",")))

# The actual timed run.
t0 <- proc.time()
prof <- profvis({
  fit <- mfsusie(X, Y, L = 10, max_iter = 30, verbose = FALSE)
})
elapsed <- (proc.time() - t0)[["elapsed"]]

cat(sprintf("wall-clock (incl profvis overhead): %.2f s\n", elapsed))

flamegraph_dir <- "/home/gw/GIT/mfsusieR/bench/profiling/flamegraphs"
dir.create(flamegraph_dir, showWarnings = FALSE, recursive = TRUE)
out_html <- file.path(flamegraph_dir, sprintf("%s.html", tag))
saveWidget(prof, out_html, selfcontained = TRUE)
cat(sprintf("flamegraph: %s\n", out_html))

# Print a tabular summary of the top hot paths so we can see them
# in the console without opening the HTML.
ev <- prof$x$message$prof
calls <- table(unlist(strsplit(
  vapply(split(ev$label, ev$time), function(stack) {
    paste(stack, collapse = " > ")
  }, character(1L)),
  " > "
)))
calls_top <- sort(calls, decreasing = TRUE)
cat("\nTop 20 functions by sample count (lower-is-callee):\n")
print(utils::head(calls_top, 20L))
cat(sprintf("\ntotal samples: %d\n", sum(calls_top)))
