## Compare mfsusieR vs mvf.susie.alpha on identical data: do they land
## at the same PIPs / signal locations / sigma2 estimates?
suppressPackageStartupMessages({
  library(mfsusieR)
  library(mvf.susie.alpha)
})

set.seed(1)
n <- 84; p <- 500; M <- 6; T_ <- 128
X <- matrix(rnorm(n*p), n, p)
true_idx <- c(11, 200, 450)
beta <- matrix(0, p, T_)
for (j in true_idx) beta[j, ] <- rnorm(T_, sd=0.4)
Y <- vector("list", M)
for (m in seq_len(M)) Y[[m]] <- X %*% beta + matrix(rnorm(n*T_, sd=0.6), n, T_)

cat(sprintf("Data: n=%d p=%d M=%d T=%d, true variables: %s\n\n",
            n, p, M, T_, paste(true_idx, collapse=",")))

## ---- mfsusieR ----
fit_mf <- mfsusieR::mfsusie(
  X = X, Y = Y, L = 10,
  prior_variance_scope    = "per_outcome",
  residual_variance_scope = "per_outcome",
  max_iter = 30, tol = 1e-4, verbose = FALSE
)

## ---- mvf ----
fit_mvf <- mvf.susie.alpha::multfsusie(
  Y = list(Y_f = Y), X = X, L = 10, prior = "mixture_normal",
  maxit = 30, tol = 1e-4, max_SNP_EM = p, cal_obj = FALSE,
  greedy = FALSE, backfit = FALSE, post_processing = "none",
  verbose = FALSE
)

cat("=== niter ===\n")
cat(sprintf("mfsusieR: %d  mvf: %d\n\n", fit_mf$niter, fit_mvf$iter))

cat("=== PIPs at the true variables ===\n")
cat(sprintf("var       mfsusieR     mvf\n"))
for (j in true_idx) {
  cat(sprintf("%4d      %.4f       %.4f\n", j, fit_mf$pip[j], fit_mvf$pip[j]))
}

cat("\n=== top-5 PIPs each ===\n")
top_mf  <- order(fit_mf$pip,  decreasing = TRUE)[1:5]
top_mvf <- order(fit_mvf$pip, decreasing = TRUE)[1:5]
cat("mfsusieR top-5 (var, pip): ",
    paste(sprintf("(%d, %.3f)", top_mf,  fit_mf$pip[top_mf]),  collapse=" "), "\n")
cat("mvf      top-5 (var, pip): ",
    paste(sprintf("(%d, %.3f)", top_mvf, fit_mvf$pip[top_mvf]), collapse=" "), "\n")

cat("\n=== whole-vector PIP correlation ===\n")
cat(sprintf("Pearson cor(mfsusieR$pip, mvf$pip) = %.4f\n",
            cor(fit_mf$pip, fit_mvf$pip)))
cat(sprintf("max |diff|                          = %.4f\n",
            max(abs(fit_mf$pip - fit_mvf$pip))))

cat("\n=== sigma2 (residual variance per outcome) ===\n")
sig_mf  <- vapply(fit_mf$sigma2,  function(s) mean(s), numeric(1))
sig_mvf <- if (is.list(fit_mvf$sigma2)) vapply(fit_mvf$sigma2, mean, numeric(1)) else fit_mvf$sigma2
cat("mfsusieR$sigma2 (mean per outcome): ", round(sig_mf, 3),  "\n")
cat("mvf$sigma2     (mean per outcome): ", round(sig_mvf, 3), "\n")

cat("\n=== top per-effect alpha[1, ] argmax ===\n")
cat(sprintf("mfsusieR: %d (alpha=%.3f)\n",
            which.max(fit_mf$alpha[1,]),  max(fit_mf$alpha[1,])))
if (!is.null(fit_mvf$alpha)) {
  cat(sprintf("mvf:      %d (alpha=%.3f)\n",
              which.max(fit_mvf$alpha[1,]), max(fit_mvf$alpha[1,])))
}

cat("\n=== credible sets ===\n")
cat("mfsusieR$sets$cs:\n")
str(fit_mf$sets$cs)
cat("mvf$cs:\n")
str(fit_mvf$cs)
