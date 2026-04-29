## For each profvis stack sample, walk up from leaf and find the
## highest-level mfsusieR function — that gives proper attribution
## of .Call / apply / <GC> leaves to the mfsusieR caller.
suppressPackageStartupMessages(library(dplyr))

mf_funcs <- c(
  "mfsusie", "create_mf_individual", "mf_dwt", "mf_prior_scale_mixture",
  "compute_residuals.mf_individual",
  "compute_ser_statistics.mf_individual",
  "loglik.mf_individual",
  "calculate_posterior_moments.mf_individual",
  "update_fitted_values.mf_individual",
  "update_variance_components.mf_individual",
  "update_derived_quantities.mf_individual",
  "update_model_variance.mf_individual",
  "refresh_em_cache.mf_individual",
  "refresh_lbf_kl.mf_individual",
  "get_objective.mfsusie",
  "Eloglik.mf_individual",
  "compute_kl.mf_individual",
  "SER_posterior_e_loglik.mf_individual",
  "neg_loglik.mf_individual",
  "mf_per_outcome_bhat_shat",
  "mf_sigma2_per_position",
  "mf_get_ER2_per_position",
  "optimize_prior_variance.mf_individual",
  "mf_em_likelihood_per_scale",
  "mf_em_m_step_per_scale",
  "mixture_log_bf_per_scale",
  "mixture_log_bf_per_scale_johnson",
  "mixture_posterior_per_scale",
  "mixsqp",
  "verify.likelihood.matrix",
  "mixobj",
  "compute_marginal_bhat_shat",
  ## susieR backbone functions we want to credit
  "single_effect_regression",
  "susie_workhorse",
  "susie_get_cs",
  "susie_get_pip"
)

attribute <- function(rds_path) {
  pf <- readRDS(rds_path)
  ## Per stack sample, walk up (depth from leaf upward = decreasing depth).
  ## We want the LOWEST-depth mfsusieR-namespace function (closest to leaf).
  ## Strip namespace prefix (`pkg::fn` → `fn`) before matching.
  pf$label_bare <- sub("^[A-Za-z0-9_.]+::", "", pf$label)
  pf$is_mf      <- pf$label_bare %in% mf_funcs
  pf <- pf[order(pf$time, -pf$depth), ]   ## leaf first per time
  ## For each `time`, take the FIRST mfsusieR row (deepest mf func).
  attr_mf <- pf %>%
    filter(is_mf) %>%
    group_by(time) %>%
    slice(1) %>%
    ungroup() %>%
    select(time, mf_func = label_bare)

  ## Identify leaf label per time
  leaf <- pf %>%
    group_by(time) %>%
    slice_max(depth, n = 1L, with_ties = FALSE) %>%
    ungroup() %>%
    select(time, leaf_label = label)

  out <- leaf %>%
    left_join(attr_mf, by = "time") %>%
    mutate(mf_func = ifelse(is.na(mf_func), "<no_mf_in_stack>", mf_func))

  ## Self-time per (mf_func, leaf_label)
  interval_ms <- 5
  ag <- out %>%
    count(mf_func, leaf_label, name = "samples") %>%
    mutate(self_ms  = samples * interval_ms,
           pct_self = round(100 * samples / sum(samples), 2)) %>%
    arrange(desc(samples))

  ## Roll up to mf_func only (= total time charged to this mfsusieR caller)
  rolled <- out %>%
    count(mf_func, name = "samples") %>%
    mutate(self_ms  = samples * interval_ms,
           pct_self = round(100 * samples / sum(samples), 1)) %>%
    arrange(desc(samples))

  list(by_pair = ag, rolled = rolled,
       total_ms = nrow(out) * interval_ms,
       n_samples = nrow(out))
}

cat("\n================ per_outcome ================\n")
po <- attribute("/tmp/profvis_per_outcome.rds")
cat(sprintf("total %.2fs\n\n", po$total_ms / 1000))
cat("--- self time charged to mfsusieR caller ---\n")
print(head(po$rolled, 20), n = 20)
cat("\n--- top (mf_func, leaf) pairs ---\n")
print(head(po$by_pair, 25), n = 25)

cat("\n================ per_scale ================\n")
ps <- attribute("/tmp/profvis_per_scale.rds")
cat(sprintf("total %.2fs\n\n", ps$total_ms / 1000))
cat("--- self time charged to mfsusieR caller ---\n")
print(head(ps$rolled, 20), n = 20)
cat("\n--- top (mf_func, leaf) pairs ---\n")
print(head(ps$by_pair, 25), n = 25)
