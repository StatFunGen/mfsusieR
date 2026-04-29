## Parse the profvis prof dumps from profvis_per_iter_p2000.R and
## report a self-time decomposition per top-of-stack function.

library(dplyr)

parse_profvis <- function(rds_path) {
  pf <- readRDS(rds_path)
  ## profvis profile is a data.frame: time, depth, label, filename,
  ## linenum, ... Each `time` index is one stack-sample; rows at that
  ## time are the call stack from outermost (depth = max) to leaf
  ## (depth = 1). We want self-time per leaf function.

  ## "self time" = number of samples whose maximum depth at that time
  ## index is on the leaf row. Group by `time`, take the deepest row.
  leaf <- pf %>%
    group_by(time) %>%
    slice_max(depth, n = 1L, with_ties = FALSE) %>%
    ungroup()

  ## Sample interval is fixed at 5ms in the run script.
  interval_ms <- 5
  total_ms    <- nrow(leaf) * interval_ms
  cat(sprintf("[%s] total self-time = %.2fs across %d samples\n",
              basename(rds_path), total_ms / 1000, nrow(leaf)))

  leaf %>%
    mutate(label = ifelse(is.na(label) | label == "", "<unknown>", label),
           where = ifelse(is.na(filename) | filename == "",
                          "<base>",
                          paste0(basename(filename), ":", linenum))) %>%
    count(label, where, name = "samples") %>%
    arrange(desc(samples)) %>%
    mutate(self_ms  = samples * interval_ms,
           pct_self = round(100 * samples / sum(samples), 1)) %>%
    select(label, where, samples, self_ms, pct_self)
}

cat("================ per_outcome (default) ================\n")
po <- parse_profvis("/tmp/profvis_per_outcome.rds")
print(head(po, 30), n = 30)

cat("\n================ per_scale ================\n")
ps <- parse_profvis("/tmp/profvis_per_scale.rds")
print(head(ps, 30), n = 30)

## Aggregate to function name only (drop file:line) for top-N readability.
roll <- function(df) {
  df %>%
    group_by(label) %>%
    summarize(samples  = sum(samples),
              self_ms  = sum(self_ms),
              pct_self = round(sum(pct_self), 1),
              .groups  = "drop") %>%
    arrange(desc(samples))
}

cat("\n--- per_outcome top-15 by function (rolled up) ---\n")
print(head(roll(po), 15), n = 15)

cat("\n--- per_scale top-15 by function (rolled up) ---\n")
print(head(roll(ps), 15), n = 15)

## Bucketed attribution: hand-classify the leaves into the 4
## investigation buckets (a/b/c/d).
classify <- function(label, where) {
  is_mixsqp_inner   <- grepl("mixsqp$|^mixsqp::", label) |
                       grepl("mixsqp.*\\.cpp", where)
  is_em_L_builder   <- label %in% c("mf_em_likelihood_per_scale",
                                    "mf_em_m_step_per_scale")
  is_mixture_dual   <- label %in% c("mixture_log_bf_per_scale",
                                    "mixture_posterior_per_scale")
  is_residual       <- label %in% c("compute_residuals.mf_individual",
                                    "compute_ser_statistics.mf_individual",
                                    "update_fitted_values.mf_individual",
                                    "update_derived_quantities.mf_individual",
                                    "calculate_posterior_moments.mf_individual",
                                    "loglik.mf_individual")
  is_em_cache_build <- label %in% c("refresh_em_cache.mf_individual",
                                    "update_variance_components.mf_individual")
  is_lbf_refresh    <- label %in% c("refresh_lbf_kl.mf_individual",
                                    "get_objective.mfsusie")
  is_cs_purity      <- grepl("get_pip|get_cs|get_purity", label)
  case_when_first <- function(...) {
    args <- list(...)
    out  <- rep("other", length(args[[1]]))
    for (k in seq(1, length(args) - 1L, by = 2L)) {
      cond <- args[[k]]
      lab  <- args[[k + 1L]]
      out[cond & out == "other"] <- lab
    }
    out
  }
  case_when_first(
    is_mixsqp_inner,   "mixsqp_inner",
    is_em_L_builder,   "em_L_builder",
    is_mixture_dual,   "mixture_log_bf+posterior",
    is_residual,       "residual+ser",
    is_em_cache_build, "em_cache_rebuild",
    is_lbf_refresh,    "lbf_kl_refresh",
    is_cs_purity,      "cs_purity"
  )
}

bucket_table <- function(df) {
  df %>%
    mutate(bucket = classify(label, where)) %>%
    group_by(bucket) %>%
    summarize(self_ms  = sum(self_ms),
              pct_self = round(sum(pct_self), 1),
              .groups  = "drop") %>%
    arrange(desc(self_ms))
}

cat("\n--- per_outcome bucket attribution ---\n")
print(bucket_table(po))

cat("\n--- per_scale bucket attribution ---\n")
print(bucket_table(ps))
