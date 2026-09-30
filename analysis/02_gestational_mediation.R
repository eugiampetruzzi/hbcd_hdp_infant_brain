#!/usr/bin/env Rscript
source("analysis/common.R")

set.seed(20260930)
n_boot <- as.integer(Sys.getenv("HBCD_BOOTSTRAP_REPS", "5000"))

mediation_fit <- function(data, modality, outcome, label) {
  covs <- modality_spec[[modality]]$covariates
  vars <- unique(c("primary_hdp", "ga_delivery", outcome, covs))
  d <- data[complete.cases(data[, vars, drop = FALSE]), vars, drop = FALSE]
  d[[outcome]] <- standardize(d[[outcome]])
  d$ga_delivery <- standardize(d$ga_delivery)
  d <- prepare_covariates(d, covs)
  mediator_formula <- reformulate(c("primary_hdp", covs), response = "ga_delivery")
  outcome_formula <- reformulate(c("primary_hdp", "ga_delivery", covs), response = outcome)
  a <- unname(coef(lm(mediator_formula, d))["primary_hdp"])
  direct_model <- lm(outcome_formula, d)
  b <- unname(coef(direct_model)["ga_delivery"])
  direct_hc3 <- hc3_fit(d, outcome, "primary_hdp", c("ga_delivery", covs))
  total_hc3 <- hc3_fit(d, outcome, "primary_hdp", covs)
  direct <- unname(direct_hc3["beta"])
  total <- unname(total_hc3["beta"])
  boot <- replicate(n_boot, {
    db <- d[sample.int(nrow(d), replace = TRUE), , drop = FALSE]
    ab <- unname(coef(lm(mediator_formula, db))["primary_hdp"])
    bb <- unname(coef(lm(outcome_formula, db))["ga_delivery"])
    ab * bb
  })
  indirect <- a * b
  p_boot <- 2 * min((sum(boot <= 0) + 1) / (n_boot + 1), (sum(boot >= 0) + 1) / (n_boot + 1))
  data.frame(modality = modality, outcome = label, n = nrow(d), total = total,
             a_hdp_to_ga = a, b_ga_to_outcome = b, indirect = indirect,
             indirect_ci_low = unname(quantile(boot, .025)),
             indirect_ci_high = unname(quantile(boot, .975)),
             indirect_p = min(1, p_boot), direct = direct,
             direct_se = unname(direct_hc3["se"]), direct_p = unname(direct_hc3["p"]))
}

rows <- list()
for (modality in names(modality_spec)) {
  data <- load_modality(modality)
  for (label in names(modality_spec[[modality]]$outcomes)) {
    rows[[length(rows) + 1]] <- mediation_fit(data, modality,
      modality_spec[[modality]]$outcomes[[label]], label)
  }
}
results <- do.call(rbind, rows)
results$indirect_q <- p.adjust(results$indirect_p, method = "BH")
write_table(results, file.path(RESULTS, "table_gestational_mediation.tsv"))
save_session("02_gestational_mediation")
