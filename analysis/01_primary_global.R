#!/usr/bin/env Rscript
source("analysis/common.R")

models <- list(
  diffusion = list(
    outcomes = c(MD = "global_median_md", FA = "global_median_fa",
                 MKT = "global_median_mkt"),
    covariates = c("candidate_age_weeks", "sex_m", "site",
                   "raw_neighbor_corr")),
  structural = list(
    outcomes = c(cortical_thickness = "global_area_weighted_thickness",
                 sulcal_depth = "global_area_weighted_absolute_sulc"),
    covariates = c("candidate_age_weeks", "child_sex", "site",
                   "brainswipes_summary_qc")),
  functional = list(
    outcomes = c(modularity = "modularity",
                 participation_coefficient = "mean_participation_coefficient"),
    covariates = c("candidate_age_weeks", "child_sex", "site",
                   "mean_fd_post_censoring", "retained_duration_minutes"))
)

rows <- list()
for (modality in names(models)) {
  data <- load_modality(modality)
  if (nrow(data) != expected_counts[[modality]]) {
    stop(modality, " count mismatch: ", nrow(data))
  }
  for (label in names(models[[modality]]$outcomes)) {
    outcome <- models[[modality]]$outcomes[[label]]
    estimate <- hc3_fit(data, outcome, exposure = "primary_hdp",
                        covariates = models[[modality]]$covariates)
    rows[[length(rows) + 1]] <- data.frame(
      modality = modality, outcome = label, t(estimate), check.names = FALSE)
  }
}
results <- do.call(rbind, rows)
results <- add_fdr(results)
write_table(results, file.path(RESULTS, "table_primary_global.tsv"))
save_session("01_primary_global")
