#!/usr/bin/env Rscript
source("analysis/common.R")

person <- read_table(file.path(HBCD_MD, "covariates", "hdp_missingness_person_level.csv"))
person <- person[, c("participant_id", "mother_age", "race_eth", "education", "income")]
sets <- list(
  maternal_age_and_race_ethnicity = c("mother_age", "race_eth"),
  plus_education = c("mother_age", "race_eth", "education"),
  plus_household_income = c("mother_age", "race_eth", "education", "income")
)
rows <- list()
for (set_name in names(sets)) {
  for (modality in names(modality_spec)) {
    d <- load_modality(modality)
    overlap <- intersect(setdiff(names(person), "participant_id"), names(d))
    if (length(overlap)) d <- d[, setdiff(names(d), overlap), drop = FALSE]
    d <- merge_one(d, person, all.x = TRUE)
    covs <- c(modality_spec[[modality]]$covariates, sets[[set_name]])
    result <- run_outcomes(d, modality, covariates = covs)
    result$sensitivity <- set_name
    rows[[length(rows) + 1]] <- result
  }
}
results <- do.call(rbind, rows)
results <- add_fdr(results, results$sensitivity)
write_table(results, file.path(RESULTS, "table_covariate_sensitivities.tsv"))
save_session("06_covariate_sensitivity")
