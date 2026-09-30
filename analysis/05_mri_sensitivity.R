#!/usr/bin/env Rscript
source("analysis/common.R")

scanner <- read_table(file.path(HBCD_MD, "covariates", "v02_scanner_metadata.tsv"))
manufacturer_fields <- c(diffusion = "dwi_scanner_manufacturer",
                         structural = "t2w_scanner_manufacturer",
                         functional = "bold_scanner_manufacturer")
rows <- list()
for (modality in names(modality_spec)) {
  d <- load_modality(modality)
  field <- manufacturer_fields[[modality]]
  d <- merge_one(d, scanner[, c("participant_id", field)], all.x = TRUE)
  covs <- setdiff(modality_spec[[modality]]$covariates, modality_spec[[modality]]$site)
  covs <- c(covs, field)
  result <- run_outcomes(d, modality, covariates = covs)
  result$sensitivity <- "scanner_manufacturer_instead_of_site"
  rows[[length(rows) + 1]] <- result
}

functional <- load_modality("functional")
for (threshold in c(7.5, 10)) {
  d <- functional[functional$retained_duration_minutes >= threshold, ]
  result <- run_outcomes(d, "functional")
  result$sensitivity <- paste0("functional_retained_minutes_", threshold)
  rows[[length(rows) + 1]] <- result
}
results <- do.call(rbind, rows)
results <- add_fdr(results, results$sensitivity)
write_table(results, file.path(RESULTS, "table_mri_sensitivities.tsv"))
save_session("05_mri_sensitivity")
