#!/usr/bin/env Rscript
source("analysis/common.R")

expected <- c(diffusion = 1114, structural = 1007, functional = 796)
results <- do.call(rbind, lapply(names(modality_spec), function(modality) {
  data <- load_modality(modality)
  data <- data[data$primary_hdp == 0 | data$other_complication == 0, , drop = FALSE]
  if (nrow(data) != expected[[modality]]) stop(modality, " sensitivity count mismatch: ", nrow(data))
  run_outcomes(data, modality)
}))
results <- add_fdr(results)
write_table(results, file.path(RESULTS, "table_no_complications.tsv"))
save_session("03_no_complications")
