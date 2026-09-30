#!/usr/bin/env Rscript
source("analysis/spatial_common.R")

base <- load_modality("structural")
files <- c(thickness = "v02_myerslabonte_thickness.tsv", sulcal_depth = "v02_myerslabonte_sulc.tsv")
rows <- lapply(names(files), function(metric) {
  x <- spatial_wide(base, outcome_path("structural", files[[metric]]), "structural", coverage = .80)
  x$metric <- metric
  x
})
results <- do.call(rbind, rows)
results <- add_fdr(results, results$metric)
write_table(results, file.path(RESULTS, "table_spatial_structural.tsv"))
save_session("08_spatial_structural")
