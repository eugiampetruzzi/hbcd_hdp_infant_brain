#!/usr/bin/env Rscript
source("analysis/spatial_common.R")

base <- load_modality("diffusion")
files <- c(
  FA = "v02_bundle_fa_masked_median.tsv",
  MD = "v02_bundle_md_masked_median.tsv",
  MKT = "v02_bundle_mkt_masked_median.tsv")
rows <- lapply(names(files), function(metric) {
  x <- spatial_wide(base, outcome_path("diffusion", files[[metric]]), "diffusion",
                    prefix = paste0(tolower(metric), "_"))
  x$metric <- metric
  x
})
results <- do.call(rbind, rows)
results <- add_fdr(results, results$metric)
write_table(results, file.path(RESULTS, "table_spatial_diffusion.tsv"))
save_session("07_spatial_diffusion")
