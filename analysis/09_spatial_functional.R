#!/usr/bin/env Rscript
source("analysis/spatial_common.R")

base <- load_modality("functional")
path <- outcome_path("functional", "v02_myerslabonte_system_participation.tsv")
systems <- read_table(path)
required <- c("participant_id", "myerslabonte_system", "participation_coefficient")
if (!all(required %in% names(systems))) stop("Functional system file must be long format with: ", paste(required, collapse = ", "))
d <- merge(base, systems[, required], by = "participant_id", all = FALSE)
rows <- lapply(sort(unique(d$myerslabonte_system)), function(system) {
  ds <- d[d$myerslabonte_system == system, ]
  fit <- hc3_fit(ds, "participation_coefficient", "primary_hdp",
                 modality_spec$functional$covariates)
  data.frame(modality = "functional", outcome = system, t(fit), check.names = FALSE)
})
results <- do.call(rbind, rows)
results <- add_fdr(results)
write_table(results, file.path(RESULTS, "table_spatial_functional.tsv"))
save_session("09_spatial_functional")
