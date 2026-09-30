#!/usr/bin/env Rscript
source("analysis/common.R")

contrasts <- list(
  preeclampsia_vs_gestational_hypertension = function(d) {
    d <- d[d$primary_hdp == 1 & (d$preeclampsia_eclampsia == 1 |
      (d$preeclampsia_eclampsia == 0 & d$gestational_hypertension == 1)), ]
    d$contrast <- as.integer(d$preeclampsia_eclampsia == 1)
    d
  },
  history_vs_no_history_within_hdp = function(d) {
    d <- d[d$primary_hdp == 1 & !is.na(d$chronic_hypertension), ]
    d$contrast <- as.integer(d$chronic_hypertension == 1)
    d
  }
)

rows <- list()
for (contrast_name in names(contrasts)) {
  for (modality in names(modality_spec)) {
    d <- contrasts[[contrast_name]](load_modality(modality))
    result <- run_outcomes(d, modality, exposure = "contrast")
    result$contrast <- contrast_name
    result$n_reference <- sum(d$contrast == 0)
    result$n_exposed <- sum(d$contrast == 1)
    rows[[length(rows) + 1]] <- result
  }
}
results <- do.call(rbind, rows)
results <- add_fdr(results, results$contrast)
write_table(results, file.path(RESULTS, "table_hdp_presentations.tsv"))
save_session("04_hdp_presentations")
