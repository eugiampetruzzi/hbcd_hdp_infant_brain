source("analysis/common.R")

outcome_path <- function(subdir, filename) first_existing(
  file.path(HBCD_MD, "analysis_datasets", filename),
  file.path(HBCD_MD, "outcomes", subdir, filename))

spatial_wide <- function(base, outcome_file, modality, prefix = NULL, coverage = 0) {
  outcomes <- read_table(outcome_file)
  id_cols <- intersect(c("participant_id", "session_id", "source_run"), names(outcomes))
  candidates <- setdiff(names(outcomes), id_cols)
  if (!is.null(prefix)) candidates <- candidates[startsWith(candidates, prefix)]
  numeric <- candidates[vapply(outcomes[candidates], is.numeric, logical(1))]
  d <- merge(base, outcomes[, c("participant_id", numeric), drop = FALSE], by = "participant_id", all = FALSE)
  numeric <- numeric[vapply(d[numeric], function(x) mean(!is.na(x)) >= coverage, logical(1))]
  covs <- modality_spec[[modality]]$covariates
  rows <- lapply(numeric, function(outcome) {
    fit <- hc3_fit(d, outcome, "primary_hdp", covs)
    data.frame(modality = modality, outcome = outcome, t(fit), check.names = FALSE)
  })
  if (!length(rows)) stop("No eligible outcomes found in ", outcome_file)
  do.call(rbind, rows)
}
