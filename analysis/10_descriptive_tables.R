#!/usr/bin/env Rscript
source("analysis/common.R")

continuous_summary <- function(x, variable, sample = "any_MRI") {
  do.call(rbind, lapply(c(0, 1), function(group) {
    values <- x[x$primary_hdp == group, variable]
    data.frame(sample = sample, variable = variable,
      group = ifelse(group == 1, "HDP", "control"), n = sum(!is.na(values)),
      mean = mean(values, na.rm = TRUE), sd = sd(values, na.rm = TRUE),
      median = median(values, na.rm = TRUE), q1 = unname(quantile(values, .25, na.rm = TRUE)),
      q3 = unname(quantile(values, .75, na.rm = TRUE)))
  }))
}

categorical_summary <- function(x, variable, sample = "any_MRI") {
  levels <- sort(unique(as.character(x[[variable]][!is.na(x[[variable]])])))
  do.call(rbind, lapply(c(0, 1), function(group) {
    denominator <- sum(x$primary_hdp == group & !is.na(x[[variable]]))
    do.call(rbind, lapply(levels, function(level) {
      n <- sum(x$primary_hdp == group & as.character(x[[variable]]) == level, na.rm = TRUE)
      data.frame(sample = sample, variable = variable, level = level,
        group = ifelse(group == 1, "HDP", "control"), n = n,
        denominator = denominator, percent = 100 * n / denominator)
    }))
  }))
}

modalities <- lapply(names(modality_spec), load_modality)
names(modalities) <- names(modality_spec)
membership <- unique(do.call(rbind, lapply(modalities, function(d) d[, c("participant_id", "primary_hdp")])))
if (nrow(membership) != 1403 || sum(membership$primary_hdp) != 393) {
  stop("Any-MRI sample mismatch: ", nrow(membership), " total; ", sum(membership$primary_hdp), " HDP")
}
person <- read_table(file.path(HBCD_MD, "covariates", "hdp_missingness_person_level.csv"))
person <- person[!duplicated(person$participant_id), ]
sample <- merge_one(membership, person, all.x = TRUE)

continuous <- do.call(rbind, lapply(intersect(c("mother_age", "ga_delivery"), names(sample)),
                                           function(v) continuous_summary(sample, v)))
categorical <- do.call(rbind, lapply(intersect(c("child_sex", "race_eth", "education", "income"), names(sample)),
                                    function(v) categorical_summary(sample, v)))
write_table(continuous, file.path(RESULTS, "table_sample_characteristics_continuous.tsv"))
write_table(categorical, file.path(RESULTS, "table_sample_characteristics_categorical.tsv"))

complications <- read_table(file.path(HBCD_MD, "covariates", "v02_pregnancy_complications.tsv"))
age_source <- do.call(rbind, lapply(modalities, function(d) {
  d[, c("participant_id", "candidate_age_weeks")]
}))
age_source <- aggregate(candidate_age_weeks ~ participant_id, age_source,
                        function(x) x[which(!is.na(x))[1]])
core <- sample[, c("participant_id", "primary_hdp", "ga_delivery")]
core <- merge_one(core, age_source[, c("participant_id", "candidate_age_weeks")], all.x = TRUE)
core <- merge_one(core, complications[, c("participant_id", "other_complication")], all.x = TRUE)
names(core)[names(core) == "candidate_age_weeks"] <- "chronological_age_at_scan"

global <- core
for (modality in names(modalities)) {
  d <- modalities[[modality]]
  cols <- unname(modality_spec[[modality]]$outcomes)
  part <- d[, c("participant_id", cols), drop = FALSE]
  global <- merge(global, part, by = "participant_id", all.x = TRUE)
}
labels <- c("primary_hdp", "ga_delivery", "chronological_age_at_scan",
            "global_median_md", "global_median_fa", "global_median_mkt",
            "global_area_weighted_thickness", "global_area_weighted_absolute_sulc",
            "modularity", "mean_participation_coefficient", "other_complication")
rows <- list()
for (i in seq_along(labels)) for (j in i:length(labels)) {
  x <- global[[labels[i]]]
  y <- global[[labels[j]]]
  keep <- is.finite(x) & is.finite(y)
  rows[[length(rows) + 1]] <- data.frame(
    outcome_1 = labels[i], outcome_2 = labels[j], n = sum(keep),
    r = if (sum(keep) > 2) cor(x[keep], y[keep]) else NA_real_)
}
correlations <- do.call(rbind, rows)
if (min(correlations$n) != 718 || max(correlations$n) != 1403) {
  stop("Table 2 pairwise N range mismatch: ", min(correlations$n), " to ", max(correlations$n))
}
write_table(correlations, file.path(RESULTS, "table_outcome_correlations.tsv"))
save_session("10_descriptive_tables")
