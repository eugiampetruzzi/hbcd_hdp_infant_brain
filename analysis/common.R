options(stringsAsFactors = FALSE)

read_env <- function(path = ".env") {
  if (!file.exists(path)) return(invisible(NULL))
  lines <- readLines(path, warn = FALSE)
  lines <- lines[grepl("^[A-Za-z_][A-Za-z0-9_]*=", lines)]
  for (line in lines) {
    pair <- strsplit(line, "=", fixed = TRUE)[[1]]
    key <- pair[1]
    value <- paste(pair[-1], collapse = "=")
    if (!nzchar(Sys.getenv(key))) Sys.setenv(structure(value, names = key))
  }
}

read_env()
required_r <- "4.5.1"
if (as.character(getRversion()) != required_r && Sys.getenv("HBCD_ALLOW_R_PATCH_MISMATCH") != "1") {
  stop("This analysis is pinned to R ", required_r, ". Set HBCD_ALLOW_R_PATCH_MISMATCH=1 only for development checks.")
}

HBCD_ROOT <- Sys.getenv("HBCD_ROOT", "/projects/iang/eanlp/eu_hbcd-release-2.1")
HBCD_MD <- Sys.getenv("HBCD_MANUSCRIPT_DATA", file.path(HBCD_ROOT, "imaging", "manuscript_data"))
HBCD_OUTPUT <- Sys.getenv("HBCD_OUTPUT", file.path(HBCD_MD, "reproducible_hdp_results"))
RESULTS <- file.path(HBCD_OUTPUT, "results")
dir.create(RESULTS, recursive = TRUE, showWarnings = FALSE)

read_table <- function(path) {
  if (!file.exists(path)) stop("Missing input: ", path)
  sep <- if (grepl("\\.csv$", path, ignore.case = TRUE)) "," else "\t"
  read.delim(path, sep = sep, check.names = FALSE, na.strings = c("", "NA", "NaN"))
}

write_table <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.table(x, path, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
}

first_existing <- function(...) {
  paths <- unlist(list(...))
  hit <- paths[file.exists(paths)]
  if (!length(hit)) stop("None of the candidate inputs exists:\n", paste(paths, collapse = "\n"))
  hit[1]
}

merge_one <- function(x, y, by = "participant_id", all.x = FALSE) {
  if (anyDuplicated(x[[by]]) || anyDuplicated(y[[by]])) stop("Non-unique merge key: ", by)
  merge(x, y, by = by, all.x = all.x, sort = FALSE)
}

restrict_controls <- function(data) {
  comp <- read_table(file.path(HBCD_MD, "covariates", "v02_pregnancy_complications.tsv"))
  keep_cols <- c("participant_id", "other_complication", "none_of_these",
                 "gestational_hypertension", "preeclampsia_eclampsia")
  missing <- setdiff(keep_cols, names(comp))
  if (length(missing)) stop("Missing complication fields: ", paste(missing, collapse = ", "))
  overlap <- intersect(setdiff(keep_cols, "participant_id"), names(data))
  if (length(overlap)) data <- data[, setdiff(names(data), overlap), drop = FALSE]
  data <- merge_one(data, comp[, keep_cols], all.x = FALSE)
  data[data$primary_hdp == 1 | (data$primary_hdp == 0 & data$none_of_these == 1), , drop = FALSE]
}

canonicalize_site <- function(data) {
  person <- read_table(file.path(HBCD_MD, "covariates", "hdp_missingness_person_level.csv"))
  person <- person[!duplicated(person$participant_id), c("participant_id", "site")]
  data <- data[, setdiff(names(data), c("site", "recruitment_site")), drop = FALSE]
  merge_one(data, person, all.x = TRUE)
}

load_modality <- function(modality) {
  generated <- file.path(HBCD_OUTPUT, "analysis_datasets", paste0("v02_", modality, "_primary.tsv"))
  if (file.exists(generated)) return(canonicalize_site(read_table(generated)))
  if (modality == "diffusion") {
    cov_path <- first_existing(
      file.path(HBCD_MD, "analysis_datasets", "v02_dwi_primary_covariates.tsv"),
      file.path(HBCD_MD, "covariates", "v02_dwi_primary_covariates.tsv"))
    x <- read_table(cov_path)
    x <- x[as.logical(x$primary_model_complete), ]
    fa_md <- read_table(file.path(HBCD_MD, "outcomes", "diffusion", "v02_global_cerebral_wm_fa_md.tsv"))
    mkt_path <- file.path(HBCD_MD, "outcomes", "diffusion", "v02_global_cerebral_wm_mkt.tsv")
    if (!("global_median_mkt" %in% names(x)) && file.exists(mkt_path)) {
      mkt <- read_table(mkt_path)
      x <- merge_one(x, mkt[, c("participant_id", "global_median_mkt")], all.x = TRUE)
    }
    x <- merge_one(x, fa_md[, c("participant_id", "global_median_md", "global_median_fa")], all.x = FALSE)
    return(canonicalize_site(restrict_controls(x)))
  }
  path <- switch(modality,
    structural = file.path(HBCD_MD, "analysis_datasets", "v02_global_cortical_morphometry_model_data.tsv"),
    functional = file.path(HBCD_MD, "analysis_datasets", "v02_global_functional_hdp_data.tsv"),
    stop("Unknown modality: ", modality))
  x <- restrict_controls(read_table(path))
  if (modality == "functional") x <- x[x$retained_duration_minutes >= 5, ]
  canonicalize_site(x)
}

modality_spec <- list(
  diffusion = list(
    outcomes = c(MD = "global_median_md", FA = "global_median_fa", MKT = "global_median_mkt"),
    covariates = c("candidate_age_weeks", "sex_m", "site", "raw_neighbor_corr"), site = "site"),
  structural = list(
    outcomes = c(cortical_thickness = "global_area_weighted_thickness", sulcal_depth = "global_area_weighted_absolute_sulc"),
    covariates = c("candidate_age_weeks", "child_sex", "site", "brainswipes_summary_qc"), site = "site"),
  functional = list(
    outcomes = c(modularity = "modularity", participation_coefficient = "mean_participation_coefficient"),
    covariates = c("candidate_age_weeks", "child_sex", "site", "mean_fd_post_censoring", "retained_duration_minutes"), site = "site")
)

standardize <- function(x) as.numeric(scale(x))

prepare_covariates <- function(d, covariates) {
  character_covs <- covariates[vapply(d[covariates], function(x) is.character(x) || is.factor(x), logical(1))]
  for (v in character_covs) {
    d[[v]] <- as.character(d[[v]])
    if (v == "site") {
      counts <- table(d[[v]])
      common <- sort(names(counts[counts > 3]))
      if (!length(common)) stop("No recruitment site has more than three observations")
      d[[v]][d[[v]] %in% names(counts[counts <= 3])] <- common[1]
    }
    d[[v]] <- factor(d[[v]])
  }
  d
}

hc3_fit <- function(data, outcome, exposure = "primary_hdp", covariates, standardize_outcome = TRUE) {
  variables <- unique(c(outcome, exposure, covariates))
  d <- data[complete.cases(data[, variables, drop = FALSE]), variables, drop = FALSE]
  d$.outcome <- if (standardize_outcome) standardize(d[[outcome]]) else d[[outcome]]
  d$.exposure <- d[[exposure]]
  d <- prepare_covariates(d, covariates)
  formula <- reformulate(c(".exposure", covariates), response = ".outcome")
  model <- lm(formula, data = d, x = TRUE, y = TRUE)
  X <- model$x
  residual <- residuals(model)
  leverage <- hatvalues(model)
  bread <- solve(crossprod(X))
  adjusted_residual <- residual / (1 - leverage)
  score <- X * as.numeric(adjusted_residual)
  vcov_hc3 <- bread %*% crossprod(score) %*% bread
  idx <- match(".exposure", names(coef(model)))
  if (is.na(idx)) stop("Exposure coefficient unavailable: ", exposure)
  beta <- unname(coef(model)[idx])
  se <- sqrt(vcov_hc3[idx, idx])
  z <- beta / se
  c(n = nrow(d), beta = beta, se = se, ci_low = beta - qnorm(.975) * se,
    ci_high = beta + qnorm(.975) * se, z = z, p = 2 * pnorm(-abs(z)))
}

run_outcomes <- function(data, modality, exposure = "primary_hdp", covariates = NULL) {
  spec <- modality_spec[[modality]]
  if (is.null(covariates)) covariates <- spec$covariates
  rows <- lapply(names(spec$outcomes), function(label) {
    fit <- hc3_fit(data, spec$outcomes[[label]], exposure, covariates)
    data.frame(modality = modality, outcome = label, t(fit), check.names = FALSE)
  })
  do.call(rbind, rows)
}

add_fdr <- function(x, family = rep("all", nrow(x))) {
  x$q <- ave(x$p, family, FUN = function(p) p.adjust(p, method = "BH"))
  x
}

save_session <- function(name) {
  path <- file.path(HBCD_OUTPUT, "provenance", paste0(name, "_sessionInfo.txt"))
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  capture.output(sessionInfo(), file = path)
}

expected_counts <- c(diffusion = 1229, structural = 1100, functional = 864)
