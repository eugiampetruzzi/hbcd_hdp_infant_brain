suppressPackageStartupMessages({
  library(mgcv)
  library(ComBatFamily)
})

read_env <- function(path = ".env") {
  if (!file.exists(path)) return(invisible())
  kv <- grep("^[A-Za-z_]+=", readLines(path), value = TRUE)
  args <- setNames(as.list(sub("^[^=]+=", "", kv)), sub("=.*", "", kv))
  args <- args[Sys.getenv(names(args)) == ""]
  if (length(args)) do.call(Sys.setenv, args)
}

read_env()
DATA <- Sys.getenv("HBCD_HDP_DATA")
if (DATA == "") stop("Set HBCD_HDP_DATA in .env (see .env.example)")
RELEASE <- Sys.getenv("HBCD_RELEASE")
OUT <- file.path(DATA, "derivatives")
dir.create(file.path(OUT, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT, "figures"), showWarnings = FALSE)

ID <- "participant_id"; HDP <- "primary_hdp"; GA <- "ga_delivery"; PMA <- "pma"
K <- 10
MIN_BATCH <- 10

MODALITIES <- list(
  "Diffusion MRI" = list(
    scanner = "dwi", qc = "raw_neighbor_corr",
    outcomes = c(global_median_md = "Mean diffusivity",
                 global_median_fa = "Fractional anisotropy",
                 global_median_mkt = "Mean kurtosis tensor")),
  "Structural MRI" = list(
    scanner = "t2w", qc = "brainswipes_summary_qc",
    outcomes = c(global_area_weighted_thickness = "Cortical thickness",
                 global_area_weighted_absolute_sulc = "Absolute sulcal depth")),
  "Functional MRI" = list(
    scanner = "bold", qc = c("mean_fd_post_censoring", "retained_duration_minutes"),
    outcomes = c(modularity = "Network modularity",
                 mean_participation_coefficient = "Participation coefficient")))

read_input <- function(rel) {
  p <- file.path(DATA, rel)
  if (grepl("\\.csv$", p)) read.csv(p, check.names = FALSE) else read.delim(p, check.names = FALSE)
}

write_input <- function(x, rel) {
  p <- file.path(DATA, rel)
  dir.create(dirname(p), recursive = TRUE, showWarnings = FALSE)
  if (grepl("\\.csv$", p)) write.csv(x, p, row.names = FALSE, na = "")
  else write.table(x, p, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
}

# Raw release tables, used only by prepare/.
read_release <- function(rel, ...) {
  if (RELEASE == "") stop("Set HBCD_RELEASE in .env (see .env.example)")
  read.delim(file.path(RELEASE, rel), check.names = FALSE, ...)
}

# One row per infant per modality with outcomes, covariates, PMA and scanner batch.
# Batch is scanner serial x software version; infants scanned on >1 scanner get NA.
load_modalities <- function() {
  site <- read_input("covariates/hdp_missingness_person_level.csv")
  site <- site[!duplicated(site[[ID]]), c(ID, "site")]
  dwi <- read_input("covariates/v02_dwi_primary_covariates.tsv")
  dwi <- dwi[dwi$primary_model_complete %in% TRUE, setdiff(names(dwi), "global_median_mkt")]
  diffusion <- read_input("outcomes/diffusion/v02_global_cerebral_wm_fa_md.tsv")
  dwi <- merge(dwi, diffusion[, c(ID, names(MODALITIES[[1]]$outcomes))], by = ID, all.x = TRUE)
  smri <- read_input("analysis_datasets/v02_global_cortical_morphometry_model_data.tsv")
  fmri <- read_input("analysis_datasets/v02_global_functional_hdp_data.tsv")
  smri$site <- NULL; fmri$site <- NULL
  smri$sex_m <- as.numeric(smri$child_sex == 1)
  fmri$sex_m <- as.numeric(fmri$child_sex == 1)
  scanners <- read_input("covariates/v02_scanner_metadata.tsv")
  frames <- list("Diffusion MRI" = dwi, "Structural MRI" = smri, "Functional MRI" = fmri)
  for (m in names(frames)) {
    pre <- MODALITIES[[m]]$scanner
    single <- scanners[[paste0(pre, "_unique_scanners")]] == 1
    batch <- ifelse(single, paste(scanners[[paste0(pre, "_scanner_serial_number")]],
                                  scanners[[paste0(pre, "_scanner_software")]], sep = " | "), NA)
    frames[[m]] <- merge(frames[[m]], site, by = ID, all.x = TRUE)
    frames[[m]] <- merge(frames[[m]], setNames(data.frame(scanners[[ID]], batch), c(ID, "batch")),
                         by = ID, all.x = TRUE)
    frames[[m]][[PMA]] <- frames[[m]][[GA]] + frames[[m]]$candidate_age_weeks
  }
  frames
}

drop_small_batches <- function(d, min_n = MIN_BATCH) {
  n <- table(d$batch)
  d[d$batch %in% names(n)[n >= min_n], ]
}

# Analytic sample for a modality: complete cases, single-scanner infants, batches >= min_n.
modality_sample <- function(frames, m, extra = character(0), min_n = MIN_BATCH) {
  spec <- MODALITIES[[m]]
  d <- frames[[m]][, unique(c(ID, HDP, GA, PMA, "sex_m", "site", spec$qc,
                              names(spec$outcomes), "batch", extra))]
  d <- d[complete.cases(d[, setdiff(names(d), "site")]), ]
  drop_small_batches(d, min_n)
}

smooth <- function(x, k = K) paste0("s(", x, ", k = ", k, ")")

# CovBat-GAM; `rhs` lists the terms preserved during harmonization and must match the
# analysis model. Image quality is left out because it is confounded with scanner.
harmonize <- function(d, ys, rhs, eb = FALSE) {
  f <- as.formula(paste("y ~", rhs))
  covars <- d[, all.vars(f)[-1], drop = FALSE]
  d[, ys] <- covfam(as.matrix(d[, ys]), factor(d$batch), covars, mgcv::gam, f, eb = eb)$dat.covbat
  d
}

zscore <- function(v) (v - mean(v)) / sqrt(mean((v - mean(v))^2))

# Coefficient on `term` with heteroskedasticity-robust (sandwich) SE.
fit_term <- function(d, y, rhs, term = HDP) {
  d$z <- zscore(d[[y]])
  g <- gam(as.formula(paste("z ~", rhs)), data = d, method = "REML")
  b <- coef(g)[[term]]
  se <- sqrt(vcov(g, sandwich = TRUE)[term, term])
  data.frame(n = nrow(d), n_hdp = sum(d[[HDP]] == 1), beta = b, se = se,
             ci_lo = b - 1.96 * se, ci_hi = b + 1.96 * se, p = 2 * pnorm(-abs(b / se)))
}

primary_rhs <- function(m, extra = character(0), ga = FALSE) {
  paste(c(HDP, "sex_m", MODALITIES[[m]]$qc, extra, smooth(PMA), if (ga) smooth(GA)), collapse = " + ")
}

harmonization_rhs <- function(extra = character(0), ga = FALSE) {
  paste(c(HDP, "sex_m", extra, smooth(PMA), if (ga) smooth(GA)), collapse = " + ")
}

write_table <- function(x, name) write.csv(x, file.path(OUT, "tables", name), row.names = FALSE)
read_table <- function(name) read.csv(file.path(OUT, "tables", name), check.names = FALSE)
