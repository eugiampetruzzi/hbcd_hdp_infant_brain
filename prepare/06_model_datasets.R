# Modality datasets: exposure, covariates, age at scan and image quality joined to the
# imaging outcomes from steps 03-05.

source("R/common.R")

exposure <- read_input("covariates/v02_pregnancy_complications.tsv")[, c(
  ID,
  HDP,
  "chronic_hypertension"
)]
person <- read_input("covariates/hdp_missingness_person_level.csv")[, c(
  ID,
  GA,
  "child_sex",
  "site"
)]
covariates <- function(d) {
  merge(
    merge(d, exposure, by = ID, all.x = TRUE),
    person,
    by = ID,
    all.x = TRUE
  )
}
weeks <- function(years) as.numeric(years) * 365.25 / 7
v02 <- function(x) x[x$session_id == "ses-V02", ]
outcomes <- lapply(MODALITIES, function(m) names(m$outcomes))
qc <- lapply(MODALITIES, `[[`, "qc")

# Diffusion: infants with valid global measures; postnatal age and neighboring-DWI
# correlation from the QSIPrep image QC table.
q <- "img_qsiprep_space-ACPC_desc-image_qc_"
qsiprep <- v02(read_release(
  "imaging/dwi_qsiprep/img_qsiprep_space-ACPC_desc-image_qc.tsv"
))
qsiprep <- qsiprep[!duplicated(qsiprep[[ID]]), ]
dwi <- data.frame(
  qsiprep[[ID]],
  weeks(qsiprep[[paste0(q, "candidate_age")]]),
  qsiprep[[paste0(q, "raw_neighbor_corr")]]
)
names(dwi) <- c(ID, "candidate_age_weeks", "raw_neighbor_corr")
diffusion <- read_input("outcomes/diffusion/v02_global_cerebral_wm_fa_md.tsv")
dwi <- merge(
  diffusion[!is.na(diffusion$global_median_mkt), ID, drop = FALSE],
  dwi,
  by = ID,
  all.x = TRUE
)
dwi <- covariates(dwi)
dwi$sex_m <- dwi$child_sex
dwi$primary_model_complete <- complete.cases(dwi[, c(
  HDP,
  "candidate_age_weeks",
  GA,
  "sex_m",
  "raw_neighbor_corr",
  "site"
)])
write_input(
  dwi[, c(
    ID,
    HDP,
    "chronic_hypertension",
    GA,
    "sex_m",
    "candidate_age_weeks",
    "raw_neighbor_corr",
    "primary_model_complete"
  )],
  "covariates/v02_dwi_primary_covariates.tsv"
)

# Structural: postnatal age and BrainSwipes score for the run used by XCP-D. Postnatal age
# from this table is also used for the functional sample.
b <- "img_brainswipes_xcpd_hash-0f306a2f+0ef9c88a_T2w_"
swipes <- v02(read_release(
  "imaging/qc_brainswipes/img_brainswipes_xcpd_hash-0f306a2f+0ef9c88a_T2w.tsv"
))
swipes <- data.frame(
  swipes[[ID]],
  as.numeric(sub("^run-", "", swipes$run_id)),
  weeks(swipes[[paste0(b, "candidate_age")]]),
  swipes[[paste0(b, "summary_QC")]]
)
names(swipes) <- c(
  ID,
  "source_run",
  "candidate_age_weeks",
  "brainswipes_summary_qc"
)
morphometry <- read_input(
  "outcomes/structural/v02_global_cortical_morphometry.tsv"
)
structural <- covariates(merge(
  morphometry,
  swipes,
  by = c(ID, "source_run"),
  all.x = TRUE
))
base <- c(
  ID,
  HDP,
  "chronic_hypertension",
  GA,
  "child_sex",
  "candidate_age_weeks"
)
smri <- c(qc[["Structural MRI"]], outcomes[["Structural MRI"]])
keep <- complete.cases(structural[, c(
  setdiff(base, "chronic_hypertension"),
  "site",
  smri
)])
write_input(
  structural[keep, c(base, smri)],
  "analysis_datasets/v02_global_cortical_morphometry_model_data.tsv"
)

# Functional: complete network metrics and at least 5 minutes of retained low-motion data.
fmri <- c(qc[["Functional MRI"]], outcomes[["Functional MRI"]])
functional <- merge(
  structural[, c(base, "site")],
  read_input("covariates/v02_functional_qc.tsv")[, c(
    ID,
    qc[["Functional MRI"]]
  )],
  by = ID
)
network <- read_input("outcomes/functional/v02_functional_network_metrics.tsv")
functional <- merge(
  functional,
  network[, c(ID, outcomes[["Functional MRI"]], "mean_network_connectivity")],
  by = ID
)
keep <- complete.cases(functional[, c(
  HDP,
  "candidate_age_weeks",
  "child_sex",
  "site",
  fmri,
  "mean_network_connectivity"
)])
write_input(
  functional[keep & functional$retained_duration_minutes >= 5, c(base, fmri)],
  "analysis_datasets/v02_global_functional_hdp_data.tsv"
)
