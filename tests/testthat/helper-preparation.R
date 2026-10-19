write_preparation_inputs <- function(a) {
  e <- a$env
  ids <- sprintf("sub-synthetic%03d", 1:8)
  exposure <- data.frame(
    participant_id = ids,
    primary_hdp = c(0, 1, NA, 0, 1, 1, 0, 1),
    chronic_hypertension = 0
  )
  person <- data.frame(
    participant_id = ids,
    ga_delivery = c(38, 39, 36, 40, 37, 37, 38, 40),
    child_sex = c(0, 1, 0, 1, NA, 1, 0, 0),
    site = "synthetic-site"
  )
  e$write_input(exposure[8:1, ], "covariates/v02_pregnancy_complications.tsv")
  e$write_input(person[8:1, ], "covariates/hdp_missingness_person_level.csv")
  write_diffusion_inputs(a, ids)
  write_structural_inputs(a, ids)
  e$write_input(
    data.frame(
      participant_id = ids,
      mean_fd_post_censoring = .1,
      retained_duration_minutes = c(5, 4.99, rep(6, 6))
    ),
    "covariates/v02_functional_qc.tsv"
  )
  e$write_input(
    data.frame(
      participant_id = ids,
      modularity = c(rep(.3, 6), NA, .3),
      mean_participation_coefficient = .4,
      mean_network_connectivity = .1
    ),
    "outcomes/functional/v02_functional_network_metrics.tsv"
  )
}

write_diffusion_inputs <- function(a, ids) {
  qc <- data.frame(participant_id = ids, session_id = "ses-V02")
  prefix <- "img_qsiprep_space-ACPC_desc-image_qc_"
  qc[[paste0(prefix, "candidate_age")]] <- (1:8) * 7 / 365.25
  qc[[paste0(prefix, "raw_neighbor_corr")]] <- .8 + (1:8) / 100
  prior_visit <- qc[1, ]
  prior_visit$session_id <- "ses-V01"
  prior_visit[[paste0(prefix, "raw_neighbor_corr")]] <- .01
  duplicate <- qc[2, ]
  duplicate[[paste0(prefix, "raw_neighbor_corr")]] <- .02
  write_release_fixture(
    a,
    "imaging/dwi_qsiprep/img_qsiprep_space-ACPC_desc-image_qc.tsv",
    rbind(prior_visit, qc, duplicate)
  )
  a$env$write_input(
    data.frame(
      participant_id = ids,
      global_median_md = .001,
      global_median_fa = .2,
      global_median_mkt = c(rep(.5, 6), NA, .5)
    ),
    "outcomes/diffusion/v02_global_cerebral_wm_fa_md.tsv"
  )
}

write_structural_inputs <- function(a, ids) {
  swipes <- data.frame(
    participant_id = ids,
    session_id = "ses-V02",
    run_id = c(rep("run-1", 3), "run-2", rep("run-1", 4))
  )
  prefix <- "img_brainswipes_xcpd_hash-0f306a2f+0ef9c88a_T2w_"
  swipes[[paste0(prefix, "candidate_age")]] <- (3:10) * 7 / 365.25
  swipes[[paste0(prefix, "summary_QC")]] <- 2
  write_release_fixture(
    a,
    "imaging/qc_brainswipes/img_brainswipes_xcpd_hash-0f306a2f+0ef9c88a_T2w.tsv",
    swipes[8:1, ]
  )
  a$env$write_input(
    data.frame(
      participant_id = ids,
      session_id = "ses-V02",
      source_run = 1,
      global_area_weighted_thickness = c(rep(2, 5), NA, 2, 2),
      global_area_weighted_absolute_sulc = 1
    ),
    "outcomes/structural/v02_global_cortical_morphometry.tsv"
  )
}

run_model_preparation <- function(a) {
  script <- file.path(a$repo, "prepare/06_model_datasets.R")
  output <- system2(
    file.path(R.home("bin"), "Rscript"),
    c("--vanilla", shQuote(script)),
    stdout = TRUE,
    stderr = TRUE
  )
  expect_null(attr(output, "status"), info = paste(output, collapse = "\n"))
}
