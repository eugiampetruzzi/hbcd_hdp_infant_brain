expand_release <- function(a, rel, kind) {
  release_frame(kind, a$env$read_release(rel))
}

test_that("competing QC fields do not alter ages, quality covariates or run matching", {
  a <- local_analysis()
  # Supply scanner identities so the prepared tables can be followed into model loading.
  write_modality_inputs(a, synthetic_frames(8))
  write_preparation_inputs(a)
  drel <- "imaging/dwi_qsiprep/img_qsiprep_space-ACPC_desc-image_qc.tsv"
  srel <- "imaging/qc_brainswipes/img_brainswipes_xcpd_hash-0f306a2f+0ef9c88a_T2w.tsv"
  d <- expand_release(a, drel, "diffusion_qc")
  d <- d[
    d$session_id == "ses-V02" &
      !duplicated(paste(d$participant_id, d$session_id)),
  ]
  q <- "img_qsiprep_space-ACPC_desc-image_qc_"
  d[[paste0(q, "candidate_age")]] <- 42 / 365.25
  # Tabulated candidate_age is years; adjusted_age is weeks (HBCD agevariables docs).
  d[[paste0(q, "adjusted_age")]] <- 2
  d[[paste0(q, "raw_neighbor_corr")]] <- .81
  d[[paste0(q, "raw_masked_neighbor_corr")]] <- .91
  d[[paste0(q, "t1_neighbor_corr")]] <- .71
  d[[paste0(q, "t1_masked_neighbor_corr")]] <- .61
  s <- expand_release(a, srel, "structural_qc")
  b <- "img_brainswipes_xcpd_hash-0f306a2f+0ef9c88a_T2w_"
  s[[paste0(b, "adjusted_age")]] <- 2
  unused <- grep("(_mean|_nrev)$", names(s), value = TRUE)
  s[, unused] <- 9
  extra <- s[s$participant_id == "sub-synthetic001", ]
  extra$run_id <- "run-2"
  extra[[paste0(b, "candidate_age")]] <- 70 / 365.25
  extra[[paste0(b, "summary_QC")]] <- 1
  s <- rbind(extra, s)
  outputs <- list()
  for (iteration in 1:2) {
    write_release_fixture(a, drel, d)
    write_release_fixture(a, srel, s)
    run_model_preparation(a)
    outputs[[iteration]] <- a$env$load_modalities()
    d <- d[nrow(d):1, ]
    d[[paste0(q, "adjusted_age")]] <- 100
    d[[paste0(q, "raw_masked_neighbor_corr")]] <- .11
    s <- s[nrow(s):1, ]
    s[[paste0(b, "adjusted_age")]] <- 100
    s[, unused] <- 99
  }
  expect_equal(outputs[[1]], outputs[[2]])
  first <- outputs[[1]][[1]]
  i <- match("sub-synthetic001", first$participant_id)
  expect_equal(first$candidate_age_weeks[i], 6)
  expect_equal(first$pma[i], 44)
  expect_true(all(first$raw_neighbor_corr == .81))
  structural <- outputs[[1]][[2]]
  expect_equal(
    structural$participant_id,
    sprintf("sub-synthetic%03d", c(1, 2, 7, 8))
  )
  expect_equal(structural$candidate_age_weeks, c(3, 4, 9, 10))
  expect_equal(structural$brainswipes_summary_qc, rep(2, 4))
})

test_that("demographic alternatives and non-tissue volumes cannot change covariates", {
  a <- local_analysis()
  write_covariate_demographics(a)
  write_covariate_exposures(a)
  write_covariate_volumes(a)
  drel <- "demographics_sed/sed_basic_demographics.tsv"
  vrel <- "imaging/struct_bibsnet/img_bibsnet_space-T2w_desc-aseg_volumes.tsv"
  d <- expand_release(a, drel, "demographics")
  d$sed_basic_demographics_mother_age_v01 <- 29
  d$sed_basic_demographics_screen_mother_race <- 6
  d$sed_basic_demographics_child_race <- 5
  d$sed_basic_demographics_recruitment_site <- "wrong-source"
  v <- expand_release(a, vrel, "brain_volumes")
  prefix <- "img_bibsnet_space-T2w_desc-aseg_volumes_"
  tissue <- c(
    "Cerebral-White-Matter",
    "Cerebral-Cortex",
    "Cerebellum-Cortex",
    "Cerebellum-White-Matter",
    "Thalamus-Proper*",
    "Caudate",
    "Putamen",
    "Pallidum",
    "Hippocampus",
    "Amygdala",
    "Accumbens-area",
    "VentralDC"
  )
  parts <- paste0(
    prefix,
    c(paste0("Left-", tissue), paste0("Right-", tissue), "Brain-Stem")
  )
  v[1, parts] <- seq_len(25)
  # BIBSNet counts discrete labels separately, including Vermis (label 172).
  v[[paste0(prefix, "Vermis")]] <- 0
  excluded <- grep(
    "CSF|Ventricle|choroid|Unknown|_age$",
    names(v),
    value = TRUE,
    ignore.case = TRUE
  )
  results <- list()
  for (iteration in 1:2) {
    v[, excluded] <- iteration * 10000
    d$sed_basic_demographics_mother_age_v01 <- 29 + iteration
    write_release_fixture(a, drel, d)
    write_release_fixture(a, vrel, v)
    output <- system2(
      file.path(R.home("bin"), "Rscript"),
      c("--vanilla", shQuote(file.path(a$repo, "prepare/01_covariates.R"))),
      stdout = TRUE,
      stderr = TRUE
    )
    expect_null(attr(output, "status"), info = paste(output, collapse = "\n"))
    results[[iteration]] <- list(
      person = a$env$read_input("covariates/hdp_missingness_person_level.csv"),
      volume = a$env$read_input("covariates/v02_bibsnet_brain_volume.tsv")
    )
  }
  expect_equal(results[[1]], results[[2]])
  person <- results[[1]]$person
  i <- match(c("sub-001", "sub-002"), person$participant_id)
  expect_equal(person$mother_age[i], c(24, 32))
  expect_equal(person$race_eth[i], c(3, 2))
  expect_equal(person$site[i], c("latest", "site-B"))
  v <- results[[1]]$volume
  expect_equal(
    v$total_brain_volume_mm3[match("sub-001", v$participant_id)],
    325
  )
  transformed <- local_sensitivity_helpers()$brain_volume_covariate(v[
    nrow(v):1,
  ])
  expect_equal(
    transformed$log_tbv[match("sub-001", transformed$participant_id)],
    log(325)
  )
  expect_true(all(is.na(transformed$log_tbv[
    transformed$participant_id != "sub-001"
  ])))
  raw <- a$env$read_release(vrel)
  raw[[paste0(prefix, "Vermis")]][1] <- 7
  write_release_fixture(a, vrel, raw)
  output <- system2(file.path(R.home("bin"), "Rscript"),
                    c("--vanilla", shQuote(file.path(a$repo, "prepare/01_covariates.R"))),
                    stdout = TRUE, stderr = TRUE)
  expect_null(attr(output, "status"), info = paste(output, collapse = "\n"))
  updated <- a$env$read_input("covariates/v02_bibsnet_brain_volume.tsv")
  expect_equal(updated$total_brain_volume_mm3[match("sub-001", updated$participant_id)], 332)
  transformed <- local_sensitivity_helpers()$brain_volume_covariate(updated)
  expect_equal(transformed$log_tbv[match("sub-001", transformed$participant_id)], log(332))
})
