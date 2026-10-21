# Invented cohort matches the release-size guard without changing production code.
write_covariate_demographics <- function(a) {
  demo <- data.frame(
    participant_id = c("sub-002", "sub-001", "sub-003", "sub-001", "sub-002"),
    session_id = c("ses-V02", "ses-V02", "ses-V01", "ses-V01", "ses-V01"),
    sex = c(1, 0, NA, 1, 9),
    mother_age_delivery = c(32, 25, NA, 24, NA),
    gestational_age_delivery = c(39, 40, 37, 38, NA),
    rc_mother_education = c(4, 3, NA, 2, NA),
    rc_mother_income = c(6, 5, NA, 4, NA),
    rc_mother_ethnoracial_aou_race_ethnicity = c(2, 1, NA, 3, NA)
  )
  names(demo)[3:8] <- paste0("sed_basic_demographics_", names(demo)[3:8])
  write_release_fixture(a, "demographics_sed/sed_basic_demographics.tsv", demo)
  visits <- data.frame(
    participant_id = c("sub-001", "sub-002", "sub-001", "sub-outsider"),
    session_id = c("ses-V01", "ses-V01", "ses-V03", "ses-V03"),
    par_visit_data_site = c("old", "site-B", "latest", "unrelated")
  )
  write_release_fixture(a, "parent_visit/par_visit_data.tsv", visits)
  # Old locations exist but contain incorrect data; using them must fail assertions.
  demo[, 3:8] <- 99
  visits$par_visit_data_site <- "wrong-folder"
  write_release_fixture(a, "demographics/sed_basic_demographics.tsv", demo)
  write_release_fixture(a, "demographics/par_visit_data.tsv", visits)
}

write_covariate_exposures <- function(a) {
  ids <- sprintf("sub-%03d", seq_len(2175))
  pregnancy <- data.frame(participant_id = ids, session_id = "ses-V02")
  prefix <- "pex_bm_healthv2_preg__compl_001___"
  for (j in 1:11) {
    pregnancy[[paste0(prefix, j)]] <- 0L
  }
  pregnancy[[paste0(prefix, 2)]][1565:1870] <- 1L
  pregnancy[[paste0(prefix, 3)]][1871:2059] <- 1L
  other <- paste0(prefix, c(1, 4:11))
  pregnancy[1, other] <- NA
  pregnancy[2, other] <- NA
  pregnancy[2, other[1]] <- 1L
  write_release_fixture(
    a,
    "pregnancy_exposure/pex_bm_healthv2_preg.tsv",
    pregnancy
  )
  chronic <- data.frame(participant_id = ids, session_id = "ses-V01")
  chronic$pex_bm_health_preg__chroncond_001___1 <- c(rep(0, 1870), rep(1, 305))
  # Visit 2 must not override the chronic history used for HDP classification.
  later <- chronic[1, ]
  later$session_id <- "ses-V02"
  later$pex_bm_health_preg__chroncond_001___1 <- 1L
  write_release_fixture(
    a,
    "pregnancy_exposure/pex_bm_health_preg__chroncond.tsv",
    rbind(later, chronic[2175:1, ])
  )
  for (name in c(
    "pex_bm_healthv2_inf",
    "pex_bm_health_preg__healthhx",
    "pex_bm_health_preg__meds",
    "pex_ch_tlfb"
  )) {
    write_release_fixture(
      a,
      paste0("pregnancy_exposure/", name, ".tsv"),
      data.frame(participant_id = "sub-001", code = "00123")
    )
  }
}

write_covariate_volumes <- function(a) {
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
  prefix <- "img_bibsnet_space-T2w_desc-aseg_volumes_"
  parts <- c(paste0("Left-", tissue), paste0("Right-", tissue), "Brain-Stem")
  d <- data.frame(
    participant_id = sprintf("sub-%03d", 1:5),
    session_id = "ses-V02"
  )
  for (part in parts) {
    d[[paste0(prefix, part)]] <- 2
  }
  d[2, paste0(prefix, parts[1])] <- NA
  d[3, paste0(prefix, parts[1])] <- -1
  d[4, paste0(prefix, parts[1:12])] <- 0
  d$session_id[5] <- "ses-V01"
  # CSF must never be counted as tissue.
  d[[paste0(prefix, "Left-Lateral-Ventricle")]] <- 10000
  write_release_fixture(
    a,
    "imaging/struct_bibsnet/img_bibsnet_space-T2w_desc-aseg_volumes.tsv",
    d
  )
}

run_covariate_fixture <- function(a) {
  write_covariate_demographics(a)
  write_covariate_exposures(a)
  write_covariate_volumes(a)
  output <- system2(
    file.path(R.home("bin"), "Rscript"),
    c("--vanilla", shQuote(file.path(a$repo, "prepare/01_covariates.R"))),
    stdout = TRUE,
    stderr = TRUE
  )
  expect_null(attr(output, "status"), info = paste(output, collapse = "\n"))
}
