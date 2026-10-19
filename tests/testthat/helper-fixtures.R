# Every source operation runs with synthetic environment paths and no local .env.
local_analysis <- function(.local_envir = parent.frame()) {
  root <- withr::local_tempdir(.local_envir = .local_envir)
  withr::local_envvar(c(HBCD_HDP_DATA = file.path(root, "data"),
                      HBCD_RELEASE = file.path(root, "release")),
                    .local_envir = .local_envir)
  withr::local_dir(root, .local_envir = .local_envir)
  repo <- getOption("hbcd.test.root")
  if (is.null(repo)) stop("Use Rscript tests/testthat.R to run these tests")
  # Preparation scripts source R/common.R relative to their working directory.
  dir.create("R")
  file.copy(file.path(repo, "R/common.R"), "R/common.R")
  env <- new.env(parent = globalenv())
  sys.source("R/common.R", envir = env)
  list(env = env, root = root, repo = repo)
}

write_release_fixture <- function(a, rel, data) {
  path <- file.path(a$env$RELEASE, rel)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.table(data, path, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
}

synthetic_frames <- function(n = 24L) {
  i <- seq_len(n)
  base <- data.frame(participant_id = sprintf("sub-synthetic%03d", i),
                     primary_hdp = rep(0:1, length.out = n),
                     ga_delivery = 35 + (i %% 6), pma = 41 + i / n,
                     candidate_age_weeks = 6 + i / n - (i %% 6),
                     sex_m = rep(c(0, 0, 1, 1), length.out = n),
                     site = "synthetic-site", batch = rep(c("A", "B"), each = n / 2))
  list(
    "Diffusion MRI" = transform(base, raw_neighbor_corr = .9,
                                global_median_md = .001 + i * .00001,
                                global_median_fa = .2 + i * .001,
                                global_median_mkt = .5 + i * .001),
    "Structural MRI" = transform(base, brainswipes_summary_qc = 2,
                                 global_area_weighted_thickness = 2 + i * .01,
                                 global_area_weighted_absolute_sulc = 1 + i * .01),
    "Functional MRI" = transform(base, mean_fd_post_censoring = .1,
                                 retained_duration_minutes = 6,
                                 modularity = .3 + i * .001,
                                 mean_participation_coefficient = .2 + i * .001)
  )
}

write_modality_inputs <- function(a, frames = synthetic_frames()) {
  e <- a$env
  ids <- frames[[1]]$participant_id
  e$write_input(data.frame(participant_id = ids, site = "synthetic-site"),
                "covariates/hdp_missingness_person_level.csv")
  dwi <- frames[[1]]
  dwi$primary_model_complete <- TRUE
  dwi$site <- dwi$batch <- dwi$pma <- NULL
  dwi$global_median_md <- dwi$global_median_fa <- NULL
  # Deliberately stale MKT must be replaced by the outcome file.
  dwi$global_median_mkt <- -999
  e$write_input(dwi, "covariates/v02_dwi_primary_covariates.tsv")
  e$write_input(frames[[1]][, c("participant_id", names(e$MODALITIES[[1]]$outcomes))],
                "outcomes/diffusion/v02_global_cerebral_wm_fa_md.tsv")
  for (j in 2:3) {
    d <- frames[[j]]
    d$child_sex <- d$sex_m
    d$sex_m <- d$batch <- d$pma <- NULL
    rel <- c("v02_global_cortical_morphometry_model_data.tsv",
             "v02_global_functional_hdp_data.tsv")[[j - 1]]
    e$write_input(d, file.path("analysis_datasets", rel))
  }
  scanners <- data.frame(participant_id = ids)
  for (prefix in c("dwi", "t2w", "bold")) {
    scanners[[paste0(prefix, "_unique_scanners")]] <- 1L
    scanners[[paste0(prefix, "_scanner_serial_number")]] <- frames[[1]]$batch
    scanners[[paste0(prefix, "_scanner_software")]] <- "v1"
  }
  e$write_input(scanners, "covariates/v02_scanner_metadata.tsv")
  invisible(scanners)
}
