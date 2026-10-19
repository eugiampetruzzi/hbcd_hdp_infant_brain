test_that("shared utilities load without running harmonization", {
  a <- local_analysis()
  expect_true(is.function(a$env$fit_term))
  expect_true(dir.exists(file.path(a$env$OUT, "tables")))
  expect_false(file.exists(file.path(a$env$OUT, "tables/global_models.csv")))
})

test_that("dotenv preserves explicit environment variables", {
  a <- local_analysis()
  withr::local_envvar(c(SYNTH_EXISTING = "explicit", SYNTH_NEW = NA))
  writeLines(c("# ignored", "SYNTH_EXISTING=file", "SYNTH_NEW=a=b", "invalid"), "settings.env")
  a$env$read_env("settings.env")
  expect_identical(Sys.getenv("SYNTH_EXISTING"), "explicit")
  expect_identical(Sys.getenv("SYNTH_NEW"), "a=b")
  expect_invisible(a$env$read_env("absent.env"))
})

test_that("CSV and TSV preserve numeric missingness, names, and participant IDs", {
  e <- local_analysis()$env
  data <- data.frame(participant_id = c("sub-001", "sub-010"),
                     value = c(1.25, NA_real_), check.names = FALSE)
  data[["region-name"]] <- c(2, 3)
  for (extension in c("csv", "tsv")) {
    rel <- paste0("nested/roundtrip.", extension)
    e$write_input(data, rel)
    expect_equal(e$read_input(rel), data)
  }
  e$write_table(data, "roundtrip.csv")
  expect_equal(e$read_table("roundtrip.csv"), data)
})

test_that("raw release reading uses the configured root and preserves names", {
  a <- local_analysis()
  data <- data.frame(participant_id = "sub-001", check.names = FALSE)
  data[["column-with-hyphens"]] <- "001"
  write_release_fixture(a, "raw.tsv", data)
  expect_equal(a$env$read_release("raw.tsv", colClasses = "character"), data)
  a$env$RELEASE <- ""
  expect_error(a$env$read_release("raw.tsv"), "HBCD_RELEASE")
})

test_that("batch threshold is inclusive and unknown batches are excluded", {
  e <- local_analysis()$env
  data <- data.frame(id = 1:20, batch = c(rep("A", 10), rep("B", 9), NA))
  expect_equal(e$drop_small_batches(data)$id, 1:10)
  expect_equal(e$drop_small_batches(data, min_n = 9)$id, 1:19)
})

test_that("complete-case filtering precedes batch-size filtering", {
  e <- local_analysis()$env
  frames <- synthetic_frames()
  frames[[1]]$raw_neighbor_corr[1:3] <- NA
  frames[[1]]$site[13] <- NA
  result <- e$modality_sample(frames, "Diffusion MRI")
  expect_equal(result$participant_id, frames[[1]]$participant_id[13:24])
  expect_equal(nrow(result), 12L)
  expect_true(is.na(result$site[1]))
})

test_that("extra covariates restrict the analytic sample", {
  e <- local_analysis()$env
  frames <- synthetic_frames()
  frames[[2]]$extra <- seq_len(24)
  frames[[2]]$extra[1] <- NA
  result <- e$modality_sample(frames, "Structural MRI", extra = "extra")
  expect_equal(nrow(result), 23L)
  expect_false(anyNA(result$extra))
})

test_that("modality joins use participant identity and current diffusion outcomes", {
  a <- local_analysis()
  expected <- synthetic_frames()
  scanners <- write_modality_inputs(a, expected)
  scanners$dwi_unique_scanners[1] <- 2L
  a$env$write_input(scanners[nrow(scanners):1, ], "covariates/v02_scanner_metadata.tsv")
  frames <- a$env$load_modalities()
  expect_identical(names(frames), names(expected))
  expect_true(is.na(frames[[1]]$batch[1]))
  expect_equal(frames[[1]]$batch[13], "B | v1")
  for (j in 1:3) {
    expect_equal(frames[[j]]$participant_id, expected[[j]]$participant_id)
    expect_equal(frames[[j]]$pma, expected[[j]]$pma)
    expect_equal(frames[[j]]$sex_m, expected[[j]]$sex_m)
    expect_true(all(frames[[j]]$site == "synthetic-site"))
  }
  expect_equal(frames[[1]]$global_median_mkt, expected[[1]]$global_median_mkt)
  expect_false("global_median_mkt.x" %in% names(frames[[1]]))
})

test_that("incomplete diffusion participants are omitted before joining", {
  a <- local_analysis()
  write_modality_inputs(a)
  rel <- "covariates/v02_dwi_primary_covariates.tsv"
  data <- a$env$read_input(rel)
  data$primary_model_complete[1:2] <- c(FALSE, NA)
  a$env$write_input(data, rel)
  expect_equal(nrow(a$env$load_modalities()[[1]]), 22L)
})

test_that("standardization uses population variance and is affine invariant", {
  e <- local_analysis()$env
  expect_equal(e$zscore(c(1, 2, 3)), c(-sqrt(1.5), 0, sqrt(1.5)))
  expect_equal(e$zscore(c(12, 14, 16)), e$zscore(c(1, 2, 3)))
})

test_that("formulas preserve exposure and age and include modality-specific QC", {
  e <- local_analysis()$env
  expect_equal(e$smooth("age", k = 5), "s(age, k = 5)")
  rhs <- e$primary_rhs("Functional MRI", extra = "bmi", ga = TRUE)
  expect_setequal(all.vars(as.formula(paste("y ~", rhs))),
                  c("y", "primary_hdp", "sex_m", "mean_fd_post_censoring",
                    "retained_duration_minutes", "bmi", "pma", "ga_delivery"))
  preserved <- e$harmonization_rhs(extra = "bmi", ga = TRUE)
  expect_setequal(all.vars(as.formula(paste("y ~", preserved))),
                  c("y", "primary_hdp", "sex_m", "bmi", "pma", "ga_delivery"))
})
