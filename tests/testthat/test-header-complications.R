test_that("None of these is not a complication and false does not establish absence", {
  a <- local_analysis()
  write_covariate_demographics(a)
  write_covariate_exposures(a)
  write_covariate_volumes(a)
  rel <- "pregnancy_exposure/pex_bm_healthv2_preg.tsv"
  d <- release_frame("pregnancy", a$env$read_release(rel))
  # User-supplied dictionary: item 12 is None of these; 0 false, 1 true.
  d$pex_bm_healthv2_preg__compl_001___12 <- 0
  d$pex_bm_healthv2_preg__compl_001___12[d$participant_id == "sub-003"] <- 1
  other <- paste0("pex_bm_healthv2_preg__compl_001___", c(1, 4:11))
  d[d$participant_id == "sub-004", other] <- NA
  d[d$participant_id == "sub-004", other[1]] <- 0
  write_release_fixture(a, rel, d)
  output <- system2(
    file.path(R.home("bin"), "Rscript"),
    c("--vanilla", shQuote(file.path(a$repo, "prepare/01_covariates.R"))),
    stdout = TRUE,
    stderr = TRUE
  )
  expect_null(attr(output, "status"), info = paste(output, collapse = "\n"))
  got <- a$env$read_input("covariates/v02_pregnancy_complications.tsv")
  i <- match(
    c("sub-001", "sub-002", "sub-003", "sub-004", "sub-1565", "sub-1871"),
    got$participant_id
  )
  expect_equal(got$other_complication[i], c(NA, 1, 0, 0, 0, 0))
  expect_equal(got$primary_hdp[i], c(0, 0, 0, 0, 1, 1))
  # A non-HDP infant can still have other complications; these are different groups.
})
