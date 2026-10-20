test_that("Release 2.1 covariates join by ID with earliest values and latest site", {
  a <- local_analysis()
  run_covariate_fixture(a)
  # write_input serializes missing values as empty CSV fields, including strings.
  d <- read.csv(
    file.path(a$env$DATA, "covariates/hdp_missingness_person_level.csv"),
    na.strings = c("", "NA")
  )
  d <- d[order(d$participant_id), ]
  expect_equal(d$participant_id, c("sub-001", "sub-002", "sub-003"))
  expect_equal(d$mother_age, c(24, 32, NA))
  expect_equal(d$ga_delivery, c(38, 39, 37))
  expect_equal(d$child_sex, c(1, 1, NA))
  expect_equal(d$education, c(2, 4, NA))
  expect_equal(d$income, c(4, 6, NA))
  expect_equal(d$race_eth, c(3, 2, NA))
  expect_equal(d$site, c("latest", "site-B", NA))
})

test_that("covariate preparation preserves exposure groups and unknown complications", {
  a <- local_analysis()
  run_covariate_fixture(a)
  d <- a$env$read_input("covariates/v02_pregnancy_complications.tsv")
  expect_equal(nrow(d), 2175L)
  expect_equal(sum(d$primary_hdp == 0, na.rm = TRUE), 1564)
  expect_equal(sum(d$primary_hdp == 1, na.rm = TRUE), 495)
  expect_equal(sum(is.na(d$primary_hdp)), 116)
  i <- match(
    c("sub-001", "sub-002", "sub-003", "sub-1565", "sub-1871", "sub-2060"),
    d$participant_id
  )
  expect_equal(d$primary_hdp[i], c(0, 0, 0, 1, 1, NA))
  expect_equal(d$other_complication[i[1:3]], c(NA, 1, 0))
  expect_equal(d$chronic_hypertension[i], c(0, 0, 0, 0, 1, 1))
})

test_that("covariate preparation excludes invalid tissue volumes and preserves source codes", {
  a <- local_analysis()
  run_covariate_fixture(a)
  d <- a$env$read_input("covariates/v02_bibsnet_brain_volume.tsv")
  expect_equal(d$participant_id, sprintf("sub-%03d", 1:4))
  expect_equal(d$total_brain_volume_mm3, c(50, NA, NA, NA))
  copied <- read.csv(
    file.path(a$env$DATA, "phenotype/pex_bm_health_preg__meds.csv"),
    colClasses = "character"
  )
  expect_equal(copied$code, "00123")
})
