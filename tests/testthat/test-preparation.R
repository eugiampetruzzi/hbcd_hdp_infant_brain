test_that("diffusion preparation filters visit, duplicate, and missing-outcome rows", {
  a <- local_analysis()
  write_preparation_inputs(a)
  run_model_preparation(a)
  d <- a$env$read_input("covariates/v02_dwi_primary_covariates.tsv")
  expect_equal(d$participant_id, sprintf("sub-synthetic%03d", c(1:6, 8)))
  expect_equal(d$raw_neighbor_corr[1:2], c(.81, .82))
  expect_equal(d$candidate_age_weeks, c(1:6, 8))
  expect_equal(
    d$primary_model_complete,
    c(TRUE, TRUE, FALSE, TRUE, FALSE, TRUE, TRUE)
  )
})

test_that("structural preparation matches run and requires complete modeling variables", {
  a <- local_analysis()
  write_preparation_inputs(a)
  run_model_preparation(a)
  d <- a$env$read_input(
    "analysis_datasets/v02_global_cortical_morphometry_model_data.tsv"
  )
  expect_equal(d$participant_id, sprintf("sub-synthetic%03d", c(1, 2, 7, 8)))
  expect_equal(d$candidate_age_weeks, c(3, 4, 9, 10))
  expect_equal(d$ga_delivery, c(38, 39, 38, 40))
  expect_equal(d$primary_hdp, c(0, 1, 0, 1))
})

test_that("functional preparation includes exactly five minutes and needs complete metrics", {
  a <- local_analysis()
  write_preparation_inputs(a)
  run_model_preparation(a)
  d <- a$env$read_input("analysis_datasets/v02_global_functional_hdp_data.tsv")
  expect_equal(d$participant_id, sprintf("sub-synthetic%03d", c(1, 6, 8)))
  expect_equal(d$retained_duration_minutes, c(5, 6, 6))
  expect_equal(d$candidate_age_weeks, c(3, 8, 10))
  expect_false(anyDuplicated(d$participant_id) > 0)
})
