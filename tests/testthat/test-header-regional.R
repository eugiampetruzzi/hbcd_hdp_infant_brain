test_that("network summaries retain PostDefault and separate thickening from thinning", {
  e <- new.env(parent = globalenv())
  sys.source(file.path(getOption("hbcd.test.root"), "R/regional_summary.R"), envir = e)
  names <- sub("_hash-.*$", "", release_columns("parcel_thickness"))
  features <- c("PostDefault", "L_PostDefault_1", "L_PostDefault_2")
  expect_true(all(features %in% names))
  results <- data.frame(metric = c(rep("Parcel cortical thickness", 3), "Parcel sulcal depth"),
                        feature = c(features, features[1]), beta = c(.3, -.2, .4, .8),
                        q = c(.01, .02, .2, .001))
  got <- e$parcel_network_summary(results)
  expect_equal(got$network, "PostDefault")
  expect_equal(got$parcels, 3)
  expect_equal(got$significant, 1)
  expect_equal(got$thinner, 1)
  expect_equal(e$parcel_network_summary(results[4:1, ]), got)
})
