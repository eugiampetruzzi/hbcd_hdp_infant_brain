test_that("group difference and robust uncertainty match a hand-computed example", {
  e <- local_analysis()$env
  d <- data.frame(primary_hdp = rep(0:1, each = 6),
                  outcome = c(4, 5, 6, 4, 5, 6, 5, 7, 9, 5, 7, 9))
  result <- e$fit_term(d, "outcome", "primary_hdp")
  expect_equal(result$n, 12L)
  expect_equal(result$n_hdp, 6L)
  expect_equal(result$beta, sqrt(1.5), tolerance = 1e-7)
  expect_equal(result$se, .5, tolerance = 1e-7)
  expect_equal(result$ci_lo, sqrt(1.5) - .98, tolerance = 1e-7)
  expect_equal(result$ci_hi, sqrt(1.5) + .98, tolerance = 1e-7)
  expect_equal(result$p, .0143058784354296, tolerance = 1e-7)
})

test_that("covariate-adjusted coefficients match independent OLS and HC1 algebra", {
  e <- local_analysis()$env
  withr::local_seed(19)
  n <- 80
  d <- data.frame(primary_hdp = rep(0:1, each = n / 2), x = rnorm(n))
  d$outcome <- 2 * d$primary_hdp + .7 * d$x + rnorm(n, sd = 1 + d$primary_hdp)
  y <- (d$outcome - mean(d$outcome)) / sqrt(mean((d$outcome - mean(d$outcome))^2))
  X <- cbind(1, d$primary_hdp, d$x)
  bread <- solve(crossprod(X))
  beta <- bread %*% crossprod(X, y)
  residual <- drop(y - X %*% beta)
  covariance <- n / (n - ncol(X)) * bread %*% crossprod(X * residual) %*% bread
  result <- e$fit_term(d, "outcome", "primary_hdp + x", term = "x")
  expect_equal(result$beta, beta[3, 1], tolerance = 1e-7)
  expect_equal(result$se, sqrt(covariance[3, 3]), tolerance = 1e-7)
})

test_that("a nonlinear age model recovers the balanced exposure contrast", {
  e <- local_analysis()$env
  d <- expand.grid(pma = seq(40, 50, length.out = 60), sex_m = 0:1, primary_hdp = 0:1)
  d$outcome <- 2 * sin(d$pma) + .7 * d$sex_m + 1.5 * d$primary_hdp + .05 * cos(7 * d$pma)
  scale <- sqrt(mean((d$outcome - mean(d$outcome))^2))
  result <- e$fit_term(d, "outcome", "primary_hdp + sex_m + s(pma, k = 10)")
  expect_equal(result$beta, 1.5 / scale, tolerance = 1e-6)
  expect_true(is.finite(result$se) && result$se > 0)
  expect_true(result$ci_lo < 1.5 / scale && result$ci_hi > 1.5 / scale)
})

test_that("real CovBat reduces a balanced scanner shift while preserving exposure", {
  skip_if_not_installed("ComBatFamily")
  e <- local_analysis()$env
  d <- expand.grid(batch = c("A", "B"), primary_hdp = 0:1, replicate = 1:30)
  noise <- sin(d$replicate) + cos(2 * d$replicate)
  d$first <- 2 * d$primary_hdp + 5 * (d$batch == "B") + noise
  d$second <- d$primary_hdp - 3 * (d$batch == "B") + cos(d$replicate)
  result <- e$harmonize(d, c("first", "second"), "primary_hdp")
  expect_equal(result[, c("batch", "primary_hdp", "replicate")],
                d[, c("batch", "primary_hdp", "replicate")])
  expect_true(all(is.finite(as.matrix(result[, c("first", "second")]))))
  expect_lt(abs(mean(result$first[d$batch == "B"]) - mean(result$first[d$batch == "A"])), .1)
  expect_equal(mean(result$first[d$primary_hdp == 1]) - mean(result$first[d$primary_hdp == 0]),
                2, tolerance = .05)
})
