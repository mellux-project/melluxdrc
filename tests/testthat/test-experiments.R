population_df <- virtual_experiment(200)

test_that("generate_two_samples generates appropriate samples", {
  lux_1 <- 100
  lux_2 <- 200
  nsample <- 20
  vals_b <- generate_two_samples(T, lux_1, lux_2, nsample, population_df)
  expect_equal(length(vals_b$vals_1), nsample)
  expect_equal(length(vals_b$vals_2), nsample)
  vals_w <- generate_two_samples(F, lux_1, lux_2, nsample, population_df)
  expect_equal(length(vals_w$vals_1), nsample)
  expect_equal(length(vals_w$vals_2), nsample)

  diff_b <- vals_b$vals_2-vals_b$vals_1
  diff_w <- vals_w$vals_2-vals_w$vals_1
  expect_true(sd(diff_w) < sd(diff_b))
})

test_that("comparison_experiment yields reasonable test results", {

  lux_1 <- 10
  lux_2 <- 30
  n <- 20

  # between
  is_between <- TRUE
  expect_true(comparison_experiment(is_between, lux_1, lux_2, n, population_df) %in% c(0, 1))

  # suppose should be easier to detect a difference between more widely spaced
  # lux values
  nreps <- 20
  results_short <- purrr::map_dbl(1:nreps, ~comparison_experiment(is_between, lux_1, lux_2, n, population_df))
  results_wide <- purrr::map_dbl(1:nreps, ~comparison_experiment(is_between, lux_1, 2000, n, population_df))
  expect_true(mean(results_short) <= mean(results_wide))

  # reverse order of luxes and same result should hold
  results_short <- purrr::map_dbl(1:nreps, ~comparison_experiment(is_between, lux_2, lux_1, n, population_df))
  results_wide <- purrr::map_dbl(1:nreps, ~comparison_experiment(is_between, 2000, lux_1, n, population_df))
  expect_true(mean(results_short) <= mean(results_wide))

  # within
  is_between <- F
  expect_true(comparison_experiment(is_between, lux_1, lux_2, n, population_df) %in% c(0, 1))

  # suppose should be easier to detect a difference between more widely spaced
  # lux values
  n <- 5
  results_short <- purrr::map_dbl(1:nreps, ~comparison_experiment(is_between, lux_1, lux_2, n, population_df))
  results_wide <- purrr::map_dbl(1:nreps, ~comparison_experiment(is_between, lux_1, 2000, n, population_df))
  expect_true(mean(results_short) < mean(results_wide))

  # reverse order of luxes and same result should hold
  results_short <- purrr::map_dbl(1:nreps, ~comparison_experiment(is_between, lux_2, lux_1, n, population_df))
  results_wide <- purrr::map_dbl(1:nreps, ~comparison_experiment(is_between, 2000, lux_1, n, population_df))
  expect_true(mean(results_short) < mean(results_wide))

  # comparing within and between
  n <- 10
  results_b <- purrr::map_dbl(1:nreps, ~comparison_experiment(T, lux_1, lux_2, n, population_df))
  results_w <- purrr::map_dbl(1:nreps, ~comparison_experiment(F, lux_1, lux_2, n, population_df))

  expect_true(mean(results_w) >= mean(results_b))
})

