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

test_that("comparison_test yields reasonable test results", {

  lux_1 <- 10
  lux_2 <- 30
  n <- 20
  nreps <- 20

  # between
  is_between <- TRUE
  results <- comparison_test(is_between, lux_1, lux_2, n, population_df, nreps=nreps)
  expect_equal(mean(results$result %in% c(0, 1)), 1)
  expect_equal(nrow(results), nreps)
  cnames <- colnames(results)
  cnames_expected <- c("replicate", "result", "p_value")
  expect_equal(mean(cnames %in% cnames_expected), 1)

  # suppose should be easier to detect a difference between more widely spaced
  # lux values
  results_short <- comparison_test(is_between, lux_1, lux_2, n, population_df, nreps=nreps)
  results_wide <- comparison_test(is_between, lux_1, 2000, n, population_df, nreps=nreps)
  expect_true(mean(results_short$result) <= mean(results_wide$result))

  # reverse order of luxes and same result should hold
  results_short <- comparison_test(is_between, lux_2, lux_1, n, population_df, nreps=nreps)
  results_wide <- comparison_test(is_between, 2000, lux_1, n, population_df, nreps=nreps)
  expect_true(mean(results_short$result) <= mean(results_wide$result))

  # within
  is_between <- F
  results <- comparison_test(is_between, lux_1, lux_2, n, population_df, nreps=nreps)
  expect_equal(mean(results$result %in% c(0, 1)), 1)

  # suppose should be easier to detect a difference between more widely spaced
  # lux values
  n <- 5
  results_short <- comparison_test(is_between, lux_1, lux_2, n, population_df, nreps=nreps)
  results_wide <- comparison_test(is_between, lux_1, 2000, n, population_df, nreps=nreps)
  expect_true(mean(results_short$result) < mean(results_wide$result))

  # reverse order of luxes and same result should hold
  results_short <- comparison_test(is_between, lux_2, lux_1, n, population_df, nreps=nreps)
  results_wide <- comparison_test(is_between, 2000, lux_1, n, population_df, nreps=nreps)
  expect_true(mean(results_short$result) < mean(results_wide$result))

  # comparing within and between
  n <- 10
  results_b <- comparison_test(T, lux_1, lux_2, n, population_df, nreps=nreps)
  results_w <- comparison_test(F, lux_1, lux_2, n, population_df, nreps=nreps)

  expect_true(mean(results_w$result) >= mean(results_b$result))
})

