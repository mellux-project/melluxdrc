population_df <- virtual_experiment(200)
population_df_treat_lowered50 <- virtual_experiment(200,
                                                    treated_ed50_multiplier=0.2) %>%
  dplyr::mutate(treated=TRUE)
population_df_treat_lowered50_1 <- virtual_experiment(200,
                                                    treated_ed50_multiplier=0.8) %>%
  dplyr::mutate(treated=TRUE)
population_df_treat_highered50 <- virtual_experiment(200,
                                                     treated_ed50_multiplier=5) %>%
  dplyr::mutate(treated=TRUE)

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
  results_short <- comparison_test(is_between, lux_1, lux_2, n, population_df, nreps=nreps) %>%
    dplyr::mutate(result=dplyr::if_else((result==1) & (p_value < 0.05), 1, 0))
  results_wide <- comparison_test(is_between, lux_1, 2000, n, population_df, nreps=nreps) %>%
    dplyr::mutate(result=dplyr::if_else((result==1) & (p_value < 0.05), 1, 0))
  expect_true(mean(results_short$result) <= mean(results_wide$result))

  # reverse order of luxes and same result should hold
  results_short <- comparison_test(is_between, lux_2, lux_1, n, population_df, nreps=nreps) %>%
    dplyr::mutate(result=dplyr::if_else((result==1) & (p_value < 0.05), 1, 0))
  results_wide <- comparison_test(is_between, 2000, lux_1, n, population_df, nreps=nreps) %>%
    dplyr::mutate(result=dplyr::if_else((result==1) & (p_value < 0.05), 1, 0))
  expect_true(mean(results_short$result) <= mean(results_wide$result))

  # within
  is_between <- F
  results <- comparison_test(is_between, lux_1, lux_2, n, population_df, nreps=nreps)
  expect_equal(mean(results$result %in% c(0, 1)), 1)

  # suppose should be easier to detect a difference between more widely spaced
  # lux values
  n <- 5
  results_short <- comparison_test(is_between, lux_1, lux_2, n, population_df, nreps=nreps) %>%
    dplyr::mutate(result=dplyr::if_else((result==1) & (p_value < 0.05), 1, 0))
  results_wide <- comparison_test(is_between, lux_1, 2000, n, population_df, nreps=nreps) %>%
    dplyr::mutate(result=dplyr::if_else((result==1) & (p_value < 0.05), 1, 0))
  expect_true(mean(results_short$result) < mean(results_wide$result))

  # reverse order of luxes and same result should hold
  results_short <- comparison_test(is_between, lux_2, lux_1, n, population_df, nreps=nreps) %>%
    dplyr::mutate(result=dplyr::if_else((result==1) & (p_value < 0.05), 1, 0))
  results_wide <- comparison_test(is_between, 2000, lux_1, n, population_df, nreps=nreps) %>%
    dplyr::mutate(result=dplyr::if_else((result==1) & (p_value < 0.05), 1, 0))
  expect_true(mean(results_short$result) < mean(results_wide$result))

  # comparing within and between
  n <- 10
  results_b <- comparison_test(T, lux_1, lux_2, n, population_df, nreps=nreps) %>%
    dplyr::mutate(result=dplyr::if_else((result==1) & (p_value < 0.05), 1, 0))
  results_w <- comparison_test(F, lux_1, lux_2, n, population_df, nreps=nreps) %>%
    dplyr::mutate(result=dplyr::if_else((result==1) & (p_value < 0.05), 1, 0))

  expect_true(mean(results_w$result) >= mean(results_b$result))
})

test_that("generate_two_samples_at_one_lux generates reasonable samples", {
  lux <- 10
  n <- 10

  is_between <- TRUE
  pop_df <- virtual_treatment_experiment(n = 50, treated_ed50_multiplier=0.1,
                                         is_between=is_between)
  vals <- generate_two_samples_at_one_lux(is_between, lux, n, pop_df)
  expect_equal(length(vals$untreated), n)
  expect_equal(length(vals$treated), n)

  is_between <- FALSE
  pop_df <- virtual_treatment_experiment(n = 50, treated_ed50_multiplier=0.1,
                                         is_between=is_between)
  vals <- generate_two_samples_at_one_lux(is_between, lux, n, pop_df)
  expect_equal(length(vals$untreated), n)
  expect_equal(length(vals$treated), n)
})

test_that("comparison_test_treatment works as desired" , {
  population_df_treat_lowered50 <- population_df_treat_lowered50 %>%
    dplyr::bind_rows(population_df %>%
                     dplyr::mutate(treated=FALSE))
  population_df_treat_lowered50_1 <- population_df_treat_lowered50_1 %>%
    dplyr::bind_rows(population_df %>%
                     dplyr::mutate(treated=FALSE))
  population_df_treat_highered50 <- population_df_treat_highered50 %>%
    dplyr::bind_rows(population_df %>%
                     dplyr::mutate(treated=FALSE))

  # when the ed50 is lower the response should be higher at a given lux
  nreps <- 20
  lux <- 10
  n <- 10
  is_between <- TRUE
  is_treated_higher <- TRUE
  results <- comparison_test_treatment(
    is_between, lux, n, population_df_treat_lowered50,
    is_treated_higher, nreps=nreps)
  expect_equal(mean(results$result %in% c(0, 1)), 1)
  expect_equal(nrow(results), nreps)
  cnames <- colnames(results)
  cnames_expected <- c("replicate", "result", "p_value")
  expect_equal(mean(cnames %in% cnames_expected), 1)

  # mistakenly flip sign and check that test results change
  results <- results %>%
    dplyr::mutate(result=dplyr::if_else((result==1) & (p_value < 0.05), 1, 0))
  is_treated_higher <- FALSE
  results_wrong <- comparison_test_treatment(
    is_between, lux, n, population_df_treat_lowered50,
    is_treated_higher, nreps=nreps) %>%
    dplyr::mutate(result=dplyr::if_else((result==1) & (p_value < 0.05), 1, 0))
  expect_true(mean(results$result) > mean(results_wrong$result))

  # same but for case where treatment is weaker
  is_treated_higher <- TRUE
  results_weak <- comparison_test_treatment(
    is_between, lux, n, population_df_treat_lowered50_1,
    is_treated_higher, nreps=nreps) %>%
    dplyr::mutate(result=dplyr::if_else((result==1) & (p_value < 0.05), 1, 0))
  expect_true(mean(results$result) > mean(results_weak$result))

  # when the ed50 is higher the response should be lower at a given lux
  is_treated_higher <- FALSE
  results_high <- comparison_test_treatment(
    is_between, lux, n, population_df_treat_highered50,
    is_treated_higher, nreps=nreps) %>%
    dplyr::mutate(result=dplyr::if_else((result==1) & (p_value < 0.05), 1, 0))
  expect_true(mean(results_high$result) > 0)
})
