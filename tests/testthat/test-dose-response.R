test_that("logistic_noise returns melatonin values between [0, 1]", {
  df <- logistic_noise(0.4, runif(1, 0, 5), runif(1, 0, 5))
  for(i in seq_along(df$lux)) {
    expect_true(df$y[i] >= 0 & df$y[i] <= 1)
  }
})

test_that("logistic_noise returns data of correct length", {
  luxes <- 10^(1:3)
  df <- logistic_noise(0.4, runif(1, 0, 5), runif(1, 0, 5),
                       lux=luxes)
  expect_equal(nrow(df), length(luxes))
})

test_that("logistic_noise shifts if ed50 changes", {
  df <- logistic_noise(0.001, 1, 5,
                       lux=10^2)
  df1 <- logistic_noise(0.001, 2, 5,
                        lux=10^2)
  expect_true(df1$y < df$y)
})

test_that("logistic_noise: noise increases as sigma increases", {
  ed50 <- 100
  p1 <- log10(ed50)
  df <- logistic_noise(0.0, p1, 5,
                       lux=10^2)
  df1 <- logistic_noise(1.0, p1, 5,
                        lux=10^2)
  expect_true(abs(df1$y - 0.5) > abs(df$y - 0.5))
})

test_that("logistic_noise: zero noise reduces to logistic curves", {
  ed50 <- 100
  p1 <- log10(ed50)
  p2 <- 5
  lux <- 10^seq(1, 3, length.out = 10)
  df <- logistic_noise(0.0, p1, p2,
                       lux=lux)
  untrue <- 0
  for(i in seq_along(df$y)) {
    untrue <- untrue + dplyr::if_else(abs(df$y[i] - logistic_2(lux[i], p1, p2)) < 0.00001, 0, 1)
  }
  expect_equal(untrue, 0)
})

test_that("treated_p1 results in proportional change in ed50", {
  fraction <- 0.75
  old_p1 <- 2
  old_ed50 <- ed(0.5, old_p1, 3)
  expect_equal(ed(0.5, treated_p1(fraction, old_p1), 3), old_ed50 * fraction)

  fraction <- 2
  old_p1 <- 2
  old_ed50 <- ed(0.5, old_p1, 3)
  expect_equal(ed(0.5, treated_p1(fraction, old_p1), 3), old_ed50 * fraction)

  fraction <- 1
  old_p1 <- 2
  old_ed50 <- ed(0.5, old_p1, 3)
  expect_equal(ed(0.5, treated_p1(fraction, old_p1), 3), old_ed50 * fraction)
})

test_that("sample_p1_p2 returns (p1, p2) distribution with/without reduced individual variation", {
  n <- 200

  weight <- 1
  p1_distribution_parameters <- list(cdf_full = chronodoseresponse::cdf,
                                     cdf_inv_full = chronodoseresponse::cdf_inv,
                                     normal_sigma = 0.05,
                                     weight = weight)

  p2_distribution_parameters <- chronodoseresponse::p1_p2_regression_draws
  p2_distribution_parameters$weight <- weight
  p1_p2_full_variation <- sample_p1_p2(n, p1_distribution_parameters, p2_distribution_parameters)
  expect_equal(nrow(p1_p2_full_variation), n)

  weight <- 0.2
  p1_distribution_parameters$weight <- weight
  p2_distribution_parameters$weight <- weight
  p1_p2_reduced <- sample_p1_p2(n, p1_distribution_parameters, p2_distribution_parameters)
  expect_true(sd(p1_p2_full_variation$p1) > sd(p1_p2_reduced$p1))
  expect_true(sd(p1_p2_full_variation$p2) > sd(p1_p2_reduced$p2))
})


test_that("valid_individual returns individual with ed25 and ed75 within thresholds", {

  p1_distribution_parameters <- list(cdf_full = chronodoseresponse::cdf,
                                    cdf_inv_full = chronodoseresponse::cdf_inv,
                                    normal_sigma = 0.05,
                                    weight = 1)

  p2_distribution_parameters <- chronodoseresponse::p1_p2_regression_draws
  p2_distribution_parameters$weight <- 0.8

  eds_25 <- estimates$ed_25
  eds_75 <- estimates$ed_75

  thresh_25 <- 0.8
  thresh_75 <- 1.5
  lower <- thresh_25 * min(eds_25)
  upper <- thresh_75 * max(eds_75)

  count_lower <- 0
  count_upper <- 0
  for(i in 1:10) {
    indiv <- valid_individual(thresh_25, thresh_75, eds_25, eds_75,
                              p1_distribution_parameters,
                              p2_distribution_parameters)
    ed_25_sim <- ed(0.25, p1 = indiv$p1[1], p2 = indiv$p2[1])
    ed_75_sim <- ed(0.75, p1 = indiv$p1[1], p2 = indiv$p2[1])
    count_lower <- count_lower + dplyr::if_else(ed_25_sim > lower, 0, 1)
    count_upper <- count_upper + dplyr::if_else(ed_75_sim < upper, 0, 1)
  }
  expect_equal(count_lower, 0)
  expect_equal(count_upper, 0)
})


test_that("virtual_population generates population of correct size", {
  n <- 10
  df <- virtual_population(n, 0.1, 2, 0.5, 0.5)
  expect_equal(nrow(df), n)
})

test_that("virtual_population generates individuals all with ed25s and ed75s within thresholds", {
  n <- 10
  lower_thres <- 0.5
  upper_thres <- 1.1
  eds_25 <- purrr::map_dbl(seq_along(estimates$p1), ~ed(0.25, estimates$p1[.], estimates$p2[.]))
  eds_75 <- purrr::map_dbl(seq_along(estimates$p1), ~ed(0.75, estimates$p1[.], estimates$p2[.]))
  lower <- lower_thres * min(eds_25)
  upper <- upper_thres * max(eds_75)

  df <- virtual_population(n, lower_thres, upper_thres, 1, 1)
  count_lower <- 0
  count_upper <- 0
  for(i in 1:nrow(df)) {
    ed_25_sim <- ed(0.25, p1 = df$p1[i], p2 = df$p2[i])
    ed_75_sim <- ed(0.75, p1 = df$p1[i], p2 = df$p2[i])
    count_lower <- count_lower + dplyr::if_else(ed_25_sim > lower, 0, 1)
    count_upper <- count_upper + dplyr::if_else(ed_75_sim < upper, 0, 1)
  }
  expect_equal(count_lower, 0)
  expect_equal(count_upper, 0)
})

test_that("sample_sigma produces reasonable sigma values", {
  n <- 2000
  sigmas <- vector(length = n)
  for(i in 1:n) {
    sigmas[i] <- sample_sigma()
  }
  a <- chronodoseresponse::sigma_fit_draws$a
  b <- chronodoseresponse::sigma_fit_draws$b
  mean_sigma <- mean(purrr::map2_dbl(a, b, ~.x/.y))
  expect_true(abs(mean(sigmas) - mean_sigma) < 0.1)
})


test_that("virtual_experiment produces populations with reasonable characteristics", {
  n <- 5
  lux <- 10^seq(1, 3, length.out = 4)
  pop_df <- virtual_experiment(n, thresh_25=0.25, thresh_75=1.5, lux=lux)
  expect_equal(dplyr::n_distinct(pop_df$id), n)
  expect_equal(dplyr::n_distinct(pop_df$lux), length(lux))
  expect_true(sum(unique(pop_df$lux) %in% lux) == length(lux))
  expect_equal(dplyr::n_distinct(pop_df$sigma), n)
  expect_equal(dplyr::n_distinct(pop_df$p1), n)
  expect_equal(dplyr::n_distinct(pop_df$p2), n)
})

test_that("virtual_experiment produces populations with reasonable suppression values", {
  # this is mainly just a test that our simulated distribution doesn't change over time
  pop_df <- virtual_experiment(n=200, lux=c(1, 10, 100, 1000)) %>%
    dplyr::group_by(lux) %>%
    dplyr::summarise(y=mean(y)) %>%
    tidyr::pivot_wider(names_from = lux,
                       values_from=y)
  expect_true(abs(pop_df$`1` - 0) < 0.05)
  expect_true(abs(pop_df$`10` - 0.23) < 0.2)
  expect_true(abs(pop_df$`100` - 0.71) < 0.2)
  expect_true(abs(pop_df$`1000` - 0.95) < 0.2)
})


test_that("virtual_experiment allows reduction of individual variation", {
  pop_df <- virtual_experiment(n=200, lux=c(100, 1000)) %>%
    dplyr::group_by(lux) %>%
    dplyr::summarise(y=sd(y)) %>%
    tidyr::pivot_wider(names_from = lux,
                       values_from=y) %>%
    dplyr::mutate(type="full")

  pop_df_reduced <- virtual_experiment(n=200, lux=c(100, 1000),
                                       individual_variation_level=0.1) %>%
    dplyr::group_by(lux) %>%
    dplyr::summarise(y=sd(y)) %>%
    tidyr::pivot_wider(names_from = lux,
                       values_from=y) %>%
    dplyr::mutate(type="reduced")
  pop_df_combined <- pop_df %>%
    dplyr::bind_rows(pop_df_reduced) %>%
    tidyr::pivot_longer(c(`100`, `1000`)) %>%
    tidyr::pivot_wider(id_cols = name,
                       names_from = type,
                       values_from=value)
  expect_equal(mean(pop_df_combined$reduced < pop_df_combined$full), 1)
})

test_that("virtual_experiment allows change in suppression from treatment", {
  pop_df <- virtual_experiment(n=200, lux=c(20, 100)) %>%
    dplyr::group_by(lux) %>%
    dplyr::summarise(y=median(y)) %>%
    tidyr::pivot_wider(names_from = lux,
                       values_from=y) %>%
    dplyr::mutate(type="full")

  pop_df_treated <- virtual_experiment(n=200, lux=c(20, 100),
                                       treated_ed50_multiplier=0.5) %>%
    dplyr::group_by(lux) %>%
    dplyr::summarise(y=median(y)) %>%
    tidyr::pivot_wider(names_from = lux,
                       values_from=y) %>%
    dplyr::mutate(type="treated")

  pop_df_combined <- pop_df %>%
    dplyr::bind_rows(pop_df_treated) %>%
    tidyr::pivot_longer(c(`20`, `100`)) %>%
    tidyr::pivot_wider(id_cols = name,
                       names_from = type,
                       values_from=value)
  expect_equal(mean(pop_df_combined$treated > pop_df_combined$full), 1)
})

test_that("virtual_experiment allows reduction of individual variation with treatment", {
  pop_df <- virtual_experiment(n=200, lux=c(100, 1000), treated_ed50_multiplier=0.5) %>%
    dplyr::group_by(lux) %>%
    dplyr::summarise(y=sd(y)) %>%
    tidyr::pivot_wider(names_from = lux,
                       values_from=y) %>%
    dplyr::mutate(type="full")

  pop_df_reduced <- virtual_experiment(n=200, lux=c(100, 1000),
                                       individual_variation_level=0.1,
                                       treated_ed50_multiplier=0.5) %>%
    dplyr::group_by(lux) %>%
    dplyr::summarise(y=sd(y)) %>%
    tidyr::pivot_wider(names_from = lux,
                       values_from=y) %>%
    dplyr::mutate(type="reduced")
  pop_df_combined <- pop_df %>%
    dplyr::bind_rows(pop_df_reduced) %>%
    tidyr::pivot_longer(c(`100`, `1000`)) %>%
    tidyr::pivot_wider(id_cols = name,
                       names_from = type,
                       values_from=value)
  expect_equal(mean(pop_df_combined$reduced < pop_df_combined$full), 1)
})

test_that("virtual_within_treatment_experiment measures before and after treatment", {
  nindiv <- 200
  test <- virtual_within_treatment_experiment(nindiv, treated_ed50_multiplier=0.5)
  n_lux <- dplyr::n_distinct(test$lux)
  expect_equal(nrow(test), 2 * nindiv * n_lux)
  expect_equal(sum(test$treated) / dplyr::n_distinct(test$lux), nindiv)

  df_summary <- test %>%
    dplyr::group_by(treated, lux) %>%
    dplyr::summarise(y=mean(y), .groups="drop") %>%
    tidyr::pivot_wider(id_cols=lux, names_from=treated,
                       values_from=y)
  expect_true(mean(df_summary$`TRUE` > df_summary$`FALSE`) > 0.7) # arbitrary cutoff
})

test_that("virtual_within_treatment_experiment works fine with individual variance reducer", {
  nindiv <- 200
  test <- virtual_within_treatment_experiment(nindiv, treated_ed50_multiplier=0.5) %>%
    dplyr::mutate(type="full")
  n_lux <- dplyr::n_distinct(test$lux)
  expect_equal(nrow(test), 2 * nindiv * n_lux)
  expect_equal(sum(test$treated) / dplyr::n_distinct(test$lux), nindiv)
  test_lower <- virtual_within_treatment_experiment(nindiv, treated_ed50_multiplier=0.5,
                                                    individual_variation_level=0.5) %>%
    dplyr::mutate(type="reduced")
  n_lux <- dplyr::n_distinct(test_lower$lux)
  expect_equal(nrow(test_lower), 2 * nindiv * n_lux)
  expect_equal(sum(test_lower$treated) / dplyr::n_distinct(test_lower$lux), nindiv)

  df_summary <- test %>%
    dplyr::bind_rows(test_lower) %>%
    dplyr::group_by(type, treated, lux) %>%
    dplyr::summarise(mu=mean(y), sigma=sd(y), .groups="drop") %>%
    tidyr::pivot_wider(id_cols=lux, names_from=c(type, treated),
                       values_from=c(mu, sigma))
  thresh <- 0.7 # arbitrary cutoff
  expect_true(mean(df_summary$mu_full_TRUE > df_summary$mu_full_FALSE) > thresh)
  expect_true(mean(df_summary$mu_reduced_TRUE > df_summary$mu_reduced_FALSE) > thresh)
  expect_true(mean(df_summary$sigma_full_TRUE > df_summary$sigma_reduced_TRUE) > thresh)
  expect_true(mean(df_summary$sigma_full_FALSE > df_summary$sigma_reduced_FALSE) > thresh)
})
