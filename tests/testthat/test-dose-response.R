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


test_that("valid_individual returns individual with ed25 and ed75 within thresholds", {
  alpha <- chronodoseresponse::p1_p2_regression_draws$alpha
  beta <- chronodoseresponse::p1_p2_regression_draws$beta
  sigma0 <- chronodoseresponse::p1_p2_regression_draws$sigma0
  sigma1 <- chronodoseresponse::p1_p2_regression_draws$sigma1

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
                                alpha, beta, sigma0, sigma1, cdf_inv)
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
  df <- virtual_population(n, 0.1, 2)
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

  df <- virtual_population(n, lower_thres, upper_thres)
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
