test_that("logistic_noise returns melatonin values between [0, 1]", {
  df <- logistic_noise(0.4, runif(1, 0, 5), runif(1, 0, 5))
  for(i in seq_along(df$lux)) {
    expect_true(df$y[i] >= 0 & df$y[i] <= 1)
  }
})

test_that("logistic_noise returns data of correct length", {
  luxes <- 1:3
  df <- logistic_noise(0.4, runif(1, 0, 5), runif(1, 0, 5),
                       log_lux=luxes)
  expect_equal(nrow(df), length(luxes))
})

test_that("logistic_noise shifts if ed50 changes", {
  df <- logistic_noise(0.001, 1, 5,
                       log_lux=2)
  df1 <- logistic_noise(0.001, 2, 5,
                        log_lux=2)
  expect_true(df1$y < df$y)
})

test_that("logistic_noise: noise increases as sigma increases", {
  ed50 <- 100
  p1 <- log10(ed50)
  df <- logistic_noise(0.0, p1, 5,
                       log_lux=2)
  df1 <- logistic_noise(1.0, p1, 5,
                        log_lux=2)
  expect_true(abs(df1$y - 0.5) > abs(df$y - 0.5))
})

test_that("logistic_noise: zero noise reduces to logistic curves", {
  ed50 <- 100
  p1 <- log10(ed50)
  p2 <- 5
  log_lux <- seq(1, 3, length.out = 10)
  df <- logistic_noise(0.0, p1, p2,
                       log_lux=log_lux)
  untrue <- 0
  for(i in seq_along(df$y)) {
    untrue <- untrue + dplyr::if_else(abs(df$y[i] - logistic_2(log_lux[i], p1, p2)) < 0.00001, 0, 1)
  }
  expect_equal(untrue, 0)
})
