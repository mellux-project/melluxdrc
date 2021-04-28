test_that("logistic function returns correct values", {
  x <- 2
  p1 <- 2
  p2 <- 3
  expect_equal(0.5, logistic_2(x, p1, p2))
})

test_that("ed function returns correct values",{
  ## relies on ed50 = 10^p1
  expect_equal(ed(0.5, 1.3, 4.2), 10^1.3)
})

test_that("Logit function returns correct values", {
  expect_equal(logit(0.2), log(0.2 / 0.8))
  expect_equal(logit(0), -Inf)
  expect_equal(logit(1), Inf)
})

test_that("Logistic function returns correct values", {
  expect_equal(logistic(0), 0.5)
  expect_equal(logistic(100), 1)
  expect_equal(logistic(-100), 0)
})

test_that("noise_logit function returns values between [0, 1]", {
  ys <- runif(100)
  count <- 0
  for(i in seq_along(ys)) {
    y_noise <- noise_logit(ys[i], 0.3)
    count <- count + dplyr::if_else(y_noise >= 0 & y_noise <= 1, 0, 1)
  }
  expect_equal(count, 0)
})

test_that("f_sample_n returns values that appear uniform with uniform cdf_inv", {
  n <- 100000
  x <- f_sample_n(n, qunif)
  expect_true(abs(mean(x) - 0.5) < 0.1)
  expect_true(abs(var(x) - 1/12) < 0.1)
})

test_that("f_sample_n returns values that appear normal with normal cdf_inv", {
  n <- 100000
  x <- f_sample_n(n, qnorm)
  expect_true(abs(mean(x) - 0) < 0.1)
  expect_true(abs(var(x) - 1) < 0.1)
})
