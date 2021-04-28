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
