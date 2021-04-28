test_that("logistic function returns correct values", {
  x <- 2
  p1 <- 2
  p2 <- 3
  expect_equal(0.5, logistic_2(x, p1, p2))
})
