#' Two-parameter dose-response logistic curve
#'
#' Yields a melatonin suppression response for a given level of log_lux
#'
#' @param log_lux a log lux value
#' @param p1 the log of the ed50 lux
#' @param p2 a positive shape parameter
#'
#' @return a melatonin suppression value between 0 and 1
#' @export
#'
#' @examples
#' # below returns 0.5 since log_lux=log_ed50
#' logistic_2(10, 10, 5)
logistic_2 <- function(log_lux, p1, p2) {
  1 - 1 / (1 + (log_lux / p1)^p2)
}

#' Calculates ed for a given quantile
#'
#' @param edx a quantile between 0 and 1
#' @param p1 the log of the ed50 lux
#' @param p2 a positive shape parameter
#'
#' @return yields a lux value
#' @export
#'
#' @examples
#' # at ed50 the returned value is 10^p1
#' ed(0.5, 1.3, 4.4) # = 10^1.3
ed <- function(edx, p1, p2) {
  10^(p1 * (-1 + 1 / (1 - edx))^(1 / p2))
}

#' Logit function
#'
#' @param y a value between 0 and 1
#'
#' @return an unbounded value
logit <- function(y) {
  return(log(y / (1 - y)))
}

#' Logistic function
#'
#' @param y an unbounded input
#'
#' @return a value between 0 and 1
logistic <- function(y) {
  return(1 / (1 + exp(-y)))
}

#' Adds noise to a data point on the logit scale
#'
#' Adds Gaussian noise to a point on the logit scale:
#' \deqn{logit(y_noise) = logit(y) + e}
#' where:
#' \deqn{e \sim normal(0, sigma)}
#'
#' Then we convert back to the normal scale via a logistic transform.
#'
#' @param y A melatonin suppression value between 0 and 1
#' @param sigma Scale of Gaussian noise
#'
#' @return A noisy version of y bound between 0 and 1
noise_logit <- function(y, sigma) {
  e <- stats::rnorm(1, 0, sigma)
  y_logit <- logit(y) + e
  logistic(y_logit)
}

#' Calculates cumulative density for a given pdf
#'
#' @param x a value
#' @param probability_density an integratable probability density function
#' @param lower a lower bound to integrate from
#'
#' @return a cdf value between 0 and 1
f_cdf <- function(x, probability_density, lower) {
  stats::integrate(probability_density, lower, x, subdivisions=2000)$value
}

#' Uses inverse transform sampling to generate an independent draw
#'
#' @param cdf_inv an inverse cumulative density function
#'
#' @return a draw from underlying probability density
f_sample_one <- function(cdf_inv) {
  u <- stats::runif(1)
  p1 <- cdf_inv(u)
  while(is.na(p1)) {
    u <- stats::runif(1)
    p1 <- cdf_inv(u)
  }
  p1
}

#' Generates n draws from underlying density using inverse transform sampling
#'
#' @param n the number of draws
#' @param cdf_inv the inverse cumulative density function
#'
#' @return a vector of length n of draws
f_sample_n <- function(n, cdf_inv) {
  purrr::map_dbl(1:n, ~f_sample_one(cdf_inv))
}

# creating cdf_inv so it can be used elsewhere
a <- KernSmooth::bkde(chronodoseresponse::estimates$p1)
f <- approxfun(a$x, a$y, method="constant")
p1_fine <- seq(0, 3.1, 0.01)
cdf <- purrr::map_dbl(p1_fine, ~f_cdf(., f, 0.0))
cdf_inv <- approxfun(cdf, p1_fine)
