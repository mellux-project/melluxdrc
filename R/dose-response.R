#' Generates noisy dose response data at a range of (log-)luxes
#'
#' Adds noise on the logit scale, resulting in more realistic (noisy) data for
#' a given patient.
#'
#' @param sigma a noise level (represents additive Gaussian noise on logit scale)
#' @param p1 the log of the ed50 lux
#' @param p2 a positive shape parameter
#' @param log_lux a vector of log lux values which defaults to log10(c(10, 30, 50, 100, 200, 400, 2000))
#' @return a tibble of paired lux and melatonin suppression values
#' @export
#' @seealso [noise_logit]
#'
#' @examples
#' logistic_noise(1, 1.3, 4.5)
logistic_noise <- function(sigma, p1, p2,
                           log_lux=log10(c(10, 30, 50, 100, 200, 400, 2000))) {
  y <- purrr::map_dbl(log_lux, ~logistic_2(., p1, p2))
  y_noise <- purrr::map_dbl(y, ~noise_logit(., sigma))
  df <- dplyr::tibble(lux=10^log_lux, y=y_noise)
  return(df)
}

#' Samples values of pairs of logistic_2 dose-response parameters using models fit to Phillips et al. (2017) estimates
#'
#' Creates a set of dose-response curves -- one per individual -- supposed to mimic participants in an experiment.
#' This is based on Bayesian models fit to parameter estimates from Phillips et al., (2017).
#'
#' @param n the number of individual dose-response curves to generate
#' @param alpha posteriors draws for the intercept in the regression of log_p2 on p1
#' @param beta posterior draws for slope parameter in the regression of log_p2 on p1
#' @param sigma0 posterior draws for the constant noise term in the regression of log_p2 on p1
#' @param sigma1 posterior draws for the heteroscedastic noise term in the regression of log_p2 on p1
#' @param cdf_inv inverse cumulative density function for distribution of p1 observed in Phillips et al., (2017)
#'
#' @return a tibble of p1 and p2 values for each individual dose-response curve
#' @importFrom magrittr "%>%"
sample_p1_p2 <- function(n, alpha, beta, sigma0, sigma1, cdf_inv) {
  ndraws <- length(alpha)
  idx <- sample(ndraws, n, replace = T)
  p1 <- vector(length = n)
  p2_log <- vector(length = n)
  for(i in 1:n) {
    a_idx <- idx[i]
    a_alpha <- alpha[a_idx]
    a_beta <- beta[a_idx]
    a_sigma0 <- sigma0[a_idx]
    a_sigma1 <- sigma1[a_idx]
    p1[i] <- f_sample_n(n=1, cdf_inv)
    p2_log[i] <- stats::rnorm(1, a_alpha + a_beta * p1[i], a_sigma0 + a_sigma1 * p1[i])
  }
  df <- dplyr::tibble(p1=p1, p2_log=p2_log) %>%
    dplyr::mutate(p2=10^p2_log) %>%
    dplyr::select(-p2_log)
  return(df)
}

#' Generates an individual with reasonable ed25 and ed75 values.
#'
#' @param thresh_25 lower bound on simulated ed25 vs observed ed25 (bound calculated as thresh_25 * observed)
#' @param thresh_75 upper bound on simulated ed75 vs observed ed75 (bound calculated as thresh_75 * observed)
#' @param eds_25
#' @param eds_75
#' @param alpha posteriors draws for the intercept in the regression of log_p2 on p1
#' @param beta posterior draws for slope parameter in the regression of log_p2 on p1
#' @param sigma0 posterior draws for the constant noise term in the regression of log_p2 on p1
#' @param sigma1 posterior draws for the heteroscedastic noise term in the regression of log_p2 on p1
#' @param cdf_inv inverse cumulative density function for distribution of p1 observed in Phillips et al., (2017)
#'
#' @return
valid_individual <- function(thresh_25, thresh_75, eds_25, eds_75, alpha, beta, sigma0, sigma1, cdf_inv) {
  lower <- thresh_25 * min(eds_25)
  upper <- thresh_75 * max(eds_75)
  ed_25_sim <- lower - 1
  ed_75_sim <- upper + 1
  while(((ed_25_sim < lower) | (ed_75_sim > upper))) {
    indiv <- sample_p1_p2(1, alpha, beta, sigma0, sigma1, cdf_inv)
    ed_25_sim <- ed(0.25, p1 = indiv$p1[1], p2 = indiv$p2[1])
    ed_75_sim <- ed(0.75, p1 = indiv$p1[1], p2 = indiv$p2[1])
  }
  return(indiv)
}

#' Generates a population of individual dose-response curves
#'
#' @param n number of individual dose-response curves to generate
#' @param thresh_25 lower bound on simulated ed25 vs observed ed25 (bound calculated as thresh_25 * observed)
#' @param thresh_75 upper bound on simulated ed75 vs observed ed75 (bound calculated as thresh_75 * observed)
#'
#' @return a tibble with n rows with pairs of logistic_2 dose-response parameters
#' @export
virtual_population_generator <- function(n, thresh_25, thresh_75) {

  # samples of parameter values
  alpha <- chronodoseresponse:::p1_p2_regression_draws$alpha
  beta <- chronodoseresponse:::p1_p2_regression_draws$beta
  sigma0 <- chronodoseresponse:::p1_p2_regression_draws$sigma0
  sigma1 <- chronodoseresponse:::p1_p2_regression_draws$sigma1

  p1 <- vector(length = n)
  p2 <- vector(length = n)
  for(i in seq_along(p1)) {
    indiv <- valid_individual(thresh_25, thresh_75,
                              chronodoseresponse:::eds$eds_25,
                              chronodoseresponse:::eds$eds_75,
                              alpha, beta, sigma0, sigma1, cdf_inv)
    p1[i] <- indiv$p1[1]
    p2[i] <- indiv$p2[1]
  }
  return(dplyr::tibble(p1=p1, p2=p2))
}
