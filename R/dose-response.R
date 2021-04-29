#' Generates noisy dose response data at a range of (log-)luxes
#'
#' Adds noise on the logit scale, resulting in more realistic (noisy) data for
#' a given patient.
#'
#' @param sigma a noise level (represents additive Gaussian noise on logit scale)
#' @param p1 the log of the ed50 lux
#' @param p2 a positive shape parameter
#' @param lux a vector of lux values which defaults to c(10, 30, 50, 100, 200, 400, 2000)
#' @return a tibble of paired lux and melatonin suppression values
#' @export
#' @seealso [noise_logit]
#'
#' @examples
#' logistic_noise(1, 1.3, 4.5)
logistic_noise <- function(sigma, p1, p2,
                           lux=c(10, 30, 50, 100, 200, 400, 2000)) {
  y <- purrr::map_dbl(lux, ~logistic_2(., p1, p2))
  y_noise <- purrr::map_dbl(y, ~noise_logit(., sigma))
  df <- dplyr::tibble(lux=lux, y=y_noise)
  df
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
  df
}

#' Generates an individual with reasonable ed25 and ed75 values.
#'
#' @param thresh_25 lower bound on simulated ed25 vs observed ed25 (bound calculated as thresh_25 * observed)
#' @param thresh_75 upper bound on simulated ed75 vs observed ed75 (bound calculated as thresh_75 * observed)
#' @param eds_25 vector of ed25s from estimates from Phillips et al., (2017)
#' @param eds_75 vector of ed75s from estimates from Phillips et al., (2017)
#' @param alpha posteriors draws for the intercept in the regression of log_p2 on p1
#' @param beta posterior draws for slope parameter in the regression of log_p2 on p1
#' @param sigma0 posterior draws for the constant noise term in the regression of log_p2 on p1
#' @param sigma1 posterior draws for the heteroscedastic noise term in the regression of log_p2 on p1
#' @param cdf_inv inverse cumulative density function for distribution of p1 observed in Phillips et al., (2017)
#'
#' @return a tibble containing dose-response parameters for an individual
valid_individual <- function(thresh_25, thresh_75, eds_25, eds_75,
                             alpha, beta, sigma0, sigma1, cdf_inv) {
  lower <- thresh_25 * min(eds_25)
  upper <- thresh_75 * max(eds_75)
  ed_25_sim <- lower - 1
  ed_75_sim <- upper + 1
  while(((ed_25_sim < lower) | (ed_75_sim > upper))) {
    indiv <- sample_p1_p2(1, alpha, beta, sigma0, sigma1, cdf_inv)
    ed_25_sim <- ed(0.25, p1 = indiv$p1[1], p2 = indiv$p2[1])
    ed_75_sim <- ed(0.75, p1 = indiv$p1[1], p2 = indiv$p2[1])
  }
  indiv
}

#' Generates a population of individual dose-response curves
#'
#' @param n number of individual dose-response curves to generate
#' @param thresh_25 lower bound on simulated ed25 vs observed ed25 (bound calculated as thresh_25 * observed)
#' @param thresh_75 upper bound on simulated ed75 vs observed ed75 (bound calculated as thresh_75 * observed)
#'
#' @return a tibble with n rows with pairs of logistic_2 dose-response parameters
#' @export
virtual_population <- function(n, thresh_25, thresh_75) {

  # samples of parameter values
  alpha <- chronodoseresponse::p1_p2_regression_draws$alpha
  beta <- chronodoseresponse::p1_p2_regression_draws$beta
  sigma0 <- chronodoseresponse::p1_p2_regression_draws$sigma0
  sigma1 <- chronodoseresponse::p1_p2_regression_draws$sigma1

  p1 <- vector(length = n)
  p2 <- vector(length = n)
  for(i in seq_along(p1)) {
    indiv <- valid_individual(thresh_25, thresh_75,
                              chronodoseresponse::estimates$ed_25,
                              chronodoseresponse::estimates$ed_75,
                              alpha, beta, sigma0, sigma1, chronodoseresponse::cdf_inv)
    p1[i] <- indiv$p1[1]
    p2[i] <- indiv$p2[1]
  }
  dplyr::tibble(p1=p1, p2=p2)
}

#' Generates noisy data from a givn virtual population
#'
#' @param sigma a noise level
#' @param simulated_pop a tibble of logistic_2 parameters
#' @param lux a vector of lux values which defaults to c(10, 30, 50, 100, 200, 400, 2000)
#'
#' @return a tibble containing simulated noisy data
experimental_data <- function(sigma, simulated_pop, lux=c(10, 30, 50, 100, 200, 400, 2000)) {
  for(i in seq_along(simulated_pop$p1)) {
    temp <- logistic_noise(sigma, simulated_pop$p1[i], simulated_pop$p2[i], lux) %>%
      dplyr::mutate(id=i)
    if(i == 1)
      big_df <- temp
    else
      big_df <- big_df %>% dplyr::bind_rows(temp)
  }
  big_df
}

#' Samples a value of sigma for the logit noise process
#'
#' The sampled value is based on fitting a gamma distribution to
#' estimates of the noise values using RMSEs from Phillips et al. (2017).
#' The gamma distribution was fit using Stan, so estimates here incorporate
#' posterior uncertainty in the parameter values.
#'
#' @return a positive value
sample_sigma <- function() {
  a <- chronodoseresponse::sigma_fit_draws$a
  b <- chronodoseresponse::sigma_fit_draws$b
  idx <- sample(1:length(a), size=1)
  a_temp <- a[idx]
  b_temp <- b[idx]
  sigma <- stats::rgamma(1, a_temp, b_temp)
  sigma
}

#' Generates data from a virtual experiment measuring individual dose-response curves
#'
#' The model used to generate these experiments comprises two elements: a model representing
#' the underlying dose-response curves (which is based on a two parameter logistic); and a
#' model of typical experimental error in these measurements. The model was fit using estimates
#' presented in Phillips et al., (2017).
#'
#' @param n the number of individual dose-response curves to generate
#' @param lux a vector of lux values which defaults to c(10, 30, 50, 100, 200, 400, 2000)
#' @param thresh_25 lower bound on simulated ed25 vs observed ed25 (bound calculated as thresh_25 * observed) which defaults to 0.5
#' @param thresh_75 upper bound on simulated ed75 vs observed ed75 (bound calculated as thresh_75 * observed) which defaults to 1.5
#'
#' @return a tibble containing 'measured' dose-response melatonin curves at each lux value for each individual
#' @export
virtual_experiment <- function(n,
                               lux=c(10, 30, 50, 100, 200, 400, 2000),
                               thresh_25=0.5, thresh_75=1.5) {

  pop_df <- virtual_population(n, thresh_25, thresh_75)
  for(i in 1:nrow(pop_df)) {
    p1_temp <- pop_df$p1[i]
    p2_temp <- pop_df$p2[i]
    sigma <- sample_sigma()
    temp <- logistic_noise(sigma, p1_temp, p2_temp, lux) %>%
      dplyr::mutate(id=i,
                    sigma=sigma,
                    p1=p1_temp,
                    p2=p2_temp)
    if(i == 1)
      big_df <- temp
    else
      big_df <- big_df %>% dplyr::bind_rows(temp)
  }
  big_df
}
