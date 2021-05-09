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

#' Constructed weighted inverse-cdf from a mixture of a cdf and a normal-cdf
#'
#' Takes the cdf of a distribution and combines it with a cdf of a N(mu, sigma)
#' in a weighted sum: cdf_mixed(x) = cdf(x) weight + normal_cdf(x|mu,sigma) (1-weight).
#'
#' @param weight a value (0<=weight<=1) specifying weight given to cdf_full
#' @param cdf_full the cdf forming first component of mixture
#' @param mu the mean of the normal distribution component
#' @param sigma the sd of the normal distribution component
#' @param p1_vals a set of p1 values used to form mixture cdf
#'
#' @return a function yielding the inverse cdf for the mixture distribution
cdf_inv_constructor <- function(weight, cdf_full, mu, sigma, p1_vals) {
  cdf_mixed <- function(x) cdf_full(x) * weight + stats::pnorm(x, mu, sigma) * (1 - weight)
  cdf_vals <- purrr::map_dbl(p1_vals, ~cdf_mixed(.))
  cdf_mixed_inv <- stats::approxfun(cdf_vals, p1_vals)
  cdf_mixed_inv
}

#' Converts p1 assuming treatment that changes ed50 by a multiplier
#'
#' @param multiplier indicates new ed50 = old ed50 * multiplier
#' @param old_p1 p1 value before treatment
#'
#' @return a new p1 value
treated_p1 <- function(multiplier, old_p1) {
  new_p1 <- log10(multiplier) + old_p1
  new_p1
}

#' Samples values of pairs of logistic_2 dose-response parameters using models fit to Phillips et al. (2017) estimates
#'
#' Creates a set of dose-response curves -- one per individual -- supposed to mimic participants in an experiment.
#' This is based on Bayesian models fit to parameter estimates from Phillips et al., (2017).
#'
#' @param n the number of individual dose-response curves to generate
#' @param p1_distribution_parameters a named list of items used to characterise the p1 distribution:
#' cdf_full, cdf_inv_full, normal_sigma, weight_full. Here cdf_inv_full and cdf_full are the cdf and inverse-cdf
#' characterising the distribution of p1 distribution from Phillips et al. (2017); normal_sigma is the standard
#' deviation of the normal distribution used as the narrow mixture component; weight_full is the weight given to the
#' p1_distribution derived from Phillips et al. (2017).
#' @param p2_distribution_parameters a named list of items used to characterise the log_p2|p1 distribution. It comprises:
#' alpha posteriors draws for the intercept in the regression of log_p2 on p1; beta posterior draws for slope
#' parameter in the regression of log_p2 on p1; sigma0 posterior draws for the constant noise term in the regression of log_p2 on p1;
#' posterior draws for the heteroscedastic noise term in the regression of log_p2 on p1.
#'
#' @return a tibble of p1 and p2 values for each individual dose-response curve
#' @importFrom magrittr "%>%"
sample_p1_p2 <- function(n, p1_distribution_parameters, p2_distribution_parameters) {

  # make p1 inverse-cdf function
  cdf_inv_full <- p1_distribution_parameters$cdf_inv_full
  cdf_full <- p1_distribution_parameters$cdf_full
  normal_sigma <- p1_distribution_parameters$normal_sigma
  weight_full <- p1_distribution_parameters$weight

  middle_p1 <- cdf_inv_full(0.5)
  p1_fine <- seq(0, 3.1, 0.01)
  cdf_inv_mixed <- cdf_inv_constructor(weight_full, cdf_full,
                                       middle_p1, normal_sigma,
                                       p1_fine)

  alpha <- p2_distribution_parameters$alpha
  beta <- p2_distribution_parameters$beta
  sigma0 <- p2_distribution_parameters$sigma0
  sigma1 <- p2_distribution_parameters$sigma1
  weight_sigma <- p2_distribution_parameters$weight

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
    p1[i] <- f_sample_n(n=1, cdf_inv_mixed)
    sigma_noise <- weight_sigma * (a_sigma0 + a_sigma1 * p1[i])
    p2_log[i] <- stats::rnorm(1, a_alpha + a_beta * p1[i], sigma_noise)
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
#' @param eds_25 vector of ed25s from estimates from Phillips et al. (2017)
#' @param eds_75 vector of ed75s from estimates from Phillips et al. (2017)
#' @inheritParams sample_p1_p2
#' @return a tibble containing dose-response parameters for an individual
valid_individual <- function(thresh_25, thresh_75, eds_25, eds_75,
                             p1_distribution_parameters, p2_distribution_parameters) {
  lower <- thresh_25 * min(eds_25)
  upper <- thresh_75 * max(eds_75)
  ed_25_sim <- lower - 1
  ed_75_sim <- upper + 1
  while(((ed_25_sim < lower) | (ed_75_sim > upper))) {
    indiv <- sample_p1_p2(1, p1_distribution_parameters, p2_distribution_parameters)
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
#' @param weight_p1 the weight (0<=weight<=1) given to the purely empirical distribution of p1 derived from Phillips et al. (2017)
#' @param weight_sigma the weight (0<=weight<=1) specifying factor to suppress the estimated standard deviation of p2_log|p1
#' @param normal_sigma the standard deviation of the normal mixture distrution used to characterise p1
#'
#' @return a tibble with n rows with pairs of logistic_2 dose-response parameters
#' @export
virtual_population <- function(n, thresh_25, thresh_75, weight_p1, weight_sigma, normal_sigma) {

  p1_distribution_parameters<- list(cdf_full = chronodoseresponse::cdf,
                                    cdf_inv_full = chronodoseresponse::cdf_inv,
                                    normal_sigma = normal_sigma,
                                    weight = weight_p1)

  p2_distribution_parameters <- chronodoseresponse::p1_p2_regression_draws
  p2_distribution_parameters$weight <- weight_sigma

  p1 <- vector(length = n)
  p2 <- vector(length = n)
  for(i in seq_along(p1)) {
    indiv <- valid_individual(thresh_25, thresh_75,
                              chronodoseresponse::estimates$ed_25,
                              chronodoseresponse::estimates$ed_75,
                              p1_distribution_parameters, p2_distribution_parameters)
    p1[i] <- indiv$p1[1]
    p2[i] <- indiv$p2[1]
  }
  dplyr::tibble(p1=p1, p2=p2)
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
#' @inheritParams virtual_population
#' @inheritParams logistic_noise
#' @param individual_variation_reduction a value (0<=value<=1) which reduces individual variability in dose-response curves
#' @param treated_ed50_multiplier a value which multiplies natural ed50 to result in a treated ed50 = natural ed50 * treated_ed50_multiplier
#' @return a tibble containing 'measured' dose-response melatonin curves at each lux value for each individual.
#' The tibble also contains the simulated natural p1 and p2 values and treated p1 value (which may be the same as the
#' natural one if the individual is untreated), as well as a Boolean indicating if a patient is treated. There is
#' also a Boolean indicating if an individual's treated p1 value was manually adjusted to: this likely happens
#' if applying a treatment that is too extreme.
#'
#' @export
virtual_experiment <- function(n,
                               lux=c(10, 30, 50, 100, 200, 400, 2000),
                               thresh_25=0.5, thresh_75=1.5,
                               individual_variation_reduction=1,
                               treated_ed50_multiplier=1) {
  if(individual_variation_reduction > 1)
    stop("individual_variation_reduction must be < 1.")
  weight_p1 <- individual_variation_reduction
  weight_sigma <- individual_variation_reduction
  normal_sigma <- 0.05 # hard code as unlikely a user would want to change
  pop_df <- virtual_population(n, thresh_25, thresh_75, weight_p1, weight_sigma, normal_sigma)
  for(i in 1:nrow(pop_df)) {
    p1_natural <- pop_df$p1[i]
    p1_temp <- treated_p1(treated_ed50_multiplier, p1_natural)
    p1_truncated <- FALSE
    if(p1_temp < 0) {
      p1_temp <- 0
      p1_truncated <- TRUE
    }
    p2_temp <- pop_df$p2[i]
    sigma <- sample_sigma()
    temp <- logistic_noise(sigma, p1_temp, p2_temp, lux) %>%
      dplyr::mutate(id=i,
                    sigma=sigma,
                    p1=p1_natural,
                    p2=p2_temp,
                    p1_treated=p1_temp,
                    p1_truncated=p1_truncated)
    if(i == 1)
      big_df <- temp
    else
      big_df <- big_df %>% dplyr::bind_rows(temp)
  }
  big_df
}


#' Generates data from a virtual experiment where individuals are treated and repeatedly measured
#'
#' The model used to generate these experiments comprises two elements: a model representing
#' the underlying dose-response curves (which is based on a two parameter logistic); and a
#' model of typical experimental error in these measurements. The model was fit using estimates
#' presented in Phillips et al., (2017).
#'
#' In this experiment setup, individuals have their dose-response relationship measured twice:
#' once before a treatment and once after it.
#'
#' @inheritParams virtual_population
#' @inheritParams logistic_noise
#' @param individual_variation_reduction a value (0<=value<=1) which reduces individual variability in dose-response curves
#' @param treated_ed50_multiplier a value which multiplies natural ed50 to result in a treated ed50 = natural ed50 * treated_ed50_multiplier
#' @return a tibble containing 'measured' dose-response melatonin curves at each lux value for each individual twice: before and after treatment.
#' The tibble also contains the simulated natural p1 and p2 values and treated p1 value (which may be the same as the
#' natural one if the individual is untreated), as well as a Boolean indicating if a patient is treated. There is
#' also a Boolean indicating if an individual's treated p1 value was manually adjusted to: this likely happens
#' if applying a treatment that is too extreme.
#' @export
virtual_within_treatment_experiment <- function(n,
                               lux=c(10, 30, 50, 100, 200, 400, 2000),
                               thresh_25=0.5, thresh_75=1.5,
                               individual_variation_reduction=1,
                               treated_ed50_multiplier=1) {
  if(individual_variation_reduction > 1)
    stop("individual_variation_reduction must be < 1.")
  weight_p1 <- individual_variation_reduction
  weight_sigma <- individual_variation_reduction
  normal_sigma <- 0.05 # hard code as unlikely a user would want to change
  pop_df <- virtual_population(n, thresh_25, thresh_75, weight_p1, weight_sigma, normal_sigma)

  measurement_count <- 1
  for(i in 1:nrow(pop_df)) {
    p1_natural <- pop_df$p1[i]
    p2_temp <- pop_df$p2[i]
    sigma <- sample_sigma()

    for(j in 1:2) {
      if(j == 1) {
        treated = FALSE
        p1_temp <- p1_natural
      } else {
        treated = TRUE
        p1_temp <- treated_p1(treated_ed50_multiplier, p1_natural)
      }
      p1_truncated <- FALSE
      if(p1_temp < 0) {
        p1_temp <- 0
        p1_truncated <- TRUE
      }
      temp <- logistic_noise(sigma, p1_temp, p2_temp, lux) %>%
        dplyr::mutate(id=i,
                      sigma=sigma,
                      p1=p1_natural,
                      p2=p2_temp,
                      p1_treated=p1_temp,
                      p1_truncated=p1_truncated,
                      treated=treated)
      if(measurement_count == 1)
        big_df <- temp
      else
        big_df <- big_df %>% dplyr::bind_rows(temp)

      measurement_count <- measurement_count + 1
    }
  }
  big_df
}
