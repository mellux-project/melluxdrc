#' Estimates of parameter values from Phillips et al. (2017).
#'
#' A dataset containing parameter estimates of the two parameter
#' logistic curves for 41 individuals. Note, the data was provided
#' via personal communication with the paper authors.
#'
#' @references Phillips et al., (2017). High sensitivity and interindividual variability
#' in the response of the human circadian system to evening light. PNAS.
#'
#' @format A data frame with 41 and 8 variables:
#' \describe{
#'   \item{id}{individual identifier}
#'   \item{p1}{log_ed50 estimate for each individual}
#'   \item{p2}{logistic shape estimate for each individual}
#'   \item{p1_l}{5% confidence interval estimate for p1}
#'   \item{p1_u}{95% confidence interval estimate for p1}
#'   \item{p2_l}{5% confidence interval estimate for p2}
#'   \item{p2_u}{95% confidence interval estimate for p2}
#'   \item{rmse}{Root mean squared error for individual}
#'   \item{ed_25}{25% ed}
#'   \item{ed_75}{75% ed}
#'   ...
#' }
"estimates"

#' Posterior parameter draws for regression model of (log-p2) on p1.
#'
#' A Stan model was used to fit a model to the Phillips et al. (2017) point
#' estimates. The model was of the form:
#' \deqn{log p2 \sim normal(alpha + beta p1, sigma0 + sigma1 * p1)}
#'
#' @references Phillips et al., (2017). High sensitivity and interindividual variability
#' in the response of the human circadian system to evening light. PNAS.
#'
#' @format A list with 4000 parameter draws for each of:
#' \describe{
#'   \item{alpha}{intercept parameter}
#'   \item{beta}{slope parameter}
#'   \item{sigma0}{noise intercept parameter}
#'   \item{sigma1}{noise slope parameter}
#'   ...
#' }
"p1_p2_regression_draws"

#' A function representing the CDF of the p1 estimates from Phillips et al. (2017).
#'
#' @references Phillips et al., (2017). High sensitivity and interindividual variability
#' in the response of the human circadian system to evening light. PNAS.
"cdf"

#' A function representing the inverse-CDF of the p1 estimates from Phillips et al. (2017).
#'
#' @references Phillips et al., (2017). High sensitivity and interindividual variability
#' in the response of the human circadian system to evening light. PNAS.
"cdf_inv"

#' Posterior parameter draws for model of sigma
#'
#' The model (a gamma distribution) was fit in Stan.
#'
#' @format A list with 800 parameter draws for each of:
#' \describe{
#'   \item{a}{gamma parameter 1 in Stan definition of distribution}
#'   \item{b}{gamma parameter 2 in Stan definition of distribution}
#' }
"sigma_fit_draws"
