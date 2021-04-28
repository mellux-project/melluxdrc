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
#'   ...
#' }
"estimates"
