#' Generates noisy dose response data at a range of (log-)luxes
#'
#' Adds noise on the logit scale, resulting in more realistic (noisy) data for
#' a given patient.
#'
#'
#' @param sigma a noise level (represents additive Gaussian noise on logit scale)
#' @param p1 the log of the ed50 lux
#' @param p2 a positive shape parameter
#' @param log_lux a vector of log lux values which defaults to log10(c(10, 30, 50, 100, 200, 400, 2000))
#'
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
