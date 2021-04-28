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
  return(1 - 1 / (1 + (log_lux / p1)^p2))
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
  return(10^(p1 * (-1 + 1 / (1 - edx))^(1 / p2)))
}
