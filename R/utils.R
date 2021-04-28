#' Two-parameter dose-response logistic curve
#'
#' Returns a melatonin suppression response for a given level of log_lux
#'
#' @param log_lux a log lux value
#' @param log_ed50 the log of the ed50 lux
#' @param p2 a positive shape parameter
#'
#' @return
#' @export
#'
#' @examples
#' # below returns 0.5 since log_lux=log_ed50
#' logistic_2(10, 10, 5)
logistic_2 <- function(log_lux, log_ed50, p2) {
  return(1 - 1 / (1 + (log_lux / log_ed50)^p2))
}
