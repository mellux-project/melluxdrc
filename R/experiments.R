#' Generates samples of measured melatonin values at two lux levels
#'
#' These samples can either be from different individuals (for a between-type
#' study) or from the same individuals measured twice (for a within-type study).
#'
#' @param is_between a Boolean indicating whether experiment is within or between type
#' @param lux_1 first measured lux value
#' @param lux_2 second measured lux value
#' @param n number of individuals in sample
#' @param population_df a tibble representing a virtual experiment
#'
#' @return a list of two sets of measured melatonin values
generate_two_samples <- function(is_between, lux_1, lux_2, n, population_df) {
  if(!lux_1 %in% population_df$lux)
    stop("lux_1 must be in set of luxes in population_df")
  if(!lux_2 %in% population_df$lux)
    stop("lux_2 must be in set of luxes in population_df")

  n_id <- dplyr::n_distinct(population_df$id)
  if(2 * n > n_id)
    stop("number of individuals in study must be less than half population size.")

  if(is_between) {
    sample_ids <- sample(1:n_id, 2 * n)
    sample_data <- population_df %>%
      dplyr::filter(.data$id %in% sample_ids)
    df_1 <- sample_data %>%
      dplyr::filter(.data$lux == lux_1)
    df_2 <- sample_data %>%
      dplyr::filter(.data$lux == lux_2)

    # pick 1st / 2nd half of ids as population 1 / 2
    vals_1 <- df_1$y[1:n]
    vals_2 <- df_2$y[(n + 1):(2 * n)]
  } else{ # within
    sample_ids <- sample(1:n_id, n)
    sample_data <- population_df %>%
      dplyr::filter(.data$id %in% sample_ids)
    df_1 <- sample_data %>%
      dplyr::filter(.data$lux == lux_1)
    df_2 <- sample_data %>%
      dplyr::filter(.data$lux == lux_2)

    # pick 1st / 2nd half of ids as population 1 / 2
    vals_1 <- df_1$y
    vals_2 <- df_2$y
  }
  list(vals_1=vals_1, vals_2=vals_2)
}

#' Determines if a t-test has found a significant difference of the correct sign
#'
#' @inheritParams generate_two_samples
#' @param vals_1 set of melatonin values at first lux
#' @param vals_2 set of melatonin values at second lux
#' @param fit a result of running t.test
#' @param p_value a statistical test size
#'
#' @return a binary value indicating test success (if=1) or failure (if=0) to detect difference of correct sign
is_comparison_successful <- function(vals_1, vals_2, lux_1, lux_2, fit, p_value) {
  res <- 0
  if(lux_2 > lux_1) {
    if(mean(vals_2) > mean(vals_1))
      if(fit$p.value < p_value)
        res <- 1
  } else if (lux_2 < lux_1){
    if(mean(vals_2) < mean(vals_1))
      if(fit$p.value < p_value)
        res <- 1
  } else { # correct detection means no difference detected
    if(fit$p.value > p_value)
      res <- 1
  }
  res
}

#' Performs between- or within-individual experiments comparing melatonin suppression at two
#' lux levels using t tests
#'
#' For a between-individual experiment, a t-test is used to compare the melatonin suppression level at two lux values for two different
#' samples of individuals. This function requires based an (ideally large) simulated population of individual dose-response data.
#'
#' For a within-individual, a paired t-test is used to compare the melatonin suppression level at two lux values for a sample of individuals
#' measured at each of those values. This function requires based an (ideally large) simulated population of individual dose-response data.
#'
#' @inheritParams generate_two_samples
#' @inheritParams is_comparison_successful
#'
#' @return a binary value indicating test success (if=1) or failure (if=0) to detect difference of correct sign
#' @export
#' @importFrom rlang .data
#'
#' @examples
#' library(chronodoseresponse)
#'
#' # generate virtual data for 200 individuals
#' population_df <- virtual_experiment(200)
#'
#' # carry out 10 replicate between-individual experiments comparing lux_1=10 with lux_2=30
#' # for a test sample size of 20
#' is_between <- TRUE
#' results <- purrr::map_dbl(1:10, ~comparison_test(is_between, 10, 30, 20, population_df))
#'
#' # calculate percentage of experiments obtaining the correct signed difference
#' mean(results)
#' # carry out 10 replicate within-individual experiments comparing lux_1=10 with lux_2=30
#' # for a test sample size of 20
#' is_between <- FALSE
#' results <- purrr::map_dbl(1:10, ~comparison_test(is_between, 10, 30, 20, population_df))
#'
#' # calculate percentage of experiments obtaining the correct signed difference
#' mean(results)
comparison_test <- function(is_between, lux_1, lux_2, n, population_df, p_value=0.05) {

  vals <- generate_two_samples(is_between, lux_1, lux_2, n, population_df)
  vals_1 <- vals$vals_1
  vals_2 <- vals$vals_2
  if(is_between)
    is_paired <- FALSE
  else
    is_paired <- TRUE

  fit <- stats::t.test(vals_1, vals_2, paired=is_paired)

  is_success <- is_comparison_successful(vals_1, vals_2, lux_1, lux_2, fit, p_value=p_value)
  is_success
}
