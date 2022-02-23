#' Plots simulated data and dose-response curves for a population
#'
#' @param experiment_df a simulated population
#' @param plot_curves a Boolean indicating whether to plot the true dose-response curves
#'
#' @return a ggplot object
#' @export
#'
#' @examples
#' # generate dose-response data for n=200 individuals measured at four lux levels
#' experimental_data <- virtual_experiment(n=200, lux=c(1, 10, 100, 1000))
#'
#' # plot data
#' plot_doseresponse(experimental_data)
#' @importFrom ggplot2 ggplot aes geom_point geom_jitter geom_line scale_x_log10 xlab ylab scale_y_continuous
#' @importFrom rlang .data
plot_doseresponse <- function(experiment_df, plot_curves=TRUE) {

  experiment_df <- experiment_df %>%
    dplyr::mutate(type="measurement")

  if(plot_curves) {
    p1s <- unique(experiment_df$p1)
    p2s <- unique(experiment_df$p2)
    lux_fine <- 10^seq(log10(1), log10(2000), length.out = 1000)
    for(i in seq_along(p1s)) {
      y <- purrr::map_dbl(lux_fine, ~logistic_2(., p1s[i], p2s[i]))
      temp <- dplyr::tibble(lux=lux_fine, y=y) %>%
        dplyr::mutate(id = i)
      if(i == 1)
        big_df <- temp
      else
        big_df <- big_df %>% dplyr::bind_rows(temp)
    }
    big_df <- big_df %>%
      dplyr::mutate(type="actual") %>%
      dplyr::bind_rows(experiment_df)
  } else {
    big_df <- experiment_df
  }

  p <- ggplot(data=big_df, aes(x=.data$lux, y=y, group=as.factor(.data$id))) +
    geom_line(data=big_df %>% dplyr::filter(.data$type=="actual")) +
    geom_jitter(data=big_df %>% dplyr::filter(.data$type=="measurement"),
                 width = 0.01, height = 0) +
    scale_x_log10() +
    scale_y_continuous(labels=scales::percent) +
    xlab("Illuminance [photopic lux]") +
    ylab("Melatonin suppression")

  if(!plot_curves)
    p <- p + geom_line(data=big_df %>% dplyr::filter(.data$type=="measurement"))
  p
}
