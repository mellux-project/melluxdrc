
<!-- README.md is generated from README.Rmd. Please edit that file -->

# melluxdrc: Power analysis for human melatonin suppression experiments (main package)

<!-- badges: start -->

[![R-CMD-check](https://github.com/mellux-project/melluxdrc/workflows/R-CMD-check/badge.svg)](https://github.com/mellux-project/melluxdrc/actions)
<!-- badges: end -->

The goal of `melluxdrc` is to allow generation of virtual dose-response
type data typical in melatonin suppression experiments.

## Installation

You can install `melluxdrc` via:

``` r
devtools::install_github("mellux-project/melluxdrc")
```

## Example

This shows how to generate simulated experimental dose-response data,
then plots it.

``` r
library("melluxdrc")
# generate dose-response data for 41 individuals measured at four lux levels
experimental_data <- virtual_experiment(n=41, lux=c(1, 10, 100, 1000))

# plot data
plot_doseresponse(experimental_data)
```

<img src="man/figures/README-example-1.png" width="100%" />
