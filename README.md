
<!-- README.md is generated from README.Rmd. Please edit that file -->
chronodoseresponse
==================

<!-- badges: start -->
[![R-CMD-check](https://github.com/ben18785/chronodoseresponse/workflows/R-CMD-check/badge.svg)](https://github.com/ben18785/chronodoseresponse/actions) <!-- badges: end -->

The goal of chronodoseresponse is to allow generation of virtual dose-response type data typical in chronobiology experiments.

Installation
------------

You can install chronodoseresponse via:

``` r
install_github("ben18785/chronodoseresponse")
```

Example
-------

This is a basic example which shows you how to solve a common problem:

``` r
library(chronodoseresponse)
# generate dose-response data for n=200 individuals measured at four lux levels
experimental_data <- virtual_experiment(n=200, lux=c(1, 10, 100, 1000))
```
