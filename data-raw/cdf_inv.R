## code to prepare `cdf_inv` dataset goes here

kern_smooth_estimates <- KernSmooth::bkde(estimates$p1)
f <- approxfun(kern_smooth_estimates$x,
               kern_smooth_estimates$y,
               method="constant")
p1_fine <- seq(0, 3.1, 0.01)
cdf <- purrr::map_dbl(p1_fine, ~f_cdf(., f, 0.0))
cdf_inv <- approxfun(cdf, p1_fine)
cdf <- approxfun(p1_fine, cdf)

usethis::use_data(cdf_inv, overwrite = TRUE)
usethis::use_data(cdf, overwrite = TRUE)
