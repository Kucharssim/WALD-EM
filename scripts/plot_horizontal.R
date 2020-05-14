library(rstan)
library(tidyverse)
source(here::here("R", "expose_helpers_stan.R"))


x_prev <- 400; y_prev <- 300
x_min <- 1; x_max <- 800; y_min <- 1; y_max <- 600
display <- tidyr::expand_grid(x = seq(x_min, x_max, by = 10), y = seq(y_min, y_max, by = 10))
display$angle           <- NA
display$radius          <- NA
display$dist_to_border  <- NA
display$log_lik         <- NA

pb <- dplyr::progress_estimated(nrow(display))
for(i in 1:nrow(display)){
  res <- calc_angle_border(display$x[i], display$y[i], x_prev, y_prev, x_min, x_max, y_min, y_max)
  
  display$angle[i]          <- res[1]
  display$radius[i]         <- res[2]
  display$dist_to_border[i] <- res[3]
  display$log_lik[i]        <- direction_bias_lpdf(res[1], res[2], res[3], rep(1/2, 2), c(0, pi), rep(8, 2))
  pb$tick()$print()
}

display$lik <- exp(display$log_lik)

ggplot2::ggplot(display, ggplot2::aes(x = x, y = y)) + 
  ggplot2::geom_tile(ggplot2::aes(fill = angle))

ggplot2::ggplot(display, ggplot2::aes(x = x, y = y)) +
  ggplot2::geom_tile(ggplot2::aes(fill = radius))

ggplot2::ggplot(display, ggplot2::aes(x = x, y = y)) + 
  ggplot2::geom_tile(ggplot2::aes(fill = dist_to_border))

ggplot2::ggplot(display, ggplot2::aes(x = x, y = y)) + 
  ggplot2::geom_tile(ggplot2::aes(fill = log_lik))

ggplot2::ggplot(display, ggplot2::aes(x = x, y = y)) + 
  ggplot2::geom_tile(ggplot2::aes(fill = lik))


## check that the pdf integrates to 1
library(pracma)
library(circular)

integrand <- function(x, y){
  res <- calc_angle_border(x, y, x_prev, y_prev, x_min, x_max, y_min, y_max)
  #r <- sqrt((x-x_prev)^2 + (y-y_prev)^2)
  
  lpdf <- direction_bias_lpdf(res[1], res[2], res[3], rep(1/2, 2), c(0, pi), rep(8, 2))
  
  exp(lpdf)
}

int <- pracma::integral2(fun = integrand, xmin = x_min, xmax = x_max, ymin = y_min, ymax = y_max, vectorized = FALSE)
int

integrand <- function(x, y, max_r){
  r <- x; theta <- y
  d <- circular::dvonmises(x = theta, mu = pi/2, kappa = 8, log = FALSE)
  d <- d * dunif(x = r, min = 0, max = max_r, log = FALSE)
  
  return(d)
}

int <- pracma::integral2(fun = integrand, xmin = 0, xmax = 10, ymin = 0, ymax = 2*pi, max_r = 10, vectorized = FALSE)
