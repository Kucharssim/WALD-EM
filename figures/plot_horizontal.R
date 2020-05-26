library(rstan)
library(tidyverse)
source(here::here("R", "expose_helpers_stan.R"))
source(here::here("R", "colours.R"))

x_prev <- 400; y_prev <- 300
x_min <- 1; x_max <- 800; y_min <- 1; y_max <- 600
resolution <- 10
display <- tidyr::expand_grid(x = seq(x_min, x_max, by = resolution), y = seq(y_min, y_max, by = resolution))
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
  display$log_lik[i]        <- direction_bias_lpdf(res[1], res[2], res[3], rep(1/2, 2), c(0, pi), rep(15, 2))
  pb$tick()$print()
}

display$lik       <- exp(display$log_lik)
display$log_joint <- display$log_lik + log(display$radius)
display$joint     <- exp(display$log_joint)

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

ggplot2::ggplot(display, ggplot2::aes(x = x, y = y)) + 
  ggplot2::geom_tile(ggplot2::aes(fill = log_joint))

ggplot2::ggplot(display, ggplot2::aes(x = x, y = y)) + 
  ggplot2::geom_tile(ggplot2::aes(fill = joint)) +
  ggplot2::theme_void() +
  ggplot2::scale_fill_continuous(name = expression(f(theta, r)))

ggplot2::ggsave("horizontal_example.jpg", path = here::here("figures"), width = 0.6*8.5, height = 0.6*6)

## check that the pdf integrates to about 1
library(pracma)

integrand <- function(x, y){
  res <- calc_angle_border(x, y, x_prev, y_prev, x_min, x_max, y_min, y_max)
  #r <- sqrt((x-x_prev)^2 + (y-y_prev)^2)
  
  lpdf <- direction_bias_lpdf(res[1], res[2], res[3], rep(1/2, 2), c(0, pi), rep(15, 2))
  
  exp(lpdf)
}

int <- pracma::integral2(fun = integrand, xmin = x_min, xmax = x_max, ymin = y_min, ymax = y_max, vectorized = FALSE)
int
