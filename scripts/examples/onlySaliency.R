# example: only saliency
library(tidyverse)
library(imager)
library(OpenImageR)
library(saliency)
library(rstan)
rstan_options(auto_write = TRUE)
library(here)
source(here::here("R", "helpers.R"))
library(matrixStats)

par(mfrow = c(3, 2))
im <- imager::load.image("https://raw.githubusercontent.com/NUS-VIP/predicting-human-gaze-beyond-pixels/master/data/stimuli/1010.jpg")
saliency::view(im)
sal <- saliency::itti_koch(im, c_scale = c(0, 1, 2, 3, 4), d_scale = c(3, 4))
sal <- saliency::gaussian_pyramid(sal, 4)
lapply(sal, saliency::view)

#sal <- saliency::downsample(sal, get_gaussian_kernel())
#saliency::view(sal)
saldf <- matrix2df(sal[[5]], "top-left")
# we downsampled by a factor of 8 so rescale x and y coordinates
resc <- function(x, factor = 2) {(factor*(x-1) + 1 + factor*x)/2}#x*(1+factor)/2}
saldf$x <- resc(saldf$x, 16)
saldf$y <- resc(saldf$y, 16)
saldf$s_norm <- saldf$s / sum(saldf$s)
saldf$log_s <- log(saldf$s_norm)

ggplot(saldf, aes(x = x, y = y, fill = s_norm)) + 
  geom_tile() + 
  scale_fill_gradient(low="black", high = "white")


simulate <- function(t_max, sigma_a, alpha){
  cum_t <- 0
  counter <- 0
  
  f <- x <- y <- d <- numeric()
  while(cum_t < t_max){
    fix <- sample.int(nrow(saldf), 1, prob = saldf$s_norm)
    dist <- (saldf$x - saldf$x[fix])^2 + (saldf$y - saldf$y[fix])^2
    dist <- dist / 2
    nu <- -matrixStats::logSumExp(saldf$log_s - dist / sigma_a)
    
    duration <- rwald(1, alpha, nu)
    
    cum_t <- cum_t + duration
    f <- c(f, fix)
    x <- c(x, saldf$x[fix])
    y <- c(y, saldf$y[fix])
    d <- c(d, duration)
  }
  
  return(data.frame(f = f, x = x, y = y, d = d, alpha = alpha, sigma_a = sigma_a))
}

sim_data <- replicate(1000, simulate(10, rgamma(1, shape = 25, rate = 0.25), alpha = rtnorm(1, 2, 2)), simplify = FALSE)

# number of points
hist(sapply(sim_data, nrow), breaks = 100)
# mean fixation duration
hist(sapply(sim_data, function(x) mean(x$d)), breaks = 100)
# std.dev of fixation duration
hist(sapply(sim_data, function(x) sd(x$d)), breaks = 100)

stan_model <- rstan::stan_model(file = here::here("stan", "example_onlySaliency.stan"))
fits <- lapply(sim_data[1:10], function(d){
  dat <- list(N_rows = nrow(d), order = 1:nrow(d), x = d$x, y = d$y, duration = d$d,
              N_pixels = nrow(saldf), id_pixel = 1:nrow(saldf),
              x_pixel = saldf$x, y_pixel = saldf$y, log_val_pixel = saldf$log_s, log_area_pixel = 16, which_pixel = d$f)
  
  fit <- rstan::sampling(stan_model, dat, chains = 4, cores = 4, warmup = 500, iter = 1500)
  
  return(fit)
})
