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
  #while(cum_t < t_max){
  while(length(f) < t_max){
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

sim_data <- replicate(100, simulate(10, rgamma(1, shape = 5, rate = 0.1), alpha = rtnorm(1, 2, 1)), simplify = FALSE)
sim_par  <- my_sapply(sim_data, function(s) c(alpha = unique(s$alpha), sigma_a = unique(s$sigma_a)))

# number of points
par(mfrow=c(1, 1))
hist(sapply(sim_data, nrow), breaks = 100)
summary(sapply(sim_data, nrow))

# central tendency of fixation duration
hist(sapply(sim_data, function(x) mean  (x$d)), breaks = 100)
hist(sapply(sim_data, function(x) median(x$d)), breaks = 100)

# std.dev of fixation duration
hist(sapply(sim_data, function(x) sd(x$d)), breaks = 100)

stan_model <- rstan::stan_model(file = here::here("stan", "example_onlySaliency.stan"))
fits <- lapply(sim_data, function(d){
  dat <- list(N_rows = nrow(d), order = 1:nrow(d), x = d$x, y = d$y, duration = d$d,
              N_pixels = nrow(saldf), id_pixel = 1:nrow(saldf),
              x_pixel = saldf$x, y_pixel = saldf$y, log_val_pixel = saldf$log_s, log_area_pixel = 16, which_pixel = d$f)
  
  fit <- rstan::sampling(stan_model, dat, chains = 4, cores = 4, warmup = 500, iter = 1500)
  
  return(fit)
})


summaries <- my_sapply(fits, function(f){
  alpha <- rstan::extract(f)$alpha
  sigma_a <- rstan::extract(f)$sigma_a
  
  c(mean_alpha = mean(alpha), var_alpha = var(alpha), mean_sigma_a = mean(sigma_a), var_sigma_a = var(sigma_a))
})

bias_variance <- data.frame(z_alpha = (summaries$mean_alpha - sim_par$alpha) / sqrt(summaries$var_alpha),
                            contraction_alpha = 1 - summaries$var_alpha / var(sim_par$alpha),
                            z_sigma_a = (summaries$mean_sigma_a - sim_par$sigma_a) / sqrt(summaries$var_sigma_a),
                            contraction_sigma_a = 1 - summaries$var_sigma_a / var(sim_par$sigma_a))
plot(bias_variance$contraction_alpha, bias_variance$z_alpha, xlim = 0:1, bty = "l", ylim = c(-2, 2),
     xlab = "Posterior contraction", ylab = "Posterior z-score (standardized bias)", pch = 20)
abline(h = 0, lty = 2)
plot(bias_variance$contraction_sigma_a, bias_variance$z_sigma_a, xlim = 0:1, bty = "l", ylim = c(-2, 2),
     xlab = "Posterior contraction", ylab = "Posterior z-score (standardized bias)", pch = 20)
abline(h = 0, lty = 2)
