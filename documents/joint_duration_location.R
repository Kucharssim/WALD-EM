library(tidyverse)
library(imager)
library(plotrix)
library(gtools)
library(rstan)
# rstan::rstan_options(auto_write=TRUE)
source(here::here("R", "load_image.R"))
source(here::here("R", "helpers.R"))
source(here::here("R", "expose_helpers_stan.R"))

image_nr <- 1001
image_name  <- paste0(image_nr, ".jpg")
objects <- read.csv(here::here("data", "object_familiarity", "objects.csv"))
objects <- subset(objects, image_name == image_nr)

img <- load_image(image_name)
dim_img <- list(min_x = 0, min_y = 0, max_x = imager::width(img), max_y = imager::height(img))
plot(img, axes = FALSE)
points(objects$x, objects$y, pch = 10, cex = 2)
for(i in 1:nrow(objects))
  plotrix::draw.ellipse(x = objects$x[i], y = objects$y[i], a = objects$width[i]/2, b = objects$height[i]/2, angle = 0, lwd = 2)


n_sim <- 1
t_max <- 2000
set.seed(2020)
n_objects <- nrow(objects)
sim_parameters_list <- replicate(n = n_sim, expr = list(pi    = gtools::rdirichlet(1, rep(2, n_objects))[1,,drop=TRUE],
                                                        delta = rexp(1, 1),
                                                        alpha = rtnorm(1, 2, 1, 0, Inf),
                                                        sigma_attention = rgamma(1, shape = 2, rate = 0.1)
),
simplify = FALSE)
sim_parameters <- lapply(sim_parameters_list, unlist)
sim_parameters <- as.data.frame(do.call(rbind, sim_parameters))

simulate <- function(t_max, pi, delta, alpha, sigma_attention, n_max = t_max * 10){
  x <- y <- duration <- nu <- numeric()
  t <- 0
  sigma_x <- delta * objects$width / 2
  sigma_y <- delta * objects$height / 2
  
  while(t <= t_max && length(x) < n_max){
    which_factor <- sample(1:n_objects, 1, FALSE, pi)
    x_now <- rtnorm(1, objects$x[which_factor], sigma_x[which_factor], dim_img$min_x, dim_img$max_x)
    y_now <- rtnorm(1, objects$y[which_factor], sigma_y[which_factor], dim_img$min_y, dim_img$max_y)
    nu_now <- -log_integral_attention_mixture_2d(x = x_now, y = y_now, weights = pi, 
                                                 mu_x = objects$x, sigma_x = sigma_x, 
                                                 mu_y = objects$y, sigma_y = sigma_y, 
                                                 width_x = sigma_attention, width_y = sigma_attention)
    duration_now <- wald_rng(alpha, nu_now)
    x <- c(x, x_now)
    y <- c(y, y_now)
    nu <- c(nu, nu_now)
    duration <- c(duration, duration_now)
    t <- t + duration_now
  }
  
  return(data.frame(x=x, y=y, duration=duration, nu=nu))
}

sim_data_list <- lapply(seq_along(sim_parameters_list), function(i) {
  pars <- sim_parameters_list[[i]]
  d <- simulate(t_max = t_max, pi = pars$pi, delta = pars$delta, alpha = pars$alpha, sigma_attention = pars$sigma_attention)
  d$sim <- i
  return(d)
})

sim_data <- do.call(rbind, sim_data_list)
hist(sim_data$duration)

get_tile <- function(fix, width = 50, min = 0, max = 800) {
  cut(fix, breaks = seq(min, max, by = width), labels = seq_len(ceiling((max-min)/width)))
}


sim_data$tile_x <- get_tile(sim_data$x, 50, 0, 800)
sim_data$tile_y <- get_tile(sim_data$y, 50, 0, 600)
sim_data$tile   <- interaction(sim_data$tile_x, sim_data$tile_y, sep = "_")

sim_data %>%
  group_by(tile) %>%
  summarise(spatial_importance = n(), mean_duration = mean(duration)) %>%
  ggplot(aes(x=spatial_importance, y = mean_duration)) +
  geom_point()

stan_model <- rstan::stan_model(here::here("stan/examples/object_oriented_behavior.stan"))
stan_data <- list(
  N_obs            = nrow(sim_data),
  
  x                = sim_data$x,
  y                = sim_data$y,
  duration         = sim_data$duration,
  
  N_objects        = n_objects,
  objects_center_x = objects$x, 
  objects_center_y = objects$y,
  objects_width    = objects$width,
  objects_height   = objects$height
)

stan_fit <- rstan::sampling(stan_model, stan_data, chains = 2, cores = 2, warmup = 1000, iter = 2000, 
                            control = list(adapt_delta = 0.95))
stan_fit

plot(summary(stan_fit, pars = "mean_duration_pred")$summary[, "mean"],
     sim_data$duration, 
     xlab = "Predicted mean of fixation duration", 
     ylab = "Fixation duration")
abline(a = 0, b = 1, lwd = 5)
