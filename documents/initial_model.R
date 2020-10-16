# Script with the initial model simulations: later it will be rewritten in rmarkdown for a nice github document
library(tidyverse)
library(imager)
library(plotrix)
library(gtools)
library(rstan)
library(here)
source(here::here("R", "load_image.R"))
source(here::here("R", "helpers.R"))
source(here::here("R", "expose_helpers_stan.R"))
log_sum_exp <- matrixStats::logSumExp

## definition of critical radius
# subset saliency in radius of 100 from the fixation position
dist_screen <- 60  # distance from the screen (cm)
theta       <- 5   # radius of foveal vision (degrees of visual angle)
width_cm    <- 51  # width of the screen (cm)
width_pix   <- 800 # number of pixels in horizontal location
size_pix    <- width_cm / width_pix # size of one pixel
radius <- tan(theta * pi / 180)*dist_screen / size_pix


image_nr <- c(1001, 1014, 1049)
image_name  <- paste0(image_nr, ".jpg")
objects <- read.csv(here::here("data", "objects.csv")) %>% subset(image %in% image_nr)
objects_in_images <- read.csv(here::here("data", "objects_in_images.csv")) %>% subset(id_img %in% objects$id_img)
saliency <- read.csv(here::here("data", "saliency.csv")) %>% subset(image %in% image_nr)


par(mfcol = c(2, 3))
for(i in seq_along(image_name)) {
  img <- load_image(image_name[i])
  dim_img <- list(min_x = 0, min_y = 0, max_x = imager::width(img), max_y = imager::height(img))
  obj <- subset(objects, image == image_nr[i])
  sal <- subset(saliency, image == image_nr[i])
  
  plot(img, axes = FALSE)
  points(obj$x, obj$y, pch = 10, cex = 2)
   for(i in 1:nrow(obj))
     plotrix::draw.ellipse(x = obj$x[i], y = obj$y[i], a = obj$width[i]/2, b = obj$height[i]/2, angle = 0, lwd = 2)
  
  plot(imager::as.cimg(sal$value, x = max(sal$row), y = max(sal$col)), axes = FALSE)
}
par(mfrow = c(1, 1))

# sal <- subset(saliency, image == image_nr[1])
# xy <- replicate(1e5, saliency_rng(sal$value_normalized, sal$x, sal$y, 20), simplify = FALSE)
# xy <- do.call(rbind, xy)
# xy[,2] <- 600-xy[,2]
# par(mfrow = c(1,2))
# plot(xy, pch = 19, cex = 0.1)
# plot(imager::as.cimg(sal$value, x = max(sal$row), y = max(sal$col)), axes = FALSE)
# par(mfrow = c(1, 1))

# draw from priors
N_sim <- 10
N_ppt <- 20
N_obj <- nrow(objects)

if(!file.exists(here::here("documents", "initial_model_saves", "true.Rdata"))) {
  # scalar parameters
  true_parameters <- data.frame(
    sigma_center = rgamma(N_sim, 2, 0.02),
    sigma_distance = rgamma(N_sim, 2, 0.02),
    scale_obj = replicate(N_sim, trunc_normal_rng(1, 0.5, 0, Inf)),
    mu_log_alpha = rnorm(N_sim, 0, 0.5),
    sigma_log_alpha = rgamma(N_sim, 2, 5),
    mu_log_sigma_attention = rnorm(N_sim, 4, 1),
    sigma_log_sigma_attention = rgamma(N_sim, 2, 5)
  )
  
  # vector parameters
  # weights of factors
  true_weights <- as.data.frame(gtools::rdirichlet(N_sim, rep(2, 4)))
  colnames(true_weights) <- sprintf("weights[%s]", 1:4)
  
  # invlogit weights objects 
  true_z_weights_obj <- as.data.frame(matrix(rnorm(N_obj * N_sim), nrow = N_sim))
  colnames(true_z_weights_obj) <- sprintf("z_weights_obj[%s]", seq_len(N_obj))
  
  # individual parameters: alpha
  true_alpha <- sapply(seq_len(N_ppt), function(p) rnorm(N_sim, true_parameters$mu_log_alpha, true_parameters$sigma_log_alpha))
  true_alpha <- exp(true_alpha) 
  # true_alpha[true_alpha > 10] <- rexp(sum(true_alpha > 10), 1) # get rid of unrealistically high alphas
  
  true_alpha <- as.data.frame(true_alpha)
  colnames(true_alpha) <- sprintf("alpha[%s]", seq_len(N_ppt))
  
  
  # individual parameters: sigma_attention
  true_sigma_attention <- sapply(seq_len(N_ppt), function(p) rnorm(N_sim, true_parameters$mu_log_sigma_attention, true_parameters$sigma_log_sigma_attention))
  true_sigma_attention <- exp(true_sigma_attention)
  
  true_sigma_attention <- as.data.frame(true_sigma_attention)  
  colnames(true_sigma_attention) <- sprintf("sigma_attention[%s]", seq_len(N_ppt))  
  
  # save generating values:
  save(true_parameters, true_weights, true_z_weights_obj, true_alpha, true_sigma_attention, 
       file = here::here("documents", "initial_model_saves", "true.Rdata"))
} else {
  load(here::here("documents", "initial_model_saves", "true.Rdata"))
}

# prepare design:
design <- expand.grid(id_ppt = seq_len(N_ppt), id_img = seq_along(image_nr), sim = seq_len(N_sim), KEEP.OUT.ATTRS = FALSE)
design$image_nr <- image_nr[design$id_img]


simulate_trial <- function(specs, t_max = 5, n_max = t_max * 10){
  # browser()
  id_ppt <- specs[['id_ppt']]
  id_img <- specs[['id_img']]
  sim    <- specs[['sim']]
  sal    <- saliency %>% subset(image %in% specs[["image_nr"]])
  obj    <- objects %>% subset(image %in% specs[["image_nr"]])
    
  sigma_center              <- true_parameters[sim, "sigma_center", drop=TRUE]
  sigma_distance            <- true_parameters[sim, "sigma_distance", drop=TRUE]
  scale_obj                 <- true_parameters[sim, "scale_obj", drop=TRUE]
  
  mu_log_sigma_attention    <- true_parameters[sim, "mu_log_sigma_attention", drop=TRUE]
  sigma_log_sigma_attention <- true_parameters[sim, "sigma_log_sigma_attention", drop=TRUE]
  
  weights                   <- true_weights[sim,] %>% as.vector()
  
  z_weights_obj             <- true_z_weights_obj[sim, objects_in_images$from[id_img]:objects_in_images$to[id_img]] %>% as.vector()
  weights_obj               <- as.matrix(exp(z_weights_obj) / sum(exp(z_weights_obj)), ncol = 1)
  
  alpha                     <- true_alpha[sim, id_ppt, drop=TRUE]
  sigma_attention           <- true_sigma_attention[sim, id_ppt, drop=TRUE]
  
  center_obj_x <- as.matrix(obj$x, ncol = 1)
  center_obj_y <- as.matrix(obj$y, ncol = 1)
  width_obj_x <- as.matrix(scale_obj * obj$width, ncol = 1)
  width_obj_y <- as.matrix(scale_obj * obj$height, ncol = 1)
  
  x <- y <- duration <- nu <- log_lik_saliency <- n_neighbors <- numeric()
  t <- 0
  m_sq_dist <- matrix(0, nrow = 0, ncol = nrow(sal))
  saliency_log <- matrix(0, nrow = 0, ncol = nrow(sal))
  
  while(t <= t_max && length(x) < n_max) {
    att_filter   <- numeric(length = 2)
    which_factor <- sample(seq_along(weights), 1, FALSE, weights)
    
    if(which_factor == 1) { # objects
      xy_now <- mixture_trunc_normals_rng(weights_obj, center_obj_x, width_obj_x, center_obj_y, width_obj_y, 0, 800, 0, 600)
      x_now <- xy_now[1]
      y_now <- xy_now[2]
      
    } else if(which_factor == 2) { # saliency
      xy_now <- saliency_rng(sal$value_normalized, sal$x, sal$y, 20)
      x_now <- xy_now[1] - 0.5
      y_now <- xy_now[2] - 0.5
      
    } else if(which_factor == 3) { # exploitation
      if(t == 0) {
        x_now <- trunc_normal_rng(400, sigma_distance, 0, 800)
        y_now <- trunc_normal_rng(300, sigma_distance, 0, 600)
      } else {
        x_now <- trunc_normal_rng(x[length(x)], sigma_distance, 0, 800)
        y_now <- trunc_normal_rng(y[length(y)], sigma_distance, 0, 600)
      }
    } else { # central bias
      x_now <- trunc_normal_rng(400, sigma_center, 0, 800)
      y_now <- trunc_normal_rng(300, sigma_center, 0, 600)
    }
    distances <- sqrt((x_now - sal$x)^2 + (y_now - sal$y)^2)
    mean_sq_distances <- distances^2 / 2
    which_closest <- which.min(distances)
    which_neighbors <- distances < radius
    
    
    att_filter[1] <- log(weights[[1]]) + log_integral_attention_mixture_2d(x_now, y_now, weights_obj, center_obj_x, width_obj_x, center_obj_y, width_obj_y, sigma_attention, sigma_attention)
    att_filter[2] <- log(weights[[2]]) + log_sum_exp(sal$saliency_log - mean_sq_distances / sigma_attention^2)
    
    saliency_log_lik_now <- sal$log_lik_saliency[which_closest]
    nu_now <- log(sum(weights[1:2])) - log_sum_exp(att_filter)
    duration_now <- wald_rng(alpha, nu_now)
    
    x <- c(x, x_now)
    y <- c(y, y_now)
    duration <- c(duration, duration_now)
    nu <- c(nu, nu_now)
    log_lik_saliency <- c(log_lik_saliency, saliency_log_lik_now)
    n_neighbors <- c(n_neighbors, sum(which_neighbors))
    
    t <- t + duration_now
    
    m_sq_dist <- rbind(m_sq_dist, sort(mean_sq_distances))
    saliency_log <- rbind(saliency_log, sal$saliency_log[order(mean_sq_distances)])
  }
  m_sq_dist <- as.data.frame(m_sq_dist)
  colnames(m_sq_dist) <- sprintf("m_sq_dist[%s]", seq_len(ncol(m_sq_dist)))
  
  saliency_log <- as.data.frame(saliency_log)
  colnames(saliency_log) <- sprintf("saliency_log[%s]", seq_len(ncol(saliency_log)))
  
  data <- data.frame(order=seq_along(x), x=x, y=y, duration=duration, nu=nu, 
                     log_lik_saliency=log_lik_saliency, n_neighbors=n_neighbors)
  data <- cbind(data, m_sq_dist)
  data <- cbind(data, saliency_log)
  
  return(data)
}

if(!file.exists(here::here("documents", "initial_model_saves", "sim_data.Rdata"))){
  sim_data <- plyr::ddply(.data = design, .variables = c("sim", "id_ppt", "id_img"), 
                          .fun = simulate_trial, .progress = "text")
  save(sim_data, file = here::here("documents", "initial_model_saves", "sim_data.Rdata"))
} else {
  load(here::here("documents", "initial_model_saves", "sim_data.Rdata"))
}

mean(sim_data$duration < 1)
hist(sim_data$duration[sim_data$duration < 1])
sim_data %>% group_by(sim, id_ppt, id_img) %>% summarise(t = n()) %>% ggplot(aes(x = t)) + geom_histogram()
sim_data %>% group_by(sim) %>% summarise(t = n())

drop <- sprintf("m_sq_dist[%s]", (max(sim_data$n_neighbors)+1):300)
sim_data <- sim_data[, !(colnames(sim_data) %in% drop)]
drop <- sprintf("saliency_log[%s]", (max(sim_data$n_neighbors)+1):300)
sim_data <- sim_data[, !(colnames(sim_data) %in% drop)]

get_stan_data <- function(data) {
  list(
    N_obs             = nrow(data),
    order             = data$order,
    x                 = data$x,
    y                 = data$y,
    duration          = data$duration,
    
    N_obj             = nrow(objects),
    obj_center_x      = objects$x,
    obj_center_y      = objects$y,
    obj_width         = objects$width,
    obj_height        = objects$height,
    N_ppt             = dplyr::n_distinct(data$id_ppt),
    id_ppt            = data$id_ppt,
    N_img             = dplyr::n_distinct(data$id_img),
    id_img            = data$id_img,
    obj_index_from    = objects_in_images$from,
    obj_index_to      = objects_in_images$to,
    N_obj_in_img      = objects_in_images$n,
    log_lik_saliency  = data$log_lik_saliency,
    max_neighbors     = length(dplyr::starts_with("m_sq_dist", vars = colnames(sim_data))),
    N_neighbors       = data$n_neighbors,
    mean_sq_distances = dplyr::select(data, dplyr::starts_with("m_sq_dist")) %>% as.matrix(),
    saliency_log      = dplyr::select(data, dplyr::starts_with("saliency_log[")) %>% as.matrix()
    
    # lb_x              = 0,
    # ub_x              = 800,
    # lb_y              = 0,
    # ub_y              = 600
  )
}

stan_model <- rstan::stan_model(here::here("stan", "objects_central_distance_saliency.stan"))

fit_sim <- function(data) {
  stan_data <- get_stan_data(data)
  
  fit <- rstan::sampling(stan_model, stan_data, cores = 2, chains = 2, iter = 750, warmup = 500, refresh = 250)

  return(fit)
}

if(!file.exists(here::here("documents", "initial_model_saves", "fits.Rdata"))) {
  fits <- plyr::dlply(.data = sim_data, .variables = c("sim"), .fun = fit_sim, .progress = "tk")
  save(fits, file = here::here("documents", "initial_model_saves", "fits.Rdata"))
} else {
  load(here::here("documents", "initial_model_saves", "fits.Rdata"))
}


get_par <- function(par, true) {
  fit_summary <- t(sapply(fits, function(fit) summary(fit, pars = par)$summary[, c("mean", "2.5%", "97.5%"), drop=TRUE]))
  fit_summary <- as.data.frame(fit_summary)
  fit_summary$true <- true[, par]
  
  fit_summary
}

plot_par <- function(par, true) {
  df <- get_par(par, true)
  lim <- range(as.matrix(df))
  plot(df$true, df$mean, pch = 19, cex = 1, main = par, 
       xlab = "True", ylab = "Estimated", 
       xlim = lim, ylim = lim)
  segments(x0 = df$true, y0 = df$`2.5%`, y1 = df$`97.5%`)
  abline(a = 0, b = 1)
}

plot_vec_par <- function(par) {
  
}

par(mfrow=c(2, 4))
for(par in colnames(true_parameters)) {
  plot_par(par, true_parameters)
}

par(mfrow=c(1, 4))
for(par in colnames(true_weights)) {
  plot_par(par, true_weights)
}

par(mfrow=c(4, 4))
for(par in colnames(true_z_weights_obj)) {
  plot_par(par, true_z_weights_obj)
}


par(mfrow=c(4, 5))
for(par in colnames(true_alpha)) {
  plot_par(par, true_alpha)
}

par(mfrow=c(4, 5))
for(par in colnames(true_sigma_attention)) {
  plot_par(par, true_sigma_attention)
}
