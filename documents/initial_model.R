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

image_nr <- c(1001, 1014, 1049)
image_name  <- paste0(image_nr, ".jpg")
objects <- read.csv(here::here("data", "objects.csv")) %>% subset(image %in% image_nr)
objects_in_images <- read.csv(here::here("data", "objects_in_images.csv")) %>% subset(id_img %in% objects$id_img)
saliency <- read.csv(here::here("data", "saliency.csv")) %>% subset(image %in% image_nr)


par(mfcol = c(2, 3))
for(i in seq_along(image_name)) {
  img <- load_image(image_name[i])
  dim_img <- list(min_x = 0, min_y = 0, max_x = imager::width(img), max_y = imager::height(img))
  obj <- subset(objects, image_name == image_nr[i])
  sal <- subset(saliency, image == image_nr[i])
  
  plot(img, axes = FALSE)
  points(obj$x, obj$y, pch = 10, cex = 2)
   for(i in 1:nrow(obj))
     plotrix::draw.ellipse(x = obj$x[i], y = obj$y[i], a = obj$width[i]/2, b = obj$height[i]/2, angle = 0, lwd = 2)
  
  plot(imager::as.cimg(sal$value, x = max(sal$row), y = max(sal$col)), axes = FALSE)
}
par(mfrow = c(1, 1))


# draw from priors
N_sim <- 2
N_ppt <- 10
N_obj <- nrow(objects)

# scalar parameters
true_parameters <- data.frame(
  sim = seq_len(N_sim),
  sigma_center = rgamma(N_sim, 2, 0.02),
  sigma_distance = rgamma(N_sim, 2, 0.02),
  scale_obj = replicate(N_sim, trunc_normal_rng(1, 0.5, 0, Inf)),
  mu_log_alpha = rnorm(N_sim, 0, 0.5),
  sigma_log_alpha = rgamma(N_sim, 2, 2),
  mu_log_sigma_attention = rnorm(N_sim, 2, 1),
  sigma_log_sigma_attention = rgamma(N_sim, 2, 2)
)

# vector parameters
# weights of factors
true_weights <- as.data.frame(gtools::rdirichlet(N_sim, rep(2, 4)))
colnames(true_weights) <- sprintf("weights[%s]", 1:4)

# invlogit weights objects 
true_z_weights_obj <- as.data.frame(matrix(rnorm(N_obj * N_sim), nrow = N_sim))
colnames(true_z_weights_obj) <- sprintf("z_weights_obj[%s]", seq_len(N_obj))

# individual parameters: alpha
true_z_log_alpha <- sapply(seq_len(N_ppt), function(p) rnorm(N_sim, true_parameters$mu_log_alpha, true_parameters$sigma_log_alpha))
true_z_log_alpha <- as.data.frame(true_z_log_alpha)  
colnames(true_z_log_alpha) <- sprintf("z_log_alpha[%s]", seq_len(N_ppt))  

# individual parameters: sigma_attention
true_z_log_sigma_attention <- sapply(seq_len(N_ppt), function(p) rnorm(N_sim, true_parameters$mu_log_sigma_attention, true_parameters$sigma_log_sigma_attention))
true_z_log_sigma_attention <- as.data.frame(true_z_log_sigma_attention)  
colnames(true_z_log_sigma_attention) <- sprintf("z_log_sigma_attention[%s]", seq_len(N_ppt))  

# merge it all together:
#true_parameters <- bind_cols(true_parameters, true_weights, true_z_weights_obj, true_z_log_alpha, true_z_log_sigma_attention)
#rm(true_weights, true_z_weights_obj, true_z_log_alpha, true_z_log_sigma_attention)

# prepare design:
design <- expand.grid(id_ppt = seq_len(N_ppt), id_img = seq_along(image_nr), sim = seq_len(N_sim), KEEP.OUT.ATTRS = FALSE)
design$image_nr <- image_nr[design$id_img]


simulate_trial <- function(specs, t_max = 10, n_max = t_max * 10){
  browser()
  id_ppt <- specs[['id_ppt']]
  id_img <- specs[['id_img']]
  sim    <- specs[['sim']]
  sal    <- saliency %>% subset(image %in% specs[["image_nr"]])
  obj    <- objects %>% subset(image %in% specs[["image_nr"]])
    
  sigma_center              <- true_parameters[sim, "sigma_center", drop=TRUE]
  sigma_distance            <- true_parameters[sim, "sigma_distance", drop=TRUE]
  scale_obj                 <- true_parameters[sim, "scale_obj", drop=TRUE]
  
  mu_log_alpha              <- true_parameters[sim, "mu_log_alpha", drop=TRUE]
  sigma_log_alpha           <- true_parameters[sim, "sigma_log_alpha", drop=TRUE]
  mu_log_sigma_attention    <- true_parameters[sim, "mu_log_sigma_attention", drop=TRUE]
  sigma_log_sigma_attention <- true_parameters[sim, "sigma_log_sigma_attention", drop=TRUE]
  
  weights               <- true_weights[sim,] %>% as.vector()
  
  z_weights_obj         <- true_z_weights_obj[sim, objects_in_images$from[id_img]:objects_in_images$to[id_img]] %>% as.vector()
  weights_obj           <- exp(z_weights_obj) / sum(exp(z_weights_obj))
  
  z_log_alpha           <- true_z_log_alpha[sim, id_ppt, drop=TRUE]
  alpha                 <- exp(mu_log_alpha + sigma_log_alpha * z_log_alpha)
  
  z_log_sigma_attention <- true_z_log_sigma_attention[sim, id_ppt, drop=TRUE]
  sigma_attention       <- exp(mu_log_sigma_attention + sigma_log_sigma_attention * z_log_sigma_attention)
  
  width_obj_x <- scale_obj * obj$width
  width_obj_y <- scale_obj * obj$height
  
  x <- y <- duration <- nu <- numeric()
  t <- 0

  while(t <= t_max && length(x) < n_max) {
    att_filter   <- numeric(length = 2)
    which_factor <- sample(seq_along(weights), 1, FALSE, weights)
    
    if(which_factor == 1) { # objects
      xy_now <- mixture_trunc_normals_rng(weights_obj, obj$x, width_obj_x, obj$y, width_obj_y, 0, 800, 0, 600)
      x_now <- xy_now[1]
      y_now <- xy_now[2]
      
    } else if(which_factor == 2) { # saliency
      xy_now <- saliency_rng(sal$value_normalized, sal$x, sal$y, 20)
      x_now <- xy_now[1]
      y_now <- xy_now[2]
      
    } else if(which_factor == 3) { # exploitation
      if(t == 0) {
        x_now <- trunc_normal_rng(400, sigma_distance, 0, 800)
        y_now <- trunc_normal_rng(300, sigma_distance, 0, 600)
      } else {
        x_now <- trunc_normal_rng(x[length(x)], sigma_distance, 0, 800)
        y_now <- trunc_normwl_rng(y[length(y)], sigma_distance, 0, 600)
      }
    } else { # central bias
      x_now <- trunc_normal_rng(400, sigma_center, 0, 800)
      y_now <- trunc_normal_rng(300, sigma_center, 0, 600)
    }
    att_filter[1] <- log_integral_attention_mixture_2d(x_now, y_now, weights_obj, obj$x, width_obj_x, obj$y, width_obj_y, sigma_attention, sigma_attention)
    att_filter[2] <- -100
    
    nu_now <- log(sum(weights[1:2])) - log(sum(exp(att_filter)))
    duration_now <- wald_rng(alpha, nu_now)
    
    
    x <- c(x, x_now)
    y <- c(y, y_now)
    duration <- c(duration, duration_now)
    nu <- c(nu, duration_now)
    t <- t + duration_now
  }
  return(data.frame(x=x, y=y, duration=duration, nu=nu))
}

sim_data <- plyr::ddply(.data = design, .variables = c("sim", "id_ppt", "id_img"), 
                        .fun = simulate_trial, .progress = "text")
