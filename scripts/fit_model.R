# Application of the Dynamic model of eye movements on the Renswoude's object familiarity data
library(tidyverse)
library(rstan)
rstan_options(auto_write = TRUE)
library(here)

# load cleaned data
load(here::here("data", "cleaned_data.Rdata"))
df_fit <- df %>% subset(train)

model <- rstan::stan_model(here::here("stan", "objects_central_distance_saliency.stan"))

stan_data <- list(
  N_obs             = nrow(df_fit),
  order             = df_fit$order,
  x                 = df_fit$x,
  y                 = df_fit$y,
  duration          = df_fit$duration,
  
  N_obj             = nrow(objects),
  obj_center_x      = objects$x,
  obj_center_y      = objects$y,
  obj_width         = objects$width,
  obj_height        = objects$height,
  N_ppt             = dplyr::n_distinct(df_fit$id_ppt),
  id_ppt            = df_fit$id_ppt,
  N_img             = dplyr::n_distinct(df_fit$id_img),
  id_img            = df_fit$id_img,
  obj_index_from    = objects_in_images$from,
  obj_index_to      = objects_in_images$to,
  N_obj_in_img      = objects_in_images$n,
  log_lik_saliency  = df_fit$log_lik_saliency,
  max_neighbors     = ncol(saliency_log),
  N_neighbors       = df_fit$n_neighbors,
  mean_sq_distances = mean_sq_distances[df$train,,drop=FALSE],
  saliency_log      = saliency_log     [df$train,,drop=FALSE]
)

fit <- rstan::sampling(model, stan_data, chains = 4, cores = 4, warmup = 1000, iter = 2000)
