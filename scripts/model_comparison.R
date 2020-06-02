library(tidyverse)
library(rstan)
library(here)
library(tidybayes)
library(patchwork)
library(imager)
library(ggforce)

ggplot2::theme_set(ggplot2::theme_classic(base_size = 14))
ggplot2::theme_update(axis.ticks.length = ggplot2::unit(6, "pt"), 
                      axis.text         = ggplot2::element_text(size = 15), 
                      axis.title        = ggplot2::element_text(size = 18))

# Load data -----
load(here::here("data", "cleaned_data.Rdata"))

# load helper functions
#source(here::here("R", "expose_helpers_stan.R"))
source(here::here("R", "colours.R"))
source(here::here("R", "load_image.R"))

# get hold out data only
df_sub <- subset(df, !train)
df_sub <- dplyr::mutate(df_sub, obs = 1:nrow(df_sub))

# create list from data to pass to Stan
stan_data <- list(
  N_obs             = nrow(df_sub),
  order             = df_sub$order,
  x                 = df_sub$x,
  y                 = df_sub$y,
  duration          = df_sub$duration,
  
  N_obj             = nrow(objects),
  obj_center_x      = objects$x,
  obj_center_y      = objects$y,
  obj_width         = objects$width,
  obj_height        = objects$height,
  N_ppt             = dplyr::n_distinct(df_sub$id_ppt),
  id_ppt            = df_sub$id_ppt,
  N_img             = dplyr::n_distinct(df_sub$id_img),
  id_img            = df_sub$id_img,
  obj_index_from    = objects_in_images$from,
  obj_index_to      = objects_in_images$to,
  N_obj_in_img      = objects_in_images$n,
  log_lik_saliency  = df_sub$log_lik_saliency,
  max_neighbors     = ncol(saliency_log),
  N_neighbors       = df_sub$n_neighbors,
  mean_sq_distances = mean_sq_distances[!df$train,,drop=FALSE],
  saliency_log      = saliency_log     [!df$train,,drop=FALSE],
  
  N_pix             = max(saliency_normalized$idx),
  half_width_pixel  = 0.5 * 800 / max(saliency_normalized$row),
  saliency_center_x = saliency_normalized$x[saliency_normalized$id_img == 1],
  saliency_center_y = saliency_normalized$y[saliency_normalized$id_img == 2],
  saliency          = lapply(unique(df_sub$id_img), function(id) {
    subset(saliency_normalized, subset = id_img == id, select = "value_normalized", drop = TRUE)
  })
)


# Log-likelihood of hold out data given the original four factor model ----
log_lik_model <- rstan::stan_model(here::here("stan", "log_lik_objects_central_distance_saliency.stan"))

load(here::here("saves", "fit_model.Rdata"))
mcmc <- as.data.frame(fit)
mcmc <- mcmc %>% dplyr::select(sigma_center, sigma_distance, scale_obj, 
                               dplyr::starts_with("weights"), 
                               dplyr::starts_with("z_weights_obj"),
                               dplyr::starts_with("log_weights"),
                               dplyr::starts_with("alpha"), 
                               dplyr::starts_with("sigma_attention"))

log_lik <- rstan::gqs(log_lik_model, data = stan_data, draws = mcmc)


# Log-likelihood of hold out data given the five factor model (additionally horizontal bias) ----
log_lik_model_horizontal <- rstan::stan_model(here::here("stan", "log_lik_objects_central_distance_saliency_horizontal.stan"))

load(here::here("saves", "fit_model_horizontal.Rdata"))
mcmc <- as.data.frame(fit)
mcmc <- mcmc %>% dplyr::select(sigma_center, sigma_distance, scale_obj, kappa,
                               dplyr::starts_with("weights"), 
                               dplyr::starts_with("z_weights_obj"),
                               dplyr::starts_with("log_weights"),
                               dplyr::starts_with("alpha"), 
                               dplyr::starts_with("sigma_attention"))

log_lik_horizontal <- rstan::gqs(log_lik_model_horizontal, data = stan_data, draws = mcmc)


# compare the two models ----
comparison <- data.frame(model_1 = rstan::extract(log_lik)$log_lik, 
                         model_2 = rstan::extract(log_lik_horizontal)$log_lik)
comparison$difference <- comparison$model_2 - comparison$model_1

hist(comparison$difference, breaks = 50, 
     main = "Out of sample comparison", 
     xlab = "Log-likelihood model 2 - Log-likelihood model 1")
summary(comparison$difference)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -115.45   15.18   45.77   45.84   77.06  216.51

boxplot(comparison$model_1, comparison$model_2)
