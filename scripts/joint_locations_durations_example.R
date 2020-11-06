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

bfcor <- function(r, n) {
  bstats::bcor.testSumStat(n, r, "greater", kappa = 0.5)
}
# load data and fitted model
load(here::here("data", "cleaned_data.Rdata"))
load(here::here("saves", "fit_model.Rdata"))
# load(here::here("saves", "stan_data.Rdata"))

summary_pars <- summary(fit)$summary
# expose stan functions
#source(here::here("R", "expose_helpers_stan.R"))
source(here::here("R", "colours.R"))
source(here::here("R", "load_image.R"))
# create list from data to pass to Stan
df_sub <- subset(df, train)
df_sub <- dplyr::mutate(df_sub, obs = 1:nrow(df_sub))

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
  mean_sq_distances = mean_sq_distances[df$train,,drop=FALSE],
  saliency_log      = saliency_log     [df$train,,drop=FALSE],
  
  N_pix             = max(saliency_normalized$idx),
  half_width_pixel  = 0.5 * 800 / max(saliency_normalized$row),
  saliency_center_x = saliency_normalized$x[saliency_normalized$id_img == 1],
  saliency_center_y = saliency_normalized$y[saliency_normalized$id_img == 2],
  saliency          = lapply(unique(df_sub$id_img), function(id) {
    subset(saliency_normalized, subset = id_img == id, select = "value_normalized", drop = TRUE)
  })
)


gqs_model <- rstan::stan_model(here::here("stan", "nu_objects_central_distance_saliency.stan"))

mcmc <- as.data.frame(fit)
#mcmc <- mcmc[, c(1:360, 363:457)]
mcmc <- mcmc %>% dplyr::select(sigma_center, sigma_distance, scale_obj,
                               dplyr::starts_with("weights"),
                               dplyr::starts_with("z_weights_obj"),
                               dplyr::starts_with("log_weights"),
                               dplyr::starts_with("alpha"),
                               dplyr::starts_with("sigma_attention"))
mcmc <- mcmc %>% dplyr::sample_n(size = 100) # generate 100 predictives for every data point

posterior_predictives <- rstan::gqs(gqs_model, data = stan_data, draws = mcmc)

# rm(fit, mcmc, stan_data, saliency_log) # unload memory a little

cor(
  summary(posterior_predictives, "duration_mean")$summary[, "mean"],
  stan_data$duration
)

data <- data.frame(
  id_ppt = as.factor(stan_data$id_ppt),
  id_img = as.factor(stan_data$id_img),
  x      = stan_data$x,
  y      = stan_data$y,
  duration = stan_data$duration,
  mean_duration_rep = summary(posterior_predictives, "duration_mean")$summary[, "mean"] %>% unname()
)


data %>% 
  group_by(id_ppt) %>% 
  summarise(cor = cor(duration, mean_duration_rep), n = n(),
            p.value = cor.test(duration, mean_duration_rep)$p.value) %>%
  summary()
  print(n = 100)

data %>%
  ggplot(aes(mean_duration_rep, duration, color = as.factor(id_ppt))) +
  geom_smooth(method = "lm", se = FALSE) +
  geom_point(alpha = 0.3, size = 0.1) +
  ylim(NA, 1) + xlim(NA, 1) +
  theme_bw() + 
  theme(legend.position = "none")

plyr::ddply(data, "id_ppt", function(d) {
  fit <- lm(duration~mean_duration_rep, d)
  broom::tidy(fit)
}) %>% 
  filter(term == "mean_duration_rep") %>%
  print(digits=3)

# fit_brm <- brms::brm(duration~(duration|id_ppt), data = data, iter=50, warmup =500)
# fit_brm


# make tiles
get_tile <- function(fix, width = 50, min = 0, max = 800) {
  cut(fix, breaks = seq(min, max, by = width), labels = seq_len(ceiling((max-min)/width)))
}
data$tile_x <- get_tile(data$x, 25, 0, 800)
data$tile_y <- get_tile(data$y, 25, 0, 600)
data$tile   <- interaction(data$tile_x, data$tile_y, sep = "_")

data_tiled <- data %>% 
  group_by(id_img, tile) %>%
  summarise(mean_duration = mean(duration),
            mean_duration_rep = mean(mean_duration_rep))

data_tiled %>%
  group_by(id_img) %>%
  summarise(cor = cor(mean_duration, mean_duration_rep), 
            n = n(),
            p.value = cor.test(mean_duration, mean_duration_rep)$p.value,
            bf = bfcor(cor, n)[["greater"]][["bf"]]) %>%
  summary()

plot(data_tiled$mean_duration_rep, data_tiled$mean_duration)
