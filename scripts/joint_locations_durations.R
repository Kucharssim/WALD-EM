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
  # one-sided bayes factor for correlations with prior: stretched beta(10, 10)<0, 1>
  if(length(n) == 1) n <- rep(n, length(r))
  out <- numeric(length = length(r))
  for(i in 1:length(r)){
    out[i] <- bstats::bcor.testSumStat(n[i], r[i], "greater", kappa = 0.1)[["greater"]][["bf"]]
  }
  
  out
}

# load data and fitted model
load(here::here("data", "cleaned_data.Rdata"))
load(here::here("saves", "fit_model.Rdata"))

# expose stan functions
#source(here::here("R", "expose_helpers_stan.R"))
source(here::here("R", "colours.R"))
source(here::here("R", "load_image.R"))

# create list from data to pass to Stan
# df_sub <- subset(df, !train)
df_sub <- df
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
  mean_sq_distances = mean_sq_distances,#[!df$train,,drop=FALSE],
  saliency_log      = saliency_log     ,#[!df$train,,drop=FALSE],
  
  N_pix             = max(saliency_normalized$idx),
  half_width_pixel  = 0.5 * 800 / max(saliency_normalized$row),
  saliency_center_x = saliency_normalized$x[saliency_normalized$id_img == 1],
  saliency_center_y = saliency_normalized$y[saliency_normalized$id_img == 2],
  saliency          = lapply(unique(df_sub$id_img), function(id) {
    subset(saliency_normalized, subset = id_img == id, select = "value_normalized", drop = TRUE)
  })
)

# recalculate drift rates for each observation
nu_model <- rstan::stan_model(here::here("stan", "nu_objects_central_distance_saliency.stan"))

mcmc <- as.data.frame(fit)
#mcmc <- mcmc[, c(1:360, 363:457)]
mcmc <- mcmc %>% dplyr::select(sigma_center, sigma_distance, scale_obj,
                               dplyr::starts_with("weights"),
                               dplyr::starts_with("z_weights_obj"),
                               dplyr::starts_with("log_weights"),
                               dplyr::starts_with("alpha"),
                               dplyr::starts_with("sigma_attention"))
mcmc <- mcmc %>% dplyr::sample_n(size = 1000) # generate 1000 predictives for every data point

nu_generated <- rstan::gqs(nu_model, data = stan_data, draws = mcmc)

data <- data.frame(
  obs    = seq_len(nrow(df_sub)),
  id_ppt = as.factor(df_sub$id_ppt),
  id_img = as.factor(df_sub$id_img),
  train  = df_sub$train,
  x      = df_sub$x,
  y      = df_sub$y,
  duration = stan_data$duration,
  mean_duration_rep = summary(nu_generated, "duration_mean")$summary[, "mean"] %>% unname()
)

# generate posterior predictives
gqs_model <- rstan::stan_model(here::here("stan", "gqs_objects_central_distance_saliency.stan"))
posterior_predictives <- rstan::gqs(gqs_model, data = stan_data, draws = mcmc)
posterior_predictives <- rstan::extract(posterior_predictives)

xy_rep <- list()
for(i in c("x_rep", "y_rep")) {
  xy_rep[[i]] <- posterior_predictives[[i]] %>% 
    as_tibble() %>% 
    pivot_longer(cols = everything(), 
                 names_to = "obs", values_to = i, 
                 names_prefix = "V") %>%
    mutate(obs = as.integer(obs)) %>%
    group_by(obs) %>%
    mutate(.iter = seq_len(n())) %>%
    ungroup() %>%
    subset(.iter <= 100)
}
xy_rep <- full_join(xy_rep[['x_rep']], xy_rep[['y_rep']])
xy_rep <- left_join(
  xy_rep,
  data %>% select(obs, id_ppt, id_img, train)
)
rm(fit, mcmc, stan_data, saliency_log, nu_generated, posterior_predictives, gqs_model) # unload memory a little

data_combined <- data %>% 
  mutate(group = ifelse(train, "In sample", "Out of sample")) %>%
  bind_rows(data %>% mutate(group = "Combined")) %>%
  mutate(group = factor(group, levels = c("In sample", "Out of sample", "Combined")))

# data %>% 
#   mutate(train = ifelse(train, "In sample", "Out of sample")) %>%
#   bind_rows(data %>% mutate(train = "Combined")) %>%
#   mutate(train = factor(train, levels = c("In sample", "Out of sample", "Combined"))) %>%
data_combined %>%
  group_by(id_ppt, group) %>% 
  summarise(cor = cor(duration, mean_duration_rep), 
            n = n(),
            p.value = cor.test(duration, mean_duration_rep)$p.value, 
            log_bf = log(bfcor(cor, n))) %>%
  ungroup() %>%
  # group_by(group) %>%
  # summarise(percent_cor_positive         = mean(cor > 0),
  #           mean_cor                     = mean(cor),
  #           percent_alt = mean(log_bf > 0),
  #           percent_alt_3 = mean(log_bf > log(3)),
  #           percent_nul_3 = mean(log_bf < log(1/3)),
  #           total_alternative_log   = sum(log_bf),
  #           total_alternative       = exp(total_alternative_log))
  # arrange(train, cor) %>% 
  # print(n = 200) %>%
  ggplot(aes(x = cor, y = log_bf)) +
    geom_abline(slope = 0, intercept = 0, size = 0.5, linetype = 2) +
    geom_abline(slope = 0, intercept = log(c(1/3, 3)), size = 0.5, linetype = 3) + 
    #geom_abline(slope = 1e6, intercept = 0, size = 0.5, linetype = 2) +
    geom_point() +
    geom_rug() +
    scale_x_continuous(breaks = seq(0, 0.5, by = 0.1), labels = gsub("0.", ".", 0:5/10)) +
    xlab("Cor(predicted vs. observed fixation duration)") + 
    ylab(expression(log (BF[+0]))) +
    facet_wrap(~group)
# train           percent_cor_positive mean_cor percent_alt percent_alt_3 percent_nul_3 total_alternative_log total_alternative
#   <fct>                        <dbl>    <dbl>       <dbl>         <dbl>         <dbl>                 <dbl>             <dbl>
# 1 In sample                    0.936    0.137       0.596         0.362        0.0426                  54.9           6.65e23
# 2 Out of sample                0.957    0.163       0.723         0.617        0.0426                  92.7           1.86e40
# 3 Combined                     1        0.148       0.830         0.638        0.0213                 142.            7.36e61
ggsave(filename = here("figures/fit_model/joint_bf.png"), width = 10, height = 6)

data_combined %>%
  ggplot(aes(x = mean_duration_rep, y = duration, color = as.factor(id_ppt))) +
  geom_point(alpha = 0.15, size = 0.15) +
  geom_smooth(method = "lm", se = FALSE, size = 0.5) +
  ylim(NA, 1) + xlim(NA, 1) +
  xlab("Posterior predictives of Fixation duration (sec)") + 
  ylab("Fixation duration (sec)") +
  theme(legend.position = "none") +
  facet_wrap(~group)
ggsave(filename = here("figures/fit_model/durations_scatter.png"), width = 10, height = 6)


plyr::ddply(data_combined, c("id_ppt","group"), function(d) {
  fit <- lm(duration~mean_duration_rep, d)
  broom::tidy(fit)
}) %>% 
  filter(term == "mean_duration_rep") %>%
  print(digits=3)



# make tiles
get_tile <- function(fix, width = 50, min = 0, max = 800) {
  cut(fix, breaks = seq(min, max, by = width), labels = seq_len(ceiling((max-min)/width)))
}
data_combined$tile_x <- get_tile(data_combined$x, 50, 0, 800)
data_combined$tile_y <- get_tile(data_combined$y, 50, 0, 600)
data_combined$tile   <- interaction(data_combined$tile_x, data_combined$tile_y, sep = "_")

xy_rep$tile_x <- get_tile(xy_rep$x_rep, 50, 0, 800)
xy_rep$tile_y <- get_tile(xy_rep$y_rep, 50, 0, 600)
xy_rep$tile   <- interaction(xy_rep$tile_x, xy_rep$tile_y, sep = "_")

# calculate probability of fixating tiles per image based on the model
probs_tile <- plyr::ddply(xy_rep, c("id_img"), 
                          function(d) {
                            # browser()
                            data.frame(tile = levels(xy_rep$tile),
                                       freq = as.vector(table(d$tile))) %>%
                              mutate(prob     = freq/nrow(d)) %>%
                              mutate(log_prob = log(prob))
                          })

data_combined_tiled <- data_combined %>% 
  group_by(id_img, group, tile) %>%
  summarise(mean_duration = mean(duration),
            mean_duration_rep = mean(mean_duration_rep)) %>%
  ungroup() %>%
  left_join(probs_tile, by = c("id_img", "tile"))

data_combined_tiled %>%
  group_by(id_img, group) %>%
  summarise(cor = cor(mean_duration, mean_duration_rep), 
            n = n(),
            p.value = cor.test(mean_duration, mean_duration_rep)$p.value,
            log_bf = log(bfcor(cor, n))) %>%
  arrange(group, cor) %>%
  print(n = 100) %>%
  ggplot(aes(x = cor, y = log_bf)) +
  geom_abline(slope = 0, intercept = 0, size = 0.5, linetype = 2) +
  geom_abline(slope = 0, intercept = log(c(1/3, 3)), size = 0.5, linetype = 3) + 
  #geom_abline(slope = 1e6, intercept = 0, size = 0.5, linetype = 2) +
  geom_point() +
  geom_rug() +
  scale_x_continuous(breaks = seq(0, 0.5, by = 0.1), labels = gsub("0.", ".", 0:5/10)) +
  xlab("Cor(predicted vs. observed mean fixation duration)") + 
  ylab(expression(log (BF[+0]))) +
  facet_wrap(~group)
ggsave(filename = here("figures/fit_model/joint_bf_tiled.png"), width = 10, height = 6)


data_combined_tiled %>%
  ggplot(aes(x = mean_duration_rep, y = mean_duration, color = as.factor(id_img))) +
  geom_point(alpha = 0.15, size = 0.15) +
  geom_smooth(method = "lm", se = FALSE, size = 0.5) +
  ylim(NA, 1) + xlim(NA, 1) +
  xlab("Posterior predictives of mean fixation duration (sec)") + 
  ylab("Mean fixation duration (sec)") +
  theme(legend.position = "none") +
  facet_wrap(~group)
ggsave(filename = here("figures/fit_model/durations_scatter_tiled.png"), width = 10, height = 6)



## probabilities of fixations of particular tiles versus mean fixation durations
data_combined_tiled %>%
  group_by(id_img, group) %>%
  summarise(cor = cor(log_prob, mean_duration), 
            n = n(),
            p.value = cor.test(log_prob, mean_duration)$p.value,
            log_bf = log(bfcor(cor, n))) %>%
  arrange(group, cor) %>%
  print(n = 100) %>%
  ggplot(aes(x = cor, y = log_bf)) +
  geom_abline(slope = 0, intercept = 0, size = 0.5, linetype = 2) +
  geom_abline(slope = 0, intercept = log(c(1/3, 3)), size = 0.5, linetype = 3) + 
  #geom_abline(slope = 1e6, intercept = 0, size = 0.5, linetype = 2) +
  geom_point() +
  geom_rug() +
  scale_x_continuous(breaks = seq(0, 0.5, by = 0.1), labels = gsub("0.", ".", 0:5/10)) +
  xlab("Cor(log probability vs. log mean fixation duration)") + 
  ylab(expression(log (BF[+0]))) +
  facet_wrap(~group)
ggsave(filename = here("figures/fit_model/durations_on_location_prob_bf.png"), width = 10, height = 6)



p1 <- data_combined_tiled %>%
  ggplot(aes(x = log_prob, y = log(mean_duration), color = as.factor(id_img))) +
  geom_point(alpha = 0.15, size = 0.15) +
  geom_smooth(method = "lm", se = FALSE, size = 0.5) +
  ylim(NA, 1) +
  ylab("") + xlab("") +
  ggtitle("Data") + 
  theme(legend.position = "none") +
  facet_wrap(~group)

p2 <- data_combined_tiled %>%
  ggplot(aes(x = log_prob, y = log(mean_duration_rep), color = as.factor(id_img))) +
  geom_point(alpha = 0.15, size = 0.15) +
  geom_smooth(method = "lm", se = FALSE, size = 0.5) +
  ylim(NA, 1) +
  ylab("") + xlab("") +
  ggtitle("Model") + 
  theme(legend.position = "none") +
  facet_wrap(~group)

ylabel <- ggplot(data.frame(x=0.5, y=0.5, text = "Log mean fixation duration (sec)"), 
                 aes(x=x,y=y,label=text)) + 
  ylim(0, 1) +
  xlim(0, 1) +
  geom_text(angle = 90, size = 7) +
  theme_void()

xlabel <- ggplot(data.frame(x=0.5, y=0.5, text = "Log probability of fixation"), 
                 aes(x=x,y=y,label=text)) + 
  ylim(0, 1) +
  xlim(0, 1) +
  geom_text(size = 7) +
  theme_void()

( (ylabel | p1 / p2)                  + patchwork::plot_layout(widths = c(1, 20)) ) / 
( (patchwork::plot_spacer() | xlabel) + patchwork::plot_layout(widths = c(1, 20)) ) +
  patchwork::plot_layout(heights = c(20, 1))
ggsave(filename = here("figures/fit_model/durations_on_location_prob_scatter.png"), width = 10, height = 8)
