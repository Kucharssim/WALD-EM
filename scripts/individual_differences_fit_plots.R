library(rstan)
library(tidyverse)
library(here)
library(patchwork)
library(progress)
load(here::here("saves/stan_data.Rdata"))
load(here::here("data", "cleaned_data.Rdata"))
load(here::here("saves/fit_model.Rdata"))
pars <- rstan::extract(fit)

# data <- data.frame(obs      = seq_len(stan_data$N_obs),
#                    id_ppt   = stan_data$id_ppt,
#                    id_img   = stan_data$id_img,
#                    duration = stan_data$duration,
#                    x        = stan_data$x,
#                    y        = stan_data$y
#                    )

### in sample predictions ----
load(here::here("saves/posterior_predictives_in_sample.Rdata"))
posterior_predictives <- rstan::extract(posterior_predictives)

data <- subset(df, train)
data <- select(data, id_ppt, id_img, x, y, duration)
data$obs <- seq_len(nrow(data))
data <- data %>% 
  group_by(id_ppt, id_img) %>%
  mutate(order = 1:n())
plots_in_sample <- list()


### fixation durations ----
pp <- list()
durations <- posterior_predictives[['duration_rep']] %>% 
  tibble::as_tibble() %>% 
  tidyr::pivot_longer(cols = dplyr::everything(), values_to = c("duration"), names_to = c("obs"), names_pattern = "V(.*)") %>%
  dplyr::mutate(obs = as.integer(obs)) %>%
  dplyr::group_by(obs) %>%
  dplyr::mutate(.iter = 1:n()) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(obs, .iter)

summary_durations <- durations %>% 
  group_by(obs) %>% 
  summarise(mean  = mean(duration), 
            lower = quantile(duration, 0.2), 
            upper = quantile(duration, 0.8)) %>%
  ungroup() %>%
  left_join(data)

pp[['data_vs_pred']] <- 
summary_durations %>% 
  group_by(id_ppt) %>% 
  summarise(mean_post = mean(mean), mean_data = mean(duration)) %>%
  ggplot(aes(y=mean_post, x=mean_data)) + 
    geom_smooth(method = "lm") + 
    geom_point() + 
    xlab("Data") + 
    ylab("Predictions") + 
    # ggtitle(label = "Fixation durations", subtitle = "Means per participant") +
    # coord_fixed() + 
    coord_equal(xlim = c(0.35, 0.55), ylim = c(0.35, 0.6)) +
    theme_bw()

pp[['fixation_duration_vs_rank']] <-
summary_durations %>% 
  group_by(id_ppt) %>% 
  summarise(mean_post = mean(mean), lower_post = mean(lower), upper_post = mean(upper),
            mean_data = mean(duration), lower_data = quantile(duration, 0.2), upper_data = quantile(duration, 0.8)) %>% 
  mutate(rank_ppt = rank(mean_post)) %>%
  pivot_longer(cols = -c(id_ppt, rank_ppt), names_pattern = "(.*)_(.*)", names_to = c("stat", "source")) %>%
  pivot_wider(names_from = stat, values_from = value) %>% 
  ggplot(aes(x=rank_ppt, y=mean, ymin=lower, ymax=upper, color=source)) + 
    geom_point(position = position_dodge(0.5), size = 1.5) + 
    geom_errorbar(position = position_dodge(0.5)) + 
    ylab("Fixation duration (sec)") + 
    xlab("Participant rank") +
    scale_y_continuous(limits = c(0, NA)) +
    scale_color_discrete(name = "", labels = c("Data", "Predictions")) +
    theme_bw()

# pp$data_vs_pred + pp$fixation_duration_vs_rank
plots_in_sample[['fixation_duration']] <- pp

### saccade amplitudes ----
x <- posterior_predictives[['x_rep']] %>%
  tibble::as_tibble() %>%
  tidyr::pivot_longer(cols = dplyr::everything(), values_to = c("x"), names_to = c("obs"), names_pattern = "V(.*)") %>%
  dplyr::mutate(obs = as.integer(obs)) %>%
  dplyr::group_by(obs) %>%
  dplyr::mutate(.iter = 1:n()) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(obs, .iter)
y <- posterior_predictives[['y_rep']] %>%
  tibble::as_tibble() %>%
  tidyr::pivot_longer(cols = dplyr::everything(), values_to = c("y"), names_to = c("obs"), names_pattern = "V(.*)") %>%
  dplyr::mutate(obs = as.integer(obs)) %>%
  dplyr::group_by(obs) %>%
  dplyr::mutate(.iter = 1:n()) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(obs, .iter)

xy <- left_join(x, y)
rm(x, y)
xy_sub <- subset(xy, .iter < 6)
xy_sub$amplitude <- NA
xy_sub$angle <- NA

pb <- progress::progress_bar$new(total = nrow(xy_sub))
for(i in seq_len(nrow(xy_sub))) {
  obs_now <- xy_sub[i, "obs", drop=TRUE]
  
  if(data[data$obs==obs_now, "order", drop=TRUE] != 1L) {
    d <- subset(data, obs == obs_now-1)
    dx <- xy_sub[i, "x", drop=TRUE] - d$x
    dy <- xy_sub[i, "y", drop=TRUE] - d$y
    xy_sub$amplitude[i] <- sqrt(dx^2 + dy^2)
    # do not forget: y axis is flipped in eye-tracking data, that's why we reverse the y components of the saccade vector
    xy_sub$angle[i] <- atan2(-dy, dx)
  }
  pb$tick()
}

saccades <- xy_sub[complete.cases(xy_sub),]

### out of sample predictions ----
load(here::here("saves/posterior_predictives_out_sample.Rdata"))
posterior_predictives <- rstan::extract(posterior_predictives)

data <- subset(df, !train)
data <- select(data, id_ppt, id_img, x, y, duration)
data$obs <- seq_len(nrow(data))

plots_out_sample <- list()


### fixation durations ----
pp <- list()
durations <- posterior_predictives[['duration_rep']] %>% 
  tibble::as_tibble() %>% 
  tidyr::pivot_longer(cols = dplyr::everything(), values_to = c("duration"), names_to = c("obs"), names_pattern = "V(.*)") %>%
  dplyr::mutate(obs = as.integer(obs)) %>%
  dplyr::group_by(obs) %>%
  dplyr::mutate(.iter = 1:n()) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(obs, .iter)

summary_durations <- durations %>% 
  group_by(obs) %>% 
  summarise(mean  = mean(duration), 
            lower = quantile(duration, 0.2), 
            upper = quantile(duration, 0.8)) %>%
  ungroup() %>%
  left_join(data)

pp[['data_vs_pred']] <- 
  summary_durations %>% 
  group_by(id_ppt) %>% 
  summarise(mean_post = mean(mean), mean_data = mean(duration)) %>%
  ggplot(aes(y=mean_post, x=mean_data)) + 
  geom_smooth(method = "lm") + 
  geom_point() + 
  xlab("Data") + 
  ylab("Predictions") + 
  # ggtitle(label = "Fixation durations", subtitle = "Means per participant") +
  # coord_fixed() + 
  coord_equal(xlim = c(0.35, 0.55), ylim = c(0.35, 0.6)) +
  theme_bw()

pp[['fixation_duration_vs_rank']] <-
  summary_durations %>% 
  group_by(id_ppt) %>% 
  summarise(mean_post = mean(mean), lower_post = mean(lower), upper_post = mean(upper),
            mean_data = mean(duration), lower_data = quantile(duration, 0.2), upper_data = quantile(duration, 0.8)) %>% 
  mutate(rank_ppt = rank(mean_post)) %>%
  pivot_longer(cols = -c(id_ppt, rank_ppt), names_pattern = "(.*)_(.*)", names_to = c("stat", "source")) %>%
  pivot_wider(names_from = stat, values_from = value) %>% 
  ggplot(aes(x=rank_ppt, y=mean, ymin=lower, ymax=upper, color=source)) + 
  geom_point(position = position_dodge(0.5), size = 1.5) + 
  geom_errorbar(position = position_dodge(0.5)) + 
  ylab("Fixation duration (sec)") + 
  xlab("Participant rank") +
  scale_y_continuous(limits = c(0, NA)) +
  scale_color_discrete(name = "", labels = c("Data", "Predictions")) +
  theme_bw()

# pp$data_vs_pred + pp$fixation_duration_vs_rank
plots_out_sample[['fixation_duration']] <- pp


plots_in_sample$fixation_duration$data_vs_pred + 
  plots_in_sample$fixation_duration$fixation_duration_vs_rank +
  plots_out_sample$fixation_duration$data_vs_pred + 
  plots_out_sample$fixation_duration$fixation_duration_vs_rank

ggsave(here::here("figures/fit_model/individual_fixation_durations.png"), width = 7, height = 4)
