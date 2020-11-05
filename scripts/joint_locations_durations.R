library(rstan)
library(tidyverse)
library(here)
library(patchwork)
library(progress)
# library(bstats)
bcor.test <- function(r, n, alternative = c("two.sided", "less", "greater"), method, kappa = 1) {
  alternative <- match.arg(alternative)
  out <- numeric(length = length(r))
  
  sapply(seq_along(r), function(i) bstats::bcor.testSumStat(n = n[i], stat = r[i], kappa = kappa, method = method)[[alternative]][['bf']])
}
load(here::here("saves/stan_data.Rdata"))
load(here::here("data", "cleaned_data.Rdata"))
#load(here::here("saves/fit_model.Rdata"))


### in sample ----
data <- subset(df, train)
data <- select(data, id_ppt, id_img, x, y, duration)
data <- data %>% 
  group_by(id_ppt, id_img) %>%
  mutate(order = 1:n()) %>% 
  ungroup() %>%
  mutate(obs = 1:length(x))

# get predictives
load(here::here("saves/posterior_predictives_in_sample.Rdata"))
posterior_predictives <- rstan::extract(posterior_predictives)
posterior_predictives <- lapply(names(posterior_predictives), function(var) {
  posterior_predictives[[var]] %>% 
    tibble::as_tibble() %>% 
    tidyr::pivot_longer(cols = dplyr::everything(), values_to = "value", names_to = c("obs"), names_pattern = "V(.*)") %>%
    dplyr::mutate(obs = as.integer(obs)) %>%
    dplyr::group_by(obs) %>%
    dplyr::mutate(.iter = 1:n()) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(obs, .iter) %>%
    dplyr::mutate(variable = gsub("_rep", "", var))
}) %>% 
  do.call(rbind, .) %>%
  pivot_wider(names_from = variable, values_from = value)

posterior_predictives$id_ppt <- data$id_ppt[match(posterior_predictives$obs, data$obs)]
posterior_predictives$id_img <- data$id_img[match(posterior_predictives$obs, data$obs)]
plots_in_sample <- list()


# cut images in 16 x 12 tiles
data$tile_x <- cut(data$x, seq(0, 800, by = 100), labels = LETTERS[1:8])
data$tile_y <- cut(data$y, seq(0, 600, by = 100), labels = LETTERS[1:6])
data$tile   <- interaction(data$tile_x, data$tile_y, sep = "_")

posterior_predictives$tile_x <- cut(posterior_predictives$x, seq(0, 800, by = 100), labels = LETTERS[1:8])
posterior_predictives$tile_y <- cut(posterior_predictives$y, seq(0, 600, by = 100), labels = LETTERS[1:6])
posterior_predictives$tile   <- interaction(posterior_predictives$tile_x, posterior_predictives$tile_y, sep = "_")


# compute average fixation duration per tile
data_agg <- data %>%
  group_by(id_img, tile, .drop=FALSE) %>%
  summarise(mean_duration_data = mean(duration), 
            n_data = n(), .groups = "keep") %>%
  group_by(id_img) %>%
  mutate(rank_duration_data = rank(mean_duration_data))

pred_agg <- posterior_predictives %>%
  group_by(id_img, tile, .drop=FALSE) %>%
  summarise(mean_duration_pred = mean(duration),
            n_pred = n(), .groups = "keep") %>%
  group_by(id_img) %>%
  mutate(rank_duration_pred = rank(mean_duration_pred))

agg <- full_join(data_agg, pred_agg)
agg %>%
  ggplot(aes(x = mean_duration_data, y = mean_duration_pred, weight = n_data, color = as.factor(id_img))) +
  geom_smooth(method = "lm", se = FALSE) +
  geom_point() +
  theme_bw() + 
  theme(legend.position = "none") + 
  facet_wrap(~id_img, scales = "free")

agg %>%
  ggplot(aes(x = rank_duration_data, y = rank_duration_pred, weight = n_data, color = as.factor(id_img))) +
  geom_smooth(method = "lm", se = FALSE) +
  geom_point() +
  theme_bw() + 
  theme(legend.position = "none") + 
  facet_wrap(~id_img, scales = "free")

agg_test <- agg %>% 
  group_by(id_img) %>%
  summarise(cor = cor(mean_duration_data, mean_duration_pred, use = "p", method = "kendall"), n = n()) %>%
  mutate(bf = bcor.test(cor, n, "greater", "kendall", 0.5))

# compute number of visits and average fixation duration to tiles for each image separately (data)
get_summaries <- function(d) {
  out <- d %>% 
    group_by(tile) %>%
    summarise(n = n(), 
              mean_duration = mean(duration), 
              sd_duration = sd(duration), 
              .groups = "keep") %>%
    ungroup()
  out$rel_freq <- out$n/sum(out$n) 
  out
}
fix_dur <- plyr::ddply(.data = data, .variables = c("id_img"), .fun = get_summaries, .progress = "text")

cor.test(fix_dur$mean_duration, fix_dur$rel_freq, method = "spearman", exact = FALSE)

ggplot(fix_dur, aes(x = rel_freq, y = mean_duration, weight = n, colour = as.factor(id_img))) +
  geom_smooth(method = "lm", se = FALSE) +
  #geom_point()
  geom_abline(intercept = coef(lm(mean_duration~rel_freq, weights = n, data = fix_dur))[1],
              slope     = coef(lm(mean_duration~rel_freq, weights = n, data = fix_dur))[2]) +
  
  theme_bw() + 
  theme(legend.position = "none")

linear_fits <- plyr::dlply(.data = fix_dur, .variables = "id_img", function(d) lm(mean_duration~rel_freq, weights = n, data = d))
linear_coef <- lapply(linear_fits, coef) %>% do.call(rbind, .) %>% as_tibble()
linear_coef$p.value <- sapply(linear_fits, function(f) summary(f)[['coefficients']]['rel_freq', 'Pr(>|t|)'])

# compute number of visits and averae fixation duration to tiles for each image separately (predictions)
# this is done for each mcmc iteration separately

fix_dur_pred <- plyr::ddply(.data = posterior_predictives, .variables = c(".iter", "id_img"), .fun = get_summaries, .progress = "text")

linear_fits_pred <- plyr::ddply(.data = fix_dur_pred, .variables = c(".iter", "id_img"),
                                function(d) {
                                  fit <- lm(mean_duration~rel_freq, weights = n, data = d)
                                  coef(fit)
                                }) 
colnames(linear_fits_pred) <- c(".iter", "id_img", "intercept", "slope")
linear_fits_pred_summary <- linear_fits_pred %>%
  group_by(id_img) %>%
  summarise(intercept = mean(intercept), slope = mean(slope), .groups = "keep")

ggplot(fix_dur, aes(x = rel_freq, y = mean_duration, weight = n, colour = as.factor(id_img))) +
  geom_smooth(method = "lm", se = FALSE) +
  geom_abline(data = linear_fits_pred, mapping = aes(intercept = intercept, slope = slope, color = as.factor(id_img))) +
  theme_bw() + 
  theme(legend.position = "none") + 
  facet_wrap(.~id_img)
