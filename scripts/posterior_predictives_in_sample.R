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


gqs_model <- rstan::stan_model(here::here("stan", "gqs_objects_central_distance_saliency.stan"))

mcmc <- as.data.frame(fit)
#mcmc <- mcmc[, c(1:360, 363:457)]
mcmc <- mcmc %>% dplyr::select(sigma_center, sigma_distance, scale_obj, 
                               dplyr::starts_with("weights"), 
                               dplyr::starts_with("z_weights_obj"),
                               dplyr::starts_with("log_weights"),
                               dplyr::starts_with("alpha"), 
                               dplyr::starts_with("sigma_attention"))
mcmc <- mcmc %>% dplyr::sample_n(size = 40) # generate 100 predictives for every data point

posterior_predictives <- rstan::gqs(gqs_model, data = stan_data, draws = mcmc)

rm(fit, mcmc, stan_data, saliency_log) # unload memory a little
save(posterior_predictives, file = here::here("saves", "posterior_predictives_in_sample.Rdata"))
load(here::here("saves", "posterior_predictives_in_sample.Rdata"))

mcmc_pred <- as.data.frame(posterior_predictives)

# Duration checks ----
duration_rep <- mcmc_pred %>% 
  dplyr::select(dplyr::starts_with("duration"))
duration_rep$iter <- 1:nrow(duration_rep) 
duration_rep <- tidyr::pivot_longer(duration_rep, cols = dplyr::starts_with("duration"), names_prefix = "duration_rep",
                                    names_to = "obs", values_to = "duration")

# there is a long tail of the predictions spanning to about 15 sec
# but the proportion of the predictions that exceed max of the data is relatively small
perc_pred_below_data <- mean(duration_rep$duration < max(df_sub$duration))

p1 <- ggplot2::ggplot(df_sub, ggplot2::aes(x = duration, y = ..density..)) +
  # plot histogram of data
  ggplot2::geom_histogram(col = cols_custom$dark_teal, fill = cols_custom$light_teal, bins = 50) + 
  ggplot2::geom_rug(mapping = ggplot2::aes(x = duration), 
                    inherit.aes = FALSE, alpha = 0.05, length = ggplot2::unit(4, "pt"), sides = "b") +
  # plot density of predictions
  ggplot2::geom_density(data = duration_rep, mapping = ggplot2::aes(x = duration, group = iter),
                        col = cols_custom$mid_trans, alpha = 0.5) +
  ggplot2::geom_density(data = duration_rep, mapping = ggplot2::aes(x = duration), 
                        col = cols_custom$dark, size = 1) +
  ggplot2::xlab("Fixation duration (sec)") +
  ggplot2::ylab("Density") + 
  ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.1), add = c(0, 0)), limits = c(0, max(df_sub$duration))) +
  ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.1), add = c(0, 0)))

p2 <- ggplot2::ggplot(df_sub, ggplot2::aes(x = duration)) +
  # plot exdf of data
  ggplot2::stat_ecdf(col = cols_custom$dark_teal, size = 1, n = 100) + 
  ggplot2::geom_rug(mapping = ggplot2::aes(x = duration), 
                    inherit.aes = FALSE, outside = FALSE, alpha = 0.05, length = ggplot2::unit(4, "pt"), sides = "b") +
  # plot exdf of predictions
  ggplot2::stat_ecdf(data = duration_rep, mapping = ggplot2::aes(x = duration, group = iter),
                     col = cols_custom$mid_trans, alpha = 0.5) +
  ggplot2::stat_ecdf(data = duration_rep, mapping = ggplot2::aes(x = duration), 
                     col = cols_custom$dark, size = 1) +
  ggplot2::xlab("Fixation duration (sec)") +
  ggplot2::ylab("Cumulative probability") + 
  ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.1), add = c(0, 0)), limits = c(0, max(df_sub$duration))) +
  ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0), add = c(0, 0)))

p1_2 <- p1 + p2
p1_2

ggplot2::ggsave(filename = "fixation_durations.jpg", path = here::here("figures", "fit_model", "in_sample"),
                plot = p1_2, width = 8, height = 4)

# X and Y coordinates checks ----
x_rep <- mcmc_pred %>% 
  dplyr::select(dplyr::starts_with("x"))
x_rep$iter <- 1:nrow(x_rep) 
x_rep <- tidyr::pivot_longer(x_rep, cols = dplyr::starts_with("x"), names_prefix = "x_rep",
                             names_to = "obs", values_to = "x")

y_rep <- mcmc_pred %>% 
  dplyr::select(dplyr::starts_with("y"))
y_rep$iter <- 1:nrow(y_rep) 
y_rep <- tidyr::pivot_longer(y_rep, cols = dplyr::starts_with("y"), names_prefix = "y_rep",
                             names_to = "obs", values_to = "y")

xy_rep <- dplyr::left_join(x_rep, y_rep)
xy_rep$obs <- gsub("\\[", "", xy_rep$obs)
xy_rep$obs <- gsub("\\]", "", xy_rep$obs)
xy_rep$obs <- as.integer(xy_rep$obs)
rm(x_rep, y_rep)

xy_rep <- dplyr::full_join(xy_rep, dplyr::select(df_sub, obs, id_ppt, id_img))

pb <- dplyr::progress_estimated(n = dplyr::n_distinct(df_sub$id_img))
for(img in unique(df_sub$id_img)){
  image_name <- paste0(image_key$image[image_key$id_img == img], ".jpg")
  image <- load_image(image_name)
  image <- as.data.frame(image, wide = "c") %>% dplyr::mutate(rgb.val = rgb(c.1, c.2, c.3))
  
  # plot image
  pp_1 <- ggplot2::ggplot(image, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_raster(ggplot2::aes(fill = rgb.val)) +
    ggplot2::scale_fill_identity() + 
    ggplot2::scale_y_continuous(trans = scales::reverse_trans(), limits = c(600, 0), expand = c(0, 0)) +
    ggplot2::scale_x_continuous(limits = c(0, 800), expand = c(0, 0)) + 
    ggplot2::theme_void() +
    ggplot2::coord_fixed() +
    ggplot2::ggtitle("Stimulus")
  
  # plot observed fixations
  pp_2 <- ggplot2::ggplot(subset(df_sub, id_img == img), ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_point(alpha = 0.5, shape = 19, col = cols_custom$dark_teal, fill = cols_custom$light_teal) + 
    ggplot2::scale_y_continuous(trans = scales::reverse_trans(), limits = c(600, 0), expand = c(0, 0)) +
    ggplot2::scale_x_continuous(limits = c(0, 800), expand = c(0, 0)) + 
    ggplot2::theme_void() + 
    ggplot2::coord_fixed() +
    ggplot2::ggtitle("Observed fixations")
  
  # plot predicted fixations
  pp_3 <- ggplot2::ggplot(subset(xy_rep, id_img == img), ggplot2::aes(x = x, y = y)) +
    #ggplot2::stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
    #ggplot2::stat_density_2d(aes(fill = ..level..),   geom = "polygon", col = cols_custom$dark) +
    ggplot2::scale_fill_gradient(low = cols_custom$light, high = cols_custom$dark) + 
    ggplot2::geom_point(alpha = 0.05, shape = 19, col = cols_custom$dark_highlight, fill = cols_custom$light) + 
    ggplot2::scale_y_continuous(trans = scales::reverse_trans(), limits = c(600, 0), expand = c(0, 0)) +
    ggplot2::scale_x_continuous(limits = c(0, 800), expand = c(0, 0)) + 
    ggplot2::theme_void() +
    ggplot2::theme(legend.position = "none") + 
    ggplot2::coord_fixed() +
    ggplot2::ggtitle("Predicted fixations")
  
  # plot objects on the scene
  pp_4 <- ggplot2::ggplot(subset(objects, id_img == img), ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_point(shape = 13) +
    ggforce::geom_ellipse(ggplot2::aes(x0 = x, y0 = y, a = width/2, b = height/2, angle = 0)) +
    ggplot2::scale_y_continuous(trans = scales::reverse_trans(), limits = c(600, 0), expand = c(0, 0)) +
    ggplot2::scale_x_continuous(limits = c(0, 800), expand = c(0, 0)) + 
    ggplot2::theme_void() +
    ggplot2::coord_fixed() +
    ggplot2::ggtitle("Objects")
  
  # plot saliency
  pp_5 <- ggplot2::ggplot(subset(saliency_normalized, id_img == img), ggplot2::aes(x = x-0.5, y = y-0.5, fill = value)) +
    ggplot2::geom_raster() +
    ggplot2::scale_fill_gradient(low = "black", high = "white") + 
    ggplot2::scale_y_continuous(trans = scales::reverse_trans(), limits = c(600, 0), expand = c(0, 0)) +
    ggplot2::scale_x_continuous(limits = c(0, 800), expand = c(0, 0)) + 
    ggplot2::theme_void() +
    ggplot2::theme(legend.position = "none") +
    ggplot2::coord_fixed() +
    ggplot2::ggtitle("Saliency")
  
  # plot exploitation
  xp <- c(200, 600)
  yp <- c(400, 200)
  dat.distance   <- tidyr::expand_grid(x = seq(0, 800, by = 5), y = seq(0, 600, by = 5))
  dat.distance$z1 <- dnorm(dat.distance$x, xp[1], summary_pars["sigma_distance", "mean"]) * dnorm(dat.distance$y, yp[1], summary_pars["sigma_distance", "mean"])
  dat.distance$z2 <- dnorm(dat.distance$x, xp[2], summary_pars["sigma_distance", "mean"]) * dnorm(dat.distance$y, yp[2], summary_pars["sigma_distance", "mean"])
  dat.distance$z  <- dat.distance$z1 + dat.distance$z2
  
  dat.mock <- data.frame(x=c(rnorm(5, xp[1], summary_pars["sigma_distance", "mean"]), rnorm(3, xp[2], summary_pars["sigma_distance", "mean"])),
                         y=c(rnorm(5, yp[1], summary_pars["sigma_distance", "mean"]), rnorm(3, yp[2], summary_pars["sigma_distance", "mean"])),
                         f=1:8)
  pp_6 <- ggplot2::ggplot(dat.distance, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_raster(ggplot2::aes(fill = z)) + 
    ggplot2::geom_path(data = dat.mock, ggplot2::aes(x = x, y = y, col = f)) +
    ggplot2::geom_point(data = dat.mock, ggplot2::aes(x = x, y = y), col = cols_custom$dark_teal) + 
    ggplot2::scale_fill_gradient(low = cols_custom$light_trans, high = cols_custom$dark_highlight) +
    ggplot2::scale_color_gradient(low = cols_custom$dark_teal, high = cols_custom$mid_teal) +
    ggplot2::scale_y_continuous(trans = scales::reverse_trans(), limits = c(600, 0), expand = c(0, 0)) +
    ggplot2::scale_x_continuous(limits = c(0, 800), expand = c(0, 0)) + 
    ggplot2::theme_void() +
    ggplot2::theme(legend.position = "none") +
    ggplot2::coord_fixed() + 
    ggplot2::ggtitle("Exploitation")
  
  # plot central bias
  dat.central <- tidyr::expand_grid(x = seq(0, 800, by = 5), y = seq(0, 600, by = 5))
  dat.central$z <- dnorm(dat.central$x, 400, summary_pars["sigma_center", "mean"]) * dnorm(dat.central$y, 300, summary_pars["sigma_center", "mean"])
  dat.central$z <- dat.central$z / max(dat.central$z)
  
  pp_7 <- ggplot2::ggplot(dat.central, ggplot2::aes(x = x, y = y, fill = z)) +
    ggplot2::geom_raster() +
    ggplot2::scale_fill_gradient(low = cols_custom$light_trans, high = cols_custom$dark_highlight) + 
    ggplot2::scale_y_continuous(trans = scales::reverse_trans(), limits = c(601, 0), expand = c(0, 0)) +
    ggplot2::scale_x_continuous(limits = c(0, 801), expand = c(0, 0)) + 
    ggplot2::theme_void() +
    ggplot2::theme(legend.position = "none") +
    ggplot2::coord_fixed() +
    ggplot2::ggtitle("Central bias")
  
  # stich if together
  pp_fac <- pp_4 + pp_5 + pp_6 + pp_7 + patchwork::plot_layout(ncol = 2)
  pp <- pp_1 + pp_fac + pp_2 + pp_3 + patchwork::plot_layout(ncol = 2)
  
  # save
  ggplot2::ggsave(image_name, pp, path = here::here("figures/fit_model/in_sample/xy"),
                  width = 20, height = 16, units = "cm")
  
  pb$tick()$print()
}


# Saccade amplitude check ----
# calculate aplitudes in data
amplitude_dat <- plyr::ddply(.data = df_sub, .variables = c("id_ppt", "id_img"), .fun = function(d){
  .prev <- 1:(nrow(d)-1)
  .next <- 2:nrow(d) 
  new_d <- data.frame(id_ppt = d$id_ppt[.prev], id_img = d$id_img[.prev],
                      distance = sqrt( (d$x[.prev] - d$x[.next])^2 + (d$y[.prev] - d$y[.next])^2 )
                      )
})

# calculate amplitudes (distances of predictions for the next fixation from the observed fixation)
# xy_rep <- subset(xy_rep, (iter %% 100) == 0)
xy_rep <- dplyr::full_join(xy_rep, dplyr::select(df_sub, obs, id_ppt, id_img))

amplitude_pred <- plyr::ddply(.data = df_sub, .variables = c("id_ppt", "id_img", "obs"), .fun = function(d){
  ppt_d <- d$id_ppt[1]
  img_d <- d$id_img[1]
  obs_d <- d$obs
  x_d   <- d$x
  y_d   <- d$y
  
  pred <- subset(xy_rep, obs == obs_d + 1 & id_ppt == ppt_d & id_img == img_d)
  
  n_row <- nrow(pred)
  if(n_row == 0){
    return(data.frame(id_ppt=integer(), id_img=integer(), distance=numeric()))
  } else{
    out <- data.frame(
      id_ppt = rep(ppt_d[1], n_row),
      id_img = rep(img_d[1], n_row),
      distance = sqrt( (pred$x - x_d)^2 + (pred$y - y_d)^2 )
    )
    
    return(out)
  }
}, .progress = "text")
amplitude_pred$iter <- 1:nrow(amplitude_pred)

p1 <- ggplot2::ggplot(amplitude_dat, ggplot2::aes(x = distance, y = ..density..)) +
  # plot histogram of data
  ggplot2::geom_histogram(col = cols_custom$dark_teal, fill = cols_custom$light_teal, bins = 50) + 
  ggplot2::geom_rug(mapping = ggplot2::aes(x = distance), 
                    inherit.aes = FALSE, alpha = 0.05, length = ggplot2::unit(4, "pt"), sides = "b") +
  # plot density of predictions
  ggplot2::geom_density(data = amplitude_pred, mapping = ggplot2::aes(x = distance), 
                        col = cols_custom$dark, size = 1.5) +
  ggplot2::xlab("Distance (pixels)") +
  ggplot2::ylab("Density") + 
  ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.1), add = c(0, 0))) +
  ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.1), add = c(0, 0)))

p2 <- ggplot2::ggplot(amplitude_dat, ggplot2::aes(x = distance)) +
  # plot exdf of data
  ggplot2::stat_ecdf(col = cols_custom$dark_teal, size = 1, n = 100) + 
  ggplot2::geom_rug(mapping = ggplot2::aes(x = distance), 
                    inherit.aes = FALSE, outside = FALSE, alpha = 0.05, length = ggplot2::unit(4, "pt"), sides = "b") +
  # plot exdf of predictions
  ggplot2::stat_ecdf(data = amplitude_pred, mapping = ggplot2::aes(x = distance), 
                     col = cols_custom$dark, size = 1.5) +
  ggplot2::xlab("Distance (pixels)") +
  ggplot2::ylab("Cumulative probability") + 
  ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.1), add = c(0, 0))) +
  ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0), add = c(0, 0)))

p1_2 <- p1 + p2
p1_2

ggplot2::ggsave(filename = "amplitude.jpg", path = here::here("figures", "fit_model", "in_sample"),
                plot = p1_2, width = 8, height = 4)

pb <- dplyr::progress_estimated(n = dplyr::n_distinct(df_sub$id_img))
for(img in unique(df_sub$id_img)){
  image_name <- paste0(image_key$image[image_key$id_img == img], ".jpg")
  amplitude_dat_sub  <- subset(amplitude_dat,  id_img == img)
  amplitude_pred_sub <- subset(amplitude_pred, id_img == img)
  
  p1 <- ggplot2::ggplot(amplitude_dat_sub, ggplot2::aes(x = distance, y = ..density..)) +
    # plot histogram of data
    ggplot2::geom_histogram(col = cols_custom$dark_teal, fill = cols_custom$light_teal, bins = 30) + 
    ggplot2::geom_rug(mapping = ggplot2::aes(x = distance), 
                      inherit.aes = FALSE, alpha = 0.05, length = ggplot2::unit(4, "pt"), sides = "b") +
    # plot density of predictions
    ggplot2::geom_density(data = amplitude_pred_sub, mapping = ggplot2::aes(x = distance), 
                          col = cols_custom$dark, size = 1.5) +
    ggplot2::xlab("Distance (pixels)") +
    ggplot2::ylab("Density") + 
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.1), add = c(0, 0))) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.1), add = c(0, 0)))
  
  p2 <- ggplot2::ggplot(amplitude_dat_sub, ggplot2::aes(x = distance)) +
    # plot exdf of data
    ggplot2::stat_ecdf(col = cols_custom$dark_teal, size = 1, n = 100) + 
    ggplot2::geom_rug(mapping = ggplot2::aes(x = distance), 
                      inherit.aes = FALSE, outside = FALSE, alpha = 0.05, length = ggplot2::unit(4, "pt"), sides = "b") +
    # plot exdf of predictions
    ggplot2::stat_ecdf(data = amplitude_pred_sub, mapping = ggplot2::aes(x = distance), 
                       col = cols_custom$dark, size = 1.5) +
    ggplot2::xlab("Distance (sec)") +
    ggplot2::ylab("Cumulative probability") + 
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.1), add = c(0, 0))) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0), add = c(0, 0)))
  
  p1_2 <- p1 + p2
  p1_2
  
  ggplot2::ggsave(filename = image_name, plot = p1_2, path = here::here("figures", "fit_model", "in_sample", "amplitude"), 
                  width = 8, height = 5)
  
  pb$tick()$print()
}



# Saccade angle check ----
# atan2: 0   pi - right
#        0.5 pi - up
#        1   pi - left
#       -0.5 pi - down
angle_dat <- plyr::ddply(.data = df_sub, .variables = c("id_ppt", "id_img"), .fun = function(d){
  .prev <- 1:(nrow(d)-1)
  .next <- 2:nrow(d) 
  x <- d$x[.next] - d$x[.prev]
  y <- d$y[.next] - d$y[.prev]
  
  # calculate angles
  # do not forget: y axis is flipped in eye-tracking data, that's why we reverse the y components of the saccade vector 
  new_d <- data.frame(id_ppt = d$id_ppt[.prev], id_img = d$id_img[.prev],
                      angle = atan2(-y, x)
                      )
  
  return(new_d)
})

angle_pred <- plyr::ddply(.data = df_sub, .variables = c("id_ppt", "id_img", "obs"), .fun = function(d){
  ppt_d <- d$id_ppt[1]
  img_d <- d$id_img[1]
  obs_d <- d$obs
  x_d   <- d$x
  y_d   <- d$y
  
  pred <- subset(xy_rep, obs == obs_d + 1 & id_ppt == ppt_d & id_img == img_d)
  
  n_row <- nrow(pred)
  if(n_row == 0){
    return(data.frame(id_ppt=integer(), id_img=integer(), angle=numeric()))
  } else{
    # calculate angles
    # do not forget: y axis is flipped in eye-tracking data, that's why we reverse the y components of the saccade vector 
    x <- pred$x - x_d
    y <- pred$y - y_d
    out <- data.frame(
      id_ppt = rep(ppt_d[1], n_row),
      id_img = rep(img_d[1], n_row),
      angle  = atan2(-y, x)
    )
    
    return(out)
  }
}, .progress = "text")

p1 <- ggplot2::ggplot(angle_dat, ggplot2::aes(x = angle, y = ..density..)) + 
  ggplot2::geom_histogram(col = cols_custom$dark_teal, fill = cols_custom$mid_teal, alpha = 0.8, bins = 32) +
  ggplot2::geom_histogram(data=angle_pred, col = cols_custom$dark, fill = cols_custom$mid, alpha = 0.8, bins = 32) +
  ggplot2::coord_polar(start = 0.5*pi, direction = -1) + 
  ggplot2::scale_x_continuous(limits = c(-pi, pi),
                              breaks = seq(-0.5, 1, by = 0.5)*pi,
                              labels = c("down", "right", "up", "left")) +
  ggplot2::scale_y_continuous(expand = c(0, 0.025)) +
  ggplot2::theme_void() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(size = 15))

ggplot2::ggsave(filename = "angle.jpg", plot = p1, path = here::here("figures", "fit_model", "in_sample"), 
                width = 5, height = 5)

pb <- dplyr::progress_estimated(n = dplyr::n_distinct(df_sub$id_img))
for(img in unique(df_sub$id_img)){
  image_name <- paste0(image_key$image[image_key$id_img == img], ".jpg")
  angle_dat_sub  <- subset(angle_dat,  id_img == img)
  angle_pred_sub <- subset(angle_pred, id_img == img)
  
  p1 <- ggplot2::ggplot(angle_dat_sub, ggplot2::aes(x = angle, y = ..density..)) + 
    ggplot2::geom_histogram(col = cols_custom$dark_teal, fill = cols_custom$mid_teal, alpha = 0.8, bins = 32) +
    ggplot2::geom_histogram(data=angle_pred_sub, col = cols_custom$dark, fill = cols_custom$mid, alpha = 0.8, bins = 32) +
    ggplot2::coord_polar(start = 0.5*pi, direction = -1) + 
    ggplot2::scale_x_continuous(limits = c(-pi, pi),
                                breaks = seq(-0.5, 1, by = 0.5)*pi,
                                labels = c("down", "right", "up", "left")) +
    ggplot2::scale_y_continuous(expand = c(0, 0.025)) +
    ggplot2::theme_void() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size = 15))
  
  ggplot2::ggsave(filename = image_name, plot = p1, path = here::here("figures", "fit_model", "in_sample", "angle"), 
                  width = 5, height = 5)
  
  pb$tick()$print()
}
