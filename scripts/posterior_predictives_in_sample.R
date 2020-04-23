library(tidyverse)
library(rstan)
library(here)
library(tidybayes)
library(patchwork)

ggplot2::theme_set(ggplot2::theme_light())

# load data and fitted model
load(here::here("data", "cleaned_data.Rdata"))
load(here::here("saves", "fit_model.Rdata"))
load(here::here("saves", "stan_data.Rdata"))

# expose stan functions
source(here::here("R", "expose_helpers_stan.R"))

gqs_model <- rstan::stan_model(here::here("stan", "gqs_objects_central_distance_saliency.stan"))

mcmc <- as.matrix(fit)
#mcmc <- mcmc[, 1:356]
mcmc <- mcmc[, c(1:360, 363:457)]
posterior_predictives <- rstan::gqs(gqs_model, data = stan_data, draws = mcmc)
