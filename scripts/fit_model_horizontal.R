# Application of the Dynamic model of eye movements on the Renswoude's object familiarity data
library(rstan)
rstan_options(auto_write = FALSE)
library(here)

# load stan data
load(here::here("saves", "stan_data.Rdata"))

model <- rstan::stan_model(here::here("stan", "objects_central_distance_saliency_horizontal.stan"))

fit <- rstan::sampling(model, stan_data, chains = 10, cores = 10, warmup = 1000, iter = 2000)

save(fit, file = here::here("saves", "fit_model_horizontal.Rdata"))
