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

# check model diagnostics
rstan::check_hmc_diagnostics(fit)
rstan::stan_diag(fit)
rstan::stan_rhat(fit)
rstan::stan_ess(fit)
rstan::stan_mcse(fit)

# time (in hours) elapsed for getting 1000 warmup and 1000 sampling iterations
rstan::get_elapsed_time(fit) / 3600


# dump diagnostic plots to figures/fit_model/par_diagnostics ----
# the plots are: density per chain, histogram, autocorrelation of the chains, traceplots, cdf per chain
mcmc <- tidybayes::tidy_draws(fit)
mcmc$.chain <- as.factor(mcmc$.chain)
par_names <- names(mcmc)[! (startsWith(names(mcmc), ".") | endsWith(names(mcmc), "__") )]

pb <- dplyr::progress_estimated(n = length(par_names), min_time = 0)
for (par in par_names) {
  p_dens <- rstan::stan_dens (fit, pars = par, separate_chains = TRUE) + 
    ggplot2::theme_light() + 
    ggplot2::theme(legend.position = "none")
  p_hist <- rstan::stan_hist (fit, pars = par) + 
    ggplot2::theme_light()
  p_ac   <- rstan::stan_ac   (fit, pars = par) + 
    ggplot2::theme_light()
  p_tr   <- rstan::stan_trace(fit, pars = par) + 
    ggplot2::theme_light() + 
    ggplot2::theme(legend.position = "none")
  plot_data <- mcmc[, c(".chain",  par)]
  names(plot_data) <- c(".chain", "par")
  p_cdf  <- ggplot2::ggplot(data = plot_data, ggplot2::aes(x = par, col = .chain)) +
    ggplot2::stat_ecdf() + 
    ggplot2::theme(legend.position = "none")
  
  p_all <- p_dens + p_hist + p_ac + 
           p_tr + p_cdf + patchwork::plot_spacer() + 
           patchwork::plot_layout(ncol = 3)
  
  pb$tick()$print()
  ggsave(filename = here::here("figures", "fit_model", "par_diagnostics", paste0(par, ".png")), plot = p_all)
}

# summary stats of parameters ----
summary(fit, pars = "weights")$summary
summary(fit, pars = c("sigma_center", "sigma_distance"))$summary
summary(fit, pars = c("mu_log_alpha",           "sigma_log_alpha",
                      "mu_log_sigma_attention", "sigma_log_sigma_attention"))$summary
summary(fit, pars = "alpha")$summary
rstan::stan_plot(fit, pars = "alpha")
summary(fit, pars = "sigma_attention")$summary
rstan::stan_plot(fit, pars = "sigma_attention")
summary(fit, pars = "z_weights_obj")$summary
