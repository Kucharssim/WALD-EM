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
# mean      se_mean          sd      2.5%       25%       50%       75%     97.5%    n_eff      Rhat
# weights[1] 0.3715439 1.058174e-04 0.012029593 0.3481540 0.3633171 0.3716291 0.3796071 0.3950507 12923.74 0.9997130
# weights[2] 0.1762092 8.423902e-05 0.009677744 0.1576765 0.1696830 0.1760493 0.1826962 0.1956321 13198.41 0.9997379
# weights[3] 0.2963110 6.000721e-05 0.007235289 0.2822240 0.2914471 0.2962750 0.3011586 0.3104434 14538.00 1.0000589
# weights[4] 0.1559359 7.404572e-05 0.009301633 0.1380191 0.1494686 0.1557860 0.1623349 0.1743298 15780.42 0.9994674

summary(fit, pars = "scale_obj")$summary
# mean      se_mean         sd      2.5%       25%       50%       75%     97.5%    n_eff      Rhat
# scale_obj 0.234143 3.459549e-05 0.00413312 0.2260164 0.2313713 0.2341202 0.2369121 0.2422198 14273.06 0.9998836

summary(fit, pars = c("sigma_center", "sigma_distance"))$summary
# mean     se_mean        sd     2.5%      25%      50%      75%    97.5%    n_eff      Rhat
# sigma_center   93.84059 0.020430479 2.6224236 88.64721 92.05886 93.85710 95.60096 98.90841 16475.88 0.9994837
# sigma_distance 34.58099 0.004929172 0.7366515 33.15060 34.07996 34.57912 35.06274 36.06412 22334.50 0.9997621

summary(fit, pars = c("mu_log_alpha",           "sigma_log_alpha",
                      "mu_log_sigma_attention", "sigma_log_sigma_attention"))$summary
summary(fit, pars = "alpha")$summary
rstan::stan_plot(fit, pars = "alpha")
summary(fit, pars = "sigma_attention")$summary
rstan::stan_plot(fit, pars = "sigma_attention")
summary(fit, pars = "z_weights_obj")$summary
