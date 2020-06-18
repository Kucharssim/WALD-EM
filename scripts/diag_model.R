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
# |           |      mean|   se_mean|        sd|      2.5%|       25%|       50%|       75%|     97.5%|    n_eff|      Rhat|
# |:----------|---------:|---------:|---------:|---------:|---------:|---------:|---------:|---------:|--------:|---------:|
# |weights[1] | 0.3631689| 0.0001015| 0.0118352| 0.3399943| 0.3552289| 0.3631090| 0.3711690| 0.3863705| 13596.73| 0.9996897|
# |weights[2] | 0.1749622| 0.0000790| 0.0098290| 0.1560483| 0.1682691| 0.1748811| 0.1817028| 0.1939706| 15463.76| 0.9997381|
# |weights[3] | 0.3336633| 0.0000641| 0.0077985| 0.3185108| 0.3284412| 0.3336866| 0.3388373| 0.3488614| 14809.48| 0.9995207|
# |weights[4] | 0.1282055| 0.0000690| 0.0091851| 0.1102108| 0.1219797| 0.1281050| 0.1344478| 0.1462483| 17733.61| 0.9994569|

summary(fit, pars = "scale_obj")$summary
# |          |      mean|  se_mean|        sd|      2.5%|       25%|       50%|       75%|     97.5%|   n_eff|      Rhat|
# |:---------|---------:|--------:|---------:|---------:|---------:|---------:|---------:|---------:|-------:|---------:|
# |scale_obj | 0.2349223| 3.32e-05| 0.0042958| 0.2265693| 0.2320579| 0.2348616| 0.2378543| 0.2433253| 16697.5| 0.9996465|

summary(fit, pars = c("sigma_center", "sigma_distance"))$summary
# |               |      mean|   se_mean|        sd|     2.5%|      25%|      50%|       75%|     97.5%|    n_eff|      Rhat|
# |:--------------|---------:|---------:|---------:|--------:|--------:|--------:|---------:|---------:|--------:|---------:|
# |sigma_center   | 100.41068| 0.0199202| 2.7578825| 95.00721| 98.55171| 100.4183| 102.28209| 105.83830| 19167.42| 0.9994467|
# |sigma_distance |  36.12298| 0.0047290| 0.7786392| 34.64187| 35.59335|  36.1099|  36.63785|  37.68548| 27109.87| 0.9993467|

summary(fit, pars = c("mu_log_alpha",           "sigma_log_alpha",
                      "mu_log_sigma_attention", "sigma_log_sigma_attention"))$summary
# |                          |      mean|   se_mean|        sd|      2.5%|       25%|       50%|       75%|     97.5%|    n_eff|      Rhat|
# |:-------------------------|---------:|---------:|---------:|---------:|---------:|---------:|---------:|---------:|--------:|---------:|
# |mu_log_alpha              | 0.1302801| 0.0001486| 0.0140506| 0.1030883| 0.1208858| 0.1304210| 0.1395109| 0.1582818| 8937.828| 0.9997598|
# |sigma_log_alpha           | 0.0727955| 0.0001893| 0.0131703| 0.0492144| 0.0637067| 0.0721529| 0.0809668| 0.1005757| 4839.276| 0.9998785|
# |mu_log_sigma_attention    | 4.0962810| 0.0006375| 0.0439878| 4.0115306| 4.0669045| 4.0959970| 4.1254156| 4.1830465| 4761.457| 1.0014200|
# |sigma_log_sigma_attention | 0.2594203| 0.0005313| 0.0355989| 0.1979296| 0.2341522| 0.2570189| 0.2808913| 0.3389441| 4489.031| 1.0010596|

summary(fit, pars = "alpha")$summary
rstan::stan_plot(fit, pars = "alpha")
summary(fit, pars = "sigma_attention")$summary
rstan::stan_plot(fit, pars = "sigma_attention")
summary(fit, pars = "z_weights_obj")$summary
