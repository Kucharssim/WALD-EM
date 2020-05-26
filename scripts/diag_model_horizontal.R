library(tidyverse)
library(rstan)
library(here)
library(tidybayes)
library(patchwork)

ggplot2::theme_set(ggplot2::theme_light())

# load data and fitted model
load(here::here("data", "cleaned_data.Rdata"))
load(here::here("saves", "fit_model_horizontal.Rdata"))
load(here::here("saves", "stan_data.Rdata"))

# check model diagnostics
rstan::check_hmc_diagnostics(fit)
rstan::stan_diag(fit)
rstan::stan_rhat(fit)
rstan::stan_ess(fit)
rstan::stan_mcse(fit)

# time (in hours) elapsed for getting 1000 warmup and 1000 sampling iterations
rstan::get_elapsed_time(fit) / 3600


# dump diagnostic plots to figures/fit_model_horizontal/par_diagnostics ----
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
  ggsave(filename = here::here("figures", "fit_model_horizontal", "par_diagnostics", paste0(par, ".png")), plot = p_all)
}

# summary stats of parameters ----
summary(fit, pars = "weights")$summary
# mean      se_mean          sd      2.5%       25%       50%       75%     97.5%    n_eff      Rhat
# weights[1] 0.3631689 1.014982e-04 0.011835202 0.3399943 0.3552289 0.3631090 0.3711690 0.3863705 13596.73 0.9996897
# weights[2] 0.1749622 7.904064e-05 0.009828970 0.1560483 0.1682691 0.1748811 0.1817028 0.1939706 15463.76 0.9997381
# weights[3] 0.3336633 6.408281e-05 0.007798506 0.3185108 0.3284412 0.3336866 0.3388373 0.3488614 14809.48 0.9995207
# weights[4] 0.1282055 6.897410e-05 0.009185117 0.1102108 0.1219797 0.1281050 0.1344478 0.1462483 17733.61 0.9994569

summary(fit, pars = "scale_obj")$summary
# mean      se_mean          sd      2.5%       25%       50%       75%     97.5%   n_eff      Rhat
# scale_obj 0.2349223 3.324471e-05 0.004295842 0.2265693 0.2320579 0.2348616 0.2378543 0.2433253 16697.5 0.9996465

summary(fit, pars = c("sigma_center", "sigma_distance"))$summary
# mean     se_mean        sd     2.5%      25%      50%       75%     97.5%    n_eff      Rhat
# sigma_center   100.41068 0.019920211 2.7578825 95.00721 98.55171 100.4183 102.28209 105.83830 19167.42 0.9994467
# sigma_distance  36.12298 0.004729035 0.7786392 34.64187 35.59335  36.1099  36.63785  37.68548 27109.87 0.9993467

summary(fit, pars = c("mu_log_alpha",           "sigma_log_alpha",
                      "mu_log_sigma_attention", "sigma_log_sigma_attention"))$summary
# mean      se_mean         sd       2.5%        25%        50%        75%     97.5%    n_eff      Rhat
# mu_log_alpha              0.13028013 0.0001486204 0.01405058 0.10308830 0.12088575 0.13042101 0.13951092 0.1582818 8937.828 0.9997598
# sigma_log_alpha           0.07279548 0.0001893233 0.01317025 0.04921442 0.06370673 0.07215292 0.08096678 0.1005757 4839.276 0.9998785
# mu_log_sigma_attention    4.09628097 0.0006374739 0.04398781 4.01153057 4.06690452 4.09599705 4.12541557 4.1830465 4761.457 1.0014200
# sigma_log_sigma_attention 0.25942030 0.0005313243 0.03559885 0.19792956 0.23415216 0.25701890 0.28089129 0.3389441 4489.031 1.0010596

summary(fit, pars = "alpha")$summary
rstan::stan_plot(fit, pars = "alpha")
summary(fit, pars = "sigma_attention")$summary
rstan::stan_plot(fit, pars = "sigma_attention")
summary(fit, pars = "z_weights_obj")$summary
