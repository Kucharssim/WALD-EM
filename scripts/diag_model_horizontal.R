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
# |           |      mean|   se_mean|        sd|      2.5%|       25%|       50%|       75%|     97.5%|    n_eff|      Rhat|
# |:----------|---------:|---------:|---------:|---------:|---------:|---------:|---------:|---------:|--------:|---------:|
# |weights[1] | 0.3488414| 0.0001031| 0.0121205| 0.3255308| 0.3405815| 0.3486657| 0.3569780| 0.3731782| 13824.76| 1.0000353|
# |weights[2] | 0.1447130| 0.0000754| 0.0094504| 0.1265281| 0.1384202| 0.1446061| 0.1510627| 0.1636202| 15712.79| 0.9996939|
# |weights[3] | 0.3078565| 0.0000676| 0.0079411| 0.2922260| 0.3025367| 0.3077695| 0.3131324| 0.3234867| 13790.37| 0.9999881|
# |weights[4] | 0.1277102| 0.0000716| 0.0090563| 0.1102800| 0.1215534| 0.1275829| 0.1337922| 0.1456998| 15989.80| 0.9996954|
# |weights[5] | 0.0708789| 0.0000484| 0.0067934| 0.0579690| 0.0661168| 0.0707720| 0.0754724| 0.0843726| 19708.92| 0.9995330|

summary(fit, pars = "scale_obj")$summary
# |          |      mean|  se_mean|        sd|      2.5%|       25%|       50%|       75%|    97.5%|    n_eff|      Rhat|
# |:---------|---------:|--------:|---------:|---------:|---------:|---------:|---------:|--------:|--------:|---------:|
# |scale_obj | 0.2330572| 3.45e-05| 0.0043261| 0.2244685| 0.2301382| 0.2330579| 0.2359779| 0.241517| 15740.23| 0.9994808|

summary(fit, pars = c("sigma_center", "sigma_distance"))$summary
# |               |     mean|   se_mean|        sd|     2.5%|      25%|      50%|       75%|     97.5%|    n_eff|      Rhat|
# |:--------------|--------:|---------:|---------:|--------:|--------:|--------:|---------:|---------:|--------:|---------:|
# |sigma_center   | 98.67865| 0.0192172| 2.6928090| 93.42119| 96.86374| 98.63371| 100.50039| 103.99425| 19634.94| 0.9997780|
# |sigma_distance | 34.27622| 0.0051063| 0.7696197| 32.79543| 33.73652| 34.28326|  34.80042|  35.78228| 22716.70| 0.9993995|

summary(fit, pars = "kappa")$summary
# |      |     mean|   se_mean|       sd|     2.5%|      25%|      50%|      75%|    97.5%|    n_eff|      Rhat|
# |:-----|--------:|---------:|--------:|--------:|--------:|--------:|--------:|--------:|--------:|---------:|
# |kappa | 18.35863| 0.0183774| 2.739641| 13.40475| 16.41876| 18.20802| 20.10992| 24.15781| 22223.89| 0.9992936|

summary(fit, pars = c("mu_log_alpha",           "sigma_log_alpha",
                      "mu_log_sigma_attention", "sigma_log_sigma_attention"))$summary
# |                          |      mean|   se_mean|        sd|      2.5%|       25%|       50%|       75%|     97.5%|    n_eff|     Rhat|
# |:-------------------------|---------:|---------:|---------:|---------:|---------:|---------:|---------:|---------:|--------:|--------:|
# |mu_log_alpha              | 0.1276052| 0.0001552| 0.0141223| 0.1002310| 0.1181214| 0.1277523| 0.1371373| 0.1549365| 8281.602| 1.000401|
# |sigma_log_alpha           | 0.0717975| 0.0001780| 0.0129484| 0.0485140| 0.0628703| 0.0709617| 0.0799575| 0.0996025| 5293.341| 1.000735|
# |mu_log_sigma_attention    | 4.0918574| 0.0006175| 0.0431854| 4.0069387| 4.0631287| 4.0921738| 4.1199279| 4.1768883| 4891.813| 1.001044|
# |sigma_log_sigma_attention | 0.2550026| 0.0004926| 0.0341719| 0.1960123| 0.2307909| 0.2526294| 0.2756957| 0.3289220| 4813.112| 1.001899|


summary(fit, pars = "alpha")$summary
rstan::stan_plot(fit, pars = "alpha")
summary(fit, pars = "sigma_attention")$summary
rstan::stan_plot(fit, pars = "sigma_attention")
summary(fit, pars = "z_weights_obj")$summary
