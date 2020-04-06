Examples
================
Simon Kucharsky
4/6/2020

## Introduction

This file serves as a supplementary to the article “Dynamic models of
eye movements” and provides step by step explanation of the building and
implementation of the modelling framework for model that describes the
joint distribution of fixation locations and fixation durations.

The models are implemented in Stan (Carpenter et al., 2017)

``` r
writeLines(readLines(here::here("stan", "example_onlySaliency.stan")), )
```

    ## functions{
    ## #include helpers/load_functions.stan
    ## }
    ## data{
    ##   int N_rows;
    ##   int order[N_rows];
    ##   real x[N_rows];
    ##   real y[N_rows];
    ##   real<lower=0> duration[N_rows];
    ##   int<lower=1> N_pixels;
    ##   int<lower=1, upper=N_pixels> id_pixel[N_pixels];
    ##   vector[N_pixels] x_pixel;
    ##   vector[N_pixels] y_pixel;
    ##   vector[N_pixels] log_val_pixel;
    ##   real log_area_pixel;
    ##   int<lower=1,upper=N_pixels> which_pixel[N_rows];
    ## }
    ## transformed data{
    ##   vector[N_pixels] sq_distances[N_rows];
    ##   
    ##   for(i in 1:N_rows){
    ##     sq_distances[i] = square(x_pixel - x[i]) + square(y_pixel - y[i]);
    ##     sq_distances[i] = sq_distances[i]/2;
    ##   }
    ## }
    ## parameters{
    ##   real<lower=0> sigma_attention;
    ##   real<lower=0> alpha;
    ## }
    ## transformed parameters{
    ##   vector[N_rows] nu;
    ##   vector[N_rows] log_lik_xy;
    ##   vector[N_rows] wald_log_lik;
    ##   
    ##   for(i in 1:N_rows){
    ##     log_lik_xy[i] = log_val_pixel[which_pixel[i]] - log_area_pixel;
    ##     
    ##     nu[i] = - log_sum_exp(log_val_pixel - sq_distances[i]/square(sigma_attention));
    ##     
    ##     wald_log_lik[i] = wald_lpdf(duration[i] | alpha, nu[i]);
    ##   }
    ## }
    ## model{
    ##   target += sum(log_lik_xy);
    ##   target += sum(wald_log_lik);
    ##   
    ##   alpha ~ normal(2, 1);
    ##   sigma_attention ~ gamma(2, 0.1);
    ## }

## References

<div id="refs" class="references">

<div id="ref-carpenter2017stan">

Carpenter, B., Gelman, A., Hoffman, M. D., Lee, D., Goodrich, B.,
Betancourt, M., Brubaker, M., Guo, J., Li, P., & Riddell, A. (2017).
Stan: A probabilistic programming language. *Journal of Statistical
Software*, *76*(1).

</div>

</div>
