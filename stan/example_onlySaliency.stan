functions{
#include helpers/load_functions.stan
}
data{
  int N_rows;
  int order[N_rows];
  real x[N_rows];
  real y[N_rows];
  real<lower=0> duration[N_rows];
  int<lower=1> N_pixels;
  int<lower=1, upper=N_pixels> id_pixel[N_pixels];
  vector<lower=1, upper=N_pixels>[N_pixels] x_pixel;
  vector<lower=1, upper=N_pixels>[N_pixels] y_pixel;
  vector[N_pixels] log_val_pixel;
  real log_area_pixel;
  int<lower=1,upper=N_pixels> which_pixel[N_rows];
}
transformed data{
  vector[N_pixels] sq_distances[N_rows];
  
  for(i in 1:N_rows){
    sq_distances[i] = square(x_pixel - x[i]) + square(y_pixel - y[i]);
    sq_distances[i] = sq_distances[i]/2;
  }
}
parameters{
  real<lower=0> sigma_attention;
  real<lower=0> alpha;
}
transformed parameters{
  vector[N_rows] nu;
  vector[N_rows] log_lik_xy;
  vector[N_rows] wald_log_lik;
  
  for(i in 1:N_rows){
    log_lik_xy[i] = log_val_pixel[which_pixel[i]] - log_area_pixel;
    
    nu[i] = - log_sum_exp(log_val_pixel - sq_distances[i]/sigma_attention);
    
    wald_log_lik[i] = wald_lpdf(duration[i] | alpha, nu[i]);
  }
}
model{
  target += sum(log_lik_xy);
  target += sum(wald_log_lik);
  
  alpha ~ normal(1, 2);
  sigma_attention ~ gamma(2, 0.02);
}