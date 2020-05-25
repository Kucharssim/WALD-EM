functions{
#include stan/helpers/load_functions.stan
}
data{
  int<lower=1> N_obs;
  real angle[N_obs];
  vector[2] mu_angle;
}

parameters{
  real<lower=0> kappa; 
}
model{
  for(i in 1:N_obs) {
    target += mixture_von_mises_lpdf(angle[i] | rep_vector(1.0/2.0, 2), mu_angle, rep_vector(kappa, 2));
  }
  
  kappa ~ gamma(10, 1);
}
