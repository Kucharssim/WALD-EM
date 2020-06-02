functions{
#include stan/helpers/load_functions.stan
}
data{
  int<lower=1> N_obs;
  real x[N_obs]; // x position of fixation
  real y[N_obs]; // y position of fixation
}
transformed data{
  real<lower=-pi(),upper=pi()> angle[N_obs];
  real<lower=0> radius[N_obs];
  real<lower=0> dist_to_border[N_obs];
  vector[3] res[N_obs]; // [angle, radius, distance to border]
  real x_prev[N_obs];
  real y_prev[N_obs];
  vector[2] mu_angle;
  mu_angle[1] = 0;
  mu_angle[2] = pi();
  
  for(i in 1:N_obs){
    // coordinates of 'previous fixations'
    
    if(i == 1){
    // if it's a first fixation, calculate angle relative to center
      x_prev[i] = 400;
      y_prev[i] = 300;
    } else{
    // otherwise calculate angle relative to the previous fixation  
      x_prev[i] = x[i-1];
      y_prev[i] = y[i-1];
    }
    
    res[i] = calc_angle_border(x[i], y[i], x_prev[i], y_prev[i], 0, 800, 0, 600);
    angle[i]          = res[i][1];
    radius[i]         = res[i][2];
    dist_to_border[i] = res[i][3];
  }
}
parameters{
  real<lower=0> kappa; 
  real<lower=0> sigma;
  simplex[2] weights;
}
model{
  for(i in 1:N_obs) {
    vector[2] log_lik = log(weights); 
    log_lik[1] += direction_bias_lpdf(angle[i] | radius[i], dist_to_border[i], rep_vector(1.0/2.0, 2), mu_angle, rep_vector(kappa, 2));
    
    log_lik[2] += normal_lpdf(x[i] | 400, sigma);
    log_lik[2] += normal_lpdf(y[i] | 300, sigma);
    
    target += log_sum_exp(log_lik);
  }
  
  kappa ~ gamma(10, 1);
  sigma ~ gamma(2, 0.02);
  weights ~ dirichlet(rep_vector(2, 2));
}


