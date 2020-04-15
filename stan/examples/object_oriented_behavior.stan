functions{
#include stan/helpers/load_functions.stan
}
data{
  int<lower=1> N_obs; // number of observations
  real x[N_obs]; // x-coordinates of fixations
  real y[N_obs]; // y coordinates of fixations
  real duration[N_obs]; // fixation durations
  
  int<lower=1> N_objects; // number of objects on a scene
  vector[N_objects] objects_center_x; // x-coordinates of object centers
  vector[N_objects] objects_center_y; // y-coordinates of object centers
  vector[N_objects] objects_width; // object widths
  vector[N_objects] objects_height; // object heights
}
transformed data{
  vector[N_objects] objects_width_2 = objects_width / 2; // object radii in x-coordinate
  vector[N_objects] objects_height_2 = objects_height / 2; // object radii in y-coordinate
}
parameters{
  simplex[N_objects] weights; // importance of different objects
  real<lower=0> alpha; // decision boundary
  real<lower=0> delta; // scaling of the objects
  real<lower=0> sigma_attention; // width of attention window
}
transformed parameters{
  vector[N_objects] sigma_x = delta * objects_width_2;
  vector[N_objects] sigma_y = delta * objects_height_2;
}
model{
  for(t in 1:N_obs){
    // compute drift rate at the specific location
    real nu = -log_integral_attention_mixture_2d(x[t], y[t], weights, objects_center_x, sigma_x, objects_center_y, sigma_y, sigma_attention, sigma_attention);
    // model for fixation locations
    target += mixture_normals(x[t], y[t], weights, objects_center_x, sigma_x, objects_center_y, sigma_y);
    // model for fixation durations
    duration[t] ~ wald(alpha, nu);
  }
  
  // priors
  alpha ~ normal(2, 1);
  sigma_attention ~ gamma(2, 0.1);
  delta ~ exponential(1);
  weights ~ dirichlet(rep_vector(2, N_objects));
}

