// hierarchical wald-mixture model
functions{
#include stan/helpers/load_functions.stan
}
data{
  int N_obs; // total number of observations
  int order[N_obs]; // number that indexed observation within ppt and img
  real x[N_obs]; // x position of fixation
  real y[N_obs]; // y position of fixation
  real<lower=0> duration[N_obs]; // fixation duration
  int N_obj; // total number of objects across all scenes
  vector[N_obj] obj_center_x; // x position of the center of objects
  vector[N_obj] obj_center_y; // y position of the center of objects
  vector[N_obj] obj_width; // width of objects
  vector[N_obj] obj_height; // height of objects
  int N_ppt; // number of participants
  int<lower=1,upper=N_ppt> id_ppt[N_obs]; // number that indexes the participant
  int N_img; // number of images (stimuli)
  int<lower=1,upper=N_img> id_img[N_obs]; // number that indexed the image (stimulus)
  int obj_index_from[N_img]; // the number that starts the indexes of objects that belong to an image
  int obj_index_to[N_img]; // the number that ends the indexes of objects that belong to an image
  int N_obj_in_img[N_img]; // the number of objects in each stimulus
  real log_lik_saliency[N_obs]; // the log likelihood of (x,y) given the saliency
  int max_neighbors;
  int<lower=1,upper=max_neighbors> N_neighbors[N_obs]; // the number of closest pixels that we take into account when calculating the drift rate
  vector[max_neighbors] mean_sq_distances[N_obs];      // distances of the neighbors to the fixation locations
  vector[max_neighbors] saliency_log[N_obs];           // log of saliences of the neighbors
}
transformed data{
  real<lower=-pi(),upper=pi()> angle[N_obs];
  real<lower=0> radius[N_obs];
  real<lower=0> dist_to_border[N_obs];
  vector[2] mu_angle;
  mu_angle[1] = 0;
  mu_angle[2] = pi();
  
  for(i in 1:N_obs){
    int current_order = order[i];
    vector[3] res; // [angle, radius, distance to border]
    // coordinates of 'previous fixations'
    real x_prev;
    real y_prev;
    
    if(current_order == 1){
    // if it's a first fixation, calculate angle relative to center
      x_prev = 400;
      y_prev = 300;
    } else{
    // otherwise calculate angle relative to the previous fixation  
      x_prev = x[i-1];
      y_prev = x[i-1];
    }
    
    res = calc_angle_border(x[i], y[i], x_prev, y_prev, 0, 800, 0, 600);
    angle[i]          = res[1];
    radius[i]         = res[2];
    dist_to_border[i] = res[3];
  }
}
parameters{
  real<lower=0> sigma_center; // width of the central bias
  real<lower=0> sigma_distance; // width of the distance bias
  real<lower=0> scale_obj; // scaling factor of the object size
  real<lower=0> kappa;
  simplex[5] weights; // objects, saliency, exploitation, central bias, horizontal bias
  vector[N_obj] z_weights_obj; // logits of weights of individual objects
  real mu_log_alpha; // population mean of the log of the decision boundary
  real sigma_log_alpha; // population sigma of the log of the decision boundary 
  vector[N_ppt] z_log_alpha; // individual offsets of boundary (on the log scale) for non-centered parametrization
  real mu_log_sigma_attention; // population mean of the log of the width of attention window
  real sigma_log_sigma_attention; // population sigma of the log ot the width of attention window 
  vector[N_ppt] z_log_sigma_attention; // individual offsets of width (on the log scale) for non-centered parametrization
}
transformed parameters{
  real log_sum_exp_log_lik_xy = 0; // this variable accumulates the log-likelihood of the fixation locations
  real wald_log_lik = 0; // this variable accumulates the log-likelihood of the fixation durations
  vector[5] log_weights = log(weights);
  vector[N_ppt] alpha = exp(mu_log_alpha + sigma_log_alpha * z_log_alpha); // individual decision boundaries
  vector[N_ppt] sigma_attention = exp(mu_log_sigma_attention + sigma_log_sigma_attention * z_log_sigma_attention); // individual widths of attention

  
  for(i in 1:N_obs){
    int current_order = order[i];
    int current_ppt   = id_ppt[i];
    int current_img   = id_img[i];
    int current_nei   = N_neighbors[i];
    int from          = obj_index_from[current_img];
    int to            = obj_index_to[current_img];
    vector[N_obj_in_img[current_img]] weights_obj = softmax(z_weights_obj[from:to]);
    vector[5] log_lik_xy = log_weights; // this vector stores the logs of (x, y) under each factor
    vector[2] att_filter = log_weights[1:2]; // only objects and saliency included in the calculation of drift
    real nu; // drift
    
    // object oriented behavior
    log_lik_xy[1] += mixture_normals(x[i], y[i], weights_obj, 
                                     obj_center_x[from:to], scale_obj * obj_width [from:to],
                                     obj_center_y[from:to], scale_obj * obj_height[from:to]);  
                                     
    att_filter[1] += log_integral_attention_mixture_2d(x[i], y[i], weights_obj, 
                                                       obj_center_x[from:to], scale_obj * obj_width [from:to],
                                                       obj_center_y[from:to], scale_obj * obj_height[from:to],
                                                       sigma_attention[current_ppt], sigma_attention[current_ppt]);
                                                       
    // saliency
    log_lik_xy[2] += log_lik_saliency[i];
    
    att_filter[2] += log_sum_exp(saliency_log[i][1:current_nei] - mean_sq_distances[i][1:current_nei] / square(sigma_attention[current_ppt]));
    
    //exploitation
    if(current_order == 1){
      log_lik_xy[3] += normal_lpdf(x[i] | 400,    sigma_distance);
      log_lik_xy[3] += normal_lpdf(y[i] | 300,    sigma_distance);
    } else{
      log_lik_xy[3] += normal_lpdf(x[i] | x[i-1], sigma_distance);
      log_lik_xy[3] += normal_lpdf(y[i] | y[i-1], sigma_distance);
    }
    
    // central bias
    log_lik_xy[4] += normal_lpdf(x[i] | 400, sigma_center);
    log_lik_xy[4] += normal_lpdf(y[i] | 300, sigma_center);


    // horizontal bias
    log_lik_xy[5] += direction_bias_lpdf(angle[i] | radius[i], dist_to_border[i], rep_vector(1.0/2.0, 2), mu_angle, rep_vector(kappa, 2));
    
    log_sum_exp_log_lik_xy += log_sum_exp(log_lik_xy);
    nu = log_sum_exp(log_weights[1:2]) - log_sum_exp(att_filter);
    
    wald_log_lik += wald_lpdf(duration[i] | alpha[current_ppt], nu);
  }
  
}
model{
  target += log_sum_exp_log_lik_xy;
  target += wald_log_lik;
  
  sigma_center   ~ gamma(2, 0.02);
  sigma_distance ~ gamma(2, 0.02);
  scale_obj      ~ normal(1, 0.5);
  kappa          ~ gamma(10, 1);
  weights        ~ dirichlet(rep_vector(2, 5));
  z_weights_obj  ~ std_normal();
  
  mu_log_alpha              ~ normal(0, 0.5);
  sigma_log_alpha           ~ gamma(2, 2);
  z_log_alpha               ~ std_normal();
  
  mu_log_sigma_attention    ~ normal(2, 1);
  sigma_log_sigma_attention ~ gamma(2, 2);
  z_log_sigma_attention     ~ std_normal();
}

