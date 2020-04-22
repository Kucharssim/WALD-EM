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
  
  /*
  real lb_x;
  real ub_x;
  real lb_y;
  real ub_y;
  */
}
parameters{
  real<lower=0> sigma_center; // width of the central bias
  real<lower=0> sigma_distance; // width of the distance bias
  real<lower=0> scale_obj; // scaling factor of the object size
  simplex[4] weights; // objects, saliency, central bias, distance
  vector[N_obj] z_weights_obj; // logits of weights of individual objects
  real mu_log_alpha; // population mean of the log of the decision boundary
  real sigma_log_alpha; // population sigma of the log of the decision boundary 
  vector[N_ppt] z_log_alpha; // individual offsets of boundary (on the log scale) for non-centered parametrization
  real mu_log_sigma_attention; // population mean of the log of the width of attention window
  real sigma_log_sigma_attention; // population sigma of the log ot the width of attention window 
  vector[N_ppt] z_log_sigma_attention; // individual offsets of width (on the log scale) for non-centered parametrization
}
transformed parameters{
  vector[4] log_weights = log(weights);
  //vector[N_obs] log_sum_exp_log_lik_xy;
  //vector[N_obs] wald_log_lik;
  real log_sum_exp_log_lik_xy = 0;
  real wald_log_lik = 0;
  //vector[N_obs] nu;
  vector[N_ppt] alpha = exp(mu_log_alpha + sigma_log_alpha * z_log_alpha);
  vector[N_ppt] sigma_attention = exp(mu_log_sigma_attention + sigma_log_sigma_attention * z_log_sigma_attention);

  
  for(i in 1:N_obs){
    int current_order = order[i];
    int current_ppt   = id_ppt[i];
    int current_img   = id_img[i];
    int current_nei   = N_neighbors[i];
    int from          = obj_index_from[current_img];
    int to            = obj_index_to[current_img];
    vector[N_obj_in_img[current_img]] weights_obj = softmax(z_weights_obj[from:to]);
    vector[4] log_lik_xy = log_weights;
    vector[2] att_filter = log_weights[1:2]; // only objects and saliency
    real nu;
    
    // object oriented behavior
    /*
    log_lik_xy[1] += mixture_trunc_normals(x[i], y[i], weights_obj, 
                                           obj_center_x[from:to], scale_obj * obj_width [from:to],
                                           obj_center_y[from:to], scale_obj * obj_height[from:to],
                                           lb_x, ub_x, lb_y, ub_y);
    */
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
    
    //distance bias
    if(current_order != 1){
      log_lik_xy[3] += normal_lpdf(x[i] | x[i-1], sigma_distance);
      log_lik_xy[3] += normal_lpdf(y[i] | y[i-1], sigma_distance);
    } else{
      log_lik_xy[3] += negative_infinity();
    }
    
    // central bias
    log_lik_xy[4] += normal_lpdf(x[i] | 400, sigma_center);
    log_lik_xy[4] += normal_lpdf(y[i] | 300, sigma_center);

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
  weights        ~ dirichlet(rep_vector(2, 4));
  z_weights_obj  ~ std_normal();
  
  mu_log_alpha              ~ normal(0, 0.5);
  sigma_log_alpha           ~ gamma(2, 2);
  z_log_alpha               ~ std_normal();
  
  mu_log_sigma_attention    ~ normal(2, 1);
  sigma_log_sigma_attention ~ gamma(2, 2);
  z_log_sigma_attention     ~ std_normal();
}
/*
generated quantities{
  vector[N_obs] mean_duration_pred;
  vector[N_obs] duration_pred;
  vector[N_obs] log_lik;
  
  for(i in 1:N_obs){
    int current_ppt = id_ppt[i];
    mean_duration_pred[i] = alpha[current_ppt] / nu[i];
    duration_pred[i] = wald_rng(alpha[current_ppt], nu[i], 0);
    log_lik[i] = wald_log_lik[i] + log_sum_exp_log_lik_xy[i];
  }
}
*/

