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
  
  int N_pix;
  real half_width_pixel;
  vector[N_pix] saliency_center_x;
  vector[N_pix] saliency_center_y;
  simplex[N_pix] saliency[N_img];
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
  vector[4] log_weights;
  vector[N_ppt] alpha;
  vector[N_ppt] sigma_attention;
}
transformed parameters{
  
}
model{
}
generated quantities{
  real<lower=0> duration_rep[N_obs];
  real<lower=0,upper=800> x_rep[N_obs];
  real<lower=0,upper=600> y_rep[N_obs];
  
  for(i in 1:N_obs){
    int current_order = order[i];
    int current_ppt   = id_ppt[i];
    int current_img   = id_img[i];
    int current_nei   = N_neighbors[i];
    int from          = obj_index_from[current_img];
    int to            = obj_index_to[current_img];
    vector[N_obj_in_img[current_img]] weights_obj = softmax(z_weights_obj[from:to]);
    vector[2] att_filter = log_weights[1:2]; // only objects and saliency
    real nu;
    int factor = categorical_rng(weights);
    int object = categorical_rng(weights_obj);
    vector[2] xy_rep = rep_vector(0, 2);
    
    // prediction fixation duration
    att_filter[1] += log_integral_attention_mixture_2d(x[i], y[i], weights_obj, 
                                                       obj_center_x[from:to], scale_obj * obj_width [from:to],
                                                       obj_center_y[from:to], scale_obj * obj_height[from:to],
                                                       sigma_attention[current_ppt], sigma_attention[current_ppt]);
    
    att_filter[2] += log_sum_exp(saliency_log[i][1:current_nei] - mean_sq_distances[i][1:current_nei] / square(sigma_attention[current_ppt]));
    
    nu = log_sum_exp(log_weights[1:2]) - log_sum_exp(att_filter);
    
    duration_rep[i] = wald_rng(alpha[current_ppt], nu);
    
    // prediction location;
    if(factor == 1){ // object
      x_rep[i] = trunc_normal_rng(obj_center_x[from:to][object], scale_obj * obj_width [from:to][object], 0, 800);
      y_rep[i] = trunc_normal_rng(obj_center_y[from:to][object], scale_obj * obj_height[from:to][object], 0, 600);
    } else if(factor == 2){ // saliency
      xy_rep   = saliency_rng(saliency[current_img], saliency_center_x, saliency_center_y, half_width_pixel);
      x_rep[i] = xy_rep[1] - 0.5;
      y_rep[i] = xy_rep[2] - 0.5;
    } else if(factor == 3){ // exploitation
      if(current_order == 1){
        x_rep[i] = trunc_normal_rng(400, sigma_distance, 0, 800);
        y_rep[i] = trunc_normal_rng(300, sigma_distance, 0, 600);
      } else{
        x_rep[i] = trunc_normal_rng(x[i-1], sigma_distance, 0, 800);
        y_rep[i] = trunc_normal_rng(y[i-1], sigma_distance, 0, 600);
      }
    } else{ // central bias
      x_rep[i] = trunc_normal_rng(400, sigma_center, 0, 800);
      y_rep[i] = trunc_normal_rng(300, sigma_center, 0, 600);
    }
  }
  
}

