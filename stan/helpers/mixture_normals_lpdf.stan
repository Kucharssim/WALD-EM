  real mixture_normals(real x, real y, vector weights, vector mu_x, vector sigma_x, vector mu_y, vector sigma_y){
    int K = num_elements(weights);
    vector[K] log_lik = log(weights);
    
    for(k in 1:K){
      log_lik[k] += normal_lpdf(x | mu_x[k], sigma_x[k]);
      log_lik[k] += normal_lpdf(y | mu_y[k], sigma_y[k]);
    }
    
    return log_sum_exp(log_lik);
  }
  
  real mixture_trunc_normals(real x, real y, vector weights, vector mu_x, vector sigma_x, vector mu_y, vector sigma_y, real lb_x, real ub_x, real lb_y, real ub_y){
    int K = num_elements(weights);
    vector[K] log_lik = log(weights);
    
    for(k in 1:K){
      log_lik[k] += trunc_normal_lpdf(x | mu_x[k], sigma_x[k], lb_x, ub_x);
      log_lik[k] += trunc_normal_lpdf(y | mu_y[k], sigma_y[k], lb_y, ub_y);
    }
    
    return log_sum_exp(log_lik);
  }

  vector mixture_normals_rng(vector weights, vector mu_x, vector sigma_x, vector mu_y, vector sigma_y){
    int K = num_elements(weights);
    vector[2] out;
    
    int k = categorical_rng(weights);
    out[1] = normal_rng(mu_x[k], sigma_x[k]);
    out[2] = normal_rng(mu_y[k], sigma_y[k]);
    
    return out;
  } 
  
  vector mixture_trunc_normals_rng(vector weights, vector mu_x, vector sigma_x, vector mu_y, vector sigma_y, real lb_x, real ub_x, real lb_y, real ub_y) {
    int K = num_elements(weights);
    vector[2] out;
    
    int k = categorical_rng(weights);
    out[1] = trunc_normal_rng(mu_x[k], sigma_x[k], lb_x, ub_x);
    out[2] = trunc_normal_rng(mu_y[k], sigma_y[k], lb_y, ub_y);
    
    return out;
  }
