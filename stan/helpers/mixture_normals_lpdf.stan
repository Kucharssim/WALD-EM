  real mixture_normals(real x, real y, vector weights, vector mu_x, vector sigma_x, vector mu_y, vector sigma_y){
    int K = num_elements(weights);
    vector[K] log_lik = log(weights);
    
    for(k in 1:K){
      log_lik[k] += normal_lpdf(x | mu_x[k], sigma_x[k]);
      log_lik[k] += normal_lpdf(y | mu_y[k], sigma_y[k]);
    }
    
    return log_sum_exp(log_lik);
  }
  
  real mixture_trunc_normals(real x, real y, vector weights, vector mu_x, vector sigma_x, vector mu_y, vector sigma_y, real max_x, real max_y){
    int K = num_elements(weights);
    vector[K] log_lik = log(weights);
    
    for(k in 1:K){
      log_lik[k] += normal_lpdf(x | mu_x[k], sigma_x[k]);
      log_lik[k] -= log_diff_exp(normal_lcdf(max_x | mu_x[k], sigma_x[k]),
                                 normal_lcdf(0     | mu_x[k], sigma_x[k])
                                 ); 
      log_lik[k] += normal_lpdf(y | mu_y[k], sigma_y[k]);
      log_lik[k] -= log_diff_exp(normal_lcdf(max_y | mu_y[k], sigma_y[k]),
                                 normal_lcdf(0     | mu_y[k], sigma_y[k])
                                 ); 
    }
    
    return log_sum_exp(log_lik);
  }

