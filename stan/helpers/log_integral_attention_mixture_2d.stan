  real log_integral_attention_mixture_2d(real x, real y, vector weights, vector mu_x, vector sigma_x, vector mu_y, vector sigma_y, real width_x, real width_y){
    int K = num_elements(weights);
    vector[K] log_int = log(weights);
    
    for(k in 1:K){
      log_int[k] += log_integral_attention_1d(x, mu_x[k], width_x, sigma_x[k]);
      log_int[k] += log_integral_attention_1d(y, mu_y[k], width_y, sigma_y[k]);
    }
    
    return log_sum_exp(log_int); 
  }

