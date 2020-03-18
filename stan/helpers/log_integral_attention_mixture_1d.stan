  real log_integral_attention_mixture_1d(real center_attention, vector centers_lambda, real width_attention, vector widths_lambda, vector weights){
    int K = num_elements(centers_lambda);
    vector[K] log_int = log(weights);
    
    for(k in 1:K){
      log_int += log_integral_attention_1d(center_attention, centers_lambda[k], width_attention, widths_lambda[k]);
    }
   
    return log_sum_exp(log_int); 
  }