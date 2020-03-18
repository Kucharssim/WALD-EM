  real mixture_normals_lpdf(real y, vector weights, vector mu, vector sigma){
    int K = num_elements(mu);
    vector[K] log_lik = log(weights);
    
    for(k in 1:K){
      log_lik[k] += normal_lpdf(y | mu[k], sigma[k]);
    }
    
    return log_sum_exp(log_lik);
  }