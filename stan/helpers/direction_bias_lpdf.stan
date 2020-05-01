  real direction_bias_lpdf(real angle, real dist, vector weights, vector mu, vector kappa){
    int K = num_elements(weights);
    vector[K] log_lik = log(weights); 
    
    for(k in 1:K){
      log_lik[k] += von_mises_lpdf(angle | mu[k], kappa[k]);
    }
    
    return log_sum_exp(log_lik) - log(dist); // normalize by distance
  }
  
  