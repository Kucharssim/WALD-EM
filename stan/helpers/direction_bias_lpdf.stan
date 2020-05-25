  // returns an angle from a mixture of von mises distributions
  real mixture_von_mises_lpdf(real angle, vector weights, vector mu, vector kappa){
    int K = num_elements(weights);
    vector[K] log_lik = log(weights); 
    
    for(k in 1:K){
      log_lik[k] += von_mises_lpdf(angle | mu[k], kappa[k]);
    }
    
    return log_sum_exp(log_lik);
  }
  
  real direction_bias_lpdf(real angle, real r, real dist, vector weights, vector mu, vector kappa){
    real log_lik_angle = mixture_von_mises_lpdf(angle | weights, mu, kappa);

    /*
    // normalize by distance to border (nomrmalizing constant of uniform between 0 and distance to the screen log(dist))
    // add jacobian to account for the transformation dx dy = r dr dtheta
    */
    return log_lik_angle - log(dist) - log(r); 
  }
  