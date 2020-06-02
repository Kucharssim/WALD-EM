  real trunc_normal_lpdf(real y, real mu, real sigma, real lb, real ub){
    real out;
    
    out  = normal_lpdf(y | mu, sigma);
    out -= log_diff_exp(normal_lcdf(ub | mu, sigma),
                        normal_lcdf(lb | mu, sigma));
                        
    return out;
  }

