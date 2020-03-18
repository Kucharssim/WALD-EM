  real wald_rng(real alpha, real nu) {
    real mu = alpha/nu;
    real lambda = alpha^2;
    real nu = normal_rng(0, 1);
    real nu_sq = nu^2;
    real x = mu + (mu^2*nu_sq)/(2*lambda) - mu/(2*lambda)*sqrt(4*mu*lambda*nu_sq + mu^2*nu_sq^2);
    real z = uniform_rng(0, 1);
    real y;
    
    if(z <= mu / (mu + x)){
      y = x;
    } else{
      y = mu^2/x;
    }
    
    return y;
  }