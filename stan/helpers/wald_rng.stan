  real wald_rng(real alpha, real nu) {
    real mu = alpha/nu;
    real lambda = alpha^2;
    real zeta = normal_rng(0, 1);
    real zeta_sq = zeta^2;
    real x = mu + (mu^2*zeta_sq)/(2*lambda) - mu/(2*lambda)*sqrt(4*mu*lambda*zeta_sq + mu^2*zeta_sq^2);
    real z = uniform_rng(0, 1);
    real y;
    
    if(z <= mu / (mu + x)){
      y = x;
    } else{
      y = mu^2/x;
    }
    
    return y;
  }

