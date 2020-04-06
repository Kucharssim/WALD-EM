  real wald_lpdf(real y, real alpha, real nu) {
    real lpdf;
    
    lpdf = (
      log(alpha) - 1.0/2.0 * log(2*pi()) - 3.0/2.0*log(y) - (alpha - nu*y)^2/(2*y)
    );
    
    return lpdf;
  }

