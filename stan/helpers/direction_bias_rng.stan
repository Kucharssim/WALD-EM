// returns an angle from a mixture of von mises distributions
  real mixture_von_mises_rng(vector weights, vector mu, vector kappa){
    int k = categorical_rng(weights); 
    real theta = von_mises_rng(mu[k], kappa[k]);
    
    if(theta >  pi()) theta -= 2*pi();
    if(theta < -pi()) theta += 2*pi();
    
    return theta;
  }
  
// calculates distance to border given theta and x and y positions  
  real distance_to_border(real theta, real x, real y, real x_min, real x_max, real y_min, real y_max){
    real dist_x;
    real dist_y;
    real dist;
    // is the saccade in left direction? Take the right screen boundary, otherwise left
    real x_border = cos(theta) > 0 ? x_max : x_min;
    // id the saccade in upper direction? Take the upper screen boundary, otherwise bottom
    real y_border = sin(theta) > 0 ? y_max : y_min;
    
    // distance to the horizontal and vertical borders, respectively
    dist_x = (x_border - x) / cos(theta);
    dist_y = (y_border - y) / sin(theta);
    
    // the smaller distance applies
    dist = dist_x < dist_y ? dist_x : dist_y;
    
    return dist;
  }  
  
// calculates the displacement of the saccade
  vector xy_move(real theta, real length){
    real x_move = length * cos(theta);
    real y_move = length * sin(theta);
    vector[2] out;
    
    out[1] = x_move;
    out[2] = y_move;
    
    return out;
  }
  
// generates new fixation given old fixation
  vector direction_bias_rng(vector weights, vector mu, vector kappa, real x, real y, real x_min, real x_max, real y_min, real y_max){
    real theta  = mixture_von_mises_rng(weights, mu, kappa);
    real dist   = distance_to_border(theta, x, y, x_min, x_max, y_min, y_max);
    real length = uniform_rng(0, dist);
    
    vector[2] out = xy_move(theta, length);
    
    out[1] += x;
    out[2] += y;
    
    return out;
  }
