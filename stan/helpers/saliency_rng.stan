  vector saliency_rng(vector saliency, vector center_x, vector center_y, real half_width){
    int which_pixel = categorical_rng(saliency);
    vector[2] out;
    
    out[1] = uniform_rng(center_x[which_pixel] - half_width, center_x[which_pixel] + half_width);
    out[2] = uniform_rng(center_y[which_pixel] - half_width, center_y[which_pixel] + half_width);
    
    return out;
  }
  
  