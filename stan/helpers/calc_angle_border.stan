  vector calc_angle_border(real x, real y, real x_prev, real y_prev, real x_min, real x_max, real y_min, real y_max){
    // 'vector' of saccade
    real x_move;
    real y_move;
    
    // 'length' of saccade
    real length;
    
    // which border is the saccade coming to
    real x_border;
    real y_border;
    
    // distances to the borders
    vector[2] distances;
    
    vector[2] scale;
    vector[3] out; //[angle, radius, distance_to_border]
    
    x_move = x - x_prev;
    y_move = y - y_prev;
    
    x_border = x_move > 0 ? x_max : x_min; // border on x coordinate in the direction of the saccade
    y_border = y_move > 0 ? y_max : y_min; // border on y coordinate in the direction of the saccade
    
    x_border = x_border - x_prev; // distance of the outgoing fixation to the x border
    y_border = y_border - y_prev; // distance of the outgoing fixation  to the y border
    
    /*
    angle:
          0 = rigtht
         pi = left
     1/2 pi = down
    -1/2 pi = up
    */
    out[1] = atan2(y_move, x_move);
    
    /*
    radius: compute the saccade amplitude
    */
    length = sqrt(square(x_move) + square(y_move));
    out[2] = length;
    
    /*
    distance:
    if is the outgoing fixation closer to the border on the x-coordinate, 
    we calculate the distance of the fixation to the edge of the display
    by scaling a triangle that represents the saccade on the x projection
    otherwise we do the same but with the y projection
    */
    scale[1] = fabs(x_border/x_move) * length;
    scale[2] = fabs(y_border/y_move) * length;
    
    out[3] = scale[1] < scale[2] ? scale[1] : scale[2];
    //out[2] = x_border > y_border ? x_border/fabs(x_move) * length : y_border/fabs(y_move) * length;
    
    return out;
  }
