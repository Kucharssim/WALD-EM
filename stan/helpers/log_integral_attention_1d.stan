  real log_integral_attention_1d(real center_attention, real center_lambda, real width_attention, real width_lambda){
    real var_a  = square(width_attention);
    real var_l  = square(width_lambda);
    real var_al = var_a + var_l;
    real diff_m_sq = square(center_attention - center_lambda);
    
    return log(width_attention) - 0.5 * log(var_al) - 0.5 * diff_m_sq / var_al;
  }

