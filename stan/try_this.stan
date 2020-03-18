functions{
#include helpers/load_functions.stan
}
data{
  int T;
  real y[T];
}
parameters{
  real mu;
}
model{
  y ~ normal(mu, 1);
}