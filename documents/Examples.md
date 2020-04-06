Examples
================
Simon Kucharsky
2020-04-07

This file serves as a supplement to the article “Dynamic models of eye
movements” and provides step by step explanation of the building and
implementation of the modelling framework of the joint distribution of
fixation locations and fixation durations.

In this folder, there are additional files that go through specific
examples of the model. These examples should provide some guidance in
building and implementing such a model. The models are implemented in
Stan (Carpenter et al., 2017) in combination with the R package rstan
(Team & others, 2016).

# Principled Bayesian modeling

Here, we summarise several criteria and approaches we used for
validating our model, following general reccomendations for Bayesian
modeling, see Schad et al. (2019).

## Prior predictives

Each of the models is parametrized with several parameters, some of
which do, and other do not have clear physical interpretation.

## Computational faithfullness

We check whether the Bayesian inference for the implemented model is
self consistent, i.e., whether the Bayesian updating of the prior leads
to correct posteriors. We check that by conducting simulation-based
calibration (Talts et al., 2018).

## Model sensitivity

We also investigate whether and to what extend are the estimates of the
paramaters biased, and to what extend are data informative to specific
parameters. To do this, we plot the posterior z-scores and posterior
contraction against each other for all focal parameters in the model.

## Posterior predictives

# Common features of the Stan files

All of the models we implemented share some common common computations,
some of which needed to be implemented by specifying custom functions in
Stan. These functions are loaded in each of the model’s file by
including the `load_functions.stan` file in the function block. This
file loads the following files:

    ## #include helpers/wald_lpdf.stan
    ## #include helpers/wald_rng.stan
    ## #include helpers/mixture_normals_lpdf.stan
    ## #include helpers/log_integral_attention_1d.stan
    ## #include helpers/log_integral_attention_mixture_2d.stan

Here, we go through all of the files to document the functionality we
added to Stan.

## Wald distribution

We implemented the Wald distribution in Stan. The first function
`wald_lpdf` provides the log of the density function for the Wald
distribution. The density of the Wald distribution is the same as the
density of the Inverse Gaussian distribution after reparametrization to
parameters ![\\alpha](https://latex.codecogs.com/png.latex?%5Calpha
"\\alpha") (boundary) and
![\\nu](https://latex.codecogs.com/png.latex?%5Cnu "\\nu") (drift rate):

    ##   real wald_lpdf(real y, real alpha, real nu) {
    ##     real lpdf;
    ##     
    ##     lpdf = (
    ##       log(alpha) - 1.0/2.0 * log(2*pi()) - 3.0/2.0*log(y) - (alpha - nu*y)^2/(2*y)
    ##     );
    ##     
    ##     return lpdf;
    ##   }

The second function `wald_rng` implements the RNG function of a Wald
distributed variable:

    ##   real wald_rng(real alpha, real nu) {
    ##     real mu = alpha/nu;
    ##     real lambda = alpha^2;
    ##     real zeta = normal_rng(0, 1);
    ##     real zeta_sq = zeta^2;
    ##     real x = mu + (mu^2*zeta_sq)/(2*lambda) - mu/(2*lambda)*sqrt(4*mu*lambda*zeta_sq + mu^2*zeta_sq^2);
    ##     real z = uniform_rng(0, 1);
    ##     real y;
    ##     
    ##     if(z <= mu / (mu + x)){
    ##       y = x;
    ##     } else{
    ##       y = mu^2/x;
    ##     }
    ##     
    ##     return y;
    ##   }

The description of both functions can be found under Wikipedia entry for
the [Inverse Gaussian
distribution](https://en.wikipedia.org/wiki/Inverse_Gaussian_distribution).

## Log density of a gaussian mixture

For modelling the distribution of fixation locations, we defined the
intensity function as a mixture density:

  
![
\\lambda(x, y) = \\sum\_{k=1}^{K} \\pi\_k \\times \\text{Normal}(x |
\\mu\_{kx}, \\sigma\_{kx}) \\times \\text{Normal}(y | \\mu\_{ky},
\\sigma\_{ky}),
](https://latex.codecogs.com/png.latex?%0A%5Clambda%28x%2C%20y%29%20%3D%20%5Csum_%7Bk%3D1%7D%5E%7BK%7D%20%5Cpi_k%20%5Ctimes%20%5Ctext%7BNormal%7D%28x%20%7C%20%5Cmu_%7Bkx%7D%2C%20%5Csigma_%7Bkx%7D%29%20%5Ctimes%20%5Ctext%7BNormal%7D%28y%20%7C%20%5Cmu_%7Bky%7D%2C%20%5Csigma_%7Bky%7D%29%2C%0A
"
\\lambda(x, y) = \\sum_{k=1}^{K} \\pi_k \\times \\text{Normal}(x | \\mu_{kx}, \\sigma_{kx}) \\times \\text{Normal}(y | \\mu_{ky}, \\sigma_{ky}),
")  

with ![\\mu\_{kx}](https://latex.codecogs.com/png.latex?%5Cmu_%7Bkx%7D
"\\mu_{kx}") and
![\\sigma\_{kx}](https://latex.codecogs.com/png.latex?%5Csigma_%7Bkx%7D
"\\sigma_{kx}") being the center and width of a factor
![k](https://latex.codecogs.com/png.latex?k "k") in the
![x](https://latex.codecogs.com/png.latex?x "x") dimension,
![\\mu\_{yx}](https://latex.codecogs.com/png.latex?%5Cmu_%7Byx%7D
"\\mu_{yx}") and
![\\sigma\_{yx}](https://latex.codecogs.com/png.latex?%5Csigma_%7Byx%7D
"\\sigma_{yx}") being the center and width of a factor
![k](https://latex.codecogs.com/png.latex?k "k") in the
![y](https://latex.codecogs.com/png.latex?y "y") dimension, and
![\\pi\_k](https://latex.codecogs.com/png.latex?%5Cpi_k "\\pi_k") the
weight of a factor ![k](https://latex.codecogs.com/png.latex?k "k").

To obtain the log of this density, we constructed the following
convenience function that implements this
    calculation

    ##   real mixture_normals(real x, real y, vector weights, vector mu_x, vector sigma_x, vector mu_y, vector sigma_y){
    ##     int K = num_elements(weights);
    ##     vector[K] log_lik = log(weights);
    ##     
    ##     for(k in 1:K){
    ##       log_lik[k] += normal_lpdf(x | mu_x[k], sigma_x[k])
    ##       log_lik[k] += normal_lpdf(y | mu_y[k], sigma_y[k]);
    ##     }
    ##     
    ##     return log_sum_exp(log_lik);
    ##   }

## Common functions to compute the drift rate ![\\nu](https://latex.codecogs.com/png.latex?%5Cnu "\\nu")

The model we use relies on integration of the intensity function with
the attention window. Recall that we define the intensity function as
and that we defined the attention window as a kernel of a Gaussian
distribution:

  
![
a(x, y | s) = \\exp \\left(-\\frac{(x - s\_x)^2}{2\\sigma\_a^2}\\right)
\\times \\exp \\left(-\\frac{(y - s\_y)^2}{2\\sigma\_a^2}\\right)
](https://latex.codecogs.com/png.latex?%0Aa%28x%2C%20y%20%7C%20s%29%20%3D%20%5Cexp%20%5Cleft%28-%5Cfrac%7B%28x%20-%20s_x%29%5E2%7D%7B2%5Csigma_a%5E2%7D%5Cright%29%20%5Ctimes%20%5Cexp%20%5Cleft%28-%5Cfrac%7B%28y%20-%20s_y%29%5E2%7D%7B2%5Csigma_a%5E2%7D%5Cright%29%0A
"
a(x, y | s) = \\exp \\left(-\\frac{(x - s_x)^2}{2\\sigma_a^2}\\right) \\times \\exp \\left(-\\frac{(y - s_y)^2}{2\\sigma_a^2}\\right)
")  

The following function in `log_integral_attention_1d.stan` implements
the integral of an intensity function in one dimension multiplied by the
attention window:

  
![
\\log \\int \\frac{1}{\\sqrt{2\\pi}\\sigma}\\exp\\left(-\\frac{(x -
\\mu\_x)^2}{2\\sigma\_x^2}\\right) \\exp\\left(-\\frac{(x -
s\_x)^2}{2\\sigma\_a^2}\\right) dx
](https://latex.codecogs.com/png.latex?%0A%5Clog%20%5Cint%20%5Cfrac%7B1%7D%7B%5Csqrt%7B2%5Cpi%7D%5Csigma%7D%5Cexp%5Cleft%28-%5Cfrac%7B%28x%20-%20%5Cmu_x%29%5E2%7D%7B2%5Csigma_x%5E2%7D%5Cright%29%20%5Cexp%5Cleft%28-%5Cfrac%7B%28x%20-%20s_x%29%5E2%7D%7B2%5Csigma_a%5E2%7D%5Cright%29%20dx%0A
"
\\log \\int \\frac{1}{\\sqrt{2\\pi}\\sigma}\\exp\\left(-\\frac{(x - \\mu_x)^2}{2\\sigma_x^2}\\right) \\exp\\left(-\\frac{(x - s_x)^2}{2\\sigma_a^2}\\right) dx
")  

    ##   real log_integral_attention_1d(real center_attention, real center_lambda, real width_attention, real width_lambda){
    ##     real var_a  = square(width_attention);
    ##     real var_l  = square(width_lambda);
    ##     real var_al = var_a + var_l;
    ##     real diff_m_sq = square(center_attention - center_lambda);
    ##     
    ##     return log(width_attention) - 0.5 * log(var_al) - 0.5 * diff_m_sq / var_al;
    ##   }

And the function in `log_integral_attention_mixture_2d.stan` implements
the integration:

  
![
\\int\\int \\lambda(x, y) \\times a(x, y | s) dx dy
](https://latex.codecogs.com/png.latex?%0A%5Cint%5Cint%20%5Clambda%28x%2C%20y%29%20%5Ctimes%20a%28x%2C%20y%20%7C%20s%29%20dx%20dy%0A
"
\\int\\int \\lambda(x, y) \\times a(x, y | s) dx dy
")  

    ##   real log_integral_attention_mixture_2d(real x, real y, vector weights, vector mu_x, vector sigma_x, vector mu_y, vector sigma_y, real width_x, real width_y){
    ##     int K = num_elements(weights);
    ##     vector[K] log_int = log(weights);
    ##     
    ##     for(k in 1:K){
    ##       log_int += log_integral_attention_1d(x, mu_x[k], width_x, sigma_x[k]);
    ##       log_int += log_integral_attention_1d(y, mu_y[k], width_y, sigma_y[k]);
    ##     }
    ##     
    ##     return log_sum_exp(log_int); 
    ##   }

## References

<div id="refs" class="references">

<div id="ref-carpenter2017stan">

Carpenter, B., Gelman, A., Hoffman, M. D., Lee, D., Goodrich, B.,
Betancourt, M., Brubaker, M., Guo, J., Li, P., & Riddell, A. (2017).
Stan: A probabilistic programming language. *Journal of Statistical
Software*, *76*(1).

</div>

<div id="ref-schad2019principled">

Schad, D. J., Betancourt, M., & Vasishth, S. (2019). Toward a principled
bayesian workflow in cognitive science. *arXiv Preprint
arXiv:1904.12765*.

</div>

<div id="ref-talts2018sbc">

Talts, S., Betancourt, M., Simpson, D., Vehtari, A., & Gelman, A.
(2018). Validating bayesian inference algorithms with simulation-based
calibration. *arXiv Preprint arXiv:1804.06788*.

</div>

<div id="ref-team2016rstan">

Team, S. D., & others. (2016). RStan: The r interface to stan. *R
Package Version*, *2*(1).

</div>

</div>
