Helper functions
================
Šimon Kucharský
2020-07-02

  - [Common features of the Stan
    files](#common-features-of-the-stan-files)
      - [Wald distribution](#wald-distribution)
      - [Truncated normal distribution](#truncated-normal-distribution)
      - [Log density of a gaussian
        mixture](#log-density-of-a-gaussian-mixture)
      - [Common functions to compute the drift rate
        ![\\nu](https://latex.codecogs.com/png.latex?%5Cnu
        "\\nu")](#common-functions-to-compute-the-drift-rate-nu)
      - [Saliency](#saliency)
      - [Direction bias](#direction-bias)
      - [References](#references)

# Common features of the Stan files

All of the models we implemented share some common common computations,
some of which needed to be implemented by specifying custom functions in
Stan (Carpenter et al., 2017). These functions are loaded in each of the
model’s file by including the `load_functions.stan` file in the function
block. Specifically, placing `#include stan/helpers/load_functions.stan`
inside of `functions{}` block in a `.stan` file makes accessible
functions defined in the following files:

``` r
writeLines(readLines(here::here("stan", "helpers", "load_functions.stan")))
```

    ## #include stan/helpers/wald_lpdf.stan
    ## #include stan/helpers/wald_rng.stan
    ## #include stan/helpers/trunc_normal_lpdf.stan
    ## #include stan/helpers/trunc_normal_rng.stan
    ## #include stan/helpers/mixture_normals_lpdf.stan
    ## #include stan/helpers/log_integral_attention_1d.stan
    ## #include stan/helpers/log_integral_attention_mixture_2d.stan
    ## #include stan/helpers/saliency_rng.stan
    ## #include stan/helpers/calc_angle_border.stan
    ## #include stan/helpers/direction_bias_lpdf.stan
    ## #include stan/helpers/direction_bias_rng.stan

Here, we go through all of the files to document the functionality we
added to Stan for fitting the models.

NB. These functions can be exposed to R by the script
`source(here::here("R", "expose_helpers_stan.R"))`.

## Wald distribution

We implemented the Wald distribution in Stan. Specifically, we
implemented the Wald log likelihood function, and the function that
generates data that are distributed under the Wald distribution. The
first function `wald_lpdf()` provides the log of the density function
for the Wald distribution. The density of the Wald distribution is the
same as the density of the Inverse Gaussian distribution after
reparametrization to parameters
![\\alpha](https://latex.codecogs.com/png.latex?%5Calpha "\\alpha")
(boundary) and ![\\nu](https://latex.codecogs.com/png.latex?%5Cnu
"\\nu") (drift rate):

``` r
writeLines(readLines(here::here("stan", "helpers", "wald_lpdf.stan")))
```

    ##   real wald_lpdf(real y, real alpha, real nu) {
    ##     real lpdf;
    ##     
    ##     lpdf = (
    ##       log(alpha) - 1.0/2.0 * log(2*pi()) - 3.0/2.0*log(y) - (alpha - nu*y)^2/(2*y)
    ##     );
    ##     
    ##     return lpdf;
    ##   }

The second function `wald_rng()` implements the RNG function of a Wald
distributed variable:

``` r
writeLines(readLines(here::here("stan", "helpers", "wald_rng.stan")))
```

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

## Truncated normal distribution

We implemented a truncated normal distributions in Stan. These are not,
however, used in the model itself. It could be argued that using
truncated normal distributions will make more sense for some of the
parts of the model, as the eye-tracking data is constrained by the
dimensions of the screen. However, we did not thoroughly test and use
the model with truncated normal distributions, and did not use it. In
general, the bias in parameters created by omitting the truncation
should be relatively small, as most of the mass of the gaussian
distributions in our model was concentrated far away from the truncation
bounds.

However, in some of the simulations, we used truncated normal
distribution to be able to generate data that are more realistic (i.e.,
do not go outside of the screen dimensions.)

The following function implements the density of the truncated normal
distribution (with `lb` the lower bound and `ub` the upper bound).

``` r
writeLines(readLines(here::here("stan", "helpers", "trunc_normal_lpdf.stan")))
```

    ##   real trunc_normal_lpdf(real y, real mu, real sigma, real lb, real ub){
    ##     real out;
    ##     
    ##     out  = normal_lpdf(y | mu, sigma);
    ##     out -= log_diff_exp(normal_lcdf(ub | mu, sigma),
    ##                         normal_lcdf(lb | mu, sigma));
    ##                         
    ##     return out;
    ##   }

And the function below implements generating data from truncated normal
distribution.

``` r
writeLines(readLines(here::here("stan", "helpers", "trunc_normal_rng.stan")))
```

    ##   real trunc_normal_rng(real mu, real sigma, real lb, real ub){
    ##     real p1 = normal_cdf(lb, mu, sigma);  // cdf with lower bound
    ##     real p2 = normal_cdf(ub, mu, sigma);  // cdf with upper bound
    ##     real u = uniform_rng(p1, p2);
    ##     return (sigma * inv_Phi(u)) + mu;  // inverse cdf
    ##   }

The functions should be self-explanatory, but for the sake of
completeness consult the related post on a Stan forum [related post on a
Stan
forum](https://discourse.mc-stan.org/t/rng-for-truncated-distributions/3122/3)
should there be any confusion.

## Log density of a gaussian mixture

For modelling the distribution of fixation locations, we used a mixture
of bivariate normal distributions for some of the “factors”:

  
![&#10;f(x, y) = \\sum\_{k=1}^{K} \\pi\_k \\times \\text{Normal}(x |
\\mu\_{kx}, \\sigma\_{kx}) \\times \\text{Normal}(y | \\mu\_{ky},
\\sigma\_{ky}),&#10;](https://latex.codecogs.com/png.latex?%0Af%28x%2C%20y%29%20%3D%20%5Csum_%7Bk%3D1%7D%5E%7BK%7D%20%5Cpi_k%20%5Ctimes%20%5Ctext%7BNormal%7D%28x%20%7C%20%5Cmu_%7Bkx%7D%2C%20%5Csigma_%7Bkx%7D%29%20%5Ctimes%20%5Ctext%7BNormal%7D%28y%20%7C%20%5Cmu_%7Bky%7D%2C%20%5Csigma_%7Bky%7D%29%2C%0A
"
f(x, y) = \\sum_{k=1}^{K} \\pi_k \\times \\text{Normal}(x | \\mu_{kx}, \\sigma_{kx}) \\times \\text{Normal}(y | \\mu_{ky}, \\sigma_{ky}),
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

To obtain the log of this density, we constructed the following helper
functions `mixture_normals()` and `mixture_trunc_normals()` that
implement this calculations. The latter function implements the
truncated mixture and requires to specify the lower and upper truncation
bounds on the x and y dimensions.

``` r
writeLines(readLines(here::here("stan", "helpers", "mixture_normals_lpdf.stan")))
```

    ##   real mixture_normals(real x, real y, vector weights, vector mu_x, vector sigma_x, vector mu_y, vector sigma_y){
    ##     int K = num_elements(weights);
    ##     vector[K] log_lik = log(weights);
    ##     
    ##     for(k in 1:K){
    ##       log_lik[k] += normal_lpdf(x | mu_x[k], sigma_x[k]);
    ##       log_lik[k] += normal_lpdf(y | mu_y[k], sigma_y[k]);
    ##     }
    ##     
    ##     return log_sum_exp(log_lik);
    ##   }
    ##   
    ##   real mixture_trunc_normals(real x, real y, vector weights, vector mu_x, vector sigma_x, vector mu_y, vector sigma_y, real lb_x, real ub_x, real lb_y, real ub_y){
    ##     int K = num_elements(weights);
    ##     vector[K] log_lik = log(weights);
    ##     
    ##     for(k in 1:K){
    ##       log_lik[k] += trunc_normal_lpdf(x | mu_x[k], sigma_x[k], lb_x, ub_x);
    ##       log_lik[k] += trunc_normal_lpdf(y | mu_y[k], sigma_y[k], lb_y, ub_y);
    ##     }
    ##     
    ##     return log_sum_exp(log_lik);
    ##   }

## Common functions to compute the drift rate ![\\nu](https://latex.codecogs.com/png.latex?%5Cnu "\\nu")

The model we use relies on integration of the intensity function with
the attention window. Recall that we define the intensity function as
and that we defined the attention window as a kernel of a Gaussian
distribution:

  
![&#10;a(x, y | s) = \\exp \\left(-\\frac{(x -
s\_x)^2}{2\\sigma\_a^2}\\right) \\times \\exp \\left(-\\frac{(y -
s\_y)^2}{2\\sigma\_a^2}\\right)&#10;](https://latex.codecogs.com/png.latex?%0Aa%28x%2C%20y%20%7C%20s%29%20%3D%20%5Cexp%20%5Cleft%28-%5Cfrac%7B%28x%20-%20s_x%29%5E2%7D%7B2%5Csigma_a%5E2%7D%5Cright%29%20%5Ctimes%20%5Cexp%20%5Cleft%28-%5Cfrac%7B%28y%20-%20s_y%29%5E2%7D%7B2%5Csigma_a%5E2%7D%5Cright%29%0A
"
a(x, y | s) = \\exp \\left(-\\frac{(x - s_x)^2}{2\\sigma_a^2}\\right) \\times \\exp \\left(-\\frac{(y - s_y)^2}{2\\sigma_a^2}\\right)
")  

The following function `log_integral_attention_1d()` implements the
integral of an intensity function in one dimension multiplied by the
attention window:

  
![&#10;\\log \\int
\\frac{1}{\\sqrt{2\\pi}\\sigma\_x}\\exp\\left(-\\frac{(x -
\\mu\_x)^2}{2\\sigma\_x^2}\\right) \\exp\\left(-\\frac{(x -
s\_x)^2}{2\\sigma\_a^2}\\right) dx =
\\frac{\\sigma\_a}{\\sqrt{\\sigma\_a^2+\\sigma\_x^2}}
\\exp\\left\[-\\frac{(\\mu\_x-s\_x)^2}{2(\\sigma\_a^2 +
\\sigma\_x^2)}\\right\],&#10;](https://latex.codecogs.com/png.latex?%0A%5Clog%20%5Cint%20%5Cfrac%7B1%7D%7B%5Csqrt%7B2%5Cpi%7D%5Csigma_x%7D%5Cexp%5Cleft%28-%5Cfrac%7B%28x%20-%20%5Cmu_x%29%5E2%7D%7B2%5Csigma_x%5E2%7D%5Cright%29%20%5Cexp%5Cleft%28-%5Cfrac%7B%28x%20-%20s_x%29%5E2%7D%7B2%5Csigma_a%5E2%7D%5Cright%29%20dx%20%3D%20%5Cfrac%7B%5Csigma_a%7D%7B%5Csqrt%7B%5Csigma_a%5E2%2B%5Csigma_x%5E2%7D%7D%20%5Cexp%5Cleft%5B-%5Cfrac%7B%28%5Cmu_x-s_x%29%5E2%7D%7B2%28%5Csigma_a%5E2%20%2B%20%5Csigma_x%5E2%29%7D%5Cright%5D%2C%0A
"
\\log \\int \\frac{1}{\\sqrt{2\\pi}\\sigma_x}\\exp\\left(-\\frac{(x - \\mu_x)^2}{2\\sigma_x^2}\\right) \\exp\\left(-\\frac{(x - s_x)^2}{2\\sigma_a^2}\\right) dx = \\frac{\\sigma_a}{\\sqrt{\\sigma_a^2+\\sigma_x^2}} \\exp\\left[-\\frac{(\\mu_x-s_x)^2}{2(\\sigma_a^2 + \\sigma_x^2)}\\right],
")  

which is the “information” available along the x dimension given the
fixation at location ![s\_x](https://latex.codecogs.com/png.latex?s_x
"s_x"), the attention window width
![\\sigma\_a](https://latex.codecogs.com/png.latex?%5Csigma_a
"\\sigma_a"), and the location and spread of some factor
![\\mu\_x](https://latex.codecogs.com/png.latex?%5Cmu_x "\\mu_x") and
![\\sigma\_x](https://latex.codecogs.com/png.latex?%5Csigma_x
"\\sigma_x").

``` r
writeLines(readLines(here::here("stan", "helpers", "log_integral_attention_1d.stan")))
```

    ##   real log_integral_attention_1d(real center_attention, real center_lambda, real width_attention, real width_lambda){
    ##     real var_a  = square(width_attention);
    ##     real var_l  = square(width_lambda);
    ##     real var_al = var_a + var_l;
    ##     real diff_m_sq = square(center_attention - center_lambda);
    ##     
    ##     return log(width_attention) - 0.5 * log(var_al) - 0.5 * diff_m_sq / var_al;
    ##   }

And the function in `log_integral_attention_mixture_2d()` implements the
integration:

  
![&#10;\\int\\int \\lambda(x, y) \\times a(x, y | s) dx
dy,&#10;](https://latex.codecogs.com/png.latex?%0A%5Cint%5Cint%20%5Clambda%28x%2C%20y%29%20%5Ctimes%20a%28x%2C%20y%20%7C%20s%29%20dx%20dy%2C%0A
"
\\int\\int \\lambda(x, y) \\times a(x, y | s) dx dy,
")  

which basically combines `log_integral_attention_1d()` for x and y
dimension separately for each “factor”, and weights them by their
`weights`.

``` r
writeLines(readLines(here::here("stan", "helpers", "log_integral_attention_mixture_2d.stan")))
```

    ##   real log_integral_attention_mixture_2d(real x, real y, vector weights, vector mu_x, vector sigma_x, vector mu_y, vector sigma_y, real width_x, real width_y){
    ##     int K = num_elements(weights);
    ##     vector[K] log_int = log(weights);
    ##     
    ##     for(k in 1:K){
    ##       log_int[k] += log_integral_attention_1d(x, mu_x[k], width_x, sigma_x[k]);
    ##       log_int[k] += log_integral_attention_1d(y, mu_y[k], width_y, sigma_y[k]);
    ##     }
    ##     
    ##     return log_sum_exp(log_int); 
    ##   }

## Saliency

``` r
writeLines(readLines(here::here("stan", "helpers", "saliency_rng.stan")))
```

    ##   vector saliency_rng(vector saliency, vector center_x, vector center_y, real half_width){
    ##     int which_pixel = categorical_rng(saliency);
    ##     vector[2] out;
    ##     
    ##     out[1] = uniform_rng(center_x[which_pixel] - half_width, center_x[which_pixel] + half_width);
    ##     out[2] = uniform_rng(center_y[which_pixel] - half_width, center_y[which_pixel] + half_width);
    ##     
    ##     return out;
    ##   }

## Direction bias

## References

<div id="refs" class="references hanging-indent">

<div id="ref-carpenter2017stan">

Carpenter, B., Gelman, A., Hoffman, M. D., Lee, D., Goodrich, B.,
Betancourt, M., Brubaker, M., Guo, J., Li, P., & Riddell, A. (2017).
Stan: A probabilistic programming language. *Journal of Statistical
Software*, *76*(1).

</div>

</div>
