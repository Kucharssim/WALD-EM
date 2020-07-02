Examples
================
Simon Kucharsky
2020-07-02

  - [Principled Bayesian modeling](#principled-bayesian-modeling)
      - [Prior predictives](#prior-predictives)
      - [Computational faithfullness](#computational-faithfullness)
      - [Model sensitivity](#model-sensitivity)
      - [Posterior predictives](#posterior-predictives)
  - [Contents of examples](#contents-of-examples)
      - [<span>direction\_bias</span>](#direction_bias)
      - [<span>direction\_central\_bias</span>](#direction_central_bias)
  - [References](#references)

This file serves as a supplement to the article “Dynamic models of eye
movements” and provides step by step explanation of the building and
implementation of the modelling framework of the joint distribution of
fixation locations and fixation durations.

In this folder, there are additional files that go through specific
examples of the model. These examples should provide some guidance in
building and implementing such a model. The models are implemented in
Stan (Carpenter et al., 2017) in combination with the R package rstan
(Guo et al., 2020).

# Principled Bayesian modeling

Here, we summarise several criteria and approaches we used while
building and validating our model, following general reccomendations for
Bayesian modeling, see Schad et al. (2019).

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

Posterior predictives are a useful tool for assessing the adequacy of
the model. Specifically, we can look at whether the model is able to
predict specific data points, or whether it captures certain
regularities in the data.

# Contents of examples

The [documents](documents/) folder contains some example files that show
how to build specific versions of the model. They also demonstrate that
various parts of the model were correctly implemented. Here, we list the
files and provide short descriptions:

## [direction\_bias](documents/direction_bias.md/)

## [direction\_central\_bias](documents/direction_central_bias.md)

# References

<div id="refs" class="references hanging-indent">

<div id="ref-carpenter2017stan">

Carpenter, B., Gelman, A., Hoffman, M. D., Lee, D., Goodrich, B.,
Betancourt, M., Brubaker, M., Guo, J., Li, P., & Riddell, A. (2017).
Stan: A probabilistic programming language. *Journal of Statistical
Software*, *76*(1).

</div>

<div id="ref-R-rstan">

Guo, J., Gabry, J., & Goodrich, B. (2020). *Rstan: R interface to stan*.
<https://CRAN.R-project.org/package=rstan>

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

</div>
