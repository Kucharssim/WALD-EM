Examples
================
Simon Kucharsky
2020-07-01

This file serves as a supplement to the article “Dynamic models of eye
movements” and provides step by step explanation of the building and
implementation of the modelling framework of the joint distribution of
fixation locations and fixation durations.

In this folder, there are additional files that go through specific
examples of the model. These examples should provide some guidance in
building and implementing such a model. The models are implemented in
Stan (Carpenter et al., 2017) in combination with the R package rstan
(<span class="citeproc-not-found" data-reference-id="team2016rstan">**???**</span>).

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

## References

<div id="refs" class="references hanging-indent">

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

</div>
