# Simulate a prior study

Simulates prior studies for a model family: data are generated from the
family likelihood and the parameters of interest are estimated,
generating one `bgm_parameters` object per study, as if obtained from a
previous study. Elicitation is a separate step. What is generated
depends on the class of `study`. See the method pages listed under
Details.

## Usage

``` r
simulate_prior_study(study, n_studies = 1L, ...)
```

## Arguments

- study:

  A `bgm_study` object, as returned by each family's constructor, for
  example
  [`ggm_study`](https://bayesian-graphical-modelling-lab.github.io/designbgm/reference/ggm_study.md).

- n_studies:

  Number of independent prior studies to simulate. Default 1.

- ...:

  Family-specific arguments passed to the method.

## Value

A list of `n_studies` `bgm_parameters` objects.

## Details

Methods are available for the following classes:

- [`ggm_study`](https://bayesian-graphical-modelling-lab.github.io/designbgm/reference/simulate_prior_study.ggm_study.md):

  Gaussian graphical models

## See also

Other prior elicitation:
[`elicit_prior()`](https://bayesian-graphical-modelling-lab.github.io/designbgm/reference/elicit_prior.md),
[`elicit_prior.ggm_parameters()`](https://bayesian-graphical-modelling-lab.github.io/designbgm/reference/elicit_prior.ggm_parameters.md),
[`ggm_study()`](https://bayesian-graphical-modelling-lab.github.io/designbgm/reference/ggm_study.md),
[`simulate_prior_study.ggm_study()`](https://bayesian-graphical-modelling-lab.github.io/designbgm/reference/simulate_prior_study.ggm_study.md)
