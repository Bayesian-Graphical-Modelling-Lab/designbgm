# Elicit the informative prior

Converts a `bgm_parameters` object into the elicited prior in the form
that the family's samplers and estimators expect. The elicitation
strategy depends on the class of `params`. See the method pages listed
under Details.

## Usage

``` r
elicit_prior(params, ...)
```

## Arguments

- params:

  A `bgm_parameters` object.

- ...:

  Family-specific arguments.

## Value

A `bgm_elicited` object describing the elicited prior. It also carries a
family-specific class (for example `ggm_elicited`), on which the
downstream methods dispatch.

## Details

Methods are available for the following classes:

- [`ggm_parameters`](https://bayesian-graphical-modelling-lab.github.io/designbgm/reference/elicit_prior.ggm_parameters.md):

  Gaussian graphical models

## See also

Other prior elicitation:
[`elicit_prior.ggm_parameters()`](https://bayesian-graphical-modelling-lab.github.io/designbgm/reference/elicit_prior.ggm_parameters.md),
[`ggm_study()`](https://bayesian-graphical-modelling-lab.github.io/designbgm/reference/ggm_study.md),
[`simulate_prior_study()`](https://bayesian-graphical-modelling-lab.github.io/designbgm/reference/simulate_prior_study.md),
[`simulate_prior_study.ggm_study()`](https://bayesian-graphical-modelling-lab.github.io/designbgm/reference/simulate_prior_study.ggm_study.md)
