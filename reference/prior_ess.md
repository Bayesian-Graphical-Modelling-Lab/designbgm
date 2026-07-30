# Prior effective sample size

Quantifies the prior effective sample size (ESS) of an elicited prior
described by a `bgm_parameters` object. The estimators available depend
on the class of `params`. See the method pages listed under Details.

## Usage

``` r
prior_ess(params, estimator = NULL, ...)
```

## Arguments

- params:

  A `bgm_elicited` object, as returned by
  [`elicit_prior`](https://bayesian-graphical-modelling-lab.github.io/designbgm/reference/elicit_prior.md).

- estimator:

  Character vector selecting the estimators, or `NULL` (default) to
  compute all available estimators for the family/prior. Several
  estimators are available `"VR"`, `"PR"`, `"MTM"`, `"PT"`, and
  `"ELIR"`.

- ...:

  Family-specific arguments passed to the method (e.g., `n_samples`, the
  Monte Carlo sample size for simulation-based estimators).

## Value

A `prior_ess` object. Its `estimates` component is a named numeric
vector with one entry per each specified estimator. `print` and
`summary` methods ara available.

## Details

Methods are available for the following classes:

- [`ggm_elicited`](https://bayesian-graphical-modelling-lab.github.io/designbgm/reference/prior_ess.ggm_elicited.md):

  Gaussian graphical models

When `estimator` is `NULL` all available estimators are computed
together using an optimized routine that shares intermediate quantities
and avoids recomputation.

## See also

Other prior ess:
[`prior_ess.ggm_elicited()`](https://bayesian-graphical-modelling-lab.github.io/designbgm/reference/prior_ess.ggm_elicited.md)
