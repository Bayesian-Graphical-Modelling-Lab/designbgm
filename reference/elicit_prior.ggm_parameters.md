# Elicit a prior for GGM parameters

Elicit a prior for GGM parameters

## Usage

``` r
# S3 method for class 'ggm_parameters'
elicit_prior(params, prior = NULL, ...)
```

## Arguments

- params:

  A `ggm_parameters` object, as returned by
  [`ggm_parameters`](https://bayesian-graphical-modelling-lab.github.io/designbgm/reference/ggm_parameters.md).
  It also inherits from `bgm_parameters`.

- prior:

  Character scalar naming the precision prior, or `NULL` (default) to
  infer it from the graph: `"wishart"` for a complete `G`, `"gwishart"`
  for a sparse `G`. Supplying `prior` overrides the default and is
  validated against the graph (e.g. a sparse-only prior cannot be placed
  on a complete graph).

- ...:

  Family-specific arguments.

## Value

A `bgm_elicited` object describing the elicited prior. It also carries a
family-specific class (for example `ggm_elicited`), on which the
downstream methods dispatch.

## Details

The elicitation currently applies the `K/nu` centering, which is exact
for Wishart, and approximate for G-Wishart.

## See also

[`elicit_prior()`](https://bayesian-graphical-modelling-lab.github.io/designbgm/reference/elicit_prior.md)
for the generic method and
[`prior_ess()`](https://bayesian-graphical-modelling-lab.github.io/designbgm/reference/prior_ess.md)
for prior ESS estimation.

Other prior elicitation:
[`elicit_prior()`](https://bayesian-graphical-modelling-lab.github.io/designbgm/reference/elicit_prior.md),
[`ggm_study()`](https://bayesian-graphical-modelling-lab.github.io/designbgm/reference/ggm_study.md),
[`simulate_prior_study()`](https://bayesian-graphical-modelling-lab.github.io/designbgm/reference/simulate_prior_study.md),
[`simulate_prior_study.ggm_study()`](https://bayesian-graphical-modelling-lab.github.io/designbgm/reference/simulate_prior_study.ggm_study.md)
