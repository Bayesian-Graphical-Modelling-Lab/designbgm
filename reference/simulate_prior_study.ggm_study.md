# Simulate prior studies for a Gaussian graphical model

Simulates prior studies from the Gaussian likelihood: for each study,
`nu` observations are drawn from a graph-respecting precision matrix and
the precision matrix is re-estimated under the same graph. What comes
back is what the study found, not the truth it was drawn from, so it
carries the estimation noise of a study that size.

## Usage

``` r
# S3 method for class 'ggm_study'
simulate_prior_study(study, n_studies = 1L, verbose = FALSE, ...)
```

## Arguments

- study:

  A `ggm_study` object, as returned by
  [`ggm_study`](https://bayesian-graphical-modelling-lab.github.io/designbgm/reference/ggm_study.md).

- n_studies:

  Number of independent studies to simulate. Default 1.

- verbose:

  Print diagnostics for each simulated precision matrix. Off by default.

- ...:

  Ignored, present for consistency with the generic.

## Value

A list of `n_studies` `ggm_parameters` objects, ready for
[`elicit_prior`](https://bayesian-graphical-modelling-lab.github.io/designbgm/reference/elicit_prior.md).

## Details

Studies are independent. If the study was specified with a fixed `G`,
every study uses it and only the estimated precision matrix varies; if
it was specified with a `structure`, a new graph is drawn for each
study. Call [`set.seed`](https://rdrr.io/r/base/Random.html) first for
reproducible output.

## See also

Other prior elicitation:
[`elicit_prior()`](https://bayesian-graphical-modelling-lab.github.io/designbgm/reference/elicit_prior.md),
[`elicit_prior.ggm_parameters()`](https://bayesian-graphical-modelling-lab.github.io/designbgm/reference/elicit_prior.ggm_parameters.md),
[`ggm_study()`](https://bayesian-graphical-modelling-lab.github.io/designbgm/reference/ggm_study.md),
[`simulate_prior_study()`](https://bayesian-graphical-modelling-lab.github.io/designbgm/reference/simulate_prior_study.md)
