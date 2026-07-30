# Validate a sample-size plan

Checks a sample-size plan by simulation: draws studies at the
recommended `n` and reports how often the criterion reaches its target
there. The checks available depend on the class of `plan` and on the
planning method used. See the method pages listed under Details.

## Usage

``` r
validate(plan, ...)
```

## Arguments

- plan:

  A `bgm_design` object, as returned by
  [`design`](https://bayesian-graphical-modelling-lab.github.io/designbgm/reference/design.md).

- ...:

  Method-specific arguments passed to the method (e.g., `n_sim`, the
  number of simulated studies).

## Value

A `bgm_design_validation` object. `print` and `summary` methods are
available.

## Details

Methods are available for the following classes:

- [`ggm_design`](https://bayesian-graphical-modelling-lab.github.io/designbgm/reference/validate.ggm_design.md):

  Gaussian graphical models

## See also

Other sample size planning:
[`bsda_control()`](https://bayesian-graphical-modelling-lab.github.io/designbgm/reference/bsda_control.md),
[`design()`](https://bayesian-graphical-modelling-lab.github.io/designbgm/reference/design.md),
[`design.ggm_elicited()`](https://bayesian-graphical-modelling-lab.github.io/designbgm/reference/design.ggm_elicited.md),
[`validate.ggm_design()`](https://bayesian-graphical-modelling-lab.github.io/designbgm/reference/validate.ggm_design.md)
