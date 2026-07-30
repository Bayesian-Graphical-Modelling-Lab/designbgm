# Sample-size planning

Recommends a sample size `n*` for a prospective study, given an elicited
prior described by a `bgm_elicited` object. The planning criterion
depends on `method`, and the methods available depend on the class of
`params`. See the method pages listed under Details.

## Usage

``` r
design(params, method = c("DPIR", "BFDA", "BSDA"), ...)
```

## Arguments

- params:

  A `bgm_elicited` object, as returned by
  [`elicit_prior`](https://bayesian-graphical-modelling-lab.github.io/designbgm/reference/elicit_prior.md).

- method:

  Which planning method to use, `"DPIR"` (default), `"BFDA"`, or
  `"BSDA"`. One method runs per call.

- ...:

  Method-specific arguments passed to the method (e.g., `max_n`, the
  largest sample size considered in the planning).

## Value

A `bgm_design` object. It records the recommended sample sizes ( two per
call, depending on `method`) together with the criterion values and the
elicited prior the planning was built from. `print` and `summary`
methods are available.

## Details

Both methods look for the smallest sample size at which a criterion
reaches a target threshold with a specified probability. They differ in
the criterion: `"DPIR"` uses a data-to-prior information ratio over the
model parameters, `"BFDA"` the Bayes factor at a representative edge.

Methods are available for the following classes:

- [`ggm_elicited`](https://bayesian-graphical-modelling-lab.github.io/designbgm/reference/design.ggm_elicited.md):

  Gaussian graphical models

## See also

Other sample size planning:
[`bsda_control()`](https://bayesian-graphical-modelling-lab.github.io/designbgm/reference/bsda_control.md),
[`design.ggm_elicited()`](https://bayesian-graphical-modelling-lab.github.io/designbgm/reference/design.ggm_elicited.md),
[`validate()`](https://bayesian-graphical-modelling-lab.github.io/designbgm/reference/validate.md),
[`validate.ggm_design()`](https://bayesian-graphical-modelling-lab.github.io/designbgm/reference/validate.ggm_design.md)
