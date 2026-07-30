# Control parameters for BSDA sample-size planning

Lists the sampler, model-fit, and search-algorithm settings for
`design(method = "BSDA")` and its validation. These are the deeper knobs
a user rarely changes; the planning targets (`measure`, `measure_value`,
`target_pow`) stay in
[`design`](https://bayesian-graphical-modelling-lab.github.io/designbgm/reference/design.md)'s
own arguments.

## Usage

``` r
bsda_control(
  gwish_sampler = "direct",
  gwish_tol = 1e-08,
  gwish_iter = 500L,
  gwish_burnin = 500L,
  edge_threshold = 0.5,
  fit_iterations = 10000L,
  fit_burnin = 5000L,
  alpha = 0.05,
  n_scout = 6L,
  n_main = 10L,
  scout_frac = 1/3,
  max_iter = 10L,
  n_boot = 5000L,
  eps = 0.001,
  tol_frac = 0.01,
  verbose = FALSE,
  seed = NULL
)
```

## Arguments

- gwish_sampler:

  G-Wishart sampler, `"direct"` or `"block"`.

- gwish_tol, gwish_iter, gwish_burnin:

  Tolerance, iterations, and burn-in for the G-Wishart sampler.

- edge_threshold:

  Posterior inclusion probability above which an edge is selected.
  Default 0.5.

- fit_iterations, fit_burnin:

  MCMC length and burn-in for fitting each simulated study.

- alpha:

  Significance level used in the fit. Default 0.05.

- n_scout, n_main:

  Grid sizes for the scout and main search passes.

- scout_frac:

  Fraction of `fit_iterations` used in the scout pass.

- max_iter:

  Maximum refine-evaluate-invert iterations.

- n_boot:

  Bootstrap resamples for the probit-inversion CI.

- eps:

  Numerical tolerance guarding the probit inversion.

- tol_frac:

  Tolerance fraction of the relative sample size difference used to stop
  the search. Default 0.01, matched to typical H = 50. Lower it if H
  raises substantially.

- verbose:

  Stream progress from the C++ routine. Default `FALSE`.

## Value

A named list of control settings.

## See also

Other sample size planning:
[`design()`](https://bayesian-graphical-modelling-lab.github.io/designbgm/reference/design.md),
[`design.ggm_elicited()`](https://bayesian-graphical-modelling-lab.github.io/designbgm/reference/design.ggm_elicited.md),
[`validate()`](https://bayesian-graphical-modelling-lab.github.io/designbgm/reference/validate.md),
[`validate.ggm_design()`](https://bayesian-graphical-modelling-lab.github.io/designbgm/reference/validate.ggm_design.md)
