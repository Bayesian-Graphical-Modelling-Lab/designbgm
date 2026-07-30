# Constructs an object for simulating Gaussian graphical model prior studies

Collects the inputs needed to simulate prior studies for a Gaussian
graphical model, then passes the result to
[`simulate_prior_study`](https://bayesian-graphical-modelling-lab.github.io/designbgm/reference/simulate_prior_study.md).

## Usage

``` r
ggm_study(p, nu, G = NULL, structure = NULL, ...)
```

## Arguments

- p:

  Number of nodes. Must be at least 3.

- nu:

  Prior study size: the number of observations each simulated study
  collects. Must exceed `p`.

- G:

  A `p` by `p` symmetric 0/1 adjacency matrix with zero diagonal. Supply
  this argument to hold the graph fixed across studies, so that only the
  estimated precision matrix varies.

- structure:

  A graph generator, as in
  [`generate_graph`](https://bayesian-graphical-modelling-lab.github.io/designbgm/reference/generate_graph.md).
  Supply this argument instead of `G` to draw a new graph for each
  study.

- ...:

  Further arguments for the generator, e.g. `prob` for
  `"structure = Bernoulli"`. Only used with `structure`.

## Value

A `ggm_study` object, which also inherits from `bgm_study`.

## See also

Other prior elicitation:
[`elicit_prior()`](https://bayesian-graphical-modelling-lab.github.io/designbgm/reference/elicit_prior.md),
[`elicit_prior.ggm_parameters()`](https://bayesian-graphical-modelling-lab.github.io/designbgm/reference/elicit_prior.ggm_parameters.md),
[`simulate_prior_study()`](https://bayesian-graphical-modelling-lab.github.io/designbgm/reference/simulate_prior_study.md),
[`simulate_prior_study.ggm_study()`](https://bayesian-graphical-modelling-lab.github.io/designbgm/reference/simulate_prior_study.ggm_study.md)
