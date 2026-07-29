# designbgm 0.3.0

## New features

* New planning method `"BSDA"` (Bayesian structural design analysis) in
  `design()`. It recommends a sample size at which edge selection reaches a
  target sensitivity or specificity with a given power, using the probit inversion
  of a simulated power curve.
* `ggm_parameters()` gains a `pip` argument (prior inclusion probabilities),
  supplied as a single probability or a matrix. Required by BSDA and carried
  through `elicit_prior()`.
* `bsda_control()` is an helper control list where the user can specify:
  the sampler, model-fit, and search settings for BSDA planning.
* `validate()` now supports BSDA plans, reporting the achieved power at the
  recommended sample size with a confidence interval.

## Minor improvements

* `print()` methods for `ggm_parameters` and `ggm_elicited` now display the
  prior inclusion probabilities (range) when supplied.

# designbgm 0.2.0

## New features

* New exported function `constrain_precision_to_graph()`, which projects a
  precision matrix onto a fixed undirected graph and returns the precision
  matrix constrained to the graph's zero pattern (where entries at non-edges 
  are set to zero).

* `design.ggm_elicited()` now includes an `nsim_bf` argument controlling the number of
  Monte Carlo simulations used in the Bayes-factor computation for the sparse
  (G-Wishart) BFDA path. Defaults to `1000L`, matching the underlying C++
  default. Finally, the value is recorded in the returned design object's `call_info`
  for sparse BFDA designs.

* The DPIR and dense BFDA design routines now report timing information, which
  in 0.1.0 was available only for the sparse BFDA path. Durations are returned
  in the design object: `duration_bisection_global` and `duration_bisection_pw`
  for DPIR, `duration_h0` and `duration_h1` for BFDA.

# designbgm 0.1.0

* First Github version.
