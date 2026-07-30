# Prior effective sample size for a Gaussian graphical model

Computes the prior ESS of an elicited prior distribution for the
parameters in a Gaussian graphical model. The matrix estimators (VR, PR)
work with matrices (Fisher information matrix and prior variance matrix)
and need to be reduced to a scalar first, the remaining estimators do
not. For the G-Wishart prior the quantities involved have no closed form
and are approximated by Monte Carlo simulations.

## Usage

``` r
# S3 method for class 'ggm_elicited'
prior_ess(
  params,
  estimator = NULL,
  aggregation = c("det", "tr", "mean"),
  sampler = c("direct", "gibbs"),
  n_samples = 1000L,
  tol = 1e-06,
  itermax = 1000L,
  burnin = 500L,
  init = NULL,
  compute_cond = FALSE,
  ...
)
```

## Arguments

- params:

  A `ggm_elicited` object, as returned by
  [`elicit_prior`](https://bayesian-graphical-modelling-lab.github.io/designbgm/reference/elicit_prior.md).

- estimator:

  Character vector selecting the estimators, or `NULL` (default) to
  compute all available estimators for the family/prior. Several
  estimators are available `"VR"`, `"PR"`, `"MTM"`, `"PT"`, and
  `"ELIR"`.

- aggregation:

  Character scalar selecting the aggregation method to reduce the
  estimator to a scalar: `"det"` (default), `"tr"`, or `"mean"`. Applies
  only to the matrix estimators (VR, PR).

- sampler:

  Monte Carlo sampler for the G-Wishart prior, either `"direct"`
  (default) or `"gibbs"` (Edge-wise Gibbs sampler). Ignored for the
  Wishart prior and the scalar estimators.

- n_samples:

  Number of Monte Carlo samples for the G-Wishart prior. Default 1000.

- tol:

  Convergence tolerance for the direct sampler. Default 1e-6.

- itermax:

  Number of maximum iterations for the direct sampler. Default 1000.
  This number caps the iterative procedure for one draw of the direct
  sampler.

- burnin:

  Number of burn-in iterations discarded by the Gibbs sampler. Default
  500.

- init:

  Starting value for the Gibbs sampler, a `p` by `p` symmetric
  positive-definite precision matrix on the same scale as the returned
  samples, which are centered on the elicited precision matrix. Entries
  where the graph has no edge must be zero: the sampler never updates
  them, so whatever is defined there is carried into every draw. `NULL`
  (default) starts from a diagonal matrix of the marginal precisions
  implied by the elicited prior
  (`diag(params$nu / diag(solve(params$scale)))`).

- compute_cond:

  Logical: whether to compute the condition number of numerator and
  denominator of VR and PR (default FALSE).

- ...:

  Ignored, present for consistency with the generic
  [`prior_ess()`](https://bayesian-graphical-modelling-lab.github.io/designbgm/reference/prior_ess.md).

## Value

A `prior_ess` object for the `ggm` family. Its `estimates` component is
a named numeric vector with one entry per estimator. `print` and
`summary` methods are available.

A `prior_ess` object for the `ggm` family. Its `estimates` component is
a named list, one element per estimator, each with a `global` value and
a `parameterwise` table (for the matrix estimators). When
`compute_cond = TRUE` those elements also carry `cond_numerator` and
`cond_denominator`. `print` and `summary` methods are available.

## Details

The sampling arguments (`sampler`, `n_samples`, `tol`, `itermax`,
`burnin`, `init`) apply to the G-Wishart prior only. Under the Wishart
prior the quantities involved are available in closed form and nothing
is simulated.

## See also

Other prior ess:
[`prior_ess()`](https://bayesian-graphical-modelling-lab.github.io/designbgm/reference/prior_ess.md)
