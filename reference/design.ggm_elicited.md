# Sample-size planning for a Gaussian graphical model

Recommends a sample size for a prospective Gaussian graphical model
study, given a prior elicited from a previous one. `"DPIR"` returns two
sizes, one for the graph as a whole and one at the parameterwise level
(considering only the off-diagonal elements). `"BFDA"` returns one per
hypothesis: under the null (edge absent) and under the alternative (edge
present). `"BSDA"` returns a single size at which edge selection reaches
a target sensitivity or specificity with a given power.

## Usage

``` r
# S3 method for class 'ggm_elicited'
design(
  params,
  method = c("DPIR", "BFDA", "BSDA"),
  H = 150L,
  J = 20L,
  n = NULL,
  threshold = NULL,
  n_tol = 1L,
  max_n = 5000L,
  target_probability = 0.95,
  rho_quantile = 0.5,
  pow0 = 0.8,
  pow1 = 0.8,
  nsim_bf = 1000L,
  measure = c("sen", "spe"),
  measure_value = 0.8,
  target_pow = 0.8,
  range_lower = 2L,
  control = bsda_control(),
  ...
)
```

## Arguments

- params:

  A `ggm_elicited` object, as returned by
  [`elicit_prior`](https://bayesian-graphical-modelling-lab.github.io/designbgm/reference/elicit_prior.md).
  It also inherits from `bgm_elicited`.

- method:

  Which planning method to use, `"DPIR"` (default), `"BFDA"`, or
  `"BSDA"`. One method runs per call.

- H, J:

  Outer and inner Monte Carlo replication counts. `H` draws from the
  prior (default 150), `J` datasets per draw (default 20).

- n:

  Candidate sample sizes to search over. `NULL` (default) picks a grid
  automatically and refines it by bisection.

- threshold:

  DPIR and BFDA only. Target threshold the criterion must reach: the
  DPIR ratio for `"DPIR"`, the Bayes factor for `"BFDA"`. `NULL`
  (default) uses 1.0 for DPIR and 10.0 for BFDA. BSDA uses
  `measure_value` instead.

- n_tol:

  Search tolerance, in observation units. The search stops once the
  bracket reaches this width. Default 1, i.e. exact (increasing this
  value speeds up the search but may reduce accuracy).

- max_n:

  Largest sample size considered. If the target is not met by `max_n`,
  the search stops and reports non-convergence. Default 5000.

- target_probability:

  DPIR only: how often the ratio must reach `threshold`. Default 0.95.

- rho_quantile:

  BFDA only: which edge to plan around. Edges are ordered by the
  magnitude of their partial correlation and the one at this quantile is
  selected. Default 0.5, the median edge. The selected edge is reported
  with its signed partial correlation.

- pow0, pow1:

  BFDA only: how often the Bayes factor must give decisive evidence for
  the hypothesis that is true. With `BF01` the Bayes factor for edge
  absence, `pow0` targets `Pr(BF01 > threshold | H0)`, correctly
  excluding the edge, and `pow1` targets `Pr(BF01 < 1/threshold | H1)`,
  correctly detecting the edge. Both default to 0.8.

- nsim_bf:

  BFDA only: number of Monte Carlo samples to estimate the Bayes factor
  for sparse graphs. Default 1000. Ignored for dense graphs, where the
  Bayes factor is available in closed form.

- measure, measure_value, target_pow:

  BSDA only: the measure to plan for (`"sen"` or `"spe"`), the target
  value of that measure, and the target power to achieve it. Defaults
  are `"sen"`, 0.8, and 0.8.

- range_lower:

  BSDA only: the lower bound of the sample size search range. Default 2.

- control:

  BSDA only: a list of control parameters, as returned by
  [`bsda_control`](https://bayesian-graphical-modelling-lab.github.io/designbgm/reference/bsda_control.md).
  Default
  [`bsda_control()`](https://bayesian-graphical-modelling-lab.github.io/designbgm/reference/bsda_control.md).
  \#' @param ... Ignored, present for consistency with the generic
  [`design()`](https://bayesian-graphical-modelling-lab.github.io/designbgm/reference/design.md).

## Value

A `ggm_design` object, which also inherits from `bgm_design`. Its
`results` component contains the recommended sizes: `n_star_global` and
`n_star_pw` for `"DPIR"`; `n_star_power_h0` and `n_star_power_h1` for
`"BFDA"`, each with a convergence flag, and a single `n_star` with a
confidence interval (`ci_lo`, `ci_hi`) for `"BSDA"`. For DPIR and BFDA,
a sample size not reached within the search is `NA` with a `FALSE`
convergence flag. For BSDA, `identified` is `FALSE` when the fit could
not be established (too few non-saturated grid points), and `converged`
is `FALSE` when the search stopped at `max_iter` without stabilizing.
`print` and `summary` methods are available.

## Details

Each method targets a criterion at a required frequency, but the
criterion and the arguments differ:

- DPIR:

  Targets a data-to-prior information ratio. Uses `threshold` (the ratio
  to surpass, default is 1) and `target_probability` (how often, default
  is 0.95).

- BFDA:

  Targets a Bayes factor at a representative edge. Uses `threshold` (the
  Bayes factor value to surpass, default is 10), `pow0`/`pow1` (the
  target power under each hypothesis), and `rho_quantile` (which edge to
  plan around).

- BSDA:

  Targets a selection criterion (sensitivity or specificity) at a given
  power. Uses `measure` (`"sen"` or `"spe"`), `measure_value` (the level
  to reach, default 0.8), `target_pow` (the target power, default 0.8),
  `range_lower` (the smallest `n` searched), and `control` (sampler and
  fit settings; see
  [`bsda_control`](https://bayesian-graphical-modelling-lab.github.io/designbgm/reference/bsda_control.md)).
  Requires posterior inclusion probabilities, supplied as `pip` to
  [`ggm_parameters`](https://bayesian-graphical-modelling-lab.github.io/designbgm/reference/ggm_parameters.md),
  and a sparse graph.

Arguments not required by a method are ignored.

Under BFDA the edge is tested with `H0` : "the edge is absent" against
`H1` : "the edge is present", and the Bayes factor reported throughout
is `BF01`, the Bayes factor for absence: `BF01 > threshold` excludes the
edge, `BF01 < 1/threshold` detects it. This is the reciprocal of the
inclusion Bayes factor `BF10`.

Under BSDA the search proceeds in two stages: a coarse exploration
locates the region where the target power is met, then repeatedly fits a
probit model to the power curve, inverts it for `n*`, and narrows the
grid, until the recommendation stabilizes. Because the recommendation
comes from a fitted curve rather than a direct threshold crossing, it
carries a confidence interval, and it can fail in two ways:
`identified = FALSE` (the fit could not be established) or
`converged = FALSE` (the loop hit `max_iter` without stabilizing).

Cost grows with `H * J` and with the width of the search: the criterion
is evaluated by simulation at every candidate `n`.

## See also

Other sample size planning:
[`bsda_control()`](https://bayesian-graphical-modelling-lab.github.io/designbgm/reference/bsda_control.md),
[`design()`](https://bayesian-graphical-modelling-lab.github.io/designbgm/reference/design.md),
[`validate()`](https://bayesian-graphical-modelling-lab.github.io/designbgm/reference/validate.md),
[`validate.ggm_design()`](https://bayesian-graphical-modelling-lab.github.io/designbgm/reference/validate.ggm_design.md)
