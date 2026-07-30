# Validate a sample-size plan for a Gaussian graphical model

Checks a design plan by simulating studies at the recommended `n` and
reporting how often the criterion actually reaches its target threshold.
Usually, the planning search uses fewer replications than this check, so
the two probabilities will not match exactly.

## Usage

``` r
# S3 method for class 'ggm_design'
validate(
  plan,
  H = 500L,
  J = 100L,
  which_n = NULL,
  scope = NULL,
  seed = NULL,
  ...
)
```

## Arguments

- plan:

  A `ggm_design` object, as returned by
  [`design`](https://bayesian-graphical-modelling-lab.github.io/designbgm/reference/design.md).
  It also inherits from `bgm_design`.

- H, J:

  Outer and inner Monte Carlo replication counts, as in
  [`design`](https://bayesian-graphical-modelling-lab.github.io/designbgm/reference/design.md).
  Defaults 500 and 100, larger than the planning defaults, since the
  check is run at a single `n`.

- which_n:

  Which of the plan's recommended sizes to validate at. `NULL` (default)
  picks `"global"` for a DPIR plan and `"h1"` for a BFDA plan. DPIR
  plans also accept `"pw"`, the parameterwise size; BFDA plans also
  accept `"h0"`; BSDA plans have only one size.

- scope:

  BFDA only: whether to check the edge the plan was built around
  (`"planning_edge"`) or every present edge in the graph
  (`"all_edges"`), the stricter guarantee check. `NULL` (default) means
  `"planning_edge"`. This argument is ignored for DPIR plans.

- seed:

  Random seed for reproducibility.

- ...:

  Ignored, present for consistency with the generic.

## Value

A `ggm_design_validation` object, which also inherits from
`bgm_design_validation`. Its `n_star` component is the size that was
checked and `results` holds the probability achieved there. `print` and
`summary` methods are available.

## See also

Other sample size planning:
[`bsda_control()`](https://bayesian-graphical-modelling-lab.github.io/designbgm/reference/bsda_control.md),
[`design()`](https://bayesian-graphical-modelling-lab.github.io/designbgm/reference/design.md),
[`design.ggm_elicited()`](https://bayesian-graphical-modelling-lab.github.io/designbgm/reference/design.ggm_elicited.md),
[`validate()`](https://bayesian-graphical-modelling-lab.github.io/designbgm/reference/validate.md)
