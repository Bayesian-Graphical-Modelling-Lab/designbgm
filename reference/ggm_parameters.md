# Construct GGM prior parameters

Takes the inputs for a Gaussian graphical model and bundles them into
one object, ready for prior elicitation with
[`elicit_prior`](https://bayesian-graphical-modelling-lab.github.io/designbgm/reference/elicit_prior.md).

## Usage

``` r
ggm_parameters(K, G, nu, mu = NULL, pip = NULL)
```

## Arguments

- K:

  A `p` by `p` symmetric positive-definite precision matrix, .

- G:

  A `p` by `p` symmetric 0/1 adjacency matrix with zero diagonal.

- nu:

  Prior study size (degrees of freedom).

- mu:

  Mean vector of length `p`. The default `NULL` corresponds to the
  centered (zero-mean) case that all current methods assume. Supplying a
  non-`NULL` `mu` is not yet supported.

- pip:

  Prior inclusion probabilities. The default `NULL` corresponds to the
  uniform prior. Can be a single probability or a `p` by `p` symmetric
  matrix with zero diagonal.

## Value

A `ggm_parameters` object, which also inherits from `bgm_parameters`.
