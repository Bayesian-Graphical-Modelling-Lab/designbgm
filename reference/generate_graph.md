# Generate a random graph with a given structure

Draws a random undirected graph on `p` nodes from one of the available
structure generators, rejecting degenerate draws (the empty graph and
the complete graph) and resampling until a usable structure is obtained.
The returned adjacency matrix is what
[`simulate_prior_study`](https://bayesian-graphical-modelling-lab.github.io/designbgm/reference/simulate_prior_study.md)
uses for simulating prior study parameters.

## Usage

``` r
generate_graph(
  p,
  structure = c("smallworld", "random", "scalefree", "Bernoulli"),
  ...,
  max_attempts = 1000L
)
```

## Arguments

- p:

  Number of nodes (variables) in the graph. Must be at least 3.

- structure:

  Which graph generator to use: `"smallworld"`, `"random"`,
  `"scalefree"` or `"Bernoulli"`.

- ...:

  Further arguments passed to the graph generator. The `"Bernoulli"`
  structure requires `prob`, the edge-inclusion probability in \\(0,
  1)\\. The other structures take no extra arguments.

- max_attempts:

  Maximum number of rejection-sampling draws before giving up with an
  error. Guards against parameter choices for which a non-degenerate
  graph is effectively unreachable.

## Value

A symmetric `p` by `p` adjacency matrix, 0/1 entries, zero diagonal.

## Details

Generation uses rejection sampling: a draw is accepted only if it has at
least one edge and is not fully connected. Since no non-degenerate graph
exists for `p < 3`, `p` must be at least 3

## Examples

``` r
set.seed(2026)
# Bernoulli (Erdos-Renyi) needs an edge probability
g <- generate_graph(p = 10, structure = "Bernoulli", prob = 0.2)
dim(g)
#> [1] 10 10

# \donttest{
g2 <- generate_graph(p = 10, structure = "smallworld")
# }
```
