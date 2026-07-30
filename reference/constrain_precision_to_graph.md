# Project a precision matrix onto a fixed graph (constrained estimate)

Given a precision matrix (`K`) and a fixed undirected graph (`G`),
returns the precision matrix constrained to the graph's structure, so
that entries at non-edges are (numerically) zero. Uses the node-wise
regression algorithm.

## Usage

``` r
constrain_precision_to_graph(K, G, tol = 1e-06, itermax = 1000)
```

## Arguments

- K:

  A \\p \times p\\ symmetric positive-definite precision matrix.

- G:

  A \\p \times p\\ symmetric adjacency matrix with 0/1 entries; the zero
  pattern defines the conditional-independence constraints.

- tol:

  Convergence tolerance for the iterative solver.

- itermax:

  Maximum number of iterations.

## Value

A \\p \times p\\ precision matrix with zeros at the non-edges of `G`.
