# Test the graph generators.

# -- Bernoulli
set.seed(1)
G <- generate_graph(p = 6, structure = "Bernoulli", prob = 0.4)

expect_true(is.matrix(G))
expect_equal(dim(G), c(6L, 6L))
expect_equal(G, t(G))
expect_equal(diag(G), rep(0, 6))
expect_true(all(G %in% c(0, 1)))

edges <- sum(G) / 2
expect_true(edges > 0)
expect_true(edges < 6 * 5 / 2)

# --- all generators: shape + non-degeneracy
set.seed(1)
for (s in c("random", "scalefree", "smallworld")) {
  G <- generate_graph(p = 8, structure = s)
  expect_equal(dim(G), c(8L, 8L))
  expect_equal(G, t(G))
  expect_equal(diag(G), rep(0, 8))
  expect_true(all(G %in% c(0, 1)))
  e <- sum(G) / 2
  expect_true(e > 0 && e < 8 * 7 / 2)
}

# -- reproducibility
set.seed(42); a <- generate_graph(p = 6, structure = "Bernoulli", prob = 0.4)
set.seed(42); b <- generate_graph(p = 6, structure = "Bernoulli", prob = 0.4)
expect_identical(a, b)

set.seed(42); a2 <- generate_graph(p = 8, structure = "scalefree")
set.seed(42); b2 <- generate_graph(p = 8, structure = "scalefree")
expect_identical(a2, b2)

# --- argument validation
expect_error(generate_graph(p = 2, structure = "Bernoulli", prob = 0.4))
expect_error(generate_graph(p = 3.5, structure = "Bernoulli", prob = 0.4))
expect_error(generate_graph(p = c(4, 5), structure = "Bernoulli", prob = 0.4))
expect_error(generate_graph(p = "six", structure = "Bernoulli", prob = 0.4))

expect_error(generate_graph(p = 6, structure = "Bernoulli"), "'prob' is required")
expect_error(generate_graph(p = 6, structure = "Bernoulli", prob = 0), "in \\(0, 1\\)")
expect_error(generate_graph(p = 6, structure = "Bernoulli", prob = 1), "in \\(0, 1\\)")
expect_error(generate_graph(p = 6, structure = "Bernoulli", prob = 1.5), "in \\(0, 1\\)")

expect_error(generate_graph(p = 6, structure = "nonsense"))
expect_error(generate_graph(p = 6, structure = c("random", "Bernoulli")))

# -- max_attempts
set.seed(1)
expect_error(generate_graph(p = 20, structure = "Bernoulli", prob = 1e-9,
                            max_attempts = 1L),
             "No non-degenerate graph")

# --- .graph_is_sparse (internal)
Gc <- matrix(1, 4, 4); diag(Gc) <- 0
expect_false(designbgm:::.graph_is_sparse(Gc))
expect_true(designbgm:::.graph_is_sparse(matrix(0, 4, 4)))
Gs <- matrix(0, 4, 4); Gs[1, 2] <- Gs[2, 1] <- 1
expect_true(designbgm:::.graph_is_sparse(Gs))