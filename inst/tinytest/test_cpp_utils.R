# Tests for the cpp_* exports that have no R wrapper.
# These are called directly by family-ggm.R, so any silent change inside them
# (indexing base, for instance) breaks callers without any R-level error.

# -- cpp_precision_to_partial_correlations
set.seed(1)
p <- 4L
A <- matrix(rnorm(p * p), p, p)
K <- crossprod(A) + diag(p)          # symmetric positive definite

R <- designbgm:::cpp_precision_to_partial_correlations(K)

expect_equal(dim(R), c(p, p))
expect_equal(R, t(R))
expect_true(all(is.finite(R)))

# off-diagonal: -K_ij / sqrt(K_ii K_jj)
d   <- sqrt(diag(K))
ref <- -K / outer(d, d)
off <- row(K) != col(K)
expect_equal(R[off], ref[off], tolerance = 1e-10)

# diagonal convention: ggm_parameters() and .bfda_planning_edge() both assume
# ones in the diagonal.
expect_equal(diag(R), rep(1, p))

# a diagonal K has zero partial correlations everywhere off-diagonal.
Rd <- designbgm:::cpp_precision_to_partial_correlations(diag(p))
expect_equal(Rd[off], rep(0, sum(off)))

# --- cpp_G_upper / cpp_n_vec
pg <- 4L
G <- matrix(0L, pg, pg)
G[1, 2] <- G[2, 1] <- 1L
G[3, 4] <- G[4, 3] <- 1L

GU <- designbgm:::cpp_G_upper(pg, G)
expect_true(is.matrix(GU))
expect_equal(dim(GU), c(pg, pg))
expect_true(all(is.finite(GU)))

# .bfda_planning_edge() reads the edges back with which(G_upper == 1, arr.ind),
# expecting 1-based (row, col) pairs with row < col.
ed <- which(GU == 1L, arr.ind = TRUE)
ed <- ed[order(ed[, 1], ed[, 2]), , drop = FALSE]
expect_equal(nrow(ed), 2L)
expect_true(all(ed[, 1] < ed[, 2]))
expect_equal(unname(ed[1, ]), c(1L, 2L))
expect_equal(unname(ed[2, ]), c(3L, 4L))

nv <- designbgm:::cpp_n_vec(pg, GU)
expect_true(is.numeric(nv))
expect_true(all(is.finite(nv)))


# --- constrain_precision_to_graph()
p <- 3L
# symmetric, positive-definite precision matrix
K <- matrix(c(1.0, 0.3, 0.2,
              0.3, 1.0, 0.1,
              0.2, 0.1, 1.0), p, p)

# single edge (1, 2); non-edges are (1, 3) and (2, 3)
G <- matrix(0, p, p)
G[1, 2] <- G[2, 1] <- 1

# --- input validation (fires before the C++ call)
expect_error(constrain_precision_to_graph(matrix(1:6, 2, 3), G), "square")
expect_error(constrain_precision_to_graph(K, matrix(0, 2, 2)), "same dimensions")

Kasym <- K
Kasym[1, 2] <- 0.9                 # break symmetry
expect_error(constrain_precision_to_graph(Kasym, G), "symmetric")

expect_error(constrain_precision_to_graph(matrix(0, p, p), G), "positive definite")

# --- valid call: result constrained to the graph
Kc <- constrain_precision_to_graph(K, G)
expect_equal(dim(Kc), c(p, p))
expect_true(isSymmetric(unname(Kc), tol = 1e-6))
expect_equal(Kc[1, 3], 0, tolerance = 1e-6)    # zero at non-edge (1, 3)
expect_equal(Kc[2, 3], 0, tolerance = 1e-6)    # zero at non-edge (2, 3)
expect_true(abs(Kc[1, 2]) > 0)                 # present edge survives

# --- invariant: complete graph is a no-op
Gfull <- matrix(1, p, p)
diag(Gfull) <- 0
expect_equal(unname(constrain_precision_to_graph(K, Gfull)), unname(K),
             tolerance = 1e-4)