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
