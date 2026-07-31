# Test the creation of the GGM parameters object.
p <- 3L
G <- matrix(0, p, p); G[1, 2] <- G[2, 1] <- 1
K <- diag(p); K[1, 2] <- K[2, 1] <- 0.3     # zero at non-edges (1,3), (2,3)

params <- ggm_parameters(K, G, nu = 10)
expect_inherits(params, "ggm_parameters")
expect_inherits(params, "bgm_parameters")
expect_equal(params$nu, 10L)
expect_true(is.integer(params$nu))
expect_equal(params$family, "ggm")
expect_stdout(print(params), "ggm_parameters")

# nu coerced from double, and must be whole
expect_equal(ggm_parameters(K, G, nu = 10.0)$nu, 10L)
expect_error(ggm_parameters(K, G, nu = 10.5), "must be an integer")

# nu > p
expect_error(ggm_parameters(K, G, nu = 3))
expect_error(ggm_parameters(K, G, nu = 2))

# K not symmetric
Kbad <- K; Kbad[1, 2] <- 0.9
expect_error(ggm_parameters(Kbad, G, nu = 10))

# K not positive definite
expect_error(ggm_parameters(matrix(0, p, p), G, nu = 10), "positive definite")

# G with nonzero diagonal
Gbad <- G; diag(Gbad) <- 1
expect_error(ggm_parameters(K, Gbad, nu = 10))

# G non-binary
Gbad2 <- G; Gbad2[1, 2] <- Gbad2[2, 1] <- 2
expect_error(ggm_parameters(K, Gbad2, nu = 10))

# dimension mismatch
expect_error(ggm_parameters(K, matrix(0, 2, 2), nu = 10))

# K incompatible with G: edge mass where G has none
Kinc <- K; Kinc[1, 3] <- Kinc[3, 1] <- 0.3
expect_error(ggm_parameters(Kinc, G, nu = 10), "no edge")

# mu not yet supported
expect_error(ggm_parameters(K, G, nu = 10, mu = rep(0, p)), "not yet supported")

# pip scalar is expanded to a p x p matrix with a zero diagonal
prm <- ggm_parameters(K, G, nu = 10, pip = 0.5)
expect_true(is.matrix(prm$pip))
expect_equal(dim(prm$pip), c(p, p))
expect_equal(diag(prm$pip), rep(0, p))
expect_equal(prm$pip[1, 2], 0.5)

# pip scalar out of [0, 1]
expect_error(ggm_parameters(K, G, nu = 10, pip = 1.5),  "\\[0, 1\\]")
expect_error(ggm_parameters(K, G, nu = 10, pip = -0.1), "\\[0, 1\\]")

# pip full symmetric matrix: diagonal zeroed, off-diagonals kept
S <- matrix(0.3, p, p); diag(S) <- 0.9
prm <- ggm_parameters(K, G, nu = 10, pip = S)
expect_true(isSymmetric(prm$pip))
expect_equal(diag(prm$pip), rep(0, p))
expect_equal(prm$pip[1, 2], 0.3)

# pip upper triangle only: mirrored into a symmetric matrix
U <- matrix(0, p, p)
U[1, 2] <- 0.7; U[1, 3] <- 0.2; U[2, 3] <- 0.5
prm <- ggm_parameters(K, G, nu = 10, pip = U)
expect_true(isSymmetric(prm$pip))
expect_equal(prm$pip[2, 1], 0.7)          # lower filled from upper
expect_equal(prm$pip[3, 1], 0.2)
expect_equal(prm$pip[3, 2], 0.5)
expect_equal(diag(prm$pip), rep(0, p))

# pip wrong dimensions
expect_error(ggm_parameters(K, G, nu = 10, pip = matrix(0.5, 2, 2)))

# pip lower triangle filled but not symmetric
A <- matrix(0, p, p); A[1, 2] <- 0.7; A[2, 1] <- 0.3
expect_error(ggm_parameters(K, G, nu = 10, pip = A), "symmetric")

# pip entries outside [0, 1]
B <- matrix(0.3, p, p); B[1, 2] <- B[2, 1] <- 1.5
expect_error(ggm_parameters(K, G, nu = 10, pip = B), "\\[0, 1\\]")

# pip neither a scalar probability nor a matrix
expect_error(ggm_parameters(K, G, nu = 10, pip = c(0.3, 0.4)),
             "single probability or a p x p matrix")