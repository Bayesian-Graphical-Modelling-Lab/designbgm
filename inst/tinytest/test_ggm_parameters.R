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
