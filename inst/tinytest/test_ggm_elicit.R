# Test elicitation of the prior for GGM.
p <- 3L
G <- matrix(0, p, p); G[1, 2] <- G[2, 1] <- 1
K <- diag(p); K[1, 2] <- K[2, 1] <- 0.3

# sparse -> gwishart, inferred
ep <- elicit_prior(ggm_parameters(K, G, nu = 10))
expect_inherits(ep, "ggm_elicited")
expect_inherits(ep, "bgm_elicited")
expect_equal(ep$prior, "gwishart")
expect_true(ep$sparse)
expect_equal(ep$p, p)
expect_equal(ep$scale, K / 10)          # K/nu centering
expect_stdout(print(ep), "gwishart")

# complete -> wishart, inferred
Gc <- matrix(1, p, p); diag(Gc) <- 0
Kc <- diag(p) + 0.2; diag(Kc) <- 1.5
epc <- elicit_prior(ggm_parameters(Kc, Gc, nu = 10))
expect_equal(epc$prior, "wishart")
expect_false(epc$sparse)

# regime mismatch both ways
expect_error(elicit_prior(ggm_parameters(K, G, nu = 10), prior = "wishart"),
             "not used for a sparse graph")
expect_error(elicit_prior(ggm_parameters(Kc, Gc, nu = 10), prior = "gwishart"),
             "not used for a complete graph")

# unknown prior
expect_error(elicit_prior(ggm_parameters(K, G, nu = 10), prior = "nonsense"))

# elicited object carries everything the downstream methods read
expect_equal(ep$G, G)
expect_equal(ep$nu, 10L)
expect_equal(ep$family, "ggm")
expect_stdout(print(epc), "wishart")
