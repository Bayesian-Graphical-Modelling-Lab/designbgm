# Test the creation and simulation of a GGM study object.
p <- 3L
G <- matrix(0, p, p); G[1, 2] <- G[2, 1] <- 1

# -- ggm_study validation (pure R) 
st <- ggm_study(p = p, nu = 20, G = G)
expect_inherits(st, "ggm_study")
expect_equal(st$inputs$p, 3L)
expect_equal(st$inputs$nu, 20L)
expect_true(is.integer(st$inputs$p))

expect_error(ggm_study(p = 2, nu = 20, G = G))
expect_error(ggm_study(p = p, nu = 3, G = G), "must exceed 'p'")
expect_error(ggm_study(p = p, nu = 20), "exactly one of")            # neither
expect_error(ggm_study(p = p, nu = 20, G = G, structure = "random"), "exactly one of")  # both
Gbad <- G; diag(Gbad) <- 1
expect_error(ggm_study(p = p, nu = 20, G = Gbad))
expect_error(ggm_study(p = p, nu = 20, G = matrix(0, 2, 2)))

st2 <- ggm_study(p = p, nu = 20, structure = "Bernoulli", prob = 0.5)
expect_null(st2$inputs$G)
expect_equal(st2$inputs$structure, "Bernoulli")
expect_equal(st2$inputs$gen_args$prob, 0.5)

# --- simulate_prior_study.ggm_study (fixed G)
set.seed(1)
out <- simulate_prior_study(st, n_studies = 2L)
expect_equal(length(out), 2L)
expect_inherits(out[[1]], "ggm_parameters")
expect_equal(out[[1]]$nu, 20L)
expect_equal(out[[1]]$G, G)               # fixed G reused
expect_equal(out[[2]]$G, G)
expect_true(!identical(out[[1]]$K, out[[2]]$K))   # only K varies

# output feeds elicit_prior
expect_inherits(elicit_prior(out[[1]]), "ggm_elicited")

# default n_studies = 1
set.seed(1)
expect_equal(length(simulate_prior_study(st)), 1L)

expect_error(simulate_prior_study(st, n_studies = 0))

# -- generated G per study
set.seed(1)
out2 <- simulate_prior_study(st2, n_studies = 2L)
expect_equal(length(out2), 2L)
expect_inherits(out2[[1]], "ggm_parameters")