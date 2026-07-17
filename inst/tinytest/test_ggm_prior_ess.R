# Test the prior effective sample size (ESS) - now for GGM, in the future for other families too.
p <- 3L
G <- matrix(0, p, p); G[1, 2] <- G[2, 1] <- 1
K <- diag(p); K[1, 2] <- K[2, 1] <- 0.3
nu <- 10L
ep <- elicit_prior(ggm_parameters(K, G, nu = nu))

# --- scalar estimators: exact, no simulation 
e <- prior_ess(ep, estimator = c("PT", "ELIR", "MTM"))
expect_inherits(e, "ggm_prior_ess")
expect_inherits(e, "bgm_prior_ess")
expect_equal(e$estimates$PT$global,   nu - p - 1)
expect_equal(e$estimates$ELIR$global, nu - p - 1)
expect_equal(e$estimates$MTM$global,  nu - 1)
expect_equal(e$nu, nu)
expect_equal(e$p, p)
expect_equal(e$prior, "gwishart")

# order follows the request
expect_equal(names(e$estimates), c("PT", "ELIR", "MTM"))
e2 <- prior_ess(ep, estimator = c("MTM", "PT"))
expect_equal(names(e2$estimates), c("MTM", "PT"))

expect_stdout(print(e), "prior_ess")
expect_stdout(summary(e), "Prior effective sample size")

expect_error(prior_ess(ep, estimator = "NOPE"))

# -- Wishart: closed form, matrix estimators are cheap 
Gc <- matrix(1, p, p); diag(Gc) <- 0
Kc <- diag(p) + 0.2; diag(Kc) <- 1.5
epc <- elicit_prior(ggm_parameters(Kc, Gc, nu = nu))
ec <- prior_ess(epc)                       # all five
expect_equal(names(ec$estimates), c("VR", "PR", "MTM", "PT", "ELIR"))
expect_true(is.finite(ec$estimates$VR$global))
expect_true(is.finite(ec$estimates$PR$global))

# --- G-Wishart: Monte Carlo, keep n_samples tiny 
set.seed(1)
eg <- prior_ess(ep, estimator = "VR", n_samples = 25L, burnin = 5L,
                sampler = "gibbs")
expect_true(is.finite(eg$estimates$VR$global))

set.seed(1); a <- prior_ess(ep, estimator = "VR", n_samples = 25L, sampler = "gibbs", burnin = 5L)
set.seed(1); b <- prior_ess(ep, estimator = "VR", n_samples = 25L, sampler = "gibbs", burnin = 5L)
expect_equal(a$estimates$VR$global, b$estimates$VR$global)   # set.seed reproducible

# -- init validation (pure R, instant) 
bad_init <- diag(p); bad_init[1, 3] <- bad_init[3, 1] <- 0.1   # nonzero at non-edge
expect_error(prior_ess(ep, estimator = "VR", init = bad_init), "no edge")

expect_error(prior_ess(ep, estimator = "VR", init = matrix(0, p, p)),
             "positive definite")
expect_error(prior_ess(ep, estimator = "VR", init = diag(2)))  # wrong size

asym <- diag(p); asym[1, 2] <- 0.1        # not symmetric, PD-ish, right size
expect_error(prior_ess(ep, estimator = "VR", init = asym), "must be symmetric")