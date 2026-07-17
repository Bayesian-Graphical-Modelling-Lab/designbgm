# design(method = "BFDA") on a sparse graph (G-Wishart prior), and validate()
# at both scopes. Two edges with different partial correlations, so
# rho_quantile can actually be used.
# The complete-graph / Wishart prior is in "test_ggm_design_wishart.R."

p <- 4L
G <- matrix(0, p, p)
G[1, 2] <- G[2, 1] <- 1
G[3, 4] <- G[4, 3] <- 1
K <- diag(p)
K[1, 2] <- K[2, 1] <- 0.4
K[3, 4] <- K[4, 3] <- 0.3
ep <- elicit_prior(ggm_parameters(K, G, nu = 12))
expect_true(ep$sparse)

set.seed(1)
plan <- design(ep, method = "BFDA", H = 5L, J = 5L,
               n = NULL, n_tol = 25L, max_n = 1000L)

# --- object shape 
expect_inherits(plan, "ggm_design")
expect_inherits(plan, "bgm_design")
expect_equal(plan$method, "BFDA")
expect_equal(plan$prior, "gwishart")
expect_false(is.null(plan$ep))
expect_equal(plan$call_info$threshold, 10.0)   # BFDA default
expect_equal(plan$call_info$pow0, 0.8)
expect_equal(plan$call_info$pow1, 0.8)

expect_equal(length(plan$results$edge), 2L)
expect_true(plan$results$edge[1] < plan$results$edge[2])   # C++ requires m < l
expect_true(is.numeric(plan$results$edge_rho))

expect_stdout(print(plan), "planning edge")
expect_stdout(summary(plan), "pow0")
expect_stdout(summary(plan), "pow1")
expect_stdout(summary(plan), "Planning edge")
expect_stdout(summary(plan), "all present edges")   # BFDA "next:" footer

# -- rho_quantile drives edge selection 
# default 0.5 keeps the stronger half and returns its weakest edge,
# 0 keeps every edge and returns the weakest overall.
set.seed(1)
plan_q <- design(ep, method = "BFDA", rho_quantile = 0,
                 H = 2L, J = 2L, n = 20L, n_tol = 20L, max_n = 60L)
expect_equal(length(plan_q$results$edge), 2L)
expect_false(identical(plan$results$edge, plan_q$results$edge))

expect_error(design(ep, method = "BFDA", rho_quantile = 1))
expect_error(design(ep, method = "BFDA", rho_quantile = -0.1))

# BFDA has nothing to plan for without an edge
Ge  <- matrix(0, p, p)
epe <- elicit_prior(ggm_parameters(diag(p), Ge, nu = 12))
expect_error(design(epe, method = "BFDA", H = 2L, J = 2L, n = 20L),
             "at least one edge")

# --- validate: planning_edge (default scope)
set.seed(1)
v <- validate(plan, H = 2L, J = 2L)
expect_inherits(v, "ggm_design_validation")
expect_equal(v$method, "BFDA")
expect_equal(v$scope, "planning_edge")
expect_equal(v$results$which_n, "h1")          # BFDA default
expect_stdout(print(v), "planning edge")
expect_stdout(summary(v), "power_h0")
expect_stdout(summary(v), "power_h1")

set.seed(1)
v0 <- validate(plan, H = 2L, J = 2L, which_n = "h0")
expect_equal(v0$results$which_n, "h0")

# -- validate: all_edges (the guarantee check)
set.seed(1)
va <- validate(plan, H = 2L, J = 2L, scope = "all_edges")
expect_equal(va$scope, "all_edges")
expect_equal(nrow(va$results$per_edge), 2L)      # one row per present edge
expect_true(all(c("m", "l", "power_h0", "power_h1") %in%
                names(va$results$per_edge)))
expect_true(all(va$results$per_edge$m < va$results$per_edge$l))
expect_stdout(print(va), "present edges")
expect_stdout(summary(va), "Guarantee check")

expect_error(validate(plan, scope = "nonsense"))

# A BFDA plan takes h0/h1, not the DPIR values
expect_error(validate(plan, which_n = "pw"))
expect_error(validate(plan, which_n = "global"))
