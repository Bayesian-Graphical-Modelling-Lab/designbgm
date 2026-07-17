# complete graph -> wishart -> the dense BFDA path
p <- 3L
Gc <- matrix(1, p, p); diag(Gc) <- 0
Kc <- diag(p) * 1.5; Kc[Kc == 0] <- 0.3; diag(Kc) <- 1.5
ep <- elicit_prior(ggm_parameters(Kc, Gc, nu = 10))
expect_equal(ep$prior, "wishart")
expect_false(ep$sparse)

set.seed(1)
plan <- design(ep, method = "BFDA", H = 2L, J = 2L,
               n = NULL, n_tol = 20L, max_n = 1000L)
expect_equal(plan$prior, "wishart")
expect_equal(length(plan$results$edge), 2L)

set.seed(1)
v <- validate(plan, H = 2L, J = 2L)              # -> cpp_design_bfda_edge_dense
expect_equal(v$method, "BFDA")

# DPIR on a complete graph too
set.seed(1)
pd <- suppressWarnings(design(ep, method = "DPIR", H = 2L, J = 2L,
                              n = c(20L, 40L), n_tol = 20L, max_n = 500L))
expect_equal(pd$prior, "wishart")