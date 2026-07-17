# For both DPIR and BFDA, the plan below is tuned to nonconvergence
p <- 3L
G <- matrix(0, p, p); G[1, 2] <- G[2, 1] <- 1
K <- diag(p); K[1, 2] <- K[2, 1] <- 0.3
ep <- elicit_prior(ggm_parameters(K, G, nu = 10))

# threshold no n can reach within max_n
set.seed(1)
expect_warning(
  bad <- design(ep, method = "DPIR", threshold = 1e6, target_probability = 0.99,
                H = 2L, J = 2L, n = c(12L, 20L), n_tol = 10L, max_n = 25L),
  "never reached")
nst <- designbgm:::.design_n_star(bad)
expect_true(is.na(nst$global))
expect_stdout(print(bad), "not reached")           
expect_stdout(summary(bad), "not reached")         
expect_stdout(print(bad), "not all targets converged")

expect_error(validate(bad), "no converged")        

set.seed(1)
badb <- suppressWarnings(
  design(ep, method = "BFDA", threshold = 1e6,
         H = 2L, J = 2L, n = c(12L, 20L), n_tol = 10L, max_n = 25L))
expect_stdout(print(badb), "not reached")   
expect_stdout(summary(badb), "not reached")     
expect_error(validate(badb), "no converged") 
