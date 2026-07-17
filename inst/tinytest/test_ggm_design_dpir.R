# design(method = "DPIR") on a sparse graph, and validate() at both n*.
# The plan below is tuned to converge (low threshold, low target) so the
# converged branches of print/summary are reachable. The "not reached"
# branches are tested in "test_ggm_design_nonconvergence.R".

p <- 3L
G <- matrix(0, p, p); G[1, 2] <- G[2, 1] <- 1
K <- diag(p); K[1, 2] <- K[2, 1] <- 0.3
ep <- elicit_prior(ggm_parameters(K, G, nu = 10))

set.seed(1)
plan <- design(ep, method = "DPIR", threshold = 0.5, target_probability = 0.5,
               H = 5L, J = 5L, n = NULL, n_tol = 25L, max_n = 1000L)

# ---- object shape
expect_inherits(plan, "ggm_design")
expect_inherits(plan, "bgm_design")
expect_equal(plan$method, "DPIR")
expect_equal(plan$family, "ggm")
expect_equal(plan$prior, "gwishart")
expect_false(is.null(plan$ep))              # validate() needs the prior back

nst <- designbgm:::.design_n_star(plan)
cv  <- designbgm:::.design_converged(plan)

expect_true(!is.na(nst$global))
expect_true(!is.na(nst$pw))
expect_true(cv[["global"]])
expect_true(cv[["pw"]])
expect_true(nst$pw >= nst$global)          # weakest parameter needs at least as much

# ---- method dispatch and argument validation
expect_error(design(ep, method = "NOPE"))

# threshold defaults to 1.0 for DPIR when not supplied
set.seed(1)
pdef <- suppressWarnings(
  design(ep, method = "DPIR", H = 2L, J = 2L, n = c(20L, 40L),
         n_tol = 20L, max_n = 60L))
expect_equal(pdef$call_info$threshold, 1.0)

# ---- converged print / summary branches
expect_stdout(print(plan), "DPIR")
expect_stdout(print(plan), "global : n\\* =")
expect_stdout(print(plan), "weakest parameter : n\\* =")
expect_stdout(summary(plan), "Sample-size plan")
expect_stdout(summary(plan), "Global determinant-ratio target")
expect_stdout(summary(plan), "Parameterwise target")
expect_stdout(summary(plan), "converged: TRUE")
expect_stdout(summary(plan), "parameterwise Pr at n\\*")
expect_stdout(summary(plan), "target = 0.50")       # call_info settings line
expect_stdout(summary(plan), "all parameters")

# ---- validate at global n* (default)
set.seed(1)
v <- validate(plan, H = 5L, J = 5L)
expect_inherits(v, "ggm_design_validation")
expect_inherits(v, "bgm_design_validation")
expect_equal(v$method, "DPIR")
expect_null(v$scope)
expect_equal(v$n_star, nst$global)
expect_equal(v$call_info$which_n, "global")
expect_true(is.finite(as.numeric(v$results$global_dpir_prob)[1]))

expect_stdout(print(v), "global Pr\\(DPIR > threshold\\)")
expect_stdout(print(v), "parameterwise Pr at n\\*")
expect_stdout(summary(v), "Design validation")
expect_stdout(summary(v), "parameters individually reach target")   # global branch
expect_stdout(summary(v), "full parameterwise coverage requires")

# ---- validate at parameterwise n*
set.seed(1)
vp <- validate(plan, H = 5L, J = 5L, which_n = "pw")
expect_equal(vp$n_star, nst$pw)
expect_equal(vp$call_info$which_n, "pw")
expect_stdout(summary(vp), "parameters reaching target")   # the non-global branch

# a DPIR plan takes global/pw, not the BFDA values
expect_error(validate(plan, which_n = "h1"))
expect_error(validate(plan, which_n = "h0"))
