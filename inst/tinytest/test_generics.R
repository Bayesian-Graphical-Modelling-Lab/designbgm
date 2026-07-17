# Tests for R/generics.R: the family-agnostic base classes, their constructors,
# the fallback print methods, and the .default methods of every generic.

# -----------------------------------------------------------------------------------
# Stand-in objects for a family with no methods
# (omrf / ising may become real families. Rename these in the future if this happens)
# -----------------------------------------------------------------------------------
fake_p <- designbgm:::new_bgm_parameters("omrf", list(thresholds = 1, interactions = 2))
fake_e <- designbgm:::new_bgm_elicited("ising", list(scale = 1))
fake_s <- designbgm:::new_bgm_study("omrf", list(p = 3L, nu = 10L))

# --- new_bgm_parameters 
expect_inherits(fake_p, "omrf_parameters")
expect_inherits(fake_p, "bgm_parameters")
expect_equal(fake_p$family, "omrf")
expect_equal(fake_p$thresholds, 1)

expect_error(designbgm:::new_bgm_parameters("ggm", list(family = 1)))  # 'family' reserved
expect_error(designbgm:::new_bgm_parameters(c("a", "b"), list()))      # not scalar
expect_error(designbgm:::new_bgm_parameters("ggm", "not a list"))

# ---- new_bgm_elicited 
expect_inherits(fake_e, "ising_elicited")
expect_inherits(fake_e, "bgm_elicited")
expect_equal(fake_e$family, "ising")

expect_error(designbgm:::new_bgm_elicited(c("a", "b"), list()))
expect_error(designbgm:::new_bgm_elicited("ggm", "not a list"))

# --- new_bgm_study 
expect_inherits(fake_s, "omrf_study")
expect_inherits(fake_s, "bgm_study")
expect_equal(fake_s$family, "omrf")
expect_equal(fake_s$inputs$p, 3L)

expect_error(designbgm:::new_bgm_study("ggm", "not a list"))
expect_error(designbgm:::new_bgm_study(c("a", "b"), list()))

# ---- new_bgm_prior_ess 
pe <- designbgm:::new_bgm_prior_ess(
  family = "omrf", estimates = list(MTM = list(global = 9)),
  prior = "someprior", nu = 10, p = 3)
expect_inherits(pe, "omrf_prior_ess")
expect_inherits(pe, "bgm_prior_ess")
expect_equal(designbgm:::.prior_ess_global(pe), c(MTM = 9))

expect_error(designbgm:::new_bgm_prior_ess("omrf", "not a list",
                                           prior = "p", nu = 1, p = 1))

# --- new_bgm_design 
d <- designbgm:::new_bgm_design(family = "omrf", prior = "someprior",
                                method = "DPIR", results = list())
expect_inherits(d, "omrf_design")
expect_inherits(d, "bgm_design")
expect_equal(d$method, "DPIR")
expect_error(designbgm:::new_bgm_design("omrf", "p", method = "NOPE",
                                        results = list()))

# ---- new_bgm_design_validation 
dv <- designbgm:::new_bgm_design_validation(
  family = "omrf", method = "DPIR", scope = NULL, prior = "someprior",
  n_star = 100L, results = list())
expect_inherits(dv, "omrf_design_validation")
expect_inherits(dv, "bgm_design_validation")
expect_equal(dv$n_star, 100L)

# ---------------------------------------------------------------------------
# Fallback print methods: used by any family without its own print method
# ---------------------------------------------------------------------------
expect_stdout(print(fake_p), "bgm_parameters")
expect_stdout(print(fake_p), "omrf")
expect_stdout(print(fake_p), "no print method")

expect_stdout(print(fake_e), "bgm_elicited")
expect_stdout(print(fake_e), "ising")
expect_stdout(print(fake_e), "no print method")

# print returns its argument invisibly (so we use withVisible() to test)
expect_identical(withVisible(print(fake_p))$value, fake_p)
expect_false(withVisible(print(fake_p))$visible)

# ---------------------------------------------------------------------------
# .default methods: every generic must fail on an unknown class
# ---------------------------------------------------------------------------
expect_error(elicit_prior(fake_p), "No elicit_prior")
expect_error(elicit_prior(list(a = 1)), "No elicit_prior")

expect_error(prior_ess(fake_p), "No method for class")
expect_error(prior_ess(fake_p), "Elicit a prior first")

expect_error(design(fake_p), "No design")
expect_error(design(fake_e), "No design")          # elicited, but wrong family

expect_error(validate(fake_p), "expects a design plan")
expect_error(validate(list(a = 1)), "expects a design plan")

expect_error(simulate_prior_study(fake_s), "No simulate_prior_study")
expect_error(simulate_prior_study(list(a = 1)), "No simulate_prior_study")
