# design(method = "BSDA") on a sparse graph (G-Wishart prior), and validate().
# BSDA plans a selection criterion (sensitivity / specificity) at a target power,
# and requires posterior inclusion probabilities (pip) and a sparse graph.
# The happy-path control is shrunk hard so the fit is fast. A probit fit on so
# few points may not identify n*, so the identified-dependent checks are guarded:
# the identified branch is checked when it happens, the "not identified" object
# otherwise. If your machine identifies reliably, the guard can be dropped and the
# non-identifying case moved to "test_ggm_design_nonconvergence.R".

p <- 4L
G <- matrix(0, p, p)
G[1, 2] <- G[2, 1] <- 1
G[3, 4] <- G[4, 3] <- 1
K <- diag(p)
K[1, 2] <- K[2, 1] <- 0.4
K[3, 4] <- K[4, 3] <- 0.3

# BSDA needs pip; supply a constant prior inclusion probability
ep <- elicit_prior(ggm_parameters(K, G, nu = 12, pip = 0.5))
expect_true(ep$sparse)
expect_false(is.null(ep$pip))

# --- guard errors (stop() precedes any C++ call)

# non-sparse graph: no absent edges to detect
Kfull <- diag(p); Kfull[upper.tri(Kfull)] <- 0.1; Kfull[lower.tri(Kfull)] <- 0.1
Gfull <- matrix(1, p, p); diag(Gfull) <- 0
ep_dense <- elicit_prior(ggm_parameters(Kfull, Gfull, nu = 12, pip = 0.5))
expect_false(ep_dense$sparse)
expect_error(design(ep_dense, method = "BSDA"), "sparse")

# sparse graph but pip missing
ep_nopip <- elicit_prior(ggm_parameters(K, G, nu = 12))
expect_true(is.null(ep_nopip$pip))
expect_error(design(ep_nopip, method = "BSDA"), "inclusion probabilities")

# --- valid run: tuned-down control for speed
ctrl <- bsda_control(fit_iterations = 300L, fit_burnin = 100L,
                     n_boot = 200L, n_scout = 4L, n_main = 5L,
                     gwish_tol = 1e-01, gwish_iter = 50L, gwish_burnin = 50L,
                     max_iter = 4L, seed = 1L)

set.seed(1)
plan <- design(ep, method = "BSDA",
               measure = "sen", measure_value = 0.8, target_pow = 0.8,
               H = 20L, J = 1L, range_lower = 2L, max_n = 400L,
               control = ctrl)

# --- object shape
expect_inherits(plan, "ggm_design")
expect_inherits(plan, "bgm_design")
expect_equal(plan$method, "BSDA")
expect_equal(plan$family, "ggm")
expect_equal(plan$prior,  "gwishart")
expect_false(is.null(plan$ep))                     # validate() needs the prior back

# --- echoed inputs (call_info)
expect_equal(plan$call_info$measure,       "sen")
expect_equal(plan$call_info$measure_value, 0.8)
expect_equal(plan$call_info$target_pow,    0.8)
expect_inherits(plan$call_info$control, "bsda_control")

# --- results: criterion echo and identifiability flag
r <- plan$results
expect_equal(r$measure,       "sen")
expect_equal(r$measure_value, 0.8)
expect_equal(r$target_pow,    0.8)
expect_true(is.logical(r$identified))

expect_stdout(print(plan),   "BSDA")             # runs for both outcomes
expect_stdout(summary(plan), "Sample-size plan")

# --- identified-dependent checks
if (isTRUE(r$identified)) {
  expect_true(is.finite(r$n_star))
  expect_true(is.finite(r$ci_lo))
  expect_true(is.finite(r$ci_hi))
  expect_true(r$ci_lo <= r$n_star)
  expect_true(r$n_star <= r$ci_hi)
  expect_stdout(print(plan),   "planned sample size")
  expect_stdout(summary(plan), "planned sample size")

  # validate() at the identified n*
  set.seed(1)
  v <- validate(plan, H = 20L, J = 1L)
  expect_inherits(v, "ggm_design_validation")
  expect_inherits(v, "bgm_design_validation")
  expect_equal(v$method, "BSDA")
  expect_equal(v$n_star, as.integer(round(r$n_star)))
  expect_equal(v$results$measure, "sen")
  expect_true(is.finite(as.numeric(v$results$power_at_n$power)[1]))
  expect_stdout(print(v),   "achieved power at n")
  expect_stdout(summary(v), "Design validation")
  expect_stdout(summary(v), "achieved power")
} else {
  # not identified: the object reports the reason and validate() refuses
  expect_stdout(print(plan), "not identified")
  expect_error(validate(plan), "identified")
}

# --- which_n rejects the DPIR/BFDA values
if (isTRUE(r$identified)) {
  expect_error(validate(plan, which_n = "global", H = 20L, J = 1L))   # which_n is a BSDA no-op
}


# bsda_control(): the control list for design(method = "BSDA"). Pure-R argument
# validation, so these checks are deterministic and use no compiled code.

ctrl <- bsda_control()

# --- object shape and defaults
expect_inherits(ctrl, "bsda_control")
expect_true(is.list(ctrl))
expect_equal(ctrl$gwish_sampler,  "direct")
expect_equal(ctrl$edge_threshold, 0.5)
expect_equal(ctrl$fit_iterations, 10000L)
expect_equal(ctrl$fit_burnin,     5000L)
expect_equal(ctrl$alpha,          0.05)
expect_equal(ctrl$n_scout,        6L)
expect_equal(ctrl$n_main,         10L)
expect_equal(ctrl$max_iter,       10L)
expect_equal(ctrl$n_boot,         5000L)
expect_false(ctrl$verbose)
expect_true(is.null(ctrl$seed))            # NULL default

# --- gwish_sampler
expect_equal(bsda_control(gwish_sampler = "block")$gwish_sampler, "block")
expect_error(bsda_control(gwish_sampler = "nonsensesampler"))   # matched against {direct, block}

# --- alpha in (0, 1), open at both ends
expect_error(bsda_control(alpha = 0))
expect_error(bsda_control(alpha = 1))
expect_error(bsda_control(alpha = 1.5))
expect_equal(bsda_control(alpha = 0.01)$alpha, 0.01)

# --- edge_threshold in [0, 1], boundaries allowed
expect_error(bsda_control(edge_threshold = -0.1))
expect_error(bsda_control(edge_threshold = 1.1))
expect_equal(bsda_control(edge_threshold = 0)$edge_threshold, 0)
expect_equal(bsda_control(edge_threshold = 1)$edge_threshold, 1)

# --- scout_frac in (0, 1], upper boundary allowed
expect_error(bsda_control(scout_frac = 0))
expect_error(bsda_control(scout_frac = 1.1))
expect_equal(bsda_control(scout_frac = 1)$scout_frac, 1)

# --- count bounds
expect_error(bsda_control(gwish_iter     = 0))    # >= 1
expect_error(bsda_control(fit_iterations = 0))    # >= 1
expect_error(bsda_control(fit_burnin     = -1))   # >= 0
expect_error(bsda_control(n_boot         = 0))    # >= 1
expect_equal(bsda_control(fit_burnin = 0)$fit_burnin, 0)   # zero burn-in allowed

# --- seed passthrough
expect_equal(bsda_control(seed = 42L)$seed, 42L)