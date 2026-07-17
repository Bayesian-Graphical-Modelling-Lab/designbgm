###############################################################################
#                          Parameters (family hub)                            #
###############################################################################

#' Construct a bgm_parameters object (called by <family>_parameters())
#' @param family Character scalar: "ggm", "omrf", "ising", "gvar" (only ggm is currently supported).
#' @param elements Named list of family-specific parameters, e.g. \code{list(nu = nu, K = K, G = G)} 
#'                 for a GGM, or \code{list(thresholds = thr, interactions = int)} for Ising/OMRF.
#' @return An object of class \code{c("<family>_parameters", "bgm_parameters")}.
#' @noRd
new_bgm_parameters <- function(family, elements) {
  stopifnot(is.character(family), length(family) == 1L, is.list(elements), !"family" %in% names(elements))
  return(structure(
    c(list(family = family), elements), class = c(paste0(family, "_parameters"), "bgm_parameters")
  ))
}

#' @export
print.bgm_parameters <- function(x, ...) {
  cat("<bgm_parameters>  family:", x$family, "\n")
  cat("  no print method for this family. See str(x) for details.\n")
  invisible(x)
}

###############################################################################
#                            Prior elicitation                                #
###############################################################################

#' @title Elicit the informative prior
#' @description Converts a \code{bgm_parameters} object into the elicited prior in the
#' form that the family's samplers and estimators expect. The elicitation strategy depends on 
#' the class of \code{params}. See the method pages listed under Details.
#' @param params A \code{bgm_parameters} object.
#' @param ... Family-specific arguments.
#' @return A \code{bgm_elicited} object describing the elicited prior. It also
#'   carries a family-specific class (for example \code{ggm_elicited}), on which
#'   the downstream methods dispatch.
#' @details
#' Methods are available for the following classes:
#' \describe{
#'  \item{\code{\link[=elicit_prior.ggm_parameters]{ggm_parameters}}}{Gaussian graphical models}
#' }
#' @family prior elicitation
#' @export
elicit_prior <- function(params, ...) {
  UseMethod("elicit_prior")
}

#' @export
elicit_prior.default <- function(params, ...) {
  stop("No elicit_prior() method for class ", paste(class(params), collapse = "/"), ".", call. = FALSE)
}

#' Construct a bgm_elicited object (called by elicit_prior.<family>_parameters()) 
#' @param family Character scalar: "ggm", "omrf", "ising", "gvar" (only ggm is currently supported).
#' @param elements Named list of family-specific elicited quantities.
#' @return An object of class \code{c("<family>_elicited", "bgm_elicited")}.
#' @noRd
new_bgm_elicited <- function(family, elements) {
  stopifnot(is.character(family), length(family) == 1L, is.list(elements))
  obj <- structure(
    c(list(family = family), elements), class = c(paste0(family, "_elicited"), "bgm_elicited")
  )
  return(obj)
}

#' @export
print.bgm_elicited <- function(x, ...) {
  cat("<bgm_elicited>  family:", x$family, "\n")
  cat("  no print method for this family. See str(x) for details.\n")
  invisible(x)
}

###############################################################################
#                        Prior effective sample size                          #
###############################################################################

#' @title Prior effective sample size
#' @description Quantifies the prior effective sample size (ESS) of an elicited prior 
#' described by a \code{bgm_parameters} object. The estimators available
#' depend on the class of \code{params}. See the method pages listed under Details.
#' @param params A \code{bgm_elicited} object, as returned by
#'   \code{\link{elicit_prior}}.
#' @param estimator Character vector selecting the estimators, or \code{NULL}
#'   (default) to compute all available estimators for the family/prior. 
#'   Several estimators are available \code{"VR"}, \code{"PR"}, \code{"MTM"}, 
#'   \code{"PT"}, and \code{"ELIR"}.
#' @param ... Family-specific arguments passed to the method (e.g., \code{n_samples},
#'   the Monte Carlo sample size for simulation-based estimators).
#' @return A \code{prior_ess} object. Its \code{estimates} component is a named
#'   numeric vector with one entry per each specified estimator. \code{print} and \code{summary}
#'   methods ara available.
#' @details
#' Methods are available for the following classes:
#' \describe{
#'  \item{\code{\link[=prior_ess.ggm_elicited]{ggm_elicited}}}{Gaussian graphical models}
#' }
#'
#' When \code{estimator} is \code{NULL} all available estimators are computed
#' together using an optimized routine that shares intermediate quantities and
#' avoids recomputation.
#' @family prior ess
#' @export
prior_ess <- function(params, estimator = NULL, ...) {
  UseMethod("prior_ess")
}

#' @export
prior_ess.default <- function(params, estimator = NULL, ...) {
  stop("No method for class ", 
  paste(class(params), collapse = "/"), ". Elicit a prior first: elicit_prior(ggm_parameters(...)).", 
  call. = FALSE)
}

#' Construct a prior ESS object (called by prior_ess.<family>_elicited())
#' @param family Character scalar: "ggm", "omrf", "ising", "gvar" (only ggm is
#'   currently supported).
#' @param estimates Named list of estimator results.
#' @param info Sampler/diagnostic information, or NULL.
#' @param prior Character scalar, e.g. "gwishart".
#' @param nu Prior study size.
#' @param p Number of nodes.
#' @param requested The estimator argument as supplied by the user, or NULL.
#' @return An object of class \code{c("<family>_prior_ess", "bgm_prior_ess")}.
#' @noRd
new_bgm_prior_ess <- function(family, estimates, info = NULL, prior, nu, p,
                              requested = NULL) {
  stopifnot(is.character(family), length(family) == 1L, is.list(estimates))
  structure(
    list(family = family, prior = prior, nu = nu, p = p,
         estimates = estimates, info = info, requested = requested),
    class = c(paste0(family, "_prior_ess"), "bgm_prior_ess")
  )
}


#' extractor: pulls the scalar global ESS per estimator 
#' @noRd
.prior_ess_global <- function(x) {
  vals <- vapply(x$estimates, function(e) as.numeric(e$global), numeric(1))
  names(vals) <- names(x$estimates)
  return(vals)
}

# Print method for prior_ess objects
#' @export
print.bgm_prior_ess <- function(x, ...) {
  cat("<prior_ess>  family:", x$family, " prior:", x$prior, "\n")
  g <- .prior_ess_global(x)
  for (nm in names(g)) {
    cat(sprintf("  %-5s : %g\n", nm, g[[nm]]))
  }
  return(invisible(x))
}

# Summary method for prior_ess objects
#' @export
summary.bgm_prior_ess <- function(object, ...) { # [note: display graph density and condition numbers (if compute_cond = TRUE)]
  cat("Prior effective sample size\n")
  cat("---------------------------\n")
  cat("  family :", object$family, "\n")
  cat("  prior  :", object$prior, "\n")
  cat("  p      :", object$p, "\n")
  cat("  nu     :", object$nu, "  (prior study size)\n")
  g <- .prior_ess_global(object)
  cat("  estimators run:", paste(names(g), collapse = ", "), "\n\n")
  for (nm in names(g)) {
    cat(sprintf("  %-5s global ESS : %g\n", nm, g[[nm]]))
  }
  return(invisible(object))
}

###############################################################################
#                            Sample-size planning                             #
###############################################################################

# Sample-size planning and validation: recommend n* for a prospective study   
# with prior elicited from  a previous study.
# Methods available: "DPIR", "BFDA".  ("BSDA" will be available in  a later version.)

#' @title Sample-size planning
#' @description Recommends a sample size \code{n*} for a prospective study, given
#' an elicited prior described by a \code{bgm_elicited} object. The planning
#' criterion depends on \code{method}, and the methods available depend on the
#' class of \code{params}. See the method pages listed under Details.
#' @param params A \code{bgm_elicited} object, as returned by
#'   \code{\link{elicit_prior}}.
#' @param method Which planning method to use, \code{"DPIR"} (default) or
#'   \code{"BFDA"}. One method runs per call.
#' @param ... Method-specific arguments passed to the method (e.g.,
#'   \code{max_n}, the largest sample size considered in the planning).
#' @return A \code{bgm_design} object. It records the recommended sample sizes (
#'   two per call, depending on \code{method}) together with the criterion
#'   values and the elicited prior the planning was built from. \code{print} and
#'   \code{summary} methods are available.
#' @details
#' Both methods look for the smallest sample size at which a criterion reaches a
#' target threshold with a specified probability. They differ in the criterion:
#' \code{"DPIR"} uses a data-to-prior information ratio over the model
#' parameters, \code{"BFDA"} the Bayes factor at a representative edge.
#' 
#' Methods are available for the following classes:
#' \describe{
#'  \item{\code{\link[=design.ggm_elicited]{ggm_elicited}}}{Gaussian graphical models}
#' }
#' @family sample size planning
#' @export
design <- function(params, method = c("DPIR", "BFDA"), ...) {
 UseMethod("design")
}

#' @export
design.default <- function(params, method = c("DPIR", "BFDA"), ...) {
  stop("No design() method for class ", 
  paste(class(params), collapse = "/"),". design() expects an elicited prior; see elicit_prior().", 
  call. = FALSE)
}

#' Construct a design plan object (called by design.<family>_elicited()) 
#' @param family Character scalar, e.g. "ggm" family.
#' @param prior Character scalar, e.g. "gwishart" prior.
#' @param method Character scalar referring to the selected planning method, 
#'    either \code{"DPIR"} or \code{"BFDA"}.
#' @param results List of method-specific results.
#' @param call_info List of call information.
#' @param ep Optional \code{bgm_elicited} object, the elicited prior used in the planning.
#' @return An object of class \code{c("<family>_design", "bgm_design")}
#' @noRd
new_bgm_design <- function(family, prior, method, results, call_info = list(), ep = NULL) {
  method <- match.arg(method, c("DPIR", "BFDA"))
  structure(
    list(family = family, prior = prior, method = method,
         results = results, call_info = call_info, ep = ep),
    class = c(paste0(family, "_design"), "bgm_design")
  )
}

# Shared method dependent helpers (internal; these stay family-agnostic) 

#' @noRd
.design_n_star <- function(x) {
  r <- x$results
  if (x$method == "DPIR") {
    list(global = if (!is.null(r$n_star_global)) as.integer(r$n_star_global) else NA_integer_,
         pw     = if (!is.null(r$n_star_pw))     as.integer(r$n_star_pw)     else NA_integer_)
  } else { # BFDA
    list(h0 = if (!is.null(r$n_star_power_h0)) as.integer(r$n_star_power_h0) else NA_integer_,
         h1 = if (!is.null(r$n_star_power_h1)) as.integer(r$n_star_power_h1) else NA_integer_)
  }
}

#' @noRd
.design_converged <- function(x) {
  r <- x$results
  if (x$method == "DPIR") c(global = isTRUE(r$converged_global), pw = isTRUE(r$converged_pw))
  else                    c(h0 = isTRUE(r$converged_h0),         h1 = isTRUE(r$converged_h1))
}

###############################################################################
#                             Plan validation                                 #
###############################################################################

#' @title Validate a sample-size plan
#' @description Checks a sample-size plan by simulation: draws studies at the
#' recommended \code{n} and reports how often the criterion reaches its target
#' there. The checks available depend on the class of \code{plan} and on the 
#' planning method used. See the method pages listed under Details.
#' @param plan A \code{bgm_design} object, as returned by \code{\link{design}}.
#' @param ... Method-specific arguments passed to the method (e.g.,
#'   \code{n_sim}, the number of simulated studies).
#' @return 
#'  A \code{bgm_design_validation} object. \code{print} and 
#'  \code{summary} methods are available.
#' @details
#' Methods are available for the following classes:
#' \describe{
#'  \item{\code{\link[=validate.ggm_design]{ggm_design}}}{Gaussian graphical models}
#' }
#' @family sample size planning
#' @export
validate <- function(plan, ...) {
  UseMethod("validate")
}

#' @export
validate.default <- function(plan, ...){
  stop("validate() expects a design plan from design().", call. = FALSE)
}

#' Construct a design validation object 
#' @param family Character scalar, e.g. "ggm" family.
#' @param method Character scalar referring to the selected planning method, 
#'    either \code{"DPIR"} or \code{"BFDA"}.
#' @param scope BFDA only: \code{"planning_edge"} (default) or \code{"all_edges"} 
#'    (the guarantee check over all present edges). Ignored for DPIR.
#' @param prior Character scalar, e.g. "gwishart" prior.
#' @param n_star The recommended sample size from the plan.
#' @param results List of method-specific results.
#' @param call_info List of call information.
#' @return An object of class \code{c("<family>_design_validation", "bgm_design_validation")}
#' @noRd
new_bgm_design_validation <- function(family, method, scope, prior,
                                      n_star, results, call_info = list()) {
  structure(
    list(family = family, method = method, scope = scope, prior = prior,
         n_star = n_star, results = results, call_info = call_info),
    class = c(paste0(family, "_design_validation"), "bgm_design_validation")
  )
}

###############################################################################
#                Prior study simulation (for simulation studies)              #
###############################################################################

#' @title Simulate a prior study
#' @description Simulates prior studies for a model family: data are generated
#' from the family likelihood and the parameters of interest are estimated,
#' generating one \code{bgm_parameters} object per study, as if obtained from a
#' previous study. Elicitation is a separate step. What is generated
#' depends on the class of \code{study}. See the method pages listed under
#' Details.
#' @param study A \code{bgm_study} object, as returned by each family's
#'   constructor, for example \code{\link{ggm_study}}.
#' @param n_studies Number of independent prior studies to simulate. Default 1.
#' @param ... Family-specific arguments passed to the method.
#' @return A list of \code{n_studies} \code{bgm_parameters} objects.
#' @details
#' Methods are available for the following classes:
#' \describe{
#'  \item{\code{\link[=simulate_prior_study.ggm_study]{ggm_study}}}{Gaussian graphical models}
#' }
#' @family prior elicitation
#' @export
simulate_prior_study <- function(study, n_studies = 1L, ...) {
  UseMethod("simulate_prior_study")
}

#' @export
simulate_prior_study.default <- function(study, n_studies = 1L, ...) {
  stop("No simulate_prior_study() method for class ",
       paste(class(study), collapse = "/"),
       ". Build a study with ggm_study(), or another family's *_study().",
       call. = FALSE)
}

#' Construct a base study object
#' @param family Character scalar: "ggm", "omrf", "ising", "gvar" (only "ggm" is currently supported).
#' @param inputs Named list of family-specific generating inputs,
#'   e.g. \code{list(p = p, nu = nu, G = G)} for a GGM.
#' @return An object of class \code{c("<family>_study", "bgm_study")}.
#' @noRd
new_bgm_study <- function(family, inputs) {
  stopifnot(is.character(family), length(family) == 1L, is.list(inputs))
  return(
    structure(list(family = family, inputs = inputs),
                    class = c(paste0(family, "_study"), "bgm_study"))
  )
}