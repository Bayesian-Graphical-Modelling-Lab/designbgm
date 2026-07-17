# ESS estimators for the GGM family
.GGM_ESS_ESTIMATORS <- c("VR", "PR", "MTM", "PT", "ELIR")

# Prior distributions for the GGM family.
#  -> regime   : indicates which graph the prior is valid for ("complete", "sparse", "both").
# Regime "both" is not used yet.
.GGM_PRIOR_DISTRIBUTIONS <- list(
  wishart  = list(regime = "complete"),
  gwishart = list(regime = "sparse")
)

###############################################################################
#                             GGM parameters (user input)                     #
###############################################################################

#' @title Construct GGM prior parameters
#' @description Takes the inputs for a Gaussian graphical model and bundles them into
#' one object, ready for prior elicitation with \code{\link{elicit_prior}}.
#' @param K A \code{p} by \code{p} symmetric positive-definite precision matrix, .
#' @param G A \code{p} by \code{p} symmetric 0/1 adjacency matrix with zero diagonal.
#' @param nu Prior study size (degrees of freedom).
#' @param mu Mean vector of length \code{p}. The default \code{NULL} corresponds
#'   to the centered (zero-mean) case that all current methods assume. Supplying
#'   a non-\code{NULL} \code{mu} is not yet supported.
#' @return A \code{ggm_parameters} object, which also inherits from \code{bgm_parameters}.
#' @export
ggm_parameters <- function(K, G, nu, mu = NULL) {
  K <- unname(as.matrix(K))
  G <- unname(as.matrix(G))
  p <- nrow(K)
  stopifnot(
    is.numeric(K), isSymmetric(K), ncol(K) == p,
    is.numeric(G), isSymmetric(G), nrow(G) == p, ncol(G) == p,
    all(G %in% c(0, 1)), all(diag(G) == 0),
    is.numeric(nu), length(nu) == 1L, nu > p # add here additional conditions on nu and p relationship if needed
  )

  # check positive definiteness of K 
  if (any(eigen(K, symmetric = TRUE, only.values = TRUE)$values <= 0))
    stop("'K' must be positive definite.", call. = FALSE)

  # check that K is compatible with G: K must be zero at non-edges (off-diagonal)
  Rho <- abs(cpp_precision_to_partial_correlations(K))
  off <- G == 0 & !diag(p)
  if (any(off)) {
    worst <- max(Rho[off])
    if (worst > 1e-3) {
      ij <- sort(which(off & Rho == worst, arr.ind = TRUE)[1, ])
      stop(sprintf("'K' implies a partial correlation of %.3g at (%d, %d), where 'G' has no edge. ",
                  worst, ij[1], ij[2]),
          "A precision matrix for this graph must be zero at non-edges; ",
          "fit 'K' under the constraint imposed by 'G'.", call. = FALSE)
    }
  }

  # mu must be NULL (centered case) for now. Future work will support non-centered mean modeling (using Normal-Wishart as conjugate prior, and other non-conjugate priors)
  if (!is.null(mu)) {
    stop("Mean modeling is not yet supported; leave 'mu' as NULL (centered case).",
         call. = FALSE)
  }

  # nu must be an integer (degrees of freedom = prior study size)
  if (abs(nu - round(nu)) > .Machine$double.eps^0.5) {
    stop("'nu' must be an integer (prior study size / degrees of freedom).", call. = FALSE)
  }
  nu <- as.integer(round(nu))
  return(new_bgm_parameters("ggm", list(K = K, G = G, nu = nu, mu = mu))) # family name + list of family elements
}

#' @export
print.ggm_parameters <- function(x, ...) {
  p <- nrow(x$G)
  cat("<ggm_parameters>\n")
  cat("  nodes :", p, "\n")
  cat("  edges :", sum(x$G[upper.tri(x$G)]), "of", choose(p, 2), "\n")
  cat("  nu    :", x$nu, " (prior study size)\n")
  invisible(x)
}

###############################################################################
#                             GGM prior elicitation                           #
###############################################################################

#' @title Elicit a prior for GGM parameters
#' @inheritParams elicit_prior
#' @param params A \code{ggm_parameters} object, as returned by \code{\link{ggm_parameters}}. 
#'               It also inherits from \code{bgm_parameters}.
#' @param prior Character scalar naming the precision prior, or \code{NULL}
#'   (default) to infer it from the graph: \code{"wishart"} for a complete
#'   \code{G}, \code{"gwishart"} for a sparse \code{G}. Supplying \code{prior}
#'   overrides the default and is validated against the graph (e.g. a sparse-only
#'   prior cannot be placed on a complete graph).
#' @inherit elicit_prior return 
#' @family prior elicitation
#' @seealso \code{elicit_prior()} for the generic method and \code{prior_ess()} 
#'          for prior ESS estimation.
#' @details The elicitation currently applies the \code{K/nu} centering, 
#'          which is exact for Wishart, and approximate for G-Wishart. 
#' @export
elicit_prior.ggm_parameters <- function(params, prior = NULL, ...) {
  K  <- params$K
  G  <- params$G
  nu <- params$nu
  p  <- nrow(G)
  sparse <- .graph_is_sparse(G)

  # check the prior: infer from structure if unspecified, validate otherwise
  if (is.null(prior)) {
    prior <- if (sparse) "gwishart" else "wishart"
  } else {
    prior <- match.arg(prior, names(.GGM_PRIOR_DISTRIBUTIONS))
  }

  # regime compatibility: the chosen prior must be valid for this graph
  regime <- .GGM_PRIOR_DISTRIBUTIONS[[prior]]$regime
  ok <- (sparse  && regime %in% c("sparse", "both")) || (!sparse && regime %in% c("complete", "both"))
  if (!ok) {
    stop(sprintf("Prior '%s' is not used for a %s graph; use '%s' instead.", prior, if (sparse) "sparse" else "complete", if (sparse) "gwishart" else "wishart"), call. = FALSE)
  }

  return(new_bgm_elicited("ggm", list(
    scale  = K / nu,   # K / nu centering: prior expectation nu * scale = K (exact for Wishart, approximate for G-Wishart)
    nu     = nu,
    G      = G,
    p      = p,
    sparse = sparse,
    prior  = prior     
  )))
}

#' @export
print.ggm_elicited <- function(x, ...) {
  cat("<ggm_elicited>  prior:", x$prior, "\n")
  cat("  nodes :", x$p, "\n")
  cat("  edges :", sum(x$G[upper.tri(x$G)]), "of", choose(x$p, 2), "\n")
  cat("  nu    :", x$nu, " (prior study size)\n")
  invisible(x)
}

###############################################################################
#                           GGM prior ESS estimators                          #
###############################################################################

# Scalar estimators: closed-form, computed in R.
#   PT, ELIR : nu - p - 1
#   MTM      : nu - 1
# Returned with the same list shape
.ess_ggm_scalar <- function(e, nu, p) {
  global <- switch(e,
    PT   = nu - p - 1,
    ELIR = nu - p - 1,
    MTM  = nu - 1
  )
  list(global = global, parameterwise = NULL,
       cond_numerator = NULL, cond_denominator = NULL)
}

# Dispatch the requested GGM estimators: scalars estimators in R, VR and PR via C++ routine.
.ess_ggm <- function(ep, estimator, aggregation, sampler, n_samples, tol, itermax, burnin, init, compute_cond) {

  requested <- if (is.null(estimator)) .GGM_ESS_ESTIMATORS
               else match.arg(estimator, .GGM_ESS_ESTIMATORS, several.ok = TRUE)

  matrix_est <- intersect(requested, c("VR", "PR"))
  scalar_est <- intersect(requested, c("PT", "ELIR", "MTM"))

  results <- list()
  info    <- NULL

  # check init matrix if provided
  if (!is.null(init)) {
    stopifnot(is.matrix(init), nrow(init) == ncol(init), nrow(init) == nrow(ep$scale))
    if (!isSymmetric(unname(init)))
      stop("'init' must be symmetric.", call. = FALSE)
    if (any(init[ep$G == 0 & !diag(nrow(ep$G))] != 0))
      stop("'init' must be zero wherever `G` has no edge. The sampler never updates those entries.", call. = FALSE)
    if (any(eigen(init, symmetric = TRUE, only.values = TRUE)$values <= 0))
      stop("'init' must be positive definite.", call. = FALSE)
  }

  # matrix estimators: C++ call with the right "which" argument
  if (length(matrix_est) > 0L) {
    which <- if (length(matrix_est) == 2L) "both" else matrix_est  # "VR" / "PR" / "both"
    res <- cpp_ggm_prior_ess_vr_pr(
      prior       = ep$prior,
      which       = which,
      nu          = ep$nu,
      K           = ep$scale,
      G           = ep$G,
      aggregation = aggregation,
      sampler     = sampler,
      n_samples   = n_samples,
      tol         = tol,
      itermax     = itermax,
      burnin      = burnin,
      init        = init,
      compute_cond = compute_cond
    )
    info <- res$info
    for (e in matrix_est) results[[e]] <- res[[e]]
  }

  # scalar estimators: simple arithmetic (in R)
  for (e in scalar_est) results[[e]] <- .ess_ggm_scalar(e, ep$nu, ep$p)

  # order results to match the input order
  results <- results[requested]

  list(estimates = results, info = info)
}


#' @title Prior effective sample size for a Gaussian graphical model 
#' @description Computes the prior ESS of an elicited prior distribution for the 
#'   parameters in a Gaussian graphical model. The matrix estimators (VR, PR) work 
#'   with matrices  (Fisher information matrix and prior variance matrix) and need to
#'   be reduced to a scalar first, the remaining estimators do not. For the
#'   G-Wishart prior the quantities involved have no closed form and are
#'   approximated by Monte Carlo simulations. 
#' @inheritParams prior_ess
#' @param params A \code{ggm_elicited} object, as returned by
#'   \code{\link{elicit_prior}}.
#' @param aggregation Character scalar selecting the aggregation method to reduce the 
#'   estimator to a scalar: \code{"det"} (default), \code{"tr"}, or \code{"mean"}.
#'   Applies only to the matrix estimators (VR, PR).
#' @param sampler Monte Carlo sampler for the G-Wishart prior, either
#'   \code{"direct"} (default) or \code{"gibbs"} (Edge-wise Gibbs sampler). 
#'   Ignored for the Wishart prior and the scalar estimators.
#' @param n_samples Number of Monte Carlo samples for the G-Wishart prior. Default 1000.
#' @param tol Convergence tolerance for the direct sampler. Default 1e-6.
#' @param itermax Number of maximum iterations for the direct sampler. Default 1000. 
#'   This number caps the iterative procedure for one draw of the direct sampler.
#' @param burnin Number of burn-in iterations discarded by the Gibbs sampler. Default 500.
#' @param init Starting value for the Gibbs sampler, a \code{p} by \code{p}
#'   symmetric positive-definite precision matrix on the same scale as the
#'   returned samples, which are centered on the elicited precision matrix.
#'   Entries where the graph has no edge must be zero: the sampler never updates
#'   them, so whatever is defined there is carried into every draw. \code{NULL}
#'   (default) starts from a diagonal matrix of the marginal precisions implied
#'   by the elicited prior (\code{diag(params$nu / diag(solve(params$scale)))}).
#' @param compute_cond Logical: whether to compute the condition number of numerator and 
#'   denominator of VR and PR (default FALSE).
#' @param ... Ignored, present for consistency with the generic \code{prior_ess()}.
#' @return A \code{prior_ess} object for the \code{ggm} family. Its \code{estimates} 
#'   component is a named numeric vector with one entry per estimator. \code{print} and \code{summary}
#'   methods are available.
#' @return A \code{prior_ess} object for the \code{ggm} family. Its
#'   \code{estimates} component is a named list, one element per estimator, each
#'   with a \code{global} value and a \code{parameterwise} table (for the matrix estimators). 
#'   When \code{compute_cond = TRUE} those elements also carry \code{cond_numerator} and 
#'   \code{cond_denominator}. \code{print} and \code{summary} methods are available.
#' @details
#' The sampling arguments (\code{sampler}, \code{n_samples}, \code{tol},
#'   \code{itermax}, \code{burnin}, \code{init}) apply to the G-Wishart prior
#'   only. Under the Wishart prior the quantities involved are available in
#'   closed form and nothing is simulated.
#' @family prior ess 
#' @export
prior_ess.ggm_elicited <- function(params,
                                   estimator    = NULL,
                                   aggregation  = c("det", "tr", "mean"),
                                   sampler      = c("direct", "gibbs"),
                                   n_samples    = 1000L,
                                   tol          = 1e-6,
                                   itermax      = 1000L,
                                   burnin       = 500L,
                                   init         = NULL,
                                   compute_cond = FALSE,
                                   ...) {
  aggregation <- match.arg(aggregation)
  sampler     <- match.arg(sampler)
  ep <- params

  res <- .ess_ggm(
    ep          = ep,
    estimator   = estimator,
    aggregation = aggregation,
    sampler     = sampler,
    n_samples   = n_samples,
    tol         = tol,
    itermax     = itermax,
    burnin      = burnin,
    init        = init,
    compute_cond = compute_cond
  )

  return(new_bgm_prior_ess(
    family    = "ggm",
    estimates = res$estimates,
    info      = res$info,
    prior     = ep$prior,
    nu        = ep$nu,
    p         = ep$p,
    requested = estimator
  ))
}

###############################################################################
#                  GGM sample-size planning (estimation)                      #
###############################################################################

# Planning-edge selection for BFDA (internal) 

# Selects the representative edge for BFDA edge-based planning. Among the present
# edges, takes |rho_ij| values, finds the rho_quantile-th quantile
# as a threshold, keeps edges at or above it, and returns the WEAKEST of those
# kept. Returns 1-based (row,col) upper-triangular indices. The downstream functions 
# will convert to 0-based for the C++ routines.
.bfda_planning_edge <- function(K, G, rho_quantile = 0.5) {
  stopifnot(is.numeric(rho_quantile), 
            length(rho_quantile) == 1L,
            rho_quantile >= 0, 
            rho_quantile < 1)

  p   <- nrow(G)
  Rho <- cpp_precision_to_partial_correlations(K)
  Rho[G == 0] <- 0
  G_upper     <- cpp_G_upper(p = p, G = G)

  all_edges <- which(G_upper == 1L, arr.ind = TRUE)   # 1-based, row < col
  if (nrow(all_edges) == 0L) {
    stop("No present edges in G to plan for (BFDA requires at least one edge).", call. = FALSE)
  }

  abs_rho <- abs(Rho[all_edges])
  rho_min <- as.numeric(stats::quantile(abs_rho, probs = rho_quantile))
  keep    <- abs_rho >= rho_min
  edges_k <- all_edges[keep, , drop = FALSE]
  rho_k   <- abs_rho[keep]
  idx     <- which.min(rho_k)                         # weakest in the upper tail

  m <- edges_k[idx, 1]
  l <- edges_k[idx, 2]
  stopifnot(m < l)                                    # C++ BF requires m < l

  return(list(m = m, l = l, rho = Rho[m, l],
              n_edges_considered = nrow(edges_k)))
}

#' @title Sample-size planning for a Gaussian graphical model
#' @description Recommends a sample size for a prospective Gaussian graphical
#'   model study, given a prior elicited from a previous one. \code{"DPIR"}
#'   returns two sizes, one for the graph as a whole and one at the parameterwise
#'   level (considering only the off-diagonal elements). \code{"BFDA"} returns 
#'   one per hypothesis: under the null (edge absent) and under the alternative (edge present).
#' @param params A \code{ggm_elicited} object, as returned by \code{\link{elicit_prior}}. 
#'               It also inherits from \code{bgm_elicited}.
#' @param method Which planning method to use, \code{"DPIR"} (default) or
#'   \code{"BFDA"}. One method runs per call.
#' @param H,J Outer and inner Monte Carlo replication counts. \code{H} draws
#'   from the prior (default 150), \code{J} datasets per draw (default 20).
#' @param n Candidate sample sizes to search over. \code{NULL} (default) picks a
#'   grid automatically and refines it by bisection.
#' @param threshold Target threshold the criterion must reach: the DPIR ratio for
#'   \code{"DPIR"}, the Bayes factor for \code{"BFDA"}. \code{NULL} (default)
#'   uses 1.0 for DPIR and 10.0 for BFDA.
#' @param target_probability DPIR only: how often the ratio must reach
#'   \code{threshold}. Default 0.95.
#' @param rho_quantile BFDA only: which edge to plan around. Edges are ordered
#'   by the magnitude of their partial correlation and the one at this quantile
#'   is selected. Default 0.5, the median edge. The selected edge is reported
#'   with its signed partial correlation.
#' @param pow0,pow1 BFDA only: how often the Bayes factor must give decisive
#'   evidence for the hypothesis that is true. With \code{BF01} the Bayes factor
#'   for edge absence, \code{pow0} targets \code{Pr(BF01 > threshold | H0)},
#'   correctly excluding the edge, and \code{pow1} targets
#'   \code{Pr(BF01 < 1/threshold | H1)}, correctly detecting the edge. Both default to
#'   0.8.
#' @param n_tol Search tolerance, in observation units. The search stops once the
#'   bracket reaches this width. Default 1, i.e. exact (increasing this value speeds
#'   up the search but may reduce accuracy).
#' @param max_n Largest sample size considered. If the target is not met by
#'   \code{max_n}, the search stops and reports non-convergence. Default 5000.
#' @param ... Ignored, present for consistency with the generic \code{design()}.
#' @return A \code{ggm_design} object, which also inherits from
#'   \code{bgm_design}. Its \code{results} component containes the recommended
#'   sizes: \code{n_star_global} and \code{n_star_pw} for \code{"DPIR"},
#'   \code{n_star_power_h0} and \code{n_star_power_h1} for \code{"BFDA"}, each
#'   with a convergence flag. A size that was not reached within the search is
#'   \code{NA} and its convergence flag is \code{FALSE}. \code{print} and \code{summary}
#'   methods are available.
#' @details
#' Under BFDA the edge is tested with \code{H0} the edge absent against
#' \code{H1} the edge present, and the Bayes factor reported throughout is
#' \code{BF01}, the Bayes factor for absence: \code{BF01 > threshold} excludes
#' the edge, \code{BF01 < 1/threshold} detects it. This is the reciprocal of the
#' inclusion Bayes factor \code{BF10}.
#'
#' Neither error rate is the complement of the corresponding power, since a
#' Bayes factor between \code{1/threshold} and \code{threshold} is inconclusive
#' and counts towards neither.
#'
#' Cost grows with \code{H * J} and with the width of the search: the criterion
#' is evaluated by simulation at every candidate \code{n}.
#' @family sample size planning
#' @export
design.ggm_elicited <- function(params, method = c("DPIR", "BFDA"),
                                H = 150L, J = 20L, n = NULL,
                                threshold = NULL,
                                target_probability = 0.95,
                                rho_quantile = 0.5,
                                pow0 = 0.8, pow1 = 0.8,
                                n_tol = 1L, max_n = 5000L, ...) {
  method <- match.arg(method)
  ep <- params

  K  <- ep$scale
  G  <- ep$G
  nu <- ep$nu

  if (is.null(n)) n <- seq.int(nu, max(nu * 4L, 200L), length.out = 20L)
  n <- as.integer(round(n))

  thr <- if (!is.null(threshold)) threshold else if (method == "DPIR") 1.0 else 10.0 # Default 1.0 for DPIR, 10.0 for BFDA

  results <- switch(method,
    DPIR = {
      prior <- ep$prior
      cpp_design_dpir(
        prior              = prior,
        K                  = K,
        nu                 = nu,
        G                  = G,
        H                  = H,
        J                  = J,
        n                  = n,
        threshold          = thr,
        optimize           = TRUE,
        target_probability = target_probability,
        n_tol              = n_tol,
        max_n              = max_n
      )
    },

    BFDA = {
      edge <- .bfda_planning_edge(K, G, rho_quantile = rho_quantile)
      m0   <- edge$m - 1L                              # 0-based for C++
      l0   <- edge$l - 1L

      res <- if (ep$sparse) {
        # G-Wishart BFDA: needs G_upper, n_vec, Rho
        G_upper <- cpp_G_upper(p = ep$p, G = G)
        n_vec   <- cpp_n_vec(p = ep$p, G_upper = G_upper)
        Rho     <- cpp_precision_to_partial_correlations(K)
        cpp_design_bfda_edge_sparse(
          K           = K, 
          nu          = nu, 
          G           = G, 
          G_upper     = G_upper, 
          n_vec       = n_vec, 
          Rho         = Rho,
          m           = m0, 
          l           = l0, 
          H           = H, 
          J           = J, 
          n           = n,
          pow0        = pow0, 
          pow1        = pow1,
          threshold   = thr, 
          optimize    = TRUE, 
          n_tol       = n_tol, 
          max_n       = max_n
        )
      } else {
        cpp_design_bfda_edge_dense(
          K           = K, 
          nu          = nu, 
          G           = G, 
          m           = m0, 
          l           = l0, 
          H           = H, 
          J           = J, 
          n           = n,
          pow0        = pow0, 
          pow1        = pow1,
          threshold   = thr, 
          optimize    = TRUE, 
          n_tol       = n_tol, 
          max_n       = max_n
        )
      }
      # carry the (1-based) planning edge for print/summary methods
      res$edge     <- c(edge$m, edge$l)
      res$edge_rho <- edge$rho
      res
    }
  )

  call_info <- if (method == "DPIR") {
    list(H = H, J = J, threshold = thr, target_probability = target_probability)
  } else {  # BFDA
    list(H = H, J = J, threshold = thr, pow0 = pow0, pow1 = pow1, rho_quantile = rho_quantile)
  }

  call_info$n_tol <- n_tol
  call_info$max_n <- max_n

  return(new_bgm_design(method = method, family = "ggm",
                    prior = ep$prior, results = results, 
                    call_info = call_info, ep = ep))
}


# Print and summary methods for ggm_design objects

#' @export
print.ggm_design <- function(x, ...) {

  cat("<design>  method:", x$method, " family:", x$family, " prior:", x$prior, "\n")

  r <- x$results 
  nst <- .design_n_star(x) 
  cv <- .design_converged(x)

  if (x$method == "DPIR") {
    cat("  planned sample size (DPIR determinant ratio):\n")
    if (!is.na(nst$global))
      cat(sprintf("    global : n* = %d   (Pr(DPIR > threshold) = %.3f)\n",
                  nst$global, r$prob_global_at_n_star))
    else cat("    global : not reached within the search range\n")
    if (!is.na(nst$pw)) cat(sprintf("    weakest parameter : n* = %d\n", nst$pw))
    else cat("    weakest parameter : not reached within the search range\n")
  } else {
    edge <- r$edge
    if (!is.null(edge)) {
      cat(sprintf("  planning edge: (%d, %d)", edge[1], edge[2]))
      if (!is.null(r$edge_rho)) cat(sprintf("   rho = %.3f", r$edge_rho))
      cat("\n")
    }
    cat("  planned sample size (Bayes factor):\n")
    if (!is.na(nst$h0)) cat(sprintf("    H0 (edge absent)  : n* = %d   (power = %.3f)\n",
                                    nst$h0, r$power_h0_at_n_star))
    else cat("    H0 (edge absent)  : not reached within the search range\n")
    if (!is.na(nst$h1)) cat(sprintf("    H1 (edge present) : n* = %d   (power = %.3f)\n",
                                    nst$h1, r$power_h1_at_n_star))
    else cat("    H1 (edge present) : not reached within the search range\n")
  }
  if (!all(cv)) cat("  note: not all targets converged. Consider a wider search / larger max_n.\n")
  invisible(x)
}

#' @export
summary.ggm_design <- function(object, ...) {
  
  x <- object
  r <- x$results

  cat("Sample-size plan\n----------------\n")
  cat("  method :", x$method, "\n  family :", x$family, "\n  prior  :", x$prior, "\n")

  ci <- x$call_info

  if (length(ci)) {
    common <- c(if (!is.null(ci$H)) sprintf("H = %d", ci$H),
                if (!is.null(ci$J)) sprintf("J = %d", ci$J),
                if (!is.null(ci$threshold)) sprintf("threshold = %g", ci$threshold))
    mbits <- if (x$method == "DPIR") {
      if (!is.null(ci$target_probability)) sprintf("target = %.2f", ci$target_probability)
    } else {
      c(if (!is.null(ci$pow0)) sprintf("pow0 = %.2f", ci$pow0),
        if (!is.null(ci$pow1)) sprintf("pow1 = %.2f", ci$pow1))
    }
    bits <- c(common, mbits)
    if (length(bits)) cat("  settings:", paste(bits, collapse = "  "), "\n")
  }
  cat("\n")
  nst <-  .design_n_star(x)
  cv <-   .design_converged(x)
  if (x$method == "DPIR") {
    cat("  Global determinant-ratio target\n")
    if (!is.na(nst$global))
      cat(sprintf("    n* = %d   Pr(DPIR > threshold) = %.3f   converged: %s\n",
                  nst$global, r$prob_global_at_n_star, cv[["global"]]))
    else cat(sprintf("    not reached (converged: %s)\n", cv[["global"]]))
    cat("  Parameterwise target (weakest off-diagonal parameter)\n")
    if (!is.na(nst$pw)) {
      cat(sprintf("    n* = %d   converged: %s\n", nst$pw, cv[["pw"]]))
      if (!is.null(r$prob_pw_at_n_star)) {
        pw <- as.numeric(r$prob_pw_at_n_star)
        cat(sprintf("    parameterwise Pr at n*: min %.3f, median %.3f, max %.3f (%d params)\n",
                    min(pw), stats::median(pw), max(pw), length(pw)))
      }
    } else cat(sprintf("    not reached (converged: %s)\n", cv[["pw"]]))
  } else {
    edge <- r$edge
    if (!is.null(edge)) cat(sprintf("  Planning edge: (%d, %d)%s\n", edge[1], edge[2],
                              if (!is.null(r$edge_rho)) sprintf("   rho = %.3f", r$edge_rho) else ""))
    cat("  H0: edge absent (power to exclude)\n")
    if (!is.na(nst$h0)) cat(sprintf("    n* = %d   power = %.3f   FPR = %.3f   converged: %s\n",
                              nst$h0, r$power_h0_at_n_star, r$fpr_at_n_star_power_h0, cv[["h0"]]))
    else cat(sprintf("    not reached (converged: %s)\n", cv[["h0"]]))
    cat("  H1: edge present (power to detect)\n")
    if (!is.na(nst$h1)) cat(sprintf("    n* = %d   power = %.3f   FNR = %.3f   converged: %s\n",
                              nst$h1, r$power_h1_at_n_star, r$fnr_at_n_star_power_h1, cv[["h1"]]))
    else cat(sprintf("    not reached (converged: %s)\n", cv[["h1"]]))
  }
  unit <- if (x$method == "BFDA") "all present edges" else "all parameters"
  cat(sprintf("\n  next: validate() refines power at n*; validate(scope='all_edges') checks %s.\n",
              unit))
  invisible(object)
}



###############################################################################
#                  GGM sample-size planning (validation)                      #
###############################################################################

# Helpers for validate.ggm_design() to evaluate the design at the planned n*.

#' @noRd
.design_eval_dpir_ggm <- function(ep, n_eval, H, J, threshold, target_probability, n_tol, max_n) {
  cpp_design_dpir(prior = ep$prior, K = ep$scale, nu = ep$nu, G = ep$G,
              H = H, J = J, n = as.integer(n_eval), threshold = threshold,
              optimize = FALSE, target_probability = target_probability, n_tol = n_tol, max_n = max_n)
}

#' @noRd
.design_eval_bfda_edge_ggm <- function(ep, m0, l0, n_eval, H, J, threshold, pow0, pow1) {
  if (ep$sparse) {
    G_upper <- cpp_G_upper(p = ep$p, G = ep$G)
    n_vec   <- cpp_n_vec(p = ep$p, G_upper = G_upper)
    Rho     <- cpp_precision_to_partial_correlations(ep$scale)
    cpp_design_bfda_edge_sparse(K = ep$scale, nu = ep$nu, G = ep$G, G_upper = G_upper,
                            n_vec = n_vec, Rho = Rho, m = m0, l = l0,
                            H = H, J = J, n = as.integer(n_eval),
                            pow0 = pow0, pow1 = pow1, threshold = threshold, optimize = FALSE)
  } else {
    cpp_design_bfda_edge_dense(K = ep$scale, nu = ep$nu, G = ep$G, m = m0, l = l0,
                           H = H, J = J, n = as.integer(n_eval),
                           pow0 = pow0, pow1 = pow1, threshold = threshold, optimize = FALSE)
  }
}

#' @noRd
.present_edges <- function(G) {
  ut <- which(upper.tri(G) & G == 1, arr.ind = TRUE)
  if (nrow(ut) == 0L) stop("No present edges in G to evaluate.", call. = FALSE)
  data.frame(m1 = ut[, 1], l1 = ut[, 2], m0 = ut[, 1] - 1L, l0 = ut[, 2] - 1L)
}

#' @title Validate a sample-size plan for a Gaussian graphical model
#' @description Checks a design plan by simulating studies at the recommended \code{n}
#'   and reporting how often the criterion actually reaches its target threshold.
#'   Usually, the planning search uses fewer replications than this check, so the two
#'   probabilities will not match exactly.
#' @param plan A \code{ggm_design} object, as returned by \code{\link{design}}.
#'   It also inherits from \code{bgm_design}.
#' @param H,J Outer and inner Monte Carlo replication counts, as in
#'   \code{\link{design}}. Defaults 500 and 100, larger than the planning
#'   defaults, since the check is run at a single \code{n}.
#' @param which_n Which of the plan's recommended sizes to validate at.
#'   \code{NULL} (default) picks \code{"global"} for a DPIR plan and
#'   \code{"h1"} for a BFDA plan. DPIR plans also accept \code{"pw"}, the
#'   parameterwise size; BFDA plans also accept \code{"h0"}.
#' @param scope BFDA only: whether to check the edge the plan was built
#'   around (\code{"planning_edge"}) or every present edge in the graph
#'   (\code{"all_edges"}), the stricter guarantee check. \code{NULL} (default) means
#'   \code{"planning_edge"}. This argument is ignored for DPIR plans.
#' @param ... Ignored, present for consistency with the generic.
#' @return A \code{ggm_design_validation} object, which also inherits from
#'   \code{bgm_design_validation}. Its \code{n_star} component is the size that
#'   was checked and \code{results} holds the probability achieved there.
#'   \code{print} and \code{summary} methods are available.
#' @family sample size planning
#' @export
validate.ggm_design <- function(plan, H = 500L, J = 100L,
                                which_n = NULL, scope = NULL, ...) {
  r <- plan$results 
  info <- plan$call_info 
  ep <- plan$ep

  if (is.null(ep))
    stop("The plan does not carry the elicited prior; design() must store ep.", call. = FALSE)

  # DPIR
  if (plan$method == "DPIR") {
    which_n <- match.arg(which_n %||% "global", c("global", "pw"))
    n_star  <- if (which_n == "global") r$n_star_global else r$n_star_pw
    if (is.null(n_star) || is.na(n_star))
      stop(sprintf("Plan has no converged %s n* to validate.",
                   if (which_n == "global") "global" else "parameterwise"),
           call. = FALSE)

    # one call returns global + parameterwise Pr at n_star (full table regardless)
    res <- .design_eval_dpir_ggm(ep, n_star, H, J,
                                 threshold = info$threshold %||% 1.0,
                                 target_probability = info$target_probability %||% 0.95,
                                 n_tol = info$n_tol %||% 1L, max_n = info$max_n %||% 5000L
                                 )
    return(new_bgm_design_validation(
      family = "ggm", method = "DPIR", scope = NULL, prior = plan$prior,
      n_star = n_star, results = res,
      call_info = list(H = H, J = J, threshold = info$threshold,
                       target_probability = info$target_probability,
                       which_n = which_n,
                       n_star_global = r$n_star_global, n_star_pw = r$n_star_pw)))
  }

  # BFDA
  which_n <- match.arg(which_n %||% "h1", c("h1", "h0"))
  scope   <- match.arg(scope   %||% "planning_edge", c("planning_edge", "all_edges"))
  n_star  <- if (which_n == "h1") r$n_star_power_h1 else r$n_star_power_h0
  if (is.null(n_star) || is.na(n_star))
    stop(sprintf("Plan has no converged n* for %s to validate.", toupper(which_n)),
         call. = FALSE)

  thr <- info$threshold %||% 10.0
  pow0 <- info$pow0 %||% 0.8
  pow1 <- info$pow1 %||% 0.8

  if (scope == "planning_edge") { # single edge: the one the plan was built around
    res <- .design_eval_bfda_edge_ggm(ep, r$edge[1] - 1L, r$edge[2] - 1L, n_star, H, J, thr, pow0, pow1)
    return(new_bgm_design_validation(
      family = "ggm", method = "BFDA", scope = "planning_edge", prior = plan$prior,
      n_star = n_star,
      results = list(planning_edge = res, edge = r$edge, edge_rho = r$edge_rho,
                     which_n = which_n),
      call_info = list(H = H, J = J, threshold = thr, pow0 = pow0, pow1 = pow1,
                       which_n = which_n)))
  }

  # scope == "all_edges": two-sided guarantee over present edges
  edges <- .present_edges(ep$G)
  tab <- do.call(rbind, lapply(seq_len(nrow(edges)), function(e) {
    res_e <- .design_eval_bfda_edge_ggm(ep, edges$m0[e], edges$l0[e],
                                        n_star, H, J, thr, pow0, pow1)
    data.frame(m = edges$m1[e], l = edges$l1[e],
               power_h0 = as.numeric(res_e$power_h0)[1],
               power_h1 = as.numeric(res_e$power_h1)[1],
               fpr_h0   = if (!is.null(res_e$fpr_h0)) as.numeric(res_e$fpr_h0)[1] else NA_real_,
               fnr_h1   = if (!is.null(res_e$fnr_h1)) as.numeric(res_e$fnr_h1)[1] else NA_real_)
  }))
  new_bgm_design_validation(
    family = "ggm", method = "BFDA", scope = "all_edges", prior = plan$prior,
    n_star = n_star,
    results = list(per_edge = tab, edge = r$edge, edge_rho = r$edge_rho,
                   which_n = which_n, pow0 = pow0, pow1 = pow1),
    call_info = list(H = H, J = J, threshold = thr, pow0 = pow0, pow1 = pow1,
                     which_n = which_n))
}

# Print and Summary methods

#' @export
print.ggm_design_validation <- function(x, ...) {
  cat("<design_validation>  method:", x$method,
      "  n* =", x$n_star, "\n")
  r <- x$results
  if (x$method == "DPIR") {
    cat(sprintf("  global Pr(DPIR > threshold) at n* : %.3f\n",
                as.numeric(r$global_dpir_prob)[1]))
    if (!is.null(r$pw_dpir_prob)) {
      pw <- as.numeric(r$pw_dpir_prob[, 1])
      cat(sprintf("  parameterwise Pr at n* : min %.3f / median %.3f / max %.3f\n",
                  min(pw), stats::median(pw), max(pw)))
    }
  } else if (x$scope == "planning_edge") {
    pe <- r$planning_edge
    cat(sprintf("  planning edge (%d, %d): power_h0 = %.3f  power_h1 = %.3f\n",
                r$edge[1], r$edge[2], as.numeric(pe$power_h0)[1], as.numeric(pe$power_h1)[1]))
  } else {
    tab  <- r$per_edge
    p1   <- sum(tab$power_h1 >= r$pow1, na.rm = TRUE)
    p0   <- sum(tab$power_h0 >= r$pow0, na.rm = TRUE)
    cat(sprintf("  present edges at n*: detect (power_h1>=%.2f) %d/%d ; exclude (power_h0>=%.2f) %d/%d\n",
                r$pow1, p1, nrow(tab), r$pow0, p0, nrow(tab))) # here: if there are NA, nrow(tab) should exclude them; in practice tab is always complete (add NA guard here if needed)
  }
  invisible(x)
}

#' @export
summary.ggm_design_validation <- function(object, ...) {
  x <- object
  r <- x$results 
  ci <- x$call_info
  cat("Design validation\n-----------------\n")
  cat("  method :", x$method, "  validated at n* =", x$n_star,
      sprintf("(%s)", ci$which_n %||% x$scope), "\n")
  cat(sprintf("  replication: H = %s, J = %s\n", ci$H %||% "?", ci$J %||% "?")); cat("\n")

  if (x$method == "DPIR") {
    cat(sprintf("  global Pr(DPIR > threshold) at n* : %.3f\n",
                as.numeric(r$global_dpir_prob)[1]))
    if (!is.null(r$pw_dpir_prob)) {
      pw  <- as.numeric(r$pw_dpir_prob[, 1])
      tgt <- ci$target_probability %||% 0.95
      reached <- sum(pw >= tgt)
      cat(sprintf("  parameterwise Pr at n*: min %.3f, median %.3f, max %.3f (%d params)\n",
                  min(pw), stats::median(pw), max(pw), length(pw)))
      if (identical(ci$which_n, "global")) {
        # context, NOT failure: at the (smaller) global n*, params need not all be at target
        cat(sprintf("  at the global n*, %d/%d parameters individually reach target;\n",
                    reached, length(pw)))
        if (!is.null(ci$n_star_pw))
          cat(sprintf("  full parameterwise coverage requires n* = %d (validate which_n = 'pw').\n",
                      ci$n_star_pw))
      } else {
        cat(sprintf("  parameters reaching target: %d / %d\n", reached, length(pw)))
      }
    }
  } else if (x$scope == "planning_edge") {
    pe <- r$planning_edge
    cat(sprintf("  planning edge (%d, %d), rho = %s\n",
                r$edge[1], r$edge[2], format(r$edge_rho %||% NA, digits = 3)))
    cat(sprintf("    power_h0 (exclude absent) = %.3f  (FPR = %s)\n",
                as.numeric(pe$power_h0)[1],
                if (!is.null(pe$fpr_h0)) sprintf("%.3f", as.numeric(pe$fpr_h0)[1]) else "NA"))
    cat(sprintf("    power_h1 (detect present) = %.3f  (FNR = %s)\n",
                as.numeric(pe$power_h1)[1],
                if (!is.null(pe$fnr_h1)) sprintf("%.3f", as.numeric(pe$fnr_h1)[1]) else "NA"))
  } else {  # all_edges
    tab <- r$per_edge
    p1  <- sum(tab$power_h1 >= r$pow1, na.rm = TRUE)
    p0  <- sum(tab$power_h0 >= r$pow0, na.rm = TRUE)
    cat("  Guarantee check at n* (present edges):\n")
    cat(sprintf("    detect presence (power_h1 >= %.2f): %d / %d edges\n",
                r$pow1, p1, nrow(tab)))
    cat(sprintf("    exclude absence (power_h0 >= %.2f): %d / %d edges\n",
                r$pow0, p0, nrow(tab)))
    short1 <- tab[tab$power_h1 < r$pow1, , drop = FALSE]
    short0 <- tab[tab$power_h0 < r$pow0, , drop = FALSE]
    if (nrow(short1)) {
      cat("    underpowered for detection (H1):\n")
      for (i in seq_len(nrow(short1)))
        cat(sprintf("      (%d, %d): power_h1 = %.3f\n",
                    short1$m[i], short1$l[i], short1$power_h1[i]))
    }
    if (nrow(short0)) {
      cat("    underpowered for exclusion (H0):\n")
      for (i in seq_len(nrow(short0)))
        cat(sprintf("      (%d, %d): power_h0 = %.3f\n",
                    short0$m[i], short0$l[i], short0$power_h0[i]))
    }
  }
  invisible(object)
}

###############################################################################
#                        GGM prior study simulation                           #
###############################################################################

#' @title Constructs an object for simulating Gaussian graphical model prior studies
#' @description Collects the inputs needed to simulate prior studies for a
#'   Gaussian graphical model, then passes the result to \code{\link{simulate_prior_study}}.
#' @param p Number of nodes. Must be at least 3.
#' @param nu Prior study size: the number of observations each simulated study
#'   collects. Must exceed \code{p}.
#' @param G A \code{p} by \code{p} symmetric 0/1 adjacency matrix with zero
#'   diagonal. Supply this argument to hold the graph fixed across studies, so that only
#'   the estimated precision matrix varies.
#' @param structure A graph generator, as in \code{\link{generate_graph}}. Supply
#'   this argument instead of \code{G} to draw a new graph for each study.
#' @param ... Further arguments for the generator, e.g. \code{prob} for
#'   \code{"structure = Bernoulli"}. Only used with \code{structure}.
#' @return A \code{ggm_study} object, which also inherits from \code{bgm_study}.
#' @family prior elicitation
#' @export
ggm_study <- function(p, nu, G = NULL, structure = NULL, ...) {
  stopifnot(is.numeric(p), length(p) == 1L, p >= 3,
            is.numeric(nu), length(nu) == 1L)
  p  <- as.integer(round(p))
  nu <- as.integer(round(nu))

  if (nu <= p)
    stop("'nu' must exceed 'p': a study of nu observations cannot identify a ",
         "p x p precision matrix.", call. = FALSE)

  if (is.null(G) == is.null(structure))
    stop("Supply exactly one of 'G' (a fixed graph) or 'structure' ",
         "(a new graph per study).", call. = FALSE)

  if (!is.null(G)) {
    G <- unname(as.matrix(G))
    stopifnot(is.numeric(G), isSymmetric(G), nrow(G) == p, ncol(G) == p,
              all(G %in% c(0, 1)), all(diag(G) == 0))
  }

  new_bgm_study("ggm", list(p = p, nu = nu, G = G,
                            structure = structure, gen_args = list(...)))
}


#' @title Simulate prior studies for a Gaussian graphical model
#' @description Simulates prior studies from the Gaussian likelihood: for each
#'   study, \code{nu} observations are drawn from a graph-respecting precision
#'   matrix and the precision matrix is re-estimated under the same graph. What
#'   comes back is what the study found, not the truth it was drawn from, so it
#'   carries the estimation noise of a study that size.
#' @param study A \code{ggm_study} object, as returned by \code{\link{ggm_study}}.
#' @param n_studies Number of independent studies to simulate. Default 1.
#' @param verbose Print diagnostics for each simulated precision matrix. Off by
#'   default.
#' @param ... Ignored, present for consistency with the generic.
#' @return A list of \code{n_studies} \code{ggm_parameters} objects, ready for
#'   \code{\link{elicit_prior}}.
#' @details
#' Studies are independent. If the study was specified with a fixed \code{G},
#' every study uses it and only the estimated precision matrix varies; if it was
#' specified with a \code{structure}, a new graph is drawn for each study. Call
#' \code{\link[base]{set.seed}} first for reproducible output.
#' @family prior elicitation
#' @export
simulate_prior_study.ggm_study <- function(study, n_studies = 1L,
                                           verbose = FALSE, ...) {
  stopifnot(is.numeric(n_studies), length(n_studies) == 1L, n_studies >= 1)
  n_studies <- as.integer(round(n_studies))

  inp <- study$inputs
  lapply(seq_len(n_studies), function(i) {
    G <- if (!is.null(inp$G)) inp$G
         else do.call(generate_graph,
                      c(list(p = inp$p, structure = inp$structure), inp$gen_args))
    K <- cpp_simulate_ggm_study_precision(p = inp$p, nu = inp$nu,
                                          G = G, verbose = verbose)
    ggm_parameters(K = K, G = G, nu = inp$nu)
  })
}