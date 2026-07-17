#' graph_is_sparse (internal)
#' @description Tests whether an adjacency matrix represents a sparse or a complete graph. 
#' @param G [matrix], a p x p symmetric 0/1 adjacency matrix.
#' @return [logical], TRUE if G is sparse (not complete), FALSE if complete.
#' @noRd
.graph_is_sparse <- function(G) {
  Gd <- G
  diag(Gd) <- 0
  p <- nrow(G)
  is_sparse <- sum(Gd) != p * (p - 1) 
  return(is_sparse)
}

#' generate_smallworld_graph (internal)
#' @description A function that generates a graph structure following the smallworld type structure, Watts-Strogatz model
#' @param p [integer], number of variables in the graph
#' @return A randomly generated smallworld structure
#' @noRd
.generate_smallworld_graph <- function(p) {
    prob    <- 0.2
    S       <- 0
    max_nei <- floor((p - 1) / 2)
    nei     <- sample(x = 1:max_nei, size = 1) # random nei (from igraph documentation: "the neighborhood within which the vertices of the lattice will be connected")
    if (p < 3) nei <- 1 # [deprecated because we only generate for p > 3] for small p's, set nei to 1 to avoid errors in igraph::sample_smallworld
    net <- NULL
    while (S < 1 && prob <= 1) {
        net <- igraph::sample_smallworld(dim = 1, size = p, nei = nei, p = prob)
        C   <- igraph::transitivity(net, type = "average")
        L   <- igraph::mean_distance(graph = net, unconnected = TRUE)
        Cr_vec <- rep(0, 100)
        Lr_vec <- rep(0, 100)

        for (i in 1:100) {
            rand_net    <- igraph::sample_gnm(igraph::vcount(net), igraph::ecount(net))
            Cr_vec[i]   <- igraph::transitivity(rand_net, type = "average")
            Lr_vec[i]   <- igraph::mean_distance(graph = rand_net, unconnected = TRUE)
        }

        Cr      <- mean(Cr_vec)
        Lr      <- mean(Lr_vec)
        S       <- if (Cr == 0 || Lr == 0) 0 else (C / Cr) / (L / Lr) 
        prob    <- prob + 0.1
    }

    if (S < 1) {
        warning(sprintf("generate_smallworld_graph: no graph reached small-worldness S >= 1 (last S = %.3f); returning closest candidate.", S), call. = FALSE)
    }

    G <- as.matrix(igraph::as_adjacency_matrix(net, sparse = FALSE))
    return(G)
}

#' generate_random_graph (internal)
#' @description A function that generates a graph structure following the Erdős-Rényi model
#' @param p [integer], number of variables in the graph
#' @return A randomly generated random structure
#' @noRd
.generate_random_graph <- function(p) {
    prob    <- stats::runif(1, 0, 1) # random inlcusion probability
    net     <- igraph::sample_gnp(n = p, p = prob)
    G       <- as.matrix(igraph::as_adjacency_matrix(net, sparse = FALSE))
    return(G)
}

#' generate_scalefree_graph (internal)
#' @description A function that generates a graph structure following the
#' Barabási-Albert scale-free model
#' @param p [integer], number of variables in the graph
#' @return A randomly generated scale-free structure
#' @noRd
.generate_scalefree_graph <- function(p) {
    max_m   <- floor(p / 2)
    m       <- sample(x = 1:max_m, size = 1) # random number of edges to attach from each new node
    net     <- igraph::sample_pa(n = p, m = m, directed = FALSE)
    G       <- as.matrix(igraph::as_adjacency_matrix(net, sparse = FALSE))
    return(G)
}

#' generate_Bernoulli_graph (internal)
#' @description A function that generates a graph structure following a Bernoulli model 
#' (similar to Erdős-Rényi model but with the possibility of fixing the edge inclusion probability)
#' @param p [integer], number of variables in the graph
#' @param prob [double], edge inclusion probability
#' @return A randomly generated Bernoulli structure
#' @noRd
.generate_Bernoulli_graph <- function(p, prob) {
    G <- matrix(0, nrow = p, ncol = p)
    G[upper.tri(G, diag = FALSE)] <- stats::rbinom(n = p * (p - 1) / 2, size = 1, prob = prob)
    G[lower.tri(G, diag = FALSE)] <- t(G)[lower.tri(G, diag = FALSE)]
    return(G)
}

#' @title Generate a random graph with a given structure
#' @description Draws a random undirected graph on \code{p} nodes from one of the available
#' structure generators, rejecting degenerate draws (the empty graph and the
#' complete graph) and resampling until a usable structure is obtained. The
#' returned adjacency matrix is what \code{\link{simulate_prior_study()}} uses for simulating
#' prior study parameters.
#' @param p Number of nodes (variables) in the graph. Must be at least 3. 
#' @param structure Which graph generator to use: \code{"smallworld"}, \code{"random"},
#'   \code{"scalefree"} or \code{"Bernoulli"}.
#' @param ... Further arguments passed to the graph generator. The
#'   \code{"Bernoulli"} structure requires \code{prob}, the edge-inclusion
#'   probability in \eqn{(0, 1)}. The other structures take no extra arguments.
#' @param max_attempts Maximum number of rejection-sampling draws before
#'   giving up with an error. Guards against parameter choices for which a
#'   non-degenerate graph is effectively unreachable.
#' @return A symmetric \code{p} by \code{p} adjacency matrix, 0/1 entries, zero
#'   diagonal.
#' @details
#' Generation uses rejection sampling: a draw is accepted only if it has at
#' least one edge and is not fully connected. Since no non-degenerate graph
#' exists for \code{p < 3}, \code{p} must be at least 3
#' @examples
#' set.seed(2026)
#' # Bernoulli (Erdos-Renyi) needs an edge probability
#' g <- generate_graph(p = 10, structure = "Bernoulli", prob = 0.2)
#' dim(g)
#'
#' \donttest{
#' g2 <- generate_graph(p = 10, structure = "smallworld")
#' }
#' @export
generate_graph <- function(p, structure = c("smallworld", "random", "scalefree", "Bernoulli"), ..., max_attempts = 1000L) {
  structure <- match.arg(structure)
  stopifnot(is.numeric(p), length(p) == 1, p >= 3, p == round(p))

  args <- list(...)
  if (structure == "Bernoulli") {
    if (is.null(args$prob)) stop("Argument 'prob' is required for 'Bernoulli'.")
    if (args$prob <= 0 || args$prob >= 1) stop("'prob' must lie in (0, 1).")
  }

  for (attempt in seq_len(max_attempts)) {
    G <- switch(structure,
      smallworld = .generate_smallworld_graph(p = p),
      random     = .generate_random_graph(p = p),
      scalefree  = .generate_scalefree_graph(p = p),
      Bernoulli  = .generate_Bernoulli_graph(p = p, prob = args$prob))
    diag(G) <- 0
    edges   <- sum(G) / 2
    if (edges > 0 && edges < p * (p - 1) / 2) return(G)
  }
  stop(sprintf("No non-degenerate graph after %d attempts; check 'p'/'prob'.",
               max_attempts))
}
