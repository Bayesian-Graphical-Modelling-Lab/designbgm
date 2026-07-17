#include <string>
#include <RcppArmadillo.h>
#include "prior_wishart.h"
#include "prior_gwishart.h"
#include "helpers.h"

// [[Rcpp::depends(RcppArmadillo)]]

// G-Wishart direct sampler comprises the following functions:
//
// find_neighbors
// create_indices_excluding_i
// hastie_adaptation_step
// hastie_adaptation
// rgwishart_direct
//
// Reference: Direct sampler of Lekonski, 2013. DOI: 10.1002/sta4.23

// Find neighbors of each node in the graph G
arma::field<arma::uvec> find_neighbors(const arma::umat& G) {
    arma::uword p = G.n_cols;
    arma::field<arma::uvec> neighbors(p);
    for (arma::uword i = 0; i < p; i++) {
        neighbors(i) = arma::find(G.col(i) == 1); // find neighbors of i (G is symmetric, we work with columns for simplicity)
    }
    return neighbors;
}

// Create indices excluding i from {0,...,p-1}, for each i
arma::field<arma::uvec> create_indices_excluding_i(const arma::uword& p) {
    arma::field<arma::uvec> indices_excluding_i(p);
    for (arma::uword i = 0; i < p; i++) {
        arma::uvec idx = arma::regspace<arma::uvec>(0, p - 1);
        idx.shed_row(i); // remove diagonal index
        indices_excluding_i(i) = idx;
    }
    return indices_excluding_i;
}

// Hastie adaptation step (based on Hastie et al. 2009) used by the G-Wishart sampler
double hastie_adaptation_step(arma::mat& W, const arma::mat& Sigma, const arma::uword& i, const arma::uword& p, const arma::field<arma::uvec>& indices_excluding_i, const arma::field<arma::uvec>& neighbors, arma::vec& beta_star_full) {
    // indices excluding i
    arma::uvec idx = indices_excluding_i(i);
    // index of i
    arma::uvec i_vec(1);
    i_vec(0) = i;
    // save the column before update
    arma::vec col_before = W(idx, i_vec);
    // algorithm
    arma::uvec N_i = neighbors(i); // neighbors of i
    if (N_i.n_elem > 0) {
        arma::mat W_i = W.submat(N_i, N_i);
        arma::mat Sigma_i = Sigma.submat(N_i, i_vec);
        arma::vec beta_star_i = arma::vec(arma::solve(W_i, Sigma_i, arma::solve_opts::fast)); // vector sizeof(N_i) x 1
        // BLAS (even faster) — only compute the necessary columns of W
        arma::vec update = W.cols(N_i) * beta_star_i;
        W(idx, i_vec) = update(idx);
        W(i_vec, idx) = update(idx).t();
    } else {
        W(i_vec, idx).zeros();
        W(idx, i_vec).zeros();
    }
    return arma::accu(arma::abs(W(idx, i_vec) - col_before));
}

// Iterative adaptation that returns an estimate of the precision matrix
arma::mat hastie_adaptation(arma::mat K, const arma::uword& p, const arma::field<arma::uvec>& indices_excluding_i, const arma::field<arma::uvec>& neighbors, double tol, arma::uword itermax, arma::vec& beta_star_full) {
    arma::mat Sigma = arma::inv_sympd(K);
    arma::mat W = Sigma;
    double mean_diff = 1.0;
    arma::uword niter = 0;
    while (mean_diff > tol && niter < itermax) {
        double sum_abs_change = 0.0;
        for (arma::uword i = 0; i < p; i++) {
            sum_abs_change += hastie_adaptation_step(W, Sigma, i, p, indices_excluding_i, neighbors, beta_star_full);
        }
        mean_diff = sum_abs_change/(p*p);
        niter++;
    }
    return arma::inv_sympd(W);
}

// Random G-Wishart direct sampler
arma::cube rgwishart_direct(const arma::uword& n, 
                            const arma::mat& K, 
                            const double& nu, 
                            const arma::umat& G, 
                            const double& tol, 
                            const arma::uword& itermax) {
    if (!G.is_symmetric()) {
        Rcpp::stop("G must be a symmetric matrix.");
    }
    // parametrization is nu = delta + |V| - 1 (delta as defined in Roverato, 2000), where |V| is the number of vertices in the graph (p)
    arma::uword p = G.n_cols;
    arma::mat Kchol = arma::chol(K, "upper"); // for rwishart_fast
    arma::cube X(p, p, n, arma::fill::zeros);
    arma::umat G_structure = G;
    G_structure.diag().zeros();
    bool is_complete = arma::accu(G_structure) == p * (p - 1);
    if(!is_complete) {
        arma::field<arma::uvec> neighbors = find_neighbors(G_structure);
        arma::field<arma::uvec> indices_excluding_i = create_indices_excluding_i(p);
        arma::vec beta_star_full(p, arma::fill::zeros);
        for (arma::uword i = 0; i < n; i++) {
            // Generate Wishart(S,nu)
            arma::mat K_i = rwishart_fast(nu, Kchol); // generate from Wishart distribution with nu degrees of freedom and scale matrix K
            // Adaptation of algorithm presented by Hastie et al. (2009)
            X.slice(i) = hastie_adaptation(K_i, p, indices_excluding_i, neighbors, tol, itermax, beta_star_full);
        }
    } else {
        for (arma::uword i = 0; i < n; i++) {
            X.slice(i) = rwishart_fast(nu, Kchol); // generate from Wishart distribution with nu degrees of freedom and scale matrix K
        }
    }
    return X;
}

// Block Gibbs (edgewise) sampler for the G-Wishart distribution W_G(b, D) , comprising the following functions:
//
// sym (helper function)
// block_gibbs_update
// build_edgewise_index_sets
// build_rests
// build_chol_scales
// rgwishart_gibbs
//
// Reference: Wang & Li, 2012, "Efficient Gaussian graphical model determination under G-Wishart prior distributions". Electronic Journal of Statistics, 6, 168-198.
//
// From Eq. 2.5 + Lemma 1 of Roverato [23]:
// Given a complete subset Ik of V, the Schur-complement satisfies
//
// Omega[Ik,Ik] - Omega[Ik,rest] * inv(Omega[rest,rest]) * Omega[rest,Ik]
//         | rest  ~  W(b, D[Ik,Ik])
//
//   So one Gibbs step for block Ik:
//     1. Draw A ~ W(b, D[Ik,Ik])  via rwishart_fast(b, chol(inv(D[Ik,Ik]), "upper"))
//     2. Set  Omega[Ik,Ik] = A + Omega[Ik,rest] * inv(Omega[rest,rest]) * Omega[rest,Ik]
//
//   The upper Cholesky factor of inv(D[Ik,Ik]) is computed once per block
//   before the sampling loop and reused at every iteration.
//
//   Sigma = inv(Omega) is updated via the rank-|Ik| Woodbury update:
//
//     Sigma_new = Sigma - Sigma[:,Ik] * inv( inv(Delta) + Sigma[Ik,Ik] ) * Sigma[Ik,:]
//
//   where Delta =  Omega[Ik,Ik]_old - Omega[Ik,Ik]_new.
//
//   NOTE for future improvement: full recompute of Sigma used as fallback if the Woodbury inner matrix is
//   numerically singular --> This is not implemented yet.


// Symmetry helper (enforce exact symmetry)
//
// Numerical symmetry: both Omega and Sigma are symmetrized at three points to prevent
//   drift from floating points accumulating across iterations:
//     (1) Omega_new_sub after the Schur complement addition
//     (2) Sigma after each Woodbury update
//     (3) Omega at the end of every iteration
static inline arma::mat sym(const arma::mat& X)
{
    return 0.5 * (X + X.t());
}

// Perform a single block Gibbs update (block ik)
static void block_gibbs_update(arma::mat&        Omega,
                               arma::mat&        Sigma,
                               const arma::uvec& ik,
                               const arma::uvec& rest,
                               double            nu,
                               const arma::mat&  chol_scale_ik)
{
    // Step 1: draw A ~ W(b, inv(D[Ik,Ik]))
    // b = nu - p + ik.n_elem, where nu is the G-Wishart degrees of freedom parameter, p is the total number of nodes in the graph (size of Omega) and ik.n_elem is the size of the current block (which is 1 for an isolated node, 2 for an edge).
    double p = Omega.n_rows;
    arma::mat A = rwishart_fast(nu - p + ik.n_elem, chol_scale_ik);

    // Step 2: Omega[Ik,Ik] = A + Schur complement
    arma::mat Omega_new_sub;
    if (rest.n_elem == 0) {
        Omega_new_sub = A;
    } else {
        arma::mat Sig_ee      = Sigma.submat(ik, ik);
        arma::mat Sig_ee_inv  = arma::inv(Sig_ee);
        arma::mat Sig_re      = Sigma.submat(rest, ik);       // (p-2) x 2
        arma::mat Om_rest_inv = Sigma.submat(rest, rest) - Sig_re * Sig_ee_inv * Sig_re.t();  // (p-2) x (p-2)
        arma::mat Om_ik_rest  = Omega.submat(ik, rest);
        Omega_new_sub = sym(A + Om_ik_rest * Om_rest_inv * Om_ik_rest.t());
    }

    // Step 3: update Sigma via paper's formula (p.177)
    arma::mat Delta    = Omega.submat(ik, ik) - Omega_new_sub;
    arma::mat Sig_cols = Sigma.cols(ik);
    arma::mat W        = arma::inv(arma::inv(Delta) - Sigma.submat(ik, ik));
    Sigma = sym(Sigma + Sig_cols * W * Sig_cols.t());
    // Write new block into Omega
    Omega.submat(ik, ik) = Omega_new_sub;
}


// Build index sets: edgewise  I = E union {isolated nodes}
static std::vector<arma::uvec> build_edgewise_index_sets(const arma::umat& G)
{
    arma::uword p = G.n_rows;
    std::vector<arma::uvec> index_sets;
    std::vector<bool> has_edge(p, false);

    for (arma::uword i = 0; i < p; ++i) {
        for (arma::uword j = i + 1; j < p; ++j) {
            if (G(i, j) == 1) {
                arma::uvec block = {i, j};
                index_sets.push_back(block);
                has_edge[i] = true;
                has_edge[j] = true;
            }
        }
    }
    for (arma::uword i = 0; i < p; ++i) {
        if (!has_edge[i]) {
            arma::uvec block = {i};
            index_sets.push_back(block);
        }
    }
    return index_sets;
}

// Build complement index vectors (done once, before the loop)
static std::vector<arma::uvec> build_rests(int p, const std::vector<arma::uvec>& index_sets)
{
    std::vector<arma::uvec> rests(index_sets.size());
    for (std::size_t k = 0; k < index_sets.size(); ++k) {
        const arma::uvec& ik = index_sets[k];
        std::vector<arma::uword> rv;
        for (int i = 0; i < p; ++i) {
            bool in_ik = false;
            for (arma::uword idx : ik)
                if ((arma::uword)i == idx) { in_ik = true; break; }
            if (!in_ik) rv.push_back(i);
        }
        rests[k] = arma::uvec(rv);
    }
    return rests;
}

// Precompute per-block Cholesky factors of inv(D[Ik,Ik])
static std::vector<arma::mat> build_chol_scales(
        const std::vector<arma::uvec>& index_sets,
        const arma::mat& D)
{
    std::vector<arma::mat> chols(index_sets.size());
    for (std::size_t k = 0; k < index_sets.size(); ++k) {
        const arma::uvec& ik = index_sets[k];
        arma::mat scale_k = arma::inv_sympd(D.submat(ik, ik));
        chols[k] = arma::chol(scale_k, "upper");
    }
    return chols;
}

// G-Wishart Block Gibbs sampler (edgewise) for the G-Wishart distribution W_G(b, D) with b = nu - p + 1
arma::cube rgwishart_gibbs( const arma::uword& n,
                            const arma::mat&  K,
                            const double& nu,
                            const arma::umat& G,
                            const arma::uword& burnin,
                            const Rcpp::Nullable<arma::mat>& init)
{
    int p = K.n_rows;
    arma::mat D = arma::inv_sympd(K);
    arma::mat Omega = init.isNull() ? arma::diagmat(nu / arma::diagvec(D)) : Rcpp::as<arma::mat>(init);
    std::vector<arma::uvec> index_sets  = build_edgewise_index_sets(G);
    std::vector<arma::uvec> rests       = build_rests(p, index_sets);
    std::vector<arma::mat>  chol_scales = build_chol_scales(index_sets, D);
    arma::mat  Sigma = arma::inv_sympd(sym(Omega));
    arma::cube samples(p, p, n);
    arma::uword total = burnin + n;
    for (arma::uword t = 0; t < total; ++t) {
        for (std::size_t k = 0; k < index_sets.size(); ++k)
            block_gibbs_update(Omega, Sigma, index_sets[k], rests[k], nu, chol_scales[k]);
        if (t >= burnin)
            samples.slice(t - burnin) = Omega;
    }
    return samples;
}

// G-Wishart dispatch function
arma::cube rgwishart(const arma::uword& n, 
                            const arma::mat& K, 
                            const double& nu, 
                            const arma::umat& G, 
                            const std::string& sampler,
                            const double& tol, 
                            const arma::uword& itermax,
                            const arma::uword& burnin,
                            const Rcpp::Nullable<arma::mat>& init)
{
    if (sampler == "direct") {
        return rgwishart_direct(n, K, nu, G, tol, itermax);
    } else if (sampler == "gibbs") {
        return rgwishart_gibbs(n, K, nu, G, burnin, init);
    } else {
        Rcpp::stop("Invalid sampler type. Use 'direct' or 'gibbs'.");
    }
}

// Variance of the G-Wishart distribution
arma::mat prior_variance_gwishart(const arma::umat &G, const arma::uvec &select_parameter, const arma::cube &samples){

    arma::uword nsim   = samples.n_slices;
    arma::uword p      = G.n_cols;
    //arma::uword n_pars = select_parameter.n_elem;

    // which lower-triangular linear positions are the free params?
    arma::uvec lower_indices = arma::trimatl_ind(arma::size(G));
    arma::uvec free_lin = lower_indices(select_parameter);   // linear indices into vec(p*p)

    // zero-copy view: (p*p) x nsim, row r = trajectory of vec-entry r
    const arma::mat S_flat(const_cast<double*>(samples.memptr()), p*p, nsim, false, true);

    // draws = nsim x n_pars: each column is the trajectory of one free param
    // S_flat.rows(free_lin) is n_pars x nsim; transpose to nsim x n_pars
    arma::mat draws = S_flat.rows(free_lin).t();

    return arma::cov(draws, 1);
}