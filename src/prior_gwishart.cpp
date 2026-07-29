#include <string>
#include <cmath>
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

// [[Rcpp::export]]
arma::mat cpp_constrain_precision_to_graph(const arma::mat& K,
                                           const arma::umat& G,
                                           double tol,
                                           int itermax) {
    arma::uword p = K.n_cols;
    arma::field<arma::uvec> neighbors = find_neighbors(G);
    arma::field<arma::uvec> idx = create_indices_excluding_i(p);
    arma::vec beta_star_full(p, arma::fill::zeros);
    return hastie_adaptation(K, p, idx, neighbors, tol, itermax, beta_star_full);
}


// =======================================================================================
// G-Wishart posterior samplers - for both non informative and informative priors
//      - Double Continuous-Time (DCT) BDMCMC sampler (all edges at once, 
//                                   Rao-Blackwellized posterior inclusion probabilities)
//      - Double Conditional Bayes Factor (DCBF) BDMCMC sampler (one edge at a iteration)
// =======================================================================================

// (log scale) Conditional Bayes factor term, eq.(6). Phi = chol(K_perm,"upper"); edge at
// 0-based (p-2, p-1). U = permuted rate matrix (Dn for posterior, D0 for prior).
static inline double cbf_logN(const arma::mat& Phi, const arma::mat& U, int p) {
    const int a = p - 2, b = p - 1;
    const double phi_aa = Phi(a, a);
    const double u_bb = U(b, b);
    const double u_ab = U(a, b);
    double s = 0.0;
    for (int l = 0; l < a; ++l) s += Phi(l, a) * Phi(l, b);
    const double inner = phi_aa * u_ab / u_bb - s / phi_aa;
    return std::log(phi_aa) + 0.5 * std::log(2.0 * M_PI / u_bb)
         + 0.5 * u_bb * inner * inner;
}

// numerically-stable logistic sigmoid
static inline double sigmoid(double x) {
    return x >= 0.0 ? 1.0 / (1.0 + std::exp(-x))
                    : std::exp(x) / (1.0 + std::exp(x));
}

// permutation placing (i,j) at the last two positions, others in original order
static inline arma::uvec perm_last(int p, int i, int j) {
    arma::uvec pm(p);
    int idx = 0;
    for (int k = 0; k < p; ++k) if (k != i && k != j) pm(idx++) = (arma::uword)k;
    pm(p - 2) = (arma::uword)i;
    pm(p - 1) = (arma::uword)j;
    return pm;
}

// =============================================================================
// Double Continuous-Time (DCT) BDMCMC — Hinne, Lenkoski, Heskes, van Gerven (2014)
// Exact-exchange birth-death sampler for GGM structure + precision with an
// informative (or non-informative) G-Wishart prior  K ~ W_G(nu, G0, K0/nu).
//   delta0 = nu-p+1 ,  D0 = nu*K0^{-1} ,  Dn = D0 + n*Sigma.
// Rates use the conditional Bayes factor N(.) (eq.6) + one auxiliary prior draw,
// so no normalizing constants and no Monte-Carlo integrals are needed.
// Reuses rgwishart (Lenkoski direct sampler) from the same file.
// Slower but more stable than the DCBF at shorter iterations (provides Rao-Blackwellized 
// posterior inclusion probabilities).
// =============================================================================

// Continuous time sampler for G-Wishart posterior with informative prior (or non-informative prior)
bdmcmc_dct_result bdmcmc_dct_sampler(
        const arma::mat D0,          // AK prior rate matrix
        const arma::mat Dn,            // AK posterior rate matrix
        const arma::mat scale_prior,    // = K0/nu   (auxiliary scale)
        const arma::mat scale_post,      // posterior scale
        int              n,            // new-data sample size
        double           nu,           // prior degrees of freedom (elicited)
        const arma::mat& prior_logodds,      // p x p prior edge-inclusion probs (m<l)
        const arma::mat& G_start,      // p x p starting adjacency (symmetric 0/1)
        int              n_iter,
        int              n_burnin,
        std::string      gwish_sampler,       // "direct" or "gibbs"
        double           gwish_tol,
        arma::uword      gwish_iter,
        arma::uword      gwish_burnin)
{
    Rcpp::RNGScope scope;
    const int p = (int)D0.n_rows;

    const double nu0 = nu;                       // std-Wishart df of the prior
    const double nun = nu + (double)n;           // std-Wishart df of the posterior
    if (nu - p + 1.0 <= 0.0)
        Rcpp::warning("nu must exceed p-1 for a proper prior.");

    // These are computed outside the loop to avoid repeated inversions of the same matrices
    // const arma::mat D0 = nu * arma::inv_sympd(K0);          // AK prior rate matrix
    // const arma::mat Dn = D0 + (double)n * Sigma;            // AK posterior rate matrix
    // const arma::mat scale_prior = arma::inv_sympd(D0);      // = K0/nu   (auxiliary scale)
    // const arma::mat scale_post  = arma::inv_sympd(Dn);      // posterior scale

    // g_prior was the matrix of prior edge-inclusion probabilities (m<l) from R, converted to log-odds here
    // arma::mat prior_logodds(p, p, arma::fill::zeros);
    // for (int i = 0; i < p; ++i)
    //     for (int j = i + 1; j < p; ++j) {
    //         double q = std::min(std::max(g_prior(i, j), 1e-12), 1.0 - 1e-12);
    //         prior_logodds(i, j) = std::log(q / (1.0 - q));
    //     }

    arma::mat G = arma::max(G_start, G_start.t());  G.diag().zeros();
    G.transform([](double v){ return v > 0.5 ? 1.0 : 0.0; });

    // stable weighted accumulators
    double    max_logw = -arma::datum::inf, w_sum = 0.0;
    arma::mat pip(p, p, arma::fill::zeros), K_hat(p, p, arma::fill::zeros);

    arma::mat log_rate(p, p, arma::fill::zeros);
    arma::mat pip_rb(p, p, arma::fill::zeros); // Rao-Blackwellized posterior inclusion probabilities
    arma::mat rb(p, p, arma::fill::zeros);

    const int total = n_burnin + n_iter;

    for (int iter = 0; iter < total; ++iter) { // (total cost = (n_burnin + n_iter) *  p*(p-1)/2 edges * cost of rgwishart)
        //Rcpp::Rcout << "iter=" << iter + 1 << "/" << total << "  n=" << n << std::endl;
        // posterior precision draw on the current graph
        const arma::umat G_u = (G > 0.5);
        arma::mat K = rgwishart((arma::uword)1, scale_post, nun, G_u, gwish_sampler, gwish_tol, gwish_iter, gwish_burnin, Rcpp::Nullable<arma::mat>()).slice(0);

        // toggle rate for every edge (loop over all edges, with Rao-Blackwellized posterior inclusion probabilities)
        double log_R_max = -arma::datum::inf;
        for (int i = 0; i < p; ++i) {
            for (int j = i + 1; j < p; ++j) {

                const bool present = (G(i, j) > 0.5);
                const arma::uvec pm = perm_last(p, i, j);

                // current K (data side): N(K, Dn)
                const arma::mat Phi   = arma::chol(K.submat(pm, pm), "upper");
                const double    log_N_K   = cbf_logN(Phi, Dn.submat(pm, pm), p);

                // toggled graph + auxiliary prior draw on it: N(K0, D0)
                arma::mat Gt = G; Gt(i, j) = Gt(j, i) = present ? 0.0 : 1.0;
                const arma::umat Gt_u = (Gt.submat(pm, pm) > 0.5);
                const arma::mat K0aux = rgwishart((arma::uword)1, scale_prior.submat(pm, pm), nu0, Gt_u, gwish_sampler, gwish_tol, gwish_iter, gwish_burnin, Rcpp::Nullable<arma::mat>()).slice(0);
                const arma::mat Phi0  = arma::chol(K0aux, "upper");
                const double    log_N_K0  = cbf_logN(Phi0, D0.submat(pm, pm), p);

                // rate (log scale); birth: N_K/N_K0 * theta/(1-theta), death inverts
                double lr;
                if (present)   // death event
                    lr = (log_N_K0 - log_N_K) - prior_logodds(i, j);
                else           // birth event
                    lr = (log_N_K - log_N_K0) + prior_logodds(i, j);

                log_rate(i, j) = lr;
                rb(i, j) = rb(j, i) = present ? sigmoid(-lr) : sigmoid(lr);
                if (lr > log_R_max) log_R_max = lr;
            }
        }

        // total rate R, expected sojourn weight w = 1/R   (log-stable)
        double sum_shift = 0.0;
        for (int i = 0; i < p; ++i)
            for (int j = i + 1; j < p; ++j)
                sum_shift += std::exp(log_rate(i, j) - log_R_max);
        const double logw = -(log_R_max + std::log(sum_shift));

        // accumulate (numerically stable weighted average)
        if (iter >= n_burnin) {
            if (logw > max_logw) {
                const double r = (max_logw == -arma::datum::inf) ? 0.0 : std::exp(max_logw - logw);
                // multiply old accumulators by r, then add the new sample with weight 1.0
                w_sum *= r; 
                pip *= r; 
                K_hat *= r; 
                pip_rb *= r;
                // update max_logw and add the new sample with weight 1.0
                max_logw = logw;
                w_sum += 1.0; 
                pip += G; 
                K_hat += K; 
                pip_rb += rb;
            } else {
                const double ww = std::exp(logw - max_logw);
                // add new sample with weight ww
                w_sum += ww; 
                pip += ww * G; 
                K_hat += ww * K; 
                pip_rb += ww * rb;
            }
        }

        // pick the event  e* ~ rate / R  and toggle
        double u = R::runif(0.0, 1.0) * sum_shift, acc = 0.0;
        int ei = -1, ej = -1;
        for (int i = 0; i < p && ei < 0; ++i)
            for (int j = i + 1; j < p; ++j) {
                acc += std::exp(log_rate(i, j) - log_R_max);
                if (u <= acc) { 
                    ei = i; 
                    ej = j; 
                    break; 
                }
            }
        if (ei < 0) { 
            ei = 0; 
            ej = 1; 
        }
        const double nv = (G(ei, ej) > 0.5) ? 0.0 : 1.0;
        G(ei, ej) = nv; 
        G(ej, ei) = nv;
    }

    if (w_sum > 0.0) { 
        pip /= w_sum; 
        K_hat /= w_sum; 
        pip_rb /= w_sum;
    }
    pip.diag().ones();
    pip_rb.diag().ones();   

    // return results as a struct
    bdmcmc_dct_result result;
    result.pip = pip;
    result.pip_rb = pip_rb;
    result.K_hat = K_hat;
    result.last_graph = G;

    return result;
}

// Double Conditional Bayes Factor sampler (faster but slower convergence, requires more iterations) - This is a Metropolis-Hastings sampler (discrete time)
bdmcmc_dcbf_result bdmcmc_dcbf_sampler(int n, 
                                        double nu,
                                        const arma::mat &D0, 
                                        const arma::mat &Dn, // D0 = nu*K0^{-1}, Dn = D0 + n*Sigma
                                        const arma::mat &scale_prior, 
                                        const arma::mat &scale_post, // scale_prior = inv(D0), scale_post = inv(Dn)
                                        const arma::mat& plo, 
                                        const arma::mat& G_start, 
                                        int n_iter, 
                                        int n_burnin,
                                        std::string gwish_sampler,
                                        double gwish_tol, 
                                        arma::uword gwish_iter,
                                        arma::uword gwish_burnin)                                       
{
    Rcpp::RNGScope scope;
    const int p = (int)D0.n_rows;
    const double nu0 = nu, nun = nu + (double)n;
    if (nu - p + 1.0 <= 0.0) Rcpp::warning("nu must exceed p-1 for a proper prior.");

    arma::umat G = (G_start > 0.5);

    const int M = p * (p - 1) / 2;
    double    w_sum = 0.0; long accepts = 0;
    arma::mat pip(p, p, arma::fill::zeros), K_hat(p, p, arma::fill::zeros);
    const int total = n_burnin + n_iter;

    for (int iter = 0; iter < total; ++iter) { // (total cost = (n_burnin + n_iter) * one edge * cost of rgwishart)

        // Gibbs: posterior precision on the current graph
        arma::mat K = rgwishart((arma::uword)1, scale_post, nun, G, gwish_sampler, gwish_tol, gwish_iter, gwish_burnin, Rcpp::Nullable<arma::mat>()).slice(0);

        // Record the consistent (G, K) pair  (equal weights — it's Metropolis)
        if (iter >= n_burnin) { w_sum += 1.0; pip += arma::conv_to<arma::mat>::from(G); K_hat += K; }
        
        // Propose one uniformly chosen edge among the M pairs (one edge at a time - possibly slow convergence)
        int k = (int)std::floor(R::runif(0.0, (double)M)); if (k >= M) k = M - 1;
        int i = 0, rem = k;
        while (rem >= p - 1 - i) { rem -= (p - 1 - i); ++i; }
        const int j = i + 1 + rem;

        const bool       present = (G(i, j) > 0.5);
        const arma::uvec pm      = perm_last(p, i, j);

        const double logN_K = cbf_logN(arma::chol(K.submat(pm, pm), "upper"),
                                       Dn.submat(pm, pm), p);

        arma::umat Gt = G; Gt(i, j) = Gt(j, i) = present ? 0 : 1;
        const arma::mat K0aux = rgwishart((arma::uword)1, scale_prior.submat(pm, pm),
                                          nu0, Gt.submat(pm, pm),
                                          gwish_sampler, gwish_tol, gwish_iter, gwish_burnin, Rcpp::Nullable<arma::mat>()).slice(0); 
        const double logN_K0 = cbf_logN(arma::chol(K0aux, "upper"),
                                        D0.submat(pm, pm), p);

        // Exchange CBF acceptance (eq.7); symmetric proposal --> no q ratio
        const double log_alpha = present
            ? (logN_K0 - logN_K - plo(i, j))     // death event
            : (logN_K  - logN_K0 + plo(i, j));   // birth event

        if (std::log(R::runif(0.0, 1.0)) < log_alpha) {
            G(i, j) = G(j, i) = present ? 0 : 1;
            ++accepts;
        }
        
    }

    pip /= w_sum; 
    pip.diag().ones();
    K_hat /= w_sum; 
   
    bdmcmc_dcbf_result result;
    result.pip = pip;
    result.K_hat = K_hat;
    result.last_graph = G;
    result.accept_rate = (double)accepts / (double)total;

    return result;
}
