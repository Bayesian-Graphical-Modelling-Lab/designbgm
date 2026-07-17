#include <string>
#include <RcppArmadillo.h>
#include "family_ggm.h"
#include "prior_wishart.h" 
#include "prior_gwishart.h" 
#include "helpers.h"

// [[Rcpp::depends(RcppArmadillo)]]

// Build free index for G-Wishart
arma::field<arma::uvec> build_free_idx(const arma::uvec& sel, const arma::uword p) {
    arma::uword n_vech = p * (p + 1) / 2;
    arma::uvec lut_row(n_vech), lut_col(n_vech);
    arma::uword v = 0;
    for (arma::uword j = 0; j < p; j++)
        for (arma::uword i = j; i < p; i++) {
            lut_row(v) = i; lut_col(v) = j; v++;
        }

    arma::uword n_free = sel.n_elem;
    arma::uvec row_idx(n_free), col_idx(n_free);
    for (arma::uword s = 0; s < n_free; s++) {
        row_idx(s) = lut_row(sel(s));
        col_idx(s) = lut_col(sel(s));
    }

    arma::field<arma::uvec> idx(2);
    idx(0) = row_idx;
    idx(1) = col_idx;
    return idx;
}

// Computes the free-parameter submatrix of D*(A⊗A)*D^T ,using the Magnus-Neudecker identity, with no need of Kronecker product. (used for hessian / fisher information kernels)
arma::mat kron_reduce_free_D(const arma::mat&               A, 
                          const arma::field<arma::uvec>& free_idx) {

  const arma::uvec& row_idx = free_idx(0);   // i-index for each free param
  const arma::uvec& col_idx = free_idx(1);   // j-index for each free param
  arma::uword n_free = row_idx.n_elem;

  arma::mat M(n_free, n_free, arma::fill::zeros);

  for (arma::uword a = 0; a < n_free; a++) {
    arma::uword i  = row_idx(a);
    arma::uword j  = col_idx(a);
    bool diag_a    = (i == j);

    for (arma::uword b = a; b < n_free; b++) {
      arma::uword k = row_idx(b);
      arma::uword l = col_idx(b);
      bool diag_b   = (k == l);

      double val;
      if      ( diag_a &&  diag_b)
        val = A(i,k) * A(i,k);
      else if ( diag_a && !diag_b)
        val = 2.0 * A(i,k) * A(i,l);
      else if (!diag_a &&  diag_b)
        val = 2.0 * A(i,k) * A(j,k);
      else
        val = 2.0 * (A(i,k)*A(j,l) + A(i,l)*A(j,k));

      M(a, b) = val;
      M(b, a) = val;
    }
  }
  return M;
}

// Computes the free-parameter submatrix of D_plus*(A⊗A)*D_plus^T, based on the kron_reduce_free_D code (used for inverse fisher information kernels)
arma::mat kron_reduce_free_Dplus(const arma::mat&               A,
                                 const arma::field<arma::uvec>& free_idx) {

  const arma::uvec& row_idx = free_idx(0);
  const arma::uvec& col_idx = free_idx(1);
  arma::uword n_free = row_idx.n_elem;

  arma::mat M(n_free, n_free, arma::fill::zeros);

  for (arma::uword a = 0; a < n_free; a++) {
    arma::uword i = row_idx(a);
    arma::uword j = col_idx(a);
    bool diag_a   = (i == j);

    for (arma::uword b = a; b < n_free; b++) {
      arma::uword k = row_idx(b);
      arma::uword l = col_idx(b);
      bool diag_b   = (k == l);

      double val;
      if      ( diag_a &&  diag_b)
        val = A(i,k) * A(i,k);                          
      else if ( diag_a && !diag_b)
        val = A(i,k) * A(i,l);                   
      else if (!diag_a &&  diag_b)
        val = A(i,k) * A(j,k);                          
      else
        val = 0.5 * (A(i,k)*A(j,l) + A(i,l)*A(j,k));           

      M(a, b) = val;
      M(b, a) = val;
    }
  }
  return M;
}

// Convert precision matrix to partial correlations
// [[Rcpp::export]]
arma::mat cpp_precision_to_partial_correlations(const arma::mat& K) {
    int p = K.n_rows;
    arma::mat Rho(p, p, arma::fill::eye); // Create an identity matrix (values are 1 on the diagonal)
    for (int i = 0; i < (p-1); ++i) {
        for (int j = (i+1); j < p; ++j) {
                Rho(i, j) = -K(i, j) / std::sqrt(K(i, i) * K(j, j)); // Compute partial correlation
                Rho(j, i) = Rho(i, j); 
        }
    }
    return Rho;
} 

// Random precision matrix from prior
arma::mat random_precision_from_prior(const std::string& prior, const arma::mat& K, const int& nu, const arma::umat& G, const arma::mat& Kchol) {
    int p = K.n_rows;
    arma::mat x(p, p, arma::fill::zeros);
    if (prior.compare("wishart") == 0) {
        x = rwishart_fast(nu, Kchol);  // Kchol is the upper triangular Cholesky factor of K
    } else if (prior.compare("gwishart") == 0) {
        x = rgwishart(1, K, static_cast<double>(nu), G, "direct", 1e-08, 500, 500, R_NilValue).slice(0); 
    } else {
        Rcpp::stop("Unknown prior '" + prior + "'. Priors available: gwishart, wishart.");
    }
    return x;
}

// Simulate a random GGM study precision matrix with structure G
// [[Rcpp::export]]
arma::mat cpp_simulate_ggm_study_precision(const int&p, const int& nu, arma::umat G, const bool& verbose) {
    arma::mat Theta(p, p, arma::fill::zeros);
    G.diag().zeros(); // Set diagonal to zero (no self loops)
    arma::umat upper_idx = arma::trimatu(arma::ones<arma::umat>(G.n_rows, G.n_cols), 1); // upper triangle, no diagonal
    arma::uvec upper_vals = G.elem(arma::find(upper_idx));

    // // generating a prior study dense precision matrix using the well-behaved method of "rescaled orthonormal basis", that is t(Theta) * sigma^(-1) * Theta
    arma::mat Q,R;
    arma::qr(Q, R, rnorm_arma(p, p)); // orthonormal basis from a random matrix (we only need Q)
    arma::vec sigma = runif_arma(p, 0.1, 5.0); // random variances for the p components 
    Theta = Q.t() * arma::diagmat(1 / sigma) * Q; 

    if (G.is_empty() || arma::all(upper_vals == 1)) { // generate dense precision matrix if G is empty or fully connected
        arma::mat X = mvnrnd_arma(arma::zeros<arma::vec>(p),  arma::inv_sympd(Theta), nu);
        arma::mat Sigma = (X.t() * X) / nu;
        Theta = arma::inv_sympd(Sigma);
    } else { // generating a prior study precision matrix with structure G
        arma::field<arma::uvec> neighbors = find_neighbors(G);
        arma::field<arma::uvec> indices_excluding_i = create_indices_excluding_i(p);
        arma::vec beta_star_full(p, arma::fill::zeros);
        Theta = hastie_adaptation(Theta, p, indices_excluding_i, neighbors, 1e-6, 1000, beta_star_full);
        arma::mat X = mvnrnd_arma(arma::zeros<arma::vec>(p),  arma::inv_sympd(Theta), nu);
        arma::mat Sigma = (X.t() * X) / nu;
        Theta = hastie_adaptation(arma::inv_sympd(Sigma), p, indices_excluding_i, neighbors, 1e-6, 1000, beta_star_full);
    }

    if (!Theta.is_sympd()) {
        Rcpp::stop("Simulated precision matrix is not positive definite. "
            "This can happen when 'nu' is close to 'p' or the graph is very dense. "
            "Try with a larger 'nu'.");
    }
    if (verbose) {
        Rcpp::Rcout << "\nPrecision matrix diagnostics:\n";
        Rcpp::Rcout << "   -> Condition number : " << arma::cond(Theta) << "\n";
        Rcpp::Rcout << "   -> Determinant      : " << std::abs(arma::det(Theta)) << "\n";
        Rcpp::Rcout << "--------------------------------------------------\n";  
    }
    return Theta;
}

// Expected (wrt the parameters) inverse Fisher information for GGMs under a Wishart prior predictive distribution
arma::mat expected_inverse_fisher_ggm_wishart(const arma::mat &K, const double &nu, const arma::mat &D_plus, const arma::mat &M){
    arma::vec vecK = arma::vectorise(K);
    arma::mat kronK = arma::kron(K,K);
    arma::mat F = D_plus * ((nu*nu)*kronK + nu*(vecK*vecK.t()) + nu*(M * kronK)) * D_plus.t(); // Eq (9) Hagedorn M. 2022 'Expected Value of Matrix Quadratic Forms with Wishart distributed Random Matrices'
    F *= 2.0;
    return F;
}

// Expected (wrt the parameters) inverse Fisher information for GGMs under a G-Wishart prior predictive distribution
// BLAS optimization.
//
// IDEA:
//   The per-sample inverse-Fisher reduces to entries of S_s (x) S_s (Kronecker).
//   Averaging over samples, E[S (x) S] = (1/nsim) * sum_s vec(S_s) vec(S_s)^T.
//   With S_flat = (p*p) x nsim whose column s is vec(S_s), that whole average is
//   [[ONE]] matrix product:  ESS = S_flat * S_flat^T / nsim 
//
//   ESS is p^2 x p^2 and its entry ((i,k),(j,l)) = E[ S(i,k) S(j,l) ], where the
//   linear index of matrix position (r,c) in column-major vec is r + (c*p). So the
//   four cases of the original loop become lookups into the ESS:
//     diag_a &  diag_b : E[ S(i,k) S(i,k) ]                          (i==j, k==l)
//     diag_a & !diag_b : E[ S(i,k) S(i,l) ]                          (i==j)
//    !diag_a &  diag_b : E[ S(i,k) S(j,k) ]                          (k==l)
//    !diag_a & !diag_b : 0.5 ( E[S(i,k)S(j,l)] + E[S(i,l)S(j,k)] )
//
// NOTE: 
//   this version computes the full p^2 x p^2 second moment, including non-free
//   entries (which can be wasteful when G is sparse). This version should be
//   benchmarked in case of very sparse G, to see if the loop (old) version is faster.
//
arma::mat expected_inverse_fisher_ggm_gwishart(const arma::cube &samples,
                                                    const arma::field<arma::uvec> &free_idx){
    arma::uword nsim   = samples.n_slices;
    arma::uword p      = samples.n_rows;
    const arma::uvec& row_idx = free_idx(0);
    const arma::uvec& col_idx = free_idx(1);
    arma::uword n_free = row_idx.n_elem;

    // zero-copy view of the cube as (p*p) x nsim: column s is vec(S_s)
    const arma::mat S_flat(const_cast<double*>(samples.memptr()), p * p, nsim, false, true);

    // full second moment E[S (x) S] = (1/nsim) S_flat S_flat^T   — ONE BLAS call
    arma::mat ESS = (S_flat * S_flat.t()) / static_cast<double>(nsim);

    // helper: column-major linear index of matrix position (r, c)
    auto idx = [p](arma::uword r, arma::uword c) -> arma::uword { return r + c * p; };

    arma::mat F(n_free, n_free, arma::fill::zeros);

    for (arma::uword a = 0; a < n_free; a++) {
        arma::uword i = row_idx(a);
        arma::uword j = col_idx(a);
        bool diag_a   = (i == j);

        for (arma::uword b = a; b < n_free; b++) {
            arma::uword k = row_idx(b);
            arma::uword l = col_idx(b);
            bool diag_b   = (k == l);

            double val;
            if      ( diag_a &&  diag_b)
                val = ESS(idx(i, k), idx(i, k));
            else if ( diag_a && !diag_b)
                val = ESS(idx(i, k), idx(i, l));
            else if (!diag_a &&  diag_b)
                val = ESS(idx(i, k), idx(j, k));
            else
                val = 0.5 * ( ESS(idx(i, k), idx(j, l)) +
                              ESS(idx(i, l), idx(j, k)) );

            F(a, b) = val;
            F(b, a) = val;
        }
    }

    F *= 2.0;
    return F;
}

// Expected (wrt the parameters) Fisher information for a GGM under a Wishart prior predictive distribution
arma::mat expected_fisher_ggm_wishart(const arma::mat &K, double nu, const arma::mat &D, const arma::mat &M){
    arma::uword p = K.n_cols;
    arma::uword n_pars = p*(p+1)/2;    
    arma::mat F(n_pars, n_pars, arma::fill::zeros);
    // calculate F using Von Rosen formula
    arma::mat Kinv = arma::inv_sympd(K);
    F = second_moment_inverse_wishart(Kinv,nu,static_cast<double>(p),D,M);
    F *= 0.5;
    return F;
}


// Expected (wrt the parameters) Fisher information for GGMs under a G-Wishart prior predictive distribution
// BLAS optimization
//
// Same idea as the inverse-Fisher version above, with the difference that
// (1) it operates on the INVERSES of the samples: build samples_inv first,
// (2) then the second moment is E[Sinv (x) Sinv] = Sinv_flat * Sinv_flat^T / nsim.
//
// Mapping (idx(r,c) = r + (c*p) is the column-major linear index of position (r,c)):
//   diag_a &  diag_b :        M(idx(i,k), idx(i,k))
//   diag_a & !diag_b : 2.0  * M(idx(i,k), idx(i,l))
//  !diag_a &  diag_b : 2.0  * M(idx(i,k), idx(j,k))
//  !diag_a & !diag_b : 2.0  * ( M(idx(i,k), idx(j,l)) + M(idx(i,l), idx(j,k)) )
//   then F *= 0.5
//
arma::mat expected_fisher_ggm_gwishart(const arma::cube &samples,
                                       const arma::field<arma::uvec> &free_idx){
    arma::uword p    = samples.n_rows;
    arma::uword nsim = samples.n_slices;
    const arma::uvec& row_idx = free_idx(0);
    const arma::uvec& col_idx = free_idx(1);
    arma::uword n_free = row_idx.n_elem;

    // invert each sample slice (this function works on the inverses)
    arma::cube samples_inv(p, p, nsim);
    for (arma::uword s = 0; s < nsim; s++) {
        samples_inv.slice(s) = arma::inv_sympd(samples.slice(s));
    }

    // zero-copy view of the INVERTED cube as (p*p) x nsim: column s is vec(Sinv_s)
    const arma::mat Sinv_flat(const_cast<double*>(samples_inv.memptr()),
                             p * p, nsim, false, true);

    // second moment E[Sinv (x) Sinv] = (1/nsim) Sinv_flat Sinv_flat^T  — one BLAS call
    arma::mat M = (Sinv_flat * Sinv_flat.t()) / static_cast<double>(nsim);

    auto idx = [p](arma::uword r, arma::uword c) -> arma::uword { return r + c * p; };

    arma::mat F(n_free, n_free, arma::fill::zeros);

    for (arma::uword a = 0; a < n_free; a++) {
        arma::uword i = row_idx(a);
        arma::uword j = col_idx(a);
        bool diag_a   = (i == j);

        for (arma::uword b = a; b < n_free; b++) {
            arma::uword k = row_idx(b);
            arma::uword l = col_idx(b);
            bool diag_b   = (k == l);

            double val;
            if      ( diag_a &&  diag_b)
                val =        M(idx(i, k), idx(i, k));
            else if ( diag_a && !diag_b)
                val = 2.0 *  M(idx(i, k), idx(i, l));
            else if (!diag_a &&  diag_b)
                val = 2.0 *  M(idx(i, k), idx(j, k));
            else
                val = 2.0 * ( M(idx(i, k), idx(j, l)) +
                              M(idx(i, l), idx(j, k)) );

            F(a, b) = val;
            F(b, a) = val;
        }
    }

    F *= 0.5;
    return F;
}


