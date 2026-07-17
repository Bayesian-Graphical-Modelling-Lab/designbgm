#ifndef FAMILY_GGM_H
#define FAMILY_GGM_H

#include <string>
#include <RcppArmadillo.h>

// Build free index for G-Wishart
arma::field<arma::uvec> build_free_idx(const arma::uvec& sel, const arma::uword p);

// Computes the free-parameter submatrix of D*(A⊗A)*D^T ,using the Magnus-Neudecker identity, with no need of Kronecker product. (used for hessian / fisher information kernels)
arma::mat kron_reduce_free_D(const arma::mat& A, const arma::field<arma::uvec>& free_idx);

// Computes the free-parameter submatrix of D_plus*(A⊗A)*D_plus^T, based on the kron_reduce_free_D code (used for inverse fisher information kernels)
arma::mat kron_reduce_free_Dplus(const arma::mat& A, const arma::field<arma::uvec>& free_idx);

// Convert precision matrix to partial correlations
arma::mat cpp_precision_to_partial_correlations(const arma::mat& K);

// Random precision matrix from prior
arma::mat random_precision_from_prior(const std::string& prior, 
                                        const arma::mat& K,  
                                        const int& nu,
                                        const arma::umat& G, 
                                        const arma::mat& Kchol);

// Simulate a random GGM study precision matrix with structure G
arma::mat cpp_simulate_ggm_study_precision(const int&p, 
                            const int& nu, 
                            arma::umat G, 
                            const bool& verbose = false);

// Expected (wrt the parameters) inverse Fisher information for GGMs under a Wishart prior predictive distribution
arma::mat expected_inverse_fisher_ggm_wishart(const arma::mat &K, const double &nu, const arma::mat &D_plus, const arma::mat &M);

// Expected (wrt the parameters) inverse Fisher information for GGMs under a G-Wishart prior predictive distribution
arma::mat expected_inverse_fisher_ggm_gwishart(const arma::cube &samples, const arma::field<arma::uvec> &free_idx);

// Expected (wrt the parameters) Fisher information for GGMs under a Wishart prior predictive distribution
arma::mat expected_fisher_ggm_wishart(const arma::mat &K, double nu, const arma::mat &D, const arma::mat &M);

// Expected (wrt the parameters) Fisher information for GGMs under a G-Wishart prior predictive distribution
arma::mat expected_fisher_ggm_gwishart(const arma::cube &samples, const arma::field<arma::uvec> &free_idx);

#endif
