#ifndef PRIOR_GWISHART_H 
#define PRIOR_GWISHART_H 

#include <string>
#include <RcppArmadillo.h>

// Find neighbors of each node in the graph G
arma::field<arma::uvec> find_neighbors(const arma::umat& G);

// Create indices excluding i from {0,...,p-1}, for each i
arma::field<arma::uvec> create_indices_excluding_i(const arma::uword& p);

// Hastie adaptation step used by the G-Wishart sampler
double hastie_adaptation_step(arma::mat &W,
							const arma::mat &Sigma,
							const arma::uword &i,
							const arma::uword &p,
							const arma::field<arma::uvec> &indices_excluding_i,
							const arma::field<arma::uvec> &neighbors,
							arma::vec &beta_star_full);

// Iterative adaptation that returns an estimate of the precision matrix
arma::mat hastie_adaptation(arma::mat K,
							const arma::uword &p,
							const arma::field<arma::uvec> &indices_excluding_i,
							const arma::field<arma::uvec> &neighbors,
							double tol,
							arma::uword itermax,
							arma::vec &beta_star_full);

// Random G-Wishart direct sampler
arma::cube rgwishart_direct(const arma::uword &n,
					const arma::mat &K,
					const double &nu,
					const arma::umat &G,
					const double& tol = 1e-08,
					const arma::uword& itermax = 500);
								
// Random G-Wishart sampler via Gibbs sampling
arma::cube rgwishart_gibbs(const arma::uword& n,
							const arma::mat&  K,
							const double& nu,
							const arma::umat& G,
							const arma::uword& burnin = 500,
							const Rcpp::Nullable<arma::mat>& init = R_NilValue);									

// G-Wishart sampler dispatch function
arma::cube rgwishart(const arma::uword& n,
					const arma::mat& K,
					const double& nu,
					const arma::umat& G,
					const std::string& sampler,
					const double& tol = 1e-08,
					const arma::uword& itermax = 500,
					const arma::uword& burnin = 500,
					const Rcpp::Nullable<arma::mat>& init = R_NilValue);

// Variance of the G-Wishart distribution
arma::mat prior_variance_gwishart(const arma::umat &G, const arma::uvec &select_parameter, const arma::cube &samples);

#endif