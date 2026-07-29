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

// Constrain a precision matrix to a graph using Hastie adaptation
arma::mat cpp_constrain_precision_to_graph(const arma::mat& K,
										   const arma::umat& G,
										   double tol = 1e-6,
										   int itermax = 1000);
										   

// Double Continuous-Time (DCT) BDMCMC sampler (all edges at once, Rao-Blackwellized posterior inclusion probabilities)
struct bdmcmc_dct_result {
	arma::mat pip;      // posterior inclusion probabilities (Rao-Blackwellized)
	arma::mat pip_rb;   // Rao-Blackwellized posterior inclusion probabilities
	arma::mat K_hat;    // posterior mean precision matrix
	arma::mat last_graph; // last sampled graph
};

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
		std::string      gwish_sampler = "direct",       // "direct" or "gibbs"
		double           gwish_tol  = 1e-6,
		arma::uword      gwish_iter = 500,
		arma::uword      gwish_burnin = 500);		
		
// Double Conditional Bayes Factor (DCBF) BDMCMC sampler (one edge at a iteration)
struct bdmcmc_dcbf_result {
	arma::mat pip;      // posterior inclusion probabilities
	arma::mat K_hat;    // posterior mean precision matrix
	arma::umat last_graph; // last sampled graph
	double accept_rate; // acceptance rate of edge toggles
};

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
										std::string gwish_sampler = "direct",
										double gwish_tol = 1e-8, 
										arma::uword gwish_iter = 500,
										arma::uword gwish_burnin = 500);

#endif