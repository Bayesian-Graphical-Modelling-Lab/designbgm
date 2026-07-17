#ifndef HELPERS_H
#define HELPERS_H

#include <string>
#include <RcppArmadillo.h>

// Generate a random normal matrix 
arma::mat rnorm_arma(int nrow, int ncol);

// Generate random uniform vector of size n
arma::vec runif_arma(int n, double min = 0.0, double max = 1.0);

// Generate n multivariate normal samples
arma::mat mvnrnd_arma(const arma::vec &mu, const arma::mat &Sigma, int n);

// Generate n multivariate normal samples using precomputed Cholesky factor of the covariance matrix
arma::mat mvnrnd_chol(const arma::vec &mu, const arma::mat &C, int n);

// Elimination matrix
arma::mat elimination_matrix(const int& p);

// Commutation matrix
arma::mat commutation_matrix(const arma::uword& p);

// Duplication matrix
arma::mat duplication_matrix(const int& p);

// Inverse duplication matrix
arma::mat inverse_duplication_matrix(const int& p);

// Select parameters for full vectorization (of a symmetric matrix) corresponding to the edges in the graph, including the diagonal
arma::uvec select_parameter_full(arma::umat G);

// Select parameters for half vectorization (of a symmetric matrix) corresponding to the edges in the graph, including the diagonal
arma::uvec select_parameter_half(arma::umat G);

// Schur complement
arma::mat schur_complement(const arma::mat& X, const int& i);

// Global ESS
double global_ess(const arma::mat &X, const arma::mat &Z, const std::string &aggregation);

// Parameter-wise ESS
arma::vec parameterwise_ess(const arma::mat &X, const arma::mat &Z);  

#endif 