#include <cmath>
#include <string>
#include <RcppArmadillo.h>
#include "helpers.h"

// [[Rcpp::depends(RcppArmadillo)]]

// Helper functions to generate random numbers, multivariate normal samples using Rcpp so that R RNG is used for reproducibility.

// Generate a random normal matrix 
arma::mat rnorm_arma(int nrow, int ncol) {
  Rcpp::NumericVector v = Rcpp::rnorm(nrow * ncol);      
  return arma::mat(v.begin(), nrow, ncol);  // fill arma::mat columnwise
}

// Generate a random uniform vector of size n
arma::vec runif_arma(int n, double min, double max) {
  Rcpp::NumericVector v = Rcpp::runif(n, min, max);     
  return arma::vec(v.begin(), n);           // copy to arma::vec
}

// Generate n multivariate normal samples (X is n x p, each row is a sample)
arma::mat mvnrnd_arma(const arma::vec &mu, const arma::mat &Sigma, int n) {
  int p = mu.n_elem;
  arma::mat Z = rnorm_arma(p, n);   // each column is a sample
  arma::mat C = arma::chol(Sigma, "lower");
  arma::mat X = C * Z; // from standard to correlated normal
  X.each_col() += mu;
  return X.t();
}

// fast version: accepts precomputed Cholesky factor
arma::mat mvnrnd_chol(const arma::vec &mu, const arma::mat &C, int n) {
    int p = mu.n_elem;
    arma::mat Z = rnorm_arma(p, n);
    arma::mat X = C * Z;
    X.each_col() += mu;
    return X.t();
}

// Elimination matrix
arma::mat elimination_matrix(const int& p) {
    arma::mat out((p * (p + 1)) / 2, p * p, arma::fill::zeros);
    for (int j = 0; j < p; ++j) {
        arma::rowvec e_j(p, arma::fill::zeros);
        e_j(j) = 1.0;
        for (int i = j; i < p; ++i) {
            arma::vec u((p * (p + 1)) / 2, arma::fill::zeros);
            u(j * p + i - ((j + 1) * j) / 2) = 1.0;
            arma::rowvec e_i(p, arma::fill::zeros);
            e_i(i) = 1.0;
            out += arma::kron(u, arma::kron(e_j, e_i));
        }
    }
    return out;
}

// Commutation matrix
arma::mat commutation_matrix(const arma::uword& p) {
    arma::mat K(p * p, p * p, arma::fill::zeros);
    arma::uword i, j;
    for (i = 0; i < p; i++) {
        for (j = 0; j < p; j++) {
            K(i * p + j, j * p + i) = 1.0;
        }
    }
    return K;
}

// Duplication matrix
arma::mat duplication_matrix(const int& p) {
    arma::mat out((p * (p + 1)) / 2, p * p, arma::fill::zeros);
    for (int j = 0; j < p; ++j) {
        for (int i = j; i < p; ++i) {
            arma::vec u((p * (p + 1)) / 2, arma::fill::zeros);
            u(j * p + i - ((j + 1) * j) / 2) = 1.0;
            arma::mat T(p, p, arma::fill::zeros);
            T(i, j) = 1.0;
            T(j, i) = 1.0;
            out += u * arma::trans(arma::vectorise(T));
        }
    }
    return out; // [[NOTE]] it was out.t() but we always transposed at R level because we premultiply by the kronecker product (p*p, p*p)
}

// Inverse duplication matrix
arma::mat inverse_duplication_matrix(const int& p) {
    arma::mat D = duplication_matrix(p);
    arma::mat A = D * D.t();
    arma::mat D_plus = arma::solve(A, D);
    return D_plus;
}

// Select parameters (useful for G-Wishart to select only the parameters corresponding to the edges in the graph, including the diagonal)

// select parameters on full matrix (indices between 0 and p*p-1): select parameters corresponding to the edges in the graph, including the diagonal
arma::uvec select_parameter_full(arma::umat G){
    G.diag().ones(); // make sure the diagonal is included 
    arma::uvec select_parameter = arma::find(arma::trimatl(G)); // arma::find(E * arma::vectorise(G));
    return select_parameter;
}

// select parameters on lower triangular part (indices between 0 and (p*(p+1))/2 - 1): select parameters corresponding to the edges in the graph, including the diagonal
arma::uvec select_parameter_half(arma::umat G){
    G.diag().ones(); // make sure the diagonal is included 
    arma::uvec lower_idx = arma::trimatl_ind(arma::size(G));
    arma::uvec g = G(lower_idx);         // lower tri values in column-major order
    arma::uvec select_parameter = arma::find(g); // arma::find(E * arma::vectorise(G));
    return select_parameter;
} 

// Schur complement
arma::mat schur_complement(const arma::mat& X, const arma::uvec& i) {
    int npars = X.n_rows;
    arma::uvec idx = arma::regspace<arma::uvec>(0, npars - 1);
    idx.shed_rows(i);

    // extracting blocks
    arma::mat X_ii = X(i, i);
    arma::mat X_i_j = X.submat(i,idx); 
    arma::mat X_j_i = X.submat(idx,i); 
    arma::mat X_j_j = X.submat(idx, idx);

    // solve system X_j_j * x = X_j_i
    arma::mat x = arma::solve(X_j_j, X_j_i);

    // Schur complement: conditional precision for parameter i
    arma::mat X_cond = X_ii - (X_i_j * x);

    return X_cond;
}

// Global ESS
double global_ess(const arma::mat &X, 
                const arma::mat &Z,
                const std::string &aggregation){
    double ess;
    if (aggregation.compare("mean") == 0) {
        ess = arma::mean(X.diag() / Z.diag());
    } else if (aggregation.compare("tr") == 0) {
        ess = arma::trace(X) / arma::trace(Z);
    } else if (aggregation.compare("det") == 0) {
        // [improvement] if X and Z are positive definite, we can compute the log-determinant using Cholesky decomposition for numerical stability
        double logdetX, logdetZ;
        double signX, signZ;
        arma::log_det(logdetX, signX, X);
        arma::log_det(logdetZ, signZ, Z);
        ess = (signX / signZ) * std::exp((logdetX - logdetZ) / static_cast<double>(X.n_rows));
    } else {
        ess = -1.0; // invalid aggregation method
        Rcpp::stop("Invalid aggregation method. Choose 'mean', 'tr', or 'det'.");
    }
    return ess;
}

// Parameter-wise ESS
arma::vec parameterwise_ess(const arma::mat &X, 
                            const arma::mat &Z){
    arma::vec ess(X.n_rows, arma::fill::zeros);
    arma::mat cholX = arma::chol(X); // Cholesky decomposition of X
    arma::mat cholZ = arma::chol(Z); // Cholesky decomposition of Z
    ess = cholX.diag() / cholZ.diag(); // element-wise division of the diagonals
    ess %= ess; // square the values
    return ess;
}