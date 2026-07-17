#include <string>
#include <RcppArmadillo.h>
#include "helpers.h"
#include "prior_wishart.h"

// [[Rcpp::depends(RcppArmadillo)]]

// Fast Wishart sampler using Bartlett decomposition
// K ~ Wishart(df, Psi) where Psi is the scale matrix
// C = chol(Psi, "upper")
arma::mat rwishart_fast(const double df, const arma::mat& C) {
    arma::uword p = C.n_rows;
    arma::mat psi(p, p, arma::fill::zeros);
    // diagonal: sqrt(chi^2(df - i)) = sqrt(Gamma((df - i)/2, 2))
    for (arma::uword i = 0; i < p; i++) {
        psi(i, i) = std::sqrt(R::rgamma((df - i) / 2.0, 2.0));
    }
    // upper triangle: standard normals
    for (arma::uword j = 1; j < p; j++) {
        for (arma::uword i = 0; i < j; i++) {
            psi(i, j) = R::rnorm(0.0, 1.0);
        }
    }
    // U = psi * C
    arma::mat U = arma::trimatu(psi) * C;
    // K = U.t() * U
    return U.t() * U;
}

// Wishart/inverse-Wishart Moments

// Variance for a Wishart prior predictive distribution (E is the elimination matrix and M is the commutation matrix)
arma::mat prior_variance_wishart(const arma::mat &K, const double &nu, const arma::mat &E, const arma::mat &M){
    arma::mat kron_K = arma::kron(K,K);
    arma::uword p  = K.n_cols * K.n_cols;
    arma::mat I_p = arma::eye<arma::mat>(p,p);   
    arma::mat V = E * (nu * ((I_p + M) * kron_K)) * E.t();
    return V;
}

// Second moment of inverse Wishart distribution (result of Von Rosen 1997) E[Theta^{-1} \otimes Theta^{-1}]
arma::mat second_moment_inverse_wishart(const arma::mat &Kinv, // inverse precision matrix
                    const double &nu,
                    const double &p,
                    const arma::mat &D, 
                    const arma::mat &M){
    double c2 = 1.0/((nu-p)*(nu-p-1.0)*(nu-p-3.0));
    double c1 = (nu-p-2.0)*c2;   
    arma::vec vec_Kinv = arma::vectorise(Kinv); 
    arma::mat kron_Kinv = arma::kron(Kinv,Kinv);   
    arma::mat V = D * (c1*kron_Kinv + c2*(vec_Kinv * vec_Kinv.t()) + c2*(M * kron_Kinv)) * D.t();
    return V;
}

