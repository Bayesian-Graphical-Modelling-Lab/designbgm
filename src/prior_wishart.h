#ifndef PRIOR_WISHART_H 
#define PRIOR_WISHART_H

#include <string>
#include <RcppArmadillo.h>

// Generate a Wishart random matrix using Bartlett decomposition
arma::mat rwishart_fast(const double df, const arma::mat& C);

// Variance of the Wishart distribution
arma::mat prior_variance_wishart(const arma::mat &K, const double &nu, const arma::mat &E, const arma::mat &M);

// Second moment of the inverse Wishart distribution (result of Von Rosen 1997) E[Theta^{-1} \otimes Theta^{-1}]
arma::mat second_moment_inverse_wishart(const arma::mat &Kinv,
                                        const double &nu,
                                        const double &p,
                                        const arma::mat &D,
                                        const arma::mat &M);

#endif
