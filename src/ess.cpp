#include <string>
#include <RcppArmadillo.h>
#include "helpers.h"
#include "prior_wishart.h"
#include "prior_gwishart.h"
#include "family_ggm.h"

// [[Rcpp::depends(RcppArmadillo)]]

// ---------------------------------------------- //
//         Prior ESS: family GGM                  //
// ---------------------------------------------- //

// Matrix-based GGM prior ESS estimators (VR, PR).
// Scalar estimators (PT, ELIR, MTM) are computed in R.
//
// which: "VR", "PR", or "both" — controls which ratios are calculated.
//        The expensive setup (samples, prior-variance) is shared regardless.
// [[Rcpp::export]]
Rcpp::List cpp_ggm_prior_ess_vr_pr(const std::string &prior,
                                   const std::string &which,
                                   const int &nu,
                                   const arma::mat &K,
                                   const arma::umat &G,
                                   const std::string &aggregation = "det",
                                   const std::string &sampler = "direct",
                                   const int &n_samples = 1000L,
                                   const double &tol = 1e-6,
                                   const arma::uword &itermax = 1000,
                                   const arma::uword &burnin = 500,
                                   const Rcpp::Nullable<arma::mat> &init = R_NilValue,
                                   const bool compute_cond = false) {
    Rcpp::List out;                                 
    arma::uword p = K.n_cols;

    const bool want_vr = (which == "VR" || which == "both");
    const bool want_pr = (which == "PR" || which == "both");

    if (!want_vr && !want_pr) {
        Rcpp::stop("Unknown 'which' value '" + which + "'. Use: VR, PR, or both.");
    }

    // utils
    arma::uvec select_half   = select_parameter_half(G);
    arma::uvec lower_indices = arma::trimatl_ind(arma::size(G));

    // partial correlations
    arma::mat Rho = cpp_precision_to_partial_correlations(K); 
    arma::vec rho_trimatl(p*(p+1)/2,arma::fill::zeros);
    rho_trimatl   = Rho(lower_indices);
    arma::vec rho = rho_trimatl(select_half);

    // precision parameters
    arma::vec pars_trimatl(p*(p+1)/2,arma::fill::zeros);
    pars_trimatl   = K(lower_indices);
    arma::vec pars = pars_trimatl(select_half);

    // is_diag: select diagonal parameters (1) and off-diagonal (free) parameters (0)
    arma::mat I_d = arma::diagmat(arma::ones<arma::vec>(p));
    arma::vec is_diag_trimatl(p*(p+1)/2,arma::fill::zeros);
    is_diag_trimatl   = I_d(lower_indices); 
    arma::vec is_diag = is_diag_trimatl(select_half); 

    // store parameters info
    arma::mat info(select_half.n_elem, 3);
    info.col(0) = is_diag; 
    info.col(1) = rho;
    info.col(2) = pars;
    out["info"] = info;

    // compute numerator and denominator of ESS estimator
    // VR:  X_VR = E[inv Fisher],      Z_VR = prior variance
    // PR:  X_PR = inv(prior variance), Z_PR = E[Fisher]
    // prior_variance is shared (Z_VR and X_PR both calculated from it).
    arma::mat X_VR, Z_VR, X_PR, Z_PR, prior_var;

    if (prior == "gwishart") {
        // build free index for G-Wishart and draw random samples
        arma::field<arma::uvec> free_idx = build_free_idx(select_half, p);
        arma::cube samples = rgwishart(n_samples, K, static_cast<double>(nu), G,
                                       sampler, tol, itermax, burnin, init);

        prior_var = prior_variance_gwishart(G, select_half, samples);   // shared

        if (want_vr) {
            X_VR = expected_inverse_fisher_ggm_gwishart(samples, free_idx);
            Z_VR = prior_var;
        }
        if (want_pr) {
            X_PR = arma::inv_sympd(prior_var);
            Z_PR = expected_fisher_ggm_gwishart(samples, free_idx);
        }
    } else if (prior == "wishart") {
        arma::mat E = elimination_matrix(p), M = commutation_matrix(p),
                  D = duplication_matrix(p), D_plus = inverse_duplication_matrix(p);

        prior_var = prior_variance_wishart(K, static_cast<double>(nu), E, M);          // shared

        if (want_vr) {
            X_VR = expected_inverse_fisher_ggm_wishart(K, static_cast<double>(nu), D_plus, M);
            Z_VR = prior_var;
        }
        if (want_pr) {
            X_PR = arma::inv_sympd(prior_var);
            Z_PR = expected_fisher_ggm_wishart(K, static_cast<double>(nu), D, M);
        }
    } else {
        Rcpp::stop("Unknown prior '" + prior + "'. Priors available: gwishart, wishart.");
    }

    // -------- assemble requested results (uniform sub-list shape) --------
    if (want_vr) {
        out["VR"] = Rcpp::List::create(
            Rcpp::Named("global")           = global_ess(X_VR, Z_VR, aggregation),
            Rcpp::Named("parameterwise")    = parameterwise_ess(X_VR, Z_VR),
            Rcpp::Named("cond_numerator")   = compute_cond ? arma::cond(X_VR) : NA_REAL,
            Rcpp::Named("cond_denominator") = compute_cond ? arma::cond(Z_VR) : NA_REAL);
    }
    if (want_pr) {
        out["PR"] = Rcpp::List::create(
            Rcpp::Named("global")           = global_ess(X_PR, Z_PR, aggregation),
            Rcpp::Named("parameterwise")    = parameterwise_ess(X_PR, Z_PR),
            Rcpp::Named("cond_numerator")   = compute_cond ? arma::cond(X_PR) : NA_REAL,
            Rcpp::Named("cond_denominator") = compute_cond ? arma::cond(Z_PR) : NA_REAL);
    }

    return out;
}
