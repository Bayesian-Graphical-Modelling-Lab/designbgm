#include <string>
#include <cmath>
#include <RcppArmadillo.h>
#include "helpers.h"
#include "family_ggm.h"
#include "prior_gwishart.h"

// [[Rcpp::depends(RcppArmadillo)]]

// Design Data-to-Prior Information Ratio (DPIR; global and parameterwise) for GGMs under a Wishart or G-Wishart prior predictive distribution
// [[Rcpp::export]]
Rcpp::List cpp_design_dpir(const std::string &prior, const arma::mat &K, const arma::uword &nu, // K is already the elicited prior precision matrix as K_prior/nu
                        const arma::umat &G, const arma::uword &H, 
                        const arma::uword &J, arma::uvec &n, const double &threshold = 1.0, 
                        const bool optimize = false, const double &target_probability = 0.95, 
                        const arma::uword &n_tol = 1, const arma::uword &max_n = 10000, 
                        const std::string &gwish_sampler = "direct", const double &gwish_tol = 1e-08, 
                        const int &gwish_iter = 500, const int &gwish_burnin = 500) {

    Rcpp::List out;
    arma::uword p = K.n_cols;
    double scale = static_cast<double>(nu - p - 1.0);
    arma::uvec select_half = select_parameter_half(G);
    arma::uword n_pars = select_half.n_elem;
    arma::uvec lower_indices = arma::trimatl_ind(arma::size(G));
    arma::mat Kchol = (prior == "wishart") ? arma::chol(K, "upper") : arma::mat(); // precompute Cholesky for Wishart prior

    arma::field<arma::uvec> free_idx = build_free_idx(select_half, p); // build free index for G-Wishart (used in kron_reduce_free_D)

    // info - diagonal elements of the prior information matrix 
    arma::mat Kinv = arma::inv_sympd(K);
    arma::mat prior_info_mat = 0.5 * scale * kron_reduce_free_D(Kinv, free_idx); 
    arma::vec prior_info = prior_info_mat.diag();

    // pars - precision parameters (vectorized lower triangular)
    arma::vec pars_trimatl(p*(p+1)/2,arma::fill::zeros);
    pars_trimatl = K(lower_indices);
    arma::vec pars = pars_trimatl(select_half);

    // is_diag : select diagonal parameters (1) and off-diagonal parameters (0)
    arma::mat I_d = arma::diagmat(arma::ones<arma::vec>(p));
    arma::uvec is_diag_trimatl(p*(p+1)/2,arma::fill::zeros);
    is_diag_trimatl = arma::conv_to<arma::uvec>::from(I_d(lower_indices));
    arma::uvec is_diag = is_diag_trimatl(select_half); 
    arma::uvec off_diag_idx = arma::find(is_diag == 0); // complementary index for off-diagonal parameters

    out["is_diag"] = is_diag;
    out["pars"] = pars;
    out["info"] = prior_info;

    if (optimize) {
        arma::uword n_grid_size = n.n_elem;
        double cur_prob_global = 0.0;
        arma::vec cur_prob_pw(n_pars, arma::fill::zeros);
        arma::vec mu(p, arma::fill::zeros);

        // lambda to evaluate Pr(DPIR > threshold) at a given n
        auto global_eval_at_n = [&](arma::uword n_val, double &prob_global) {
            prob_global = 0.0;
            const double inv_HJ = 1.0 / static_cast<double>(H * J);
            for (arma::uword h = 0; h < H; h++) {
                arma::mat K_h = random_precision_from_prior(prior, K, nu, G, Kchol, gwish_sampler, gwish_tol, gwish_iter, gwish_burnin);
                arma::mat S_h = arma::inv_sympd(K_h);
                arma::mat cholS_h = arma::chol(S_h, "lower");
                // arma::mat prior_info_full = 0.5 * scale * (D * arma::kron(S_h, S_h) * D.t());
                // arma::mat prior_info = prior_info_full(select_half, select_half);
                arma::mat prior_info = 0.5 * scale * kron_reduce_free_D(S_h, free_idx);
                arma::mat chol_prior = arma::chol(prior_info, "lower");
                arma::vec chol_prior_diag_sq = arma::square(chol_prior.diag());
                double logdet_prior = arma::sum(arma::log(chol_prior_diag_sq));
                for (arma::uword j = 0; j < J; j++) {
                    arma::mat X_j = mvnrnd_chol(mu, cholS_h, n_val);
                    arma::mat S_j = (X_j.t() * X_j) / static_cast<double>(n_val);
                    //arma::mat data_info_full  = 0.5 * static_cast<double>(n_val) * (D * arma::kron(S_j, S_j) * D.t()); // on all parameters
                    //arma::mat data_info = data_info_full(select_half, select_half); // only on free parameters
                    arma::mat data_info = 0.5 * static_cast<double>(n_val) * kron_reduce_free_D(S_j, free_idx); 
                    arma::mat chol_data      = arma::chol(data_info, "lower");
                    arma::vec chol_data_diag_sq = arma::square(chol_data.diag());
                    double logdet_data  = arma::sum(arma::log(chol_data_diag_sq)); // log det from Cholesky
                    double determinant_ratio_value = std::exp((logdet_data - logdet_prior) / static_cast<double>(data_info.n_rows));
                    if (determinant_ratio_value > threshold) prob_global += inv_HJ; // prob_global update
                }
            }
        };

        // lambda to evaluate Pr(DPIR > threshold) at a given n
        auto pw_eval_at_n = [&](arma::uword n_val, arma::vec &prob_pw) {
            prob_pw.zeros();
            const double inv_HJ = 1.0 / static_cast<double>(H * J);
            for (arma::uword h = 0; h < H; h++) {
                arma::mat K_h = random_precision_from_prior(prior, K, nu, G, Kchol, gwish_sampler, gwish_tol, gwish_iter, gwish_burnin);
                arma::mat S_h = arma::inv_sympd(K_h);
                arma::mat cholS_h = arma::chol(S_h, "lower");
                // arma::mat prior_info_full = 0.5 * scale * (D * arma::kron(S_h, S_h) * D.t());
                // arma::mat prior_info = prior_info_full(select_half, select_half);
                arma::mat prior_info = 0.5 * scale * kron_reduce_free_D(S_h, free_idx);
                arma::mat chol_prior = arma::chol(prior_info, "lower");
                arma::vec chol_prior_diag_sq = arma::square(chol_prior.diag());
                for (arma::uword j = 0; j < J; j++) {
                    arma::mat X_j = mvnrnd_chol(mu, cholS_h, n_val);
                    arma::mat S_j = (X_j.t() * X_j) / static_cast<double>(n_val);
                    //arma::mat data_info_full  = 0.5 * static_cast<double>(n_val) * (D * arma::kron(S_j, S_j) * D.t()); // on all parameters
                    //arma::mat data_info = data_info_full(select_half, select_half); // only on free parameters
                    arma::mat data_info = 0.5 * static_cast<double>(n_val) * kron_reduce_free_D(S_j, free_idx); 
                    arma::mat chol_data      = arma::chol(data_info, "lower");
                    arma::vec chol_data_diag_sq = arma::square(chol_data.diag());
                    arma::vec ratio = chol_data_diag_sq / chol_prior_diag_sq;
                    for (arma::uword k = 0; k < n_pars; k++) {
                        if (ratio(k) > threshold) prob_pw(k) += inv_HJ; // prob_pw update
                    }
                }
            }
        };

        // lambda to evaluate Pr(DPIR > threshold) at a given n
        auto both_eval_at_n = [&](arma::uword n_val, double &prob_global, arma::vec &prob_pw) {
            prob_global = 0.0;
            prob_pw.zeros();
            const double inv_HJ = 1.0 / static_cast<double>(H * J);
            for (arma::uword h = 0; h < H; h++) {
                arma::mat K_h = random_precision_from_prior(prior, K, nu, G, Kchol, gwish_sampler, gwish_tol, gwish_iter, gwish_burnin);
                arma::mat S_h = arma::inv_sympd(K_h);
                arma::mat cholS_h = arma::chol(S_h, "lower");
                // arma::mat prior_info_full = 0.5 * scale * (D * arma::kron(S_h, S_h) * D.t());
                // arma::mat prior_info = prior_info_full(select_half, select_half);
                arma::mat prior_info = 0.5 * scale * kron_reduce_free_D(S_h, free_idx);
                arma::mat chol_prior = arma::chol(prior_info, "lower");
                arma::vec chol_prior_diag_sq = arma::square(chol_prior.diag());
                double logdet_prior = arma::sum(arma::log(chol_prior_diag_sq));
                for (arma::uword j = 0; j < J; j++) {
                    arma::mat X_j = mvnrnd_chol(mu, cholS_h, n_val);
                    arma::mat S_j = (X_j.t() * X_j) / static_cast<double>(n_val);
                    //arma::mat data_info_full  = 0.5 * static_cast<double>(n_val) * (D * arma::kron(S_j, S_j) * D.t()); // on all parameters
                    //arma::mat data_info = data_info_full(select_half, select_half); // only on free parameters
                    arma::mat data_info = 0.5 * static_cast<double>(n_val) * kron_reduce_free_D(S_j, free_idx); 
                    arma::mat chol_data      = arma::chol(data_info, "lower");
                    arma::vec chol_data_diag_sq = arma::square(chol_data.diag());
                    double logdet_data  = arma::sum(arma::log(chol_data_diag_sq)); // log det from Cholesky
                    double determinant_ratio_value = std::exp((logdet_data - logdet_prior) / static_cast<double>(data_info.n_rows));
                    if (determinant_ratio_value > threshold) prob_global += inv_HJ; // prob_global update
                    arma::vec ratio = chol_data_diag_sq / chol_prior_diag_sq;
                    for (arma::uword k = 0; k < n_pars; k++) {
                        if (ratio(k) > threshold) prob_pw(k) += inv_HJ; // prob_pw update
                    }
                }
            }
        };

        auto find_upper_bracket = [&](arma::uword n_start, double target_probability, arma::uword &n_left_global, 
                                        arma::uword &n_right_global, arma::uword &n_left_pw, arma::uword &n_right_pw, 
                                        bool &global_bracketed, bool &pw_bracketed, const arma::uvec off_diag_idx) {
            arma::uword n_try = n_start;
            double p_gl = 0.0;
            arma::vec p_pw(n_pars, arma::fill::zeros);
            for (int attempt = 0; attempt < 15; attempt++) { // capping max doublings at 15
                if (n_try > max_n) {
                    Rcpp::warning("Upper bracket exceeds max_n");
                    break;
                }
                both_eval_at_n(n_try, p_gl, p_pw);
                if (p_gl >= target_probability - 1e-6){
                    n_right_global = n_try;
                    global_bracketed = true;
                }
                else{
                    n_left_global = n_try;
                }
                if (arma::min(p_pw(off_diag_idx)) >= target_probability - 1e-6){
                    n_right_pw = n_try;
                    pw_bracketed = true;
                }
                else{
                    n_left_pw = n_try;
                }
                if(global_bracketed && pw_bracketed) break; // stop if both brackets found
                else n_try *= 2;
            }
            if (!global_bracketed || !pw_bracketed) {
                Rcpp::warning("Upper bracket not found after 15 doublings");
            }
        };


        // --- phase 1: check current max n grid value: if both pw and global are bracketed continue, otherwise find upper bracket by doubling ---
        arma::uword n_left_global = n(0);
        arma::uword n_right_global = n(n_grid_size - 1);
        arma::uword n_left_pw = n(0);
        arma::uword n_right_pw = n(n_grid_size - 1);
        bool global_bracketed = false;
        bool pw_bracketed = false;

        // evaluate dpir at the right supplied n_grid value
        both_eval_at_n(n(n_grid_size - 1), cur_prob_global, cur_prob_pw);
        if(cur_prob_global >= target_probability - 1e-6) global_bracketed = true;
        if(arma::min(cur_prob_pw(off_diag_idx)) >= target_probability - 1e-6) pw_bracketed = true;

        // if target not reached at max n, find upper bracket by doubling from max n grid value
        if (!global_bracketed || !pw_bracketed) {
            find_upper_bracket(n(n_grid_size - 1), target_probability, n_left_global, n_right_global, n_left_pw, n_right_pw, global_bracketed, pw_bracketed, off_diag_idx);
        }
        
        // warnings if target never reached
        if (!global_bracketed) {
            Rcpp::warning("target_probability never reached for global DPIR within the provided n grid");
        }
        if (!pw_bracketed) {
            Rcpp::warning("target_probability never reached for the weakest off-diagonal parameter within the provided n grid");
        }

        // initialize timer to measure time for each bisection
        arma::wall_clock timer; 

        // --- phase 2: bisection within bracket ---
        double prob_global_at_n_star = 0.0;  
        arma::vec prob_pw_at_n_star(n_pars, arma::fill::zeros); 

        // (1) Find global bisection 
        timer.tic();
        if (global_bracketed) {
            
            while (n_right_global - n_left_global > n_tol) {
                if (n_left_global > max_n) break; // early stopping if left bracket exceeds max_n
                arma::uword n_mid = (n_left_global + n_right_global) / 2;
                global_eval_at_n(n_mid, cur_prob_global);
                if (cur_prob_global >= target_probability - 1e-6) {
                    n_right_global = n_mid;
                    prob_global_at_n_star = cur_prob_global; // store probability at current n_mid for output
                }
                else n_left_global  = n_mid;
            }
            
        }
        double global_bisection_time = timer.toc();

        // (2) Find parameterwise bisection (based on min probability criterion, this ensures all parameters meet the target probability and we find n* as max n* across off-diag parameters)
        timer.tic();
        if(pw_bracketed) {
            while(n_right_pw - n_left_pw > n_tol) {
                if (n_left_global > max_n) break; // early stopping if left bracket exceeds max_n
                arma::uword n_mid = (n_left_pw + n_right_pw) / 2;
                pw_eval_at_n(n_mid, cur_prob_pw);
                if(arma::min(cur_prob_pw(off_diag_idx)) >= target_probability - 1e-6) {
                    n_right_pw = n_mid;
                    prob_pw_at_n_star = cur_prob_pw; // store probability at current n_mid for output
                }
                else n_left_pw  = n_mid;
            }  
        }
        double pw_bisection_time = timer.toc();
        
        // return n* star and prob at the corresponding n* (n_star_global for global and n_star_pw for pw)
        bool converged_global = global_bracketed && (prob_global_at_n_star >= target_probability - 1e-6) && (n_right_global <= max_n);
        if(converged_global){
            out["n_star_global"]            = n_right_global;
            out["prob_global_at_n_star"]    = prob_global_at_n_star;
        }

        bool converged_pw = pw_bracketed && (arma::min(prob_pw_at_n_star(off_diag_idx)) >= target_probability - 1e-6) && (n_right_pw <= max_n);
        if(converged_pw){
            out["n_star_pw"]     = n_right_pw;
            out["prob_pw_at_n_star"]        = prob_pw_at_n_star;
        }

        out["duration_global_bisection"] = global_bisection_time;
        out["duration_pw_bisection"]     = pw_bisection_time;
        out["converged_global"]          = converged_global;
        out["converged_pw"]              = converged_pw;
    } 
    else {
        arma::vec global_dpir_prob(n.n_elem, arma::fill::zeros);
        arma::mat pw_dpir_prob(n_pars, n.n_elem, arma::fill::zeros);
        arma::vec mu(p, arma::fill::zeros);
        double inv_HJ = 1.0 / static_cast<double>(H * J);
        // write here grid n evalutation of dpir over H x J replications
        arma::mat K_h(p, p, arma::fill::zeros), 
                  S_h(p, p, arma::fill::zeros),
                  S_j(p, p, arma::fill::zeros),
                  cholS_h(p, p, arma::fill::zeros),
                  chol_prior(n_pars, n_pars, arma::fill::zeros),
                  chol_data(n_pars, n_pars, arma::fill::zeros),
                  data_info(n_pars, n_pars, arma::fill::zeros),
                  prior_info(n_pars, n_pars, arma::fill::zeros);
                  //prior_info_full(p*(p+1)/2, p*(p+1)/2, arma::fill::zeros),
                  //data_info_full(p*(p+1)/2, p*(p+1)/2, arma::fill::zeros);
        arma::vec chol_prior_diag_sq(n_pars, arma::fill::zeros),
                  chol_data_diag_sq(n_pars, arma::fill::zeros),
                  pw_ratio(n_pars, arma::fill::zeros);
        double logdet_prior, 
                logdet_data, 
                determinant_ratio_value;

        for (arma::uword h = 0; h < H; h++) {
                K_h = random_precision_from_prior(prior, K, nu, G, Kchol, gwish_sampler, gwish_tol, gwish_iter, gwish_burnin); // random draw from the prior distribution
                S_h = arma::inv_sympd(K_h); // convert random precision draw to random covariance matrix
                cholS_h = arma::chol(S_h, "lower");
                //prior_info_full = 0.5 * scale * (D * arma::kron(S_h,S_h) * D.t()); // full Fisher information of the prior (all parameters)
                //prior_info = prior_info_full(select_half, select_half); // Fisher information of the prior (only free parameters)
                prior_info = 0.5 * scale * kron_reduce_free_D(S_h, free_idx);
                chol_prior = arma::chol(prior_info, "lower");
                chol_prior_diag_sq = arma::square(chol_prior.diag());
                logdet_prior = arma::sum(arma::log(chol_prior_diag_sq));
            for (arma::uword j = 0; j < J; j++) {
                for(arma::uword i = 0; i < n.n_elem; i++){
                    arma::mat X_j = mvnrnd_chol(mu, cholS_h, n(i)); // simulate data from the random covariance matrix
                    S_j = (X_j.t() * X_j) / static_cast<double>(n(i)); // sample covariance matrix for the first n(i) samples of the j-th replication
                    //data_info_full =  0.5 * static_cast<double>(n(i)) * (D * arma::kron(S_j,S_j) * D.t()); // Fisher information of the data on all parameters
                    //data_info = data_info_full(select_half, select_half); // Fisher information of the data (only free parameters)
                    data_info = 0.5 * static_cast<double>(n(i)) * kron_reduce_free_D(S_j, free_idx);
                    chol_data         = arma::chol(data_info, "lower");
                    chol_data_diag_sq = arma::square(chol_data.diag());
                    logdet_data       = arma::sum(arma::log(chol_data_diag_sq)); // log det from Cholesky
                    determinant_ratio_value  = std::exp((logdet_data - logdet_prior) / static_cast<double>(data_info.n_rows));
                    if (determinant_ratio_value > threshold) global_dpir_prob(i) += inv_HJ; // global DPIR probability update for n(i)
                    pw_ratio = chol_data_diag_sq / chol_prior_diag_sq;
                    for (arma::uword k = 0; k < n_pars; k++) {
                        if (pw_ratio(k) > threshold) pw_dpir_prob(k,i) += inv_HJ; // parameterwise DPIR probability update for n(i)
                    }  
                }     
            }
        }

        // store results in output list
        out["global_dpir_prob"] = global_dpir_prob;
        out["pw_dpir_prob"] = pw_dpir_prob;
    
    }

    return out;
}


// Utility functions for BFDA -- Complete graphs 

// constant factor for conditional Bayes factor (used in compute_conditional_bf)
arma::vec get_const_conditional_bf(double nu, const arma::vec &n) {
    arma::vec x = std::lgamma(nu / 2.0) +
                  std::lgamma((nu - 1.0) / 2.0) +
                  2.0 * arma::lgamma((n + nu + 1.0) / 2.0) -
                  arma::lgamma((n + nu) / 2.0) -
                  arma::lgamma((n + nu - 1.0) / 2.0) -
                  2.0 * std::lgamma((nu + 1.0) / 2.0);
    return arma::exp(x);
}

// compute conditional Bayes factor for edge {i,j} given the rest of the graph (Giudici, 1995)
double compute_conditional_bf(const arma::mat &nS, const arma::mat &K, const arma::uword &p, double nu, 
                      arma::uword i, arma::uword j, double n, double cons) {

    // indices for {i,j} and complement b
    arma::uvec a = {i, j};
    arma::uvec b = arma::regspace<arma::uvec>(0, p - 1);
    b.shed_rows(a); 

    // (1) compute T_ij = (K[a,a])^{-1}
    arma::mat Tij = arma::inv_sympd(K.submat(a, a));

    // (2) compute sample partial deviance matrix SA_ij
    arma::mat nS_aa = nS.submat(a, a);
    arma::mat nS_ab = nS.submat(a, b);
    arma::mat nS_bb = nS.submat(b, b);
    arma::mat nS_bb_reg = nS_bb + 1e-8 * arma::eye(arma::size(nS_bb)); // temporary, for numerical stability in case nS_bb is near-singular
    arma::mat SAij  = nS_aa - nS_ab * arma::solve(nS_bb_reg, nS_ab.t()); // arma::solve used for numerical stability instead of arma::inv_sympd(nS_bb)

    // (3) compute Bayes Factor
    double rp  = Tij(0,1) / std::sqrt(Tij(0,0) * Tij(1,1)); // prior partial correlation
    double rf  = (Tij(0,1) + SAij(0,1)) / std::sqrt((Tij(0,0) + SAij(0,0)) * (Tij(1,1) + SAij(1,1))); // posterior partial correlation
    double num   = std::pow(1.0 - rf * rf, (nu + n) / 2.0);
    double denom = std::pow(1.0 - rp * rp, nu / 2.0);
    double BF = cons * (num / denom) * std::sqrt((Tij(0,0) * Tij(1,1)) / 
                ((SAij(0,0) + Tij(0,0)) * (SAij(1,1) + Tij(1,1))));

    return BF;
}

// Frequentist warm start: sample size for partial correlation using Fisher z ("Statistical Power Analysis for the Behavioral Sciences", 2nd ed., Cohen,1988)
arma::uword warm_start_n(double rho, arma::uword p, 
                          double alpha = 0.05, double power = 0.80) {
    if (std::abs(rho) < 1e-6) return 200; // fallback for zero partial correlation
    double z_rho   = std::atanh(std::abs(rho));
    double z_alpha = R::qnorm(1.0 - alpha / 2.0, 0, 1, 1, 0); // two-tailed
    double z_beta  = R::qnorm(power, 0, 1, 1, 0);
    double q       = static_cast<double>(p - 2);  // variables partialled out
    return static_cast<arma::uword>(
        std::ceil(std::pow((z_alpha + z_beta) / z_rho, 2.0) + q + 3.0)
    );
}

// BFDA for complete graphs: sample size planning for a single edge in a complete graph using conditional Bayes factor
// [[Rcpp::export]]
Rcpp::List cpp_design_bfda_edge_dense(const arma::mat &K, const arma::uword &nu, arma::umat &G, // K is already the elicited prior precision matrix as K_prior/nu
                                    const arma::uword &m, const arma::uword &l,
                                    const arma::uword &H, const arma::uword &J, arma::uvec &n,
                                    const double &pow0 = 0.8, const double &pow1 = 0.8,
                                    const double &threshold = 10.0, const bool optimize = false, 
                                    const arma::uword &n_tol = 1, const arma::uword &max_n = 10000, 
                                    const std::string &gwish_sampler = "direct", const double &gwish_tol = 1e-08, 
                                    const arma::uword &gwish_iter = 500, const arma::uword &gwish_burnin = 500) {

    Rcpp::List out;
    arma::uword p = K.n_cols;
    arma::vec mu(p, arma::fill::zeros);
    arma::mat Rho = cpp_precision_to_partial_correlations(K);
    std::string prior_h0 = "gwishart"; // prior for hypothesis 0 (edge absent)
    double density_G = arma::accu(G) / static_cast<double>(p * (p - 1)); // density of the input graph G
    std::string prior_h1 = (density_G < 1) ? "gwishart" : "wishart"; //std::string prior_h1 = "wishart";   // prior for hypothesis 1
 
    arma::mat Kchol = arma::chol(K, "upper"); // precompute Cholesky for Wishart prior

    if (optimize) {
        if (G(m,l) == 0) {
            Rcpp::stop("optimize = TRUE is only valid for prior present edges (G(m,l) = 1). Use optimize = FALSE for absent edges.");
        }

        auto eval_edge_at_n = [&](arma::uword n_val, arma::uword hypothesis,
                                double &cur_power, double &cur_error) {
            arma::vec n_vec = {static_cast<double>(n_val)};
            double cons = get_const_conditional_bf(static_cast<double>(nu), n_vec)(0);    
            arma::vec bf_vals(H * J, arma::fill::zeros);
            arma::umat G_temp = G;
            std::string prior = (hypothesis == 0) ? prior_h0 : prior_h1;
            if (hypothesis == 0) {
                G_temp(l,m) = G_temp(m,l) = 0;
            }
            arma::uword idx = 0;
            for (arma::uword h = 0; h < H; h++) {
                arma::mat K_h     = random_precision_from_prior(prior, K, nu, G_temp, Kchol, gwish_sampler, gwish_tol, gwish_iter, gwish_burnin);
                arma::mat S_h     = arma::inv_sympd(K_h);
                arma::mat cholS_h = arma::chol(S_h, "lower");
                for (arma::uword j = 0; j < J; j++) {
                    arma::mat X_j  = mvnrnd_chol(mu, cholS_h, n_val);
                    arma::mat nS_j = X_j.t() * X_j;
                    double bf_01   = compute_conditional_bf(nS_j, K, p, static_cast<double>(nu), m, l, static_cast<double>(n_val), cons);
                    bf_vals(idx++) = bf_01;
                }
            }
            double hj = static_cast<double>(H * J);
            if (hypothesis == 0) {
                cur_power = arma::accu(bf_vals > threshold)       / hj;
                cur_error = arma::accu(bf_vals < 1.0 / threshold) / hj;
            } else {
                cur_power = arma::accu(bf_vals < 1.0 / threshold) / hj;
                cur_error = arma::accu(bf_vals > threshold)        / hj;
            }
        };

        double rho_ml    = Rho(m, l);
        arma::uword n_min = static_cast<arma::uword>(p + 2);
        arma::wall_clock timer; // to measure time for each hypothesis

        
        // --- hypothesis 0: bisect power_h0 ---
        double cur_power = 0.0, 
                cur_error = 1.0, 
                pow0_at_n_star = 0.0, 
                pow1_at_n_star = 0.0, 
                fnr_at_n_star = 1.0, 
                fpr_at_n_star = 1.0;
        arma::uword nl_pow0      = n_min; // because as the prior becomes more informative, the n* decreases and can also be lower than nu
        arma::uword nr_pow0      = warm_start_n(0.01, p, 0.05, pow0); // define upper bracket as the n* for a partial correlation as low as 0.01

        timer.tic();
        if (nr_pow0 > nl_pow0) {
            while (nr_pow0 - nl_pow0 > n_tol) {
                if (nl_pow0 > max_n) break;
                arma::uword n_mid = (nl_pow0 + nr_pow0) / 2;
                eval_edge_at_n(n_mid, 0, cur_power, cur_error);
                if (cur_power >= pow0) {
                    nr_pow0 = n_mid;
                    pow0_at_n_star = cur_power; 
                    fpr_at_n_star = cur_error; 
                }
                else nl_pow0 = n_mid;
            }
        }
        double t_h0 = timer.toc();

        bool converged_h0 = (pow0_at_n_star >= pow0) && (nr_pow0 <= max_n);
        if(converged_h0) {
            out["n_star_power_h0"] = nr_pow0;
            out["fpr_at_n_star_power_h0"] = fpr_at_n_star;
            out["power_h0_at_n_star"] = pow0_at_n_star;
        }
        out["duration_h0"]  = t_h0;
        out["converged_h0"] = converged_h0;
        

        // --- hypothesis 1: bisect power_h1 ---
        cur_power = 0.0; cur_error = 1.0;
        arma::uword nl_pow1      = n_min; // because as the prior becomes more informative, the n* decreases and can also be lower than nu
        arma::uword nr_pow1      = warm_start_n(rho_ml, p, 0.05, pow1) * 10; // define upper bracket as the n* for the selected partial correlation and inflate it by 10 
        
        timer.tic();
        if (nr_pow1 > nl_pow1) {
            while (nr_pow1 - nl_pow1 > n_tol) {
                if (nl_pow1 > max_n) break;
                arma::uword n_mid = (nl_pow1 + nr_pow1) / 2;
                eval_edge_at_n(n_mid, 1, cur_power, cur_error);
                if (cur_power >= pow1) {
                    nr_pow1 = n_mid;
                    pow1_at_n_star = cur_power; 
                    fnr_at_n_star = cur_error; 
                }
                else nl_pow1 = n_mid;
            }
           
        }
        double t_h1 = timer.toc();

        // check convergence and store results 
        bool converged_h1 = (pow1_at_n_star >= pow1) && (nr_pow1 <= max_n);
        if(converged_h1) {
            out["n_star_power_h1"] = nr_pow1;
            out["fnr_at_n_star_power_h1"] = fnr_at_n_star;
            out["power_h1_at_n_star"] = pow1_at_n_star;
        }
        out["duration_h1"]  = t_h1;
        out["converged_h1"] = converged_h1;

    } else {
        arma::vec const_bf = get_const_conditional_bf(nu, arma::conv_to<arma::vec>::from(n));
        arma::mat power_h0(n.n_elem, 1, arma::fill::zeros),
                  fpr_h0(n.n_elem,  1, arma::fill::zeros),
                  power_h1(n.n_elem, 1, arma::fill::zeros),
                  fnr_h1(n.n_elem,  1, arma::fill::zeros);

        for (arma::uword hypothesis = 0; hypothesis < 2; hypothesis++) {
            arma::umat G_temp = G;
            if (hypothesis == 0) {
                G_temp(l,m) = G_temp(m,l) = 0; // H0: edge absent
            }
            std::string prior = (hypothesis == 0) ? prior_h0 : prior_h1;
            arma::mat bf_vals(n.n_elem, H * J, arma::fill::zeros);
            arma::uword counter_iter = 0;
            for (arma::uword h = 0; h < H; h++) {
                arma::mat K_h    = random_precision_from_prior(prior, K, nu, G_temp, Kchol, gwish_sampler, gwish_tol, gwish_iter, gwish_burnin); //rgwishart(1, K, nu, G_temp, 1e-08, 500).slice(0);
                arma::mat S_h    = arma::inv_sympd(K_h);
                arma::mat cholSh = arma::chol(S_h, "lower");
                for (arma::uword j = 0; j < J; j++) {
                    for (arma::uword i = 0; i < n.n_elem; i++) {
                        arma::mat X_j  = mvnrnd_chol(mu, cholSh, n(i));
                        arma::mat nS_j = X_j.t() * X_j;
                        double bf_01   = compute_conditional_bf(nS_j, K, p, static_cast<double>(nu), m, l, static_cast<double>(n(i)), const_bf(i));
                        bf_vals(i, counter_iter) = bf_01; // G(m,l) is 1 by default as we are working with a full precision matrix K (all edges are present) and testing for exclusion of edge (m,l) 
                    }
                    counter_iter++;
                }
            }
            double hj = static_cast<double>(H * J);
            if (hypothesis == 0) {
                for (arma::uword i = 0; i < n.n_elem; i++) {
                    power_h0(i, 0) = arma::accu(bf_vals.row(i) > threshold)       / hj;
                    fpr_h0(i, 0)   = arma::accu(bf_vals.row(i) < 1.0 / threshold) / hj;
                }
                out["bf_h0"] = bf_vals;
            } else {
                for (arma::uword i = 0; i < n.n_elem; i++) {
                    power_h1(i, 0) = arma::accu(bf_vals.row(i) < 1.0 / threshold) / hj;
                    fnr_h1(i, 0)   = arma::accu(bf_vals.row(i) > threshold)       / hj;
                }
                out["bf_h1"] = bf_vals;
            }
        }
        out["power_h0"] = power_h0;
        out["fpr_h0"]   = fpr_h0;
        out["power_h1"] = power_h1;
        out["fnr_h1"]   = fnr_h1;
    }

    out["m"] = m;
    out["l"] = l;

    return out;
}


// Atay-Kayis computation of the normalizing constant for G-Wishart distribution with graph G and scale matrix K, using Monte Carlo integration
// Reference: Atay-Kayis and Massam (2005), "A Monte Carlo method for computing the marginal likelihood in nondecomposable Gaussian graphical models", Biometrika, 92(2):317-335
//
// Notation follows Atay-Kayis & Massam (2005):
//   d       : degrees of freedom of W_G(d, D)
//   T       : upper triangular Cholesky of D^{-1}, D^{-1} = T'T
//   t_{js}] : T(j,s)/T(s,s)  [equation 28]  stored in t_norm(j,s)
//   n_i     : upper-triangular degree of node i in G  (rowSums of upper-tri adj)
//   V'      : non-free pairs {(i,j): i<j, (i,j) not in E of G}
//   y_{ii}  ~ sqrt(chi^2_{d+n_i})   [equation 41]
//   y_{ij}  ~ N(0,1) for (i,j) in E [equation 42]
//   y_{ij}  : non-free, computed via equations (31)-(32)
//
// Key point:
//   G and G_{-e} differ in the distribution of FREE elements:
//     Under G:     y_{mm} ~ sqrt(chi^2_{d + n_m}),   y_{ml} ~ N(0,1)
//     Under G_{-e}: y_{mm} ~ sqrt(chi^2_{d + n_m-1}), y_{ml} non-free
//  Therefore delta_E requires 2 separate Monte Carlo runs with their own draws.

// logsumexp: numerically stable log(sum(exp(x)))
static double logsumexp(const arma::vec& x)
{
    arma::vec x_finite = x.elem(arma::find_finite(x));
    if(x_finite.n_elem == 0) return -arma::datum::inf;
    double x_max = x_finite.max();
    // Sum only finite terms but divide by full S (non-finite contribute 0)
    return x_max + std::log(arma::sum(arma::exp(x_finite - x_max)));
    // Caller subtracts log(S) for the full sample size
}

// Upper triangular adjacency as arma::imat
// [[Rcpp::export]]
arma::imat cpp_G_upper(arma::uword p, const arma::umat& G) {
    arma::imat G_upper(p, p, arma::fill::zeros);
    for(arma::uword i = 0; i < p; i++)
        for(arma::uword j = i + 1; j < p; j++)
            G_upper(i, j) = static_cast<int>(G(i, j));
    return G_upper;
}

// Upper-triangular degree for each node in G
// [[Rcpp::export]]
arma::ivec cpp_n_vec(arma::uword p, const arma::imat& G_upper) {
    arma::ivec n_vec(p, arma::fill::zeros);
    for(arma::uword i = 0; i < p; i++)
        for(arma::uword j = i + 1; j < p; j++)
            n_vec(i) += G_upper(i, j);
    return n_vec;
}

// constraint_value: compute y_{r,s} via equations (31)-(32): returns the value y_{r,s} WOULD take if (r,s) were non-free (from Atay-Kayis & Massam, 2005)
static double constraint_value(const arma::mat& y,
                                const arma::mat& t_norm,
                                int r, int s)
{
    double y_rs = 0.0;

    if(r == 0) {
        // Equation (32): y_{0,s} = -sum_{k=0}^{s-1} y_{0,k} * t_norm(k,s)
        for(int k = 0; k < s; k++)
            y_rs -= y(0, k) * t_norm(k, s);

    } else {
        // Equation (31)

        // First term: -sum_{k=r}^{s-1} y_{r,k} * t_norm(k,s)
        for(int k = r; k < s; k++)
            y_rs -= y(r, k) * t_norm(k, s);

        // Second term: -sum_{q=0}^{r-1} [bracket_r * bracket_s]
        for(int q = 0; q < r; q++) {

            // bracket_r = (y_{q,r} + sum_{k=q}^{r-1} y_{q,k}*t_norm(k,r)) / y_{r,r}
            double bracket_r = y(q, r);
            for(int k = q; k < r; k++)
                bracket_r += y(q, k) * t_norm(k, r);
            bracket_r /= y(r, r);

            // bracket_s = y_{q,s} + sum_{k=q}^{s-1} y_{q,k}*t_norm(k,s)
            double bracket_s = y(q, s);
            for(int k = q; k < s; k++)
                bracket_s += y(q, k) * t_norm(k, s);

            y_rs -= bracket_r * bracket_s;
        }
    }

    return y_rs;
}

// mc_log_E: one Monte Carlo run computing log E[exp(-f/2)] (from Atay-Kayis & Massam, 2005) 
//
// Draws free elements under a given graph specification and computes
// the sum of squared non-free elements f^{(s)} for each sample.
//
// Arguments:
//   G_upper   : upper triangular adjacency of the graph being sampled (imat)
//   t_norm    : normalized T matrix, t_norm(j,s) = T(j,s)/T(s,s)
//   d         : degrees of freedom for this graph
//   n_vec     : upper-triangular degree of each node for this graph
//   p         : dimension
//   S         : number of MC samples
//
// Returns: log(mean(exp(-f/2))) = logsumexp(-f/2) - log(S)
//
static double mc_log_E(const arma::imat& G_upper,
                        const arma::mat&  t_norm,
                        double             d,
                        const arma::ivec& n_vec,
                        int                p,
                        int                S)
{
    arma::vec f_s(S, arma::fill::zeros);
    arma::mat y(p, p, arma::fill::zeros);

    for(int s = 0; s < S; s++) {

        y.zeros();

        // Diagonal: y_{ii} ~ sqrt(chi^2_{d + n_i})  [equation 41]
        for(int i = 0; i < p; i++) {
            double df_i = d + static_cast<double>(n_vec(i));
            y(i, i) = std::sqrt(R::rchisq(df_i));
        }

        // Free off-diagonal: y_{ij} ~ N(0,1) for (i,j) in E  [equation 42]
        for(int i = 0; i < p; i++)
            for(int j = i + 1; j < p; j++)
                if(G_upper(i, j) == 1)
                    y(i, j) = R::rnorm(0.0, 1.0);

        // Non-free elements: compute via constraint, accumulate into f_s
        for(int j = 1; j < p; j++)
            for(int i = 0; i < j; i++)
                if(G_upper(i, j) == 0) {
                    double y_ij = constraint_value(y, t_norm, i, j);
                    y(i, j)  = y_ij;
                    f_s(s)  += y_ij * y_ij;
                }
    }

    // log(mean(exp(-f/2))) = logsumexp(-f/2) - log(S)
    arma::vec lse = -0.5 * f_s;
    return logsumexp(lse) - std::log(static_cast<double>(S));
}

// delta_log_C (from methodological paper Arena et al. 2026)
//
// It is the closed-form part of log I_{G_{-e}}(d,D) - log I_G(d,D):
//
//   delta_log_C = - 0.5 * log(4 * pi)
//               + lgamma((d + n_m - 1) / 2)
//               - lgamma((d + n_m) / 2)
//               - log(T_{mm})
//               - log(T_{ll})
//
// where n_m = upper-triangular degree of node m in G,
//       and T_{mm}, T_{ll} = diagonal of chol(D^{-1})
//
double delta_log_C(double d,
                   int    n_m,
                   double T_mm,
                   double T_ll)
{
    double val = -0.5 * std::log(4.0 * M_PI)
               + R::lgammafn((d + n_m - 1.0) / 2.0)
               - R::lgammafn((d + n_m)        / 2.0)
               - std::log(T_mm)
               - std::log(T_ll);
    return val;
}

// delta_E (from methodological paper Arena et al. 2026)
//
// It is a Monte Carlo estimate of:
//   delta_E = log E_{G_{-e}}[exp(-f_{G_{-e}}/2)] - log E_G[exp(-f_G/2)]
//
// It uses 2 separate Monte Carlo runs because the sampling distributions differ:
//
//   Under G:
//     y_{ii}  ~ sqrt(chi^2_{d + n_i})     for all i
//     y_{ml}  ~ N(0,1)                    (free edge)
//     f_G     = sum of y_{ij}^2 for (i,j) in V'_G  (non-edges of G)
//
//   Under G_{-e}:
//     y_{mm}  ~ sqrt(chi^2_{d + n_m - 1}) (degree decreases by 1)
//     y_{ii}  ~ sqrt(chi^2_{d + n_i})     for i != m  (unchanged)
//     y_{ml}  is non-free, computed via constraint
//     f_{G-e} = sum of y_{ij}^2 for (i,j) in V'_{G-e} = V'_G union {(m,l)}
//
// Arguments:
//   G_upper : upper triangular adjacency of G, INTEGER matrix (arma::imat)
//   T       : upper triangular chol(D^{-1})
//   d       : degrees of freedom
//   n_vec   : upper-triangular degree of each node in G
//   m   : 0-based row of tested edge (m < l, G_upper(m,l)=1)
//   l   : 0-based col of tested edge
//   S       : MC samples per run
//
double delta_E(const arma::imat& G_upper,
               const arma::mat&  T,
               double             d,
               const arma::ivec& n_vec,
               arma::uword                m,
               arma::uword                l,
               int                S)
{

    int p = G_upper.n_rows;

    // Precompute t_norm(j,s) = T(j,s)/T(s,s)  [equation 28]
    arma::mat t_norm(p, p, arma::fill::zeros);
    for(int s = 0; s < p; s++)
        for(int j = 0; j <= s; j++)
            t_norm(j, s) = T(j, s) / T(s, s);

    // Run 1: MC under G
    double log_E_G = mc_log_E(G_upper, t_norm, d, n_vec, p, S);

    // Run 2: MC under G_{-e}
    // Build G_{-e} adjacency: remove edge (m, l)
    arma::imat G_upper_e = G_upper;
    G_upper_e(m, l) = 0;

    // Update n_vec: node m loses one upper-triangular neighbor
    arma::ivec n_vec_e = n_vec;
    n_vec_e(m) -= 1;

    double log_E_Ge = mc_log_E(G_upper_e, t_norm, d, n_vec_e, p, S);

    return log_E_Ge - log_E_G;
}


// Bayes Factor (BF_{01}) calculated as the exp of log_bf01
double gwishart_bf01(const arma::mat& nS,
                         const arma::mat& D, // inverse precision matrix (K^{-1})
                         const arma::mat& T_prior, // cholesky of K (precision matrix)
                         const arma::imat &G_upper, // upper triangular adjacency of G (imat)
                         const arma::ivec &n_vec, // upper-triangular degree of each node in G
                         double           d,
                         arma::uword              n,
                         arma::uword              m,
                         arma::uword              l,
                         int              S = 1000)
{
    Rcpp::RNGScope rng_scope;
    

    int n_m = n_vec(m);

    // Prior Cholesky: arma::mat T_prior = arma::chol(K), T = chol(D^{-1}) = chol(K) where K = D^{-1} is the prior precision
    double T_mm_prior = T_prior(m, m);
    double T_ll_prior = T_prior(l, l);

    // Posterior Cholesky: T* = chol((D + nS)^{-1})
    arma::mat T_post  = arma::chol(arma::inv_sympd(D + nS));
    double T_mm_post  = T_post(m, m);
    double T_ll_post  = T_post(l, l);

    // Closed-form differences
    double dC_prior = delta_log_C(d,                          n_m, T_mm_prior, T_ll_prior);
    double dC_post  = delta_log_C(d + static_cast<double>(n), n_m, T_mm_post,  T_ll_post);

    // MC differences: two separate runs each (under G and G_{-e})
    double dE_prior = delta_E(G_upper, T_prior, d,                          n_vec, m, l, S);
    double dE_post  = delta_E(G_upper, T_post,  d + static_cast<double>(n), n_vec, m, l, S);

    // log BF_{01}
    return std::exp((dC_post + dE_post) - (dC_prior + dE_prior));
}


// Bayes Factor (BF_{01}) calculated as the exp of log_bf01 but using precomputed prior terms
// (dC_prior + dE_prior) to avoid recomputing them for every (h,j,i) call.
// Prior terms are fixed for a given edge, graph, and prior hyperparameters,
// they do not depend on the data nS or sample size n.
//
double gwishart_bf01_post_only(const arma::mat& nS,
                                const arma::mat& D,
                                const arma::imat& G_upper,
                                const arma::ivec& n_vec,
                                double            d,
                                arma::uword       n,
                                arma::uword       m,
                                arma::uword       l,
                                int               S,
                                double            dC_prior,
                                double            dE_prior)
{
    Rcpp::RNGScope rng_scope;

    // Posterior Cholesky: T* = chol((D + nS)^{-1})
    arma::mat T_post = arma::chol(arma::inv_sympd(D + nS));

    // Posterior closed-form difference
    double dC_post = delta_log_C(d + static_cast<double>(n),
                                  n_vec(m),
                                  T_post(m, m),
                                  T_post(l, l));

    // Posterior MC difference
    double dE_post = delta_E(G_upper, T_post,
                              d + static_cast<double>(n),
                              n_vec, m, l, S);

    // log BF_{01} using precomputed prior terms
    return std::exp((dC_post + dE_post) - (dC_prior + dE_prior));
}


// BFDA for sparse graphs: sample size planning for a single edge in a sparse graph using Bayes factor with Atay-Kayis computation of the normalizing constant for G-Wishart distribution
// [[Rcpp::export]]
Rcpp::List cpp_design_bfda_edge_sparse(const arma::mat &K,  // K is already elicited as K_prior/nu
                                    const arma::uword &nu, 
                                    arma::umat &G,
                                    const arma::imat &G_upper, // upper triangular adjacency of G (imat)
                                    const arma::ivec &n_vec, // upper-triangular degree of each node in G
                                    const arma::mat &Rho, // partial correlation matrix corresponding to K
                                    const arma::uword &m, const arma::uword &l, // m < l for the AK-BF01
                                    const arma::uword &H, const arma::uword &J, arma::uvec &n,
                                    const double &pow0 = 0.8, const double &pow1 = 0.8,
                                    const double &threshold = 10.0, const bool optimize = false, 
                                    const int nsim_bf = 1000, const arma::uword &n_tol = 1, 
                                    const arma::uword &max_n = 10000, const std::string &gwish_sampler = "direct", 
                                    const double &gwish_tol = 1e-08, const arma::uword &gwish_iter = 500, const arma::uword &gwish_burnin = 500) {

    Rcpp::List out;
    arma::uword p = K.n_cols;
    arma::vec mu(p, arma::fill::zeros);

    arma::mat D = arma::inv_sympd(K);
    arma::mat cholK = arma::chol(K);
    // compute prior terms once — fixed for this edge, graph and prior
    int   n_m_prior   = n_vec(m); // upper-triangular degree of node m in G for the prior 
    double d_shape    = static_cast<double>(nu) - static_cast<double>(p) + 1.0;
    double dC_prior   = delta_log_C(d_shape,
                                    n_m_prior,
                                    cholK(m, m),
                                    cholK(l, l));
    double dE_prior   = delta_E(G_upper, cholK, d_shape,
                                n_vec,
                                m, l,
                                nsim_bf);

    if (optimize) {

        // check: optimize=TRUE only valid for prior present edges
        if (G(m,l) == 0) {
            Rcpp::stop("optimize=TRUE is only valid for prior present edges (G(m,l)=1). Use optimize=FALSE for absent edges.");
        }

        auto eval_edge_at_n = [&](arma::uword n_val, arma::uword hypothesis,
                                double &cur_power, double &cur_error) {
                                
            arma::vec bf_vals(H * J, arma::fill::zeros);
            arma::umat G_temp = G;
            if (hypothesis == 0) {
                G_temp(l,m) = G_temp(m,l) = 0; // H0: edge absent
            }

            // hypothesis 1: G_temp keeps original G (edge present)
            arma::uword idx = 0;
            for (arma::uword h = 0; h < H; h++) {
                arma::mat K_h     = rgwishart(1, K, nu, G_temp, gwish_sampler, gwish_tol, gwish_iter, gwish_burnin, R_NilValue).slice(0);  
                arma::mat S_h     = arma::inv_sympd(K_h);
                arma::mat cholS_h = arma::chol(S_h, "lower");
                for (arma::uword j = 0; j < J; j++) {
                    arma::mat X_j  = mvnrnd_chol(mu, cholS_h, static_cast<int>(n_val));
                    arma::mat nS_j = X_j.t() * X_j;
                    double bf_01   = gwishart_bf01_post_only(nS_j, D, G_upper, n_vec, d_shape, n_val, m, l, nsim_bf, dC_prior, dE_prior); // we provide the full G because it is set up internally to compute the BF01 as long as G(m,l)==1. The G_temp is only used to generate the data under the correct hypothesis.
                    bf_vals(idx++) = bf_01; // G(m,l)==1 guaranteed, always use BF01
                }
            }
            double hj = static_cast<double>(H * J);
            if (hypothesis == 0) {
                cur_power = arma::accu(bf_vals > threshold)       / hj; // Pr(BF01 > threshold | H0)
                cur_error = arma::accu(bf_vals < 1.0 / threshold) / hj; // FNR
            } else {
                cur_power = arma::accu(bf_vals < 1.0 / threshold) / hj; // Pr(BF01 < 1/threshold | H1)
                cur_error = arma::accu(bf_vals > threshold)        / hj; // FPR
            }
        };

        double rho_ml     = Rho(m, l);
        arma::uword n_min = static_cast<arma::uword>(p + 2);

        // --- hypothesis 0: bisect power_h0 ---
        arma::wall_clock timer; // to measure time for each hypothesis
        
        double cur_power = 0.0, 
                cur_error = 1.0, 
                pow0_at_n_star = 0.0, 
                pow1_at_n_star = 0.0, 
                fnr_at_n_star = 1.0, 
                fpr_at_n_star = 1.0;
        arma::uword nl_pow0      = n_min; // because as the prior becomes more informative, the n* decreases and can also be lower than nu
        arma::uword nr_pow0      = warm_start_n(0.01, p, 0.05, pow0); // define upper bracket as the n* for a partial correlation as low as 0.01

        timer.tic();
        if (nr_pow0 > nl_pow0) {
            while (nr_pow0 - nl_pow0 > n_tol) {
                if (nl_pow0 > max_n) break;
                arma::uword n_mid = (nl_pow0 + nr_pow0) / 2;
                eval_edge_at_n(n_mid, 0, cur_power, cur_error);
                if (cur_power >= pow0) {
                    nr_pow0 = n_mid;
                    pow0_at_n_star = cur_power; 
                    fpr_at_n_star = cur_error; 
                }
                else nl_pow0 = n_mid;
            }
        }
        double t_h0 = timer.toc();


        bool converged_h0 = (pow0_at_n_star >= pow0) && (nr_pow0 <= max_n);
        if(converged_h0) {
            out["n_star_power_h0"] = nr_pow0;
            out["fpr_at_n_star_power_h0"] = fpr_at_n_star;
            out["power_h0_at_n_star"] = pow0_at_n_star;
        }
        out["converged_h0"] = converged_h0;
        out["duration_h0"] = t_h0;

        


        // --- hypothesis 1: bisect power_h1 ---
        cur_power = 0.0; cur_error = 1.0;
        arma::uword nl_pow1      = n_min; // because as the prior becomes more informative, the n* decreases and can also be lower than nu
        arma::uword nr_pow1      = warm_start_n(rho_ml, p, 0.05, pow1) * 10; // define upper bracket as the n* for the selected partial correlation and inflate it by 10 

        timer.tic();
        if (nr_pow1 > nl_pow1) {
            while (nr_pow1 - nl_pow1 > n_tol) {
                if (nl_pow1 > max_n) break;
                arma::uword n_mid = (nl_pow1 + nr_pow1) / 2;
                eval_edge_at_n(n_mid, 1, cur_power, cur_error);
                if (cur_power >= pow1) {
                    nr_pow1 = n_mid;
                    pow1_at_n_star = cur_power; 
                    fnr_at_n_star = cur_error; 
                }
                else nl_pow1 = n_mid;
            }
        }
        double t_h1 = timer.toc();

        bool converged_h1 = (pow1_at_n_star >= pow1) && (nr_pow1 <= max_n);
        if(converged_h1) {
            out["n_star_power_h1"] = nr_pow1;
            out["fnr_at_n_star_power_h1"] = fnr_at_n_star;
            out["power_h1_at_n_star"] = pow1_at_n_star;
        }
        out["converged_h1"] = converged_h1;
        out["duration_h1"] = t_h1;

    } else {

        if(G(m,l) == 1){ // for prior present edges we quantify the BF01 directly. We handle only H0, to exclude the edge.
            arma::mat power_h0(n.n_elem, 1, arma::fill::zeros),
            fpr_h0(n.n_elem,  1, arma::fill::zeros),
            power_h1(n.n_elem, 1, arma::fill::zeros),
            fnr_h1(n.n_elem,  1, arma::fill::zeros);
            for (arma::uword hypothesis = 0; hypothesis < 2; hypothesis++) {
                arma::umat G_temp = G;
                if (hypothesis == 0) {
                    G_temp(l,m) = G_temp(m,l) = 0; // H0: edge absent (impose edge absent)
                }
                arma::mat bf_vals(n.n_elem, H * J, arma::fill::zeros);
                arma::uword counter_iter = 0;
                for (arma::uword h = 0; h < H; h++) {
                    arma::mat K_h    = rgwishart(1, K, nu, G_temp, gwish_sampler, gwish_tol, gwish_iter, gwish_burnin, R_NilValue).slice(0); 
                    arma::mat S_h    = arma::inv_sympd(K_h);
                    arma::mat cholSh = arma::chol(S_h, "lower");
                    for (arma::uword j = 0; j < J; j++) {
                        for (arma::uword i = 0; i < n.n_elem; i++) {
                            arma::mat X_j  = mvnrnd_chol(mu, cholSh, static_cast<int>(n(i)));
                            arma::mat nS_j = X_j.t() * X_j;
                            double bf_01   = gwishart_bf01_post_only(nS_j, D, G_upper, n_vec, d_shape, n(i), m, l, nsim_bf, dC_prior, dE_prior);
                            bf_vals(i, counter_iter) = bf_01;
                        }
                        counter_iter++;
                    }
                }
                double hj = static_cast<double>(H * J);
                if (hypothesis == 0) {
                    for (arma::uword i = 0; i < n.n_elem; i++) {
                        power_h0(i, 0) = arma::accu(bf_vals.row(i) > threshold)       / hj; // under H0: edge absent, power is Pr(BF01 > threshold)
                        fpr_h0(i, 0)   = arma::accu(bf_vals.row(i) < 1.0 / threshold) / hj; // under H0: edge absent, FPR (false positive rate) is Pr(BF01 < 1/threshold), because we detect an edge when actually absent
                    }
                    out["bf_h0"] = bf_vals;
                } else {
                    for (arma::uword i = 0; i < n.n_elem; i++) {
                        power_h1(i, 0) = arma::accu(bf_vals.row(i) < 1.0 / threshold) / hj; // under H1: edge present, power is Pr(BF01 < 1/threshold)
                        fnr_h1(i, 0)   = arma::accu(bf_vals.row(i) > threshold)       / hj; // under H1: edge present, FNR (false negative rate) is Pr(BF01 > threshold), because we fail to detect an edge (exclude it) when actually present
                    }
                    out["bf_h1"] = bf_vals;
                }
            }
            out["power_h0"] = power_h0;
            out["fpr_h0"]   = fpr_h0;
            out["power_h1"] = power_h1;
            out["fnr_h1"]   = fnr_h1;
        }

        if(G(m,l) == 0){ // if prior edge is absent we need to handle differently the hypothesis H1 for inclusion, with the final inversion of the BF 
            arma::mat power_h0(n.n_elem, 1, arma::fill::zeros),
            fnr_h0(n.n_elem,  1, arma::fill::zeros),
            power_h1(n.n_elem, 1, arma::fill::zeros),
            fpr_h1(n.n_elem,  1, arma::fill::zeros);
            for (arma::uword hypothesis = 0; hypothesis < 2; hypothesis++) {
                arma::umat G_temp = G;
                if (hypothesis == 1) {
                    G_temp(l,m) = G_temp(m,l) = 1; // H1: edge present
                }
                arma::mat bf_vals(n.n_elem, H * J, arma::fill::zeros);
                arma::uword counter_iter = 0;
                for (arma::uword h = 0; h < H; h++) {
                    arma::mat K_h    = rgwishart(1, K, nu, G_temp, gwish_sampler, gwish_tol, gwish_iter, gwish_burnin, R_NilValue).slice(0);  
                    arma::mat S_h    = arma::inv_sympd(K_h);
                    arma::mat cholSh = arma::chol(S_h, "lower");
                    for (arma::uword j = 0; j < J; j++) {
                        for (arma::uword i = 0; i < n.n_elem; i++) {
                            arma::mat X_j  = mvnrnd_chol(mu, cholSh, static_cast<int>(n(i)));
                            arma::mat nS_j = X_j.t() * X_j;
                            double bf_01   = gwishart_bf01_post_only(nS_j, D, G_upper, n_vec, d_shape, n(i), m, l, nsim_bf, dC_prior, dE_prior);
                            bf_vals(i, counter_iter) = 1.0 / bf_01; // we always compute the BF01 where the H0 and H1 correspond to those for prior present edges, then we invert the BF01 for prior absent edges (using encompassing prior approach)
                        }
                        counter_iter++;
                    }
                }
                double hj = static_cast<double>(H * J);
                if (hypothesis == 0) { // which is H1: rho = 0 for prior absent edges
                    for (arma::uword i = 0; i < n.n_elem; i++) {
                        power_h1(i, 0) = arma::accu(bf_vals.row(i) < 1.0 / threshold) / hj; // under H1: edge absent, power is Pr(1/BF01 < 1/threshold)
                        fpr_h1(i, 0)   = arma::accu(bf_vals.row(i) > threshold) / hj; // under H1: edge absent, FPR (false positive rate) is Pr(1/BF01 > threshold), because we detect an edge when actually absent
                    }
                    out["bf_h1"] = bf_vals;
                } else { // which is the H0 : rho \in (-1,1) for prior absent edges
                    for (arma::uword i = 0; i < n.n_elem; i++) {
                        power_h0(i, 0) = arma::accu(bf_vals.row(i) > threshold) / hj; // under H0: edge present, power is Pr(1/BF01 > threshold)
                        fnr_h0(i, 0)   = arma::accu(bf_vals.row(i) < 1.0 / threshold)       / hj; // under H0: edge present, FNR (false negative rate) is Pr(1/BF01 < 1/threshold), because we fail to detect an edge (exclude it) when actually present
                    }
                    out["bf_h0"] = bf_vals;
                }
            }
            out["power_h0"] = power_h0;
            out["fnr_h0"]   = fnr_h0;
            out["power_h1"] = power_h1;
            out["fpr_h1"]   = fpr_h1;
        }



    }

    out["m"] = m;
    out["l"] = l;
    return out;
}


// =============================================================================
// BSDA: Bayesian Structural Design Analysis via Probit inversion
// =============================================================================

// Evaluation of Pr(measure > measure_value), that is the power at a single n.
// [[Rcpp::export]]
Rcpp::List power_at_n(
        const arma::mat&  K,            // elicited precision scale from a ggm_elicited object via elicit_prior.ggm_parameters()
        const arma::mat&  G,            // true adjacency, symmetric 0/1
        const arma::mat&  pip,          // elicited PIPs (prior inclusion probabilities)
        int               p,
        int               n,           // sample size to evaluate at
        int               H,           // number of outer (K_h) draws 
        int               J,           // datasets per draw (use 1)
        int               nu,
        std::string       gwish_sampler, // "direct" or "block"
        double            gwish_tol,
        arma::uword       gwish_iter,
        arma::uword       gwish_burnin,
        double            edge_selection_threshold,
        std::string       measure,     // "sen" or "spe"
        double            measure_value,
        std::string       init          = "empty",   // "prior" or "empty" or "truth"
        int               fit_iterations = 10000,
        int               fit_burnin     = 5000,
        double            alpha          = 0.05,
        Rcpp::Nullable<int> seed         = R_NilValue) // set the same seed at every n for common random numbers
{
    Rcpp::RNGScope scope;

    // common random numbers: identical seed -> identical K_h sequence across n,
    // which cancels the between-draw variance from differences along the curve.
    if (seed.isNotNull()) {
        Rcpp::Environment base("package:base");
        Rcpp::Function set_seed = base["set.seed"];
        set_seed(Rcpp::as<int>(seed));
    }

    const bool sen = (measure == "sen");
    if (!sen && measure != "spe") Rcpp::stop("measure must be \"sen\" or \"spe\"");

    const arma::uvec ut  = arma::trimatu_ind(arma::size(G), 1);
    const arma::vec  gtr = G.elem(ut);
    const double P1 = arma::accu(gtr > 0.5);
    const double P0 = (double)gtr.n_elem - P1;
    if (sen  && P1 == 0.0) Rcpp::stop("sensitivity undefined: no true edges");
    if (!sen && P0 == 0.0) Rcpp::stop("specificity undefined: no true non-edges");


    // define the measure function to compute sensitivity or specificity from a given posterior edge probability vector
    const double thr = edge_selection_threshold;
    auto score = [&](const arma::vec& gh) {
        if (sen) { double TP = arma::accu((gh > thr) % (gtr > 0.5)); return TP / P1; }
        else     { double TN = arma::accu((gh < thr) % (gtr < 0.5)); return TN / P0; }
    };

    const arma::umat G_u = (G > 0.5);
    const arma::vec  mu(p, arma::fill::zeros);

    // precompute the prior scale matrix for the G-Wishart sampler
    // K is the elicited prior scale matrix = K_prior / nu  (nu = nu = prior df)
    const arma::mat scale_prior = K;                  // = K_prior / nu
    const arma::mat D0          = arma::inv_sympd(K); // = nu * K_prior^{-1}

    // precompute the log-odds of the elicited edge prior probabilities for use in the sampler
    arma::mat plo(p, p, arma::fill::zeros);
    for (int i = 0; i < p; ++i)
        for (int j = i + 1; j < p; ++j) {
            double q = std::min(std::max(pip(i, j), 1e-12), 1.0 - 1e-12);
            plo(i, j) = std::log(q / (1.0 - q));
        }

    // analysis-side start graph -- this is not the truth unless explicitly requested.
    // "prior"  = prior-median model 1[pip>0.5]   
    // "empty"  = no edges                        (overdispersed start)
    // "truth"  = the true G                      (optimistic)
    arma::mat G_start;
    if (init != "empty" && init != "prior" && init != "truth" && init != "random")
    Rcpp::stop("init must be one of \"empty\", \"prior\", \"truth\", \"random\"");
    if (init == "truth")      G_start = G;
    else if (init == "empty") G_start = arma::zeros<arma::mat>(p, p);
    else if (init == "random") { // make it an input in a future version, but for now just a random graph with the elicited pip
        G_start = arma::zeros<arma::mat>(p, p);
        for (arma::uword a = 0; a < (arma::uword)p; ++a)
            for (arma::uword b = a + 1; b < (arma::uword)p; ++b)
                if (R::runif(0.0, 1.0) < pip(a, b)) { G_start(a, b) = G_start(b, a) = 1.0; }
    }
    else {
        arma::mat pm = arma::zeros<arma::mat>(p, p);
        for (arma::uword a = 0; a < (arma::uword)p; ++a)
            for (arma::uword b = a + 1; b < (arma::uword)p; ++b)
                if (pip(a, b) > 0.5) { pm(a, b) = pm(b, a) = 1.0; }
        G_start = pm;
    }
    
    // fixed-H loop, Welford accumulators, NO early stopping 
    double mean_pow = 0.0, M2 = 0.0, mean_measure = 0.0;
    arma::vec pow_batches(H, arma::fill::zeros);
    arma::vec meas_batches(H, arma::fill::zeros);

    for (int h = 1; h <= H; ++h) {
        //Rcpp::Rcout << "h=" << h << "/" << H << "  n=" << n << std::endl;
        arma::mat K_h = rgwishart((arma::uword)1, scale_prior, nu, G_u, gwish_sampler,
                                  gwish_tol, gwish_iter, 0,
                                  Rcpp::Nullable<arma::mat>()).slice(0);
        const arma::mat Sigma_h = arma::inv_sympd(K_h);
        const arma::mat cholS_h = arma::chol(Sigma_h);      // upper

        double pow_j = 0.0, val_j = 0.0;
        for (int j = 0; j < J; ++j) {
            const arma::mat X_j = mvnrnd_chol(mu, cholS_h, n);
            const arma::mat XXt = X_j.t() * X_j;
            const arma::mat Dn = D0 + XXt;
            const arma::mat scale_post  = arma::inv_sympd(Dn);   
            bdmcmc_dcbf_result fit = bdmcmc_dcbf_sampler(n, nu, D0, Dn, scale_prior, scale_post, plo, G_start, fit_iterations, fit_burnin, gwish_sampler, gwish_tol, gwish_iter, gwish_burnin);
            arma::mat Ghat  = fit.pip;
            // continuous time sampler (need to include a switch function to select between the two samplers)
            // Rcpp::List fit = bdmcmc_dct_sampler(D0, Dn, scale_prior, scale_post, n, nu, plo, G_start, fit_iterations, fit_burnin, gwish_tol, gwish_iter, gwish_burnin);
            // arma::mat Ghat  = Rcpp::as<arma::mat>(fit["pip_rb"]);
            const arma::vec gh = Ghat.elem(ut);
            const double val = score(gh);
            pow_j += (val > measure_value) ? 1.0 : 0.0;
            val_j += val;
        }
        pow_j /= (double)J;
        val_j /= (double)J;
        pow_batches(h - 1)  = pow_j;
        meas_batches(h - 1) = val_j;

        const double d1 = pow_j - mean_pow;
        mean_pow     += d1 / (double)h;
        M2           += d1 * (pow_j - mean_pow);
        mean_measure += (val_j - mean_measure) / (double)h;
    }

    // batch-means CI half-width (fixed H, so this is a clean CLT interval)
    double halfwidth = arma::datum::nan;
    if (H >= 2) {
        const double var_b = M2 / (double)(H - 1);
        const double se_b  = std::sqrt(var_b / (double)H);
        const double tq    = R::qt(1.0 - alpha / 2.0, (double)(H - 1), 1, 0);
        halfwidth = tq * se_b;
    }

    return Rcpp::List::create(
        Rcpp::Named("n")            = n,
        Rcpp::Named("power")        = mean_pow,
        Rcpp::Named("halfwidth")    = halfwidth,
        Rcpp::Named("measure")      = mean_measure,
        Rcpp::Named("H")            = H,
        Rcpp::Named("J")            = J,
        Rcpp::Named("pow_batches")  = pow_batches,   // per-batch, for CI / bootstrap
        Rcpp::Named("meas_batches") = meas_batches);
}

// integer log-spaced grid on [lo, hi] with m points (unique, sorted)
static arma::ivec bsda_log_grid(int lo, int hi, int m) {
    if (lo < 2) lo = 2;
    arma::vec lg = arma::linspace(std::log((double)lo), std::log((double)hi), m);
    arma::ivec g = arma::conv_to<arma::ivec>::from(arma::round(arma::exp(lg)));
    return arma::unique(g);                        // unique() sorts ascending
}

// struct to hold grid evaluation results
struct GridEval { 
    arma::ivec n; 
    arma::vec power, hw, meas, se; 
}; 

static GridEval bsda_eval_grid(
        const arma::ivec& grid,
        const arma::mat& K, 
        const arma::mat& G, 
        const arma::mat& pip, 
        int p,
        int nu, 
        int H, 
        int J, 
        std::string gwish_sampler,
        double gwish_tol, 
        arma::uword gwish_iter,
        arma::uword gwish_burnin,
        double thr, 
        const std::string& measure, 
        double measure_value,
        const std::string& init, 
        int fit_iterations, 
        int fit_burnin, 
        double alpha,
        const char* tag, 
        bool verbose) {

    const int m = (int)grid.n_elem;

    // allocate output struct
    GridEval e; 
    e.n = grid;
    e.power.set_size(m); 
    e.hw.set_size(m); 
    e.meas.set_size(m); 
    e.se.set_size(m);

    const double tq = R::qt(1.0 - alpha / 2.0, (double)(H - 1), 1, 0);

    for (int i = 0; i < m; ++i) {
        Rcpp::List r = power_at_n(K, G, pip, p, (int)grid(i), H, J, nu,
                                  gwish_sampler, gwish_tol, gwish_iter, gwish_burnin, thr, measure, measure_value,
                                  init, fit_iterations, fit_burnin, alpha, R_NilValue);
        e.power(i) = Rcpp::as<double>(r["power"]);
        e.hw(i)    = Rcpp::as<double>(r["halfwidth"]);
        e.meas(i)  = Rcpp::as<double>(r["measure"]);
        e.se(i)    = e.hw(i) / tq;
        if (verbose) Rcpp::Rcout << "  [" << tag << "] n=" << grid(i) << "  power="
                                 << std::fixed << std::setprecision(3) << e.power(i)
                                 << "  meas=" << e.meas(i) << "\n";
    }
    return e;
}

// Refined grid from the scout in 5 steps:
// (1) find interior (moving) region padded one point each
// (2) if interior is empty, use the full range [range_lower, max_n] otherwise use the interior range
// (3) find padded range via lo_idx -=1 and hi_idx += 1, then convert to lo and hi sample sizes 
// (4) apply bsda_log_grid to the padded range to get the refined grid
// (5) fallback to bsda_log_grid(range_lower, max_n, n_main) if the refined grid has < 3 points
static arma::ivec bsda_refine_grid(const GridEval& s, 
                                    int range_lower, 
                                    int max_n,                     
                                    int n_main, 
                                    double eps) {
    arma::uvec interior = arma::find((s.power > eps) % (s.power < 1.0 - eps));
    int lo, hi;
    if (interior.n_elem == 0) {
        lo = range_lower; hi = max_n;
    } else {
        arma::uword lo_idx = interior.min();
        arma::uword hi_idx = interior.max();
        if (lo_idx > 0)              lo_idx -= 1;
        if (hi_idx + 1 < s.n.n_elem) hi_idx += 1;
        lo = (int)s.n(lo_idx); hi = (int)s.n(hi_idx);
    }
    arma::ivec g = bsda_log_grid(lo, hi, n_main);
    if ((int)g.n_elem < 3) g = bsda_log_grid(range_lower, max_n, n_main);
    return g;
}

// struct to hold the probit inversion results
struct ProbitResult {
    double a, b, n_star, ci_lo, ci_hi, keep_frac;
    int n_nonsat; 
    bool mismatch; 
    std::string slope, direction;
};

static ProbitResult bsda_probit_invert(const arma::ivec& n, 
                                        const arma::vec& power, 
                                        const arma::vec& se,
                                        const std::string& measure, 
                                        double target_pow, double eps, 
                                        int n_boot){
    ProbitResult R0;
    const int m = (int)n.n_elem;
    const bool expected_inc = (measure == "sen"); // sensitivity is expected to increase with n, specificity is expected to decrease with n
    R0.mismatch = false;

    // probit transform + delta-method sd (scalar qnorm/dnorm; arma for the rest)
    arma::vec y = arma::clamp(power, eps, 1.0 - eps);
    arma::vec z(m), sz(m), x = arma::log(arma::conv_to<arma::vec>::from(n));
    for (int i = 0; i < m; ++i) {
        z(i) = R::qnorm(y(i), 0.0, 1.0, 1, 0);
        double dz = R::dnorm(z(i), 0.0, 1.0, 0);
        if (dz < 1e-12) dz = 1e-12;
        sz(i) = se(i) / dz;
    }
    sz = arma::clamp(sz, 1e-3, arma::datum::inf);        // floor (was std::max)
    arma::vec w = 1.0 / (sz % sz);

    // fit only NON-SATURATED points
    arma::uvec keep = arma::find((power > eps) % (power < 1.0 - eps));
    R0.n_nonsat = (int)keep.n_elem;
    if (R0.n_nonsat < 3) {
        R0.a = R0.b = R0.n_star = R0.ci_lo = R0.ci_hi = R0.keep_frac = NA_REAL;
        R0.slope = "undetermined";
        R0.direction = expected_inc ? "lower_crossing" : "upper_bound";
        return R0;
    }

    arma::vec zf = z.elem(keep), wf = w.elem(keep), xf = x.elem(keep);
    arma::mat X(R0.n_nonsat, 2); X.col(0).ones(); X.col(1) = xf;
    arma::mat XtWX = X.t() * arma::diagmat(wf) * X;
    arma::vec beta = arma::solve(XtWX, X.t() * (wf % zf));
    R0.a = beta(0); R0.b = beta(1);

    arma::vec resid = zf - X * beta;
    const double sigma2 = arma::accu(wf % resid % resid) / (double)(R0.n_nonsat - 2);
    arma::mat Vb = sigma2 * arma::inv_sympd(XtWX);

    const bool obs_inc = (R0.b > 0.0);
    R0.mismatch  = (obs_inc != expected_inc);
    R0.slope     = obs_inc ? "increasing" : "decreasing";
    R0.direction = obs_inc ? "lower_crossing" : "upper_bound";

    const double zt = R::qnorm(target_pow, 0.0, 1.0, 1, 0);
    R0.n_star = (std::abs(R0.b) < 1e-8) ? NA_REAL : std::exp((zt - R0.a) / R0.b);

    // vectorized parametric bootstrap (no std::vector / push_back)
    arma::mat Lc = arma::chol(Vb);                       // Vb = Lc' Lc
    arma::mat Zr(n_boot, 2); Zr.imbue([]() { return R::rnorm(0.0, 1.0); });
    arma::mat ab = Zr * Lc; ab.col(0) += R0.a; ab.col(1) += R0.b;

    arma::vec bcol = ab.col(1); // [NOTE] if near-zero draws --> Inf in nall --> possible Inf in ci_hi --> should add a guard here
    arma::vec nall = arma::exp((zt - ab.col(0)) / bcol);         // all draws

    // observed-sign side
    arma::uvec ok;
    if (obs_inc) 
        ok = arma::find(bcol > 0.0); 
    else         
        ok = arma::find(bcol < 0.0);

    arma::vec kept = nall.elem(ok);
    R0.ci_lo = arma_quantile(kept, 0.05);
    R0.ci_hi = arma_quantile(kept, 0.95);
    R0.keep_frac = (double)ok.n_elem / (double)n_boot;
    return R0;
}


// This is the version with the point accumulating strategy, which is more efficient and robust.
// Each refined grid only evaluates n values not already present, and the probit
// is fitted to the whole accumulated set. Saturated points (power 0 or 1) are
// kept in the store but excluded from the fit, since the probit link is
// undefined there. This stops discarding the expensive scout evaluations and
// avoids re-simulating n values already seen.
// [[Rcpp::export]]
Rcpp::List cpp_bsda_probit(
        const arma::mat&  K,
        const arma::mat&  G,
        const arma::mat&  pip,
        int               p,
        int               nu,
        std::string       measure,          // "sen" or "spe"
        double            measure_value,
        double            target_pow,
        int               range_lower,
        int               max_n,
        int               n_scout        = 6,
        int               n_main         = 10,
        int               H              = 50,
        int               H_scout        = 20,
        int               J              = 1,
        double            scout_frac     = 1.0/3.0,
        std::string       gwish_sampler  = "direct",
        double            gwish_tol      = 1e-08,
        arma::uword       gwish_iter     = 500,
        arma::uword       gwish_burnin   = 500,
        double            edge_selection_threshold = 0.5,
        std::string       init           = "empty",
        int               fit_iterations = 10000,
        int               fit_burnin     = 5000,
        double            alpha          = 0.05,
        int               n_boot         = 5000,
        double            eps            = 1e-3,
        int               tol            = 2,
        double            tol_frac       = 0.01, // relative convergence: |n* - prev| / n* < tol_frac
        int               max_iter       = 10,
        bool              verbose        = true) {

    Rcpp::RNGScope scope;
    arma::wall_clock timer;
    timer.tic();

    if (range_lower < 2)      range_lower = 2;
    if (max_n <= range_lower) Rcpp::stop("max_n must exceed range_lower");
    const double thr = edge_selection_threshold;

    int iter_scout = (int)std::round(fit_iterations * scout_frac);
    int burn_scout = (int)std::round(fit_burnin     * scout_frac);
    if (iter_scout < 1) iter_scout = 1;
    if (burn_scout < 1) burn_scout = 1;

    // Persistent store of every point evaluated, across scout and all main
    // passes. Kept sorted implicitly by re-sorting before each fit.
    std::vector<int>    store_n;
    std::vector<double> store_pow, store_se, store_meas;

    // utility function to check if n has already been evaluated (exact integer match)
    auto already_have = [&](int nq) {
        for (int ns : store_n) if (ns == nq) return true;
        return false;
    };

    // utility function to append a GridEval's points to the store, skipping n already present
    auto append_to_store = [&](const GridEval& ev) {
        for (arma::uword i = 0; i < ev.n.n_elem; ++i) {
            int ni = (int)ev.n(i);
            if (already_have(ni)) continue;
            store_n.push_back(ni);
            store_pow.push_back(ev.power(i));
            store_se.push_back(ev.se(i));
            store_meas.push_back(ev.meas(i));
        }
    };

    // Step 1: explore coarse grid -- seed the store with these (they pin the tails the
    // narrow refined grids can never reach). Lower precision (H_scout), so
    // their standard error is larger and the fit downweights them accordingly.
    arma::ivec scout_grid = bsda_log_grid(range_lower, max_n, n_scout);
    if (verbose) Rcpp::Rcout << "[scout] " << scout_grid.n_elem << " pts, H=" << H_scout
                             << " iters=" << iter_scout << "\n";
    GridEval scout = bsda_eval_grid(scout_grid, K, G, pip, p, nu, H_scout, J, gwish_sampler,
                                    gwish_tol, gwish_iter, gwish_burnin, thr, measure, measure_value,
                                    init, iter_scout, burn_scout, alpha, "scout", verbose);
    append_to_store(scout);

    // initial range comes from the scout
    arma::ivec cur_grid = bsda_refine_grid(scout, range_lower, max_n, n_main, eps);

    // Steps 2-4 loop until |n* change| <= tolerance (or max_iter)
    double  n_star_prev = arma::datum::inf;
    ProbitResult fit;     // holds the probit inversion results
    bool converged = false;
    int  iter_counter = 0;

    for (int it = 0; it < max_iter; ++it) {

        iter_counter++;

        // Step 3: evaluate only the points of cur_grid not already in the store
        std::vector<int> fresh;
        for (arma::uword i = 0; i < cur_grid.n_elem; ++i) {
            int ni = (int)cur_grid(i);
            if (!already_have(ni)) fresh.push_back(ni);
        }

        if (!fresh.empty()) {
            arma::ivec fresh_grid = arma::conv_to<arma::ivec>::from(fresh);
            GridEval ev = bsda_eval_grid(fresh_grid, K, G, pip, p, nu, H, J, gwish_sampler,
                                gwish_tol, gwish_iter, gwish_burnin, thr, measure, measure_value,
                                init, fit_iterations, fit_burnin, alpha, "main", verbose);
            append_to_store(ev);
        }

        if (verbose) {
            Rcpp::Rcout << "[current] grid:";
            for (arma::uword i = 0; i < cur_grid.n_elem; ++i)
                Rcpp::Rcout << " " << cur_grid(i);
            Rcpp::Rcout << "  (" << fresh.size() << " new, " << store_n.size()
                        << " total in store)\n";
        }

        // -- assemble the accumulated set, sorted ascending by n --------------
        arma::uword M = store_n.size();
        arma::ivec all_n(M);
        arma::vec all_pow(M), all_se(M), all_meas(M);
        for (arma::uword k = 0; k < M; ++k) {
            all_n(k)    = (double)store_n[k];
            all_pow(k)  = store_pow[k];
            all_se(k)   = store_se[k];
            all_meas(k) = store_meas[k];
        }
        arma::uvec ord = arma::sort_index(all_n);
        all_n = all_n(ord); 
        all_pow = all_pow(ord);
        all_se = all_se(ord); 
        all_meas = all_meas(ord);

        // Step 4: probit inversion on the non-saturated accumulated points.
        // Saturated (power 0 or 1) stay in the store but are dropped here --
        // the probit link is undefined at the boundaries.
        arma::uvec keep = arma::find(all_pow > 0.0 && all_pow < 1.0);
        fit = bsda_probit_invert(all_n(keep), all_pow(keep), all_se(keep), measure,
                                 target_pow, eps, n_boot);

        if (verbose) Rcpp::Rcout << "a=" << fit.a << " b=" << fit.b << " [" << fit.slope
                                 << " -> " << fit.direction << "] n_nonsat=" << fit.n_nonsat
                                 << " n*=" << fit.n_star << " CI[" << fit.ci_lo << ", "
                                 << fit.ci_hi << "] keep=" << fit.keep_frac << "\n";

        // converged if EITHER the absolute or the relative change is within tolerance
        double delta = std::abs(fit.n_star - n_star_prev);
        bool within_abs = delta <= (double)tol;
        bool within_rel = std::isfinite(fit.n_star) && fit.n_star > 0 &&
                        delta / fit.n_star <= tol_frac;
        if (verbose)
            Rcpp::Rcout << "  [conv] delta=" << delta << " tol=" << tol
                        << " rel=" << (fit.n_star > 0 ? delta / fit.n_star : NA_REAL)
                        << " tol_frac=" << tol_frac << "\n";
        if (std::isfinite(fit.n_star) && (within_abs || within_rel)) {
            converged = true;
            break;
        }
        n_star_prev = fit.n_star;

        // Step 2 (re-refine): narrow the grid around the current CI.
        // If every point of the new grid is already in the store, the next
        // pass evaluates nothing, the fit is unchanged, and the loop converges
        // on the following iteration (|n* change| == 0).
        if (fit.n_nonsat >= 3 && std::isfinite(fit.ci_lo) && std::isfinite(fit.ci_hi)) {
            int lo = std::max(range_lower, (int)std::floor(fit.ci_lo * 0.9));
            int hi = std::min(max_n,       (int)std::ceil (fit.ci_hi * 1.1));
            cur_grid = bsda_log_grid(lo, hi, n_main);
            if ((int)cur_grid.n_elem < 3) break;   // window collapsed --> stop, report last fit
        } else {
            break;   // fit not identified --> stop, report last fit (converged stays FALSE)
        }
    }

    bool identified = (fit.n_nonsat >= 3) && std::isfinite(fit.n_star);
    double duration = timer.toc();

    // Return the full accumulated curve as main_*, sorted by n. This is a
    // change from the previous behaviour (last grid only): main_n now holds
    // every non-scout... in fact every point evaluated, scout included, which
    // is the complete power curve the fit was built on.
    arma::uword M = store_n.size();
    arma::vec out_n(M), out_pow(M), out_se(M), out_meas(M);
    for (arma::uword k = 0; k < M; ++k) {
        out_n(k)    = (double)store_n[k];
        out_pow(k)  = store_pow[k];
        out_se(k)   = store_se[k];
        out_meas(k) = store_meas[k];
    }
    arma::uvec ford = arma::sort_index(out_n);
    out_n = out_n(ford); out_pow = out_pow(ford);
    out_se = out_se(ford); out_meas = out_meas(ford);

    return Rcpp::List::create(
        Rcpp::Named("n_star")             = fit.n_star,
        Rcpp::Named("ci_lo")              = fit.ci_lo,
        Rcpp::Named("ci_hi")              = fit.ci_hi,
        Rcpp::Named("direction")          = fit.direction,
        Rcpp::Named("n_nonsat")           = fit.n_nonsat,
        Rcpp::Named("boot_keep_frac")     = fit.keep_frac,
        Rcpp::Named("direction_mismatch") = fit.mismatch,
        Rcpp::Named("a")                  = fit.a,
        Rcpp::Named("b")                  = fit.b,
        Rcpp::Named("slope")              = fit.slope,
        Rcpp::Named("main_n")             = out_n,
        Rcpp::Named("main_power")         = out_pow,
        Rcpp::Named("main_se")            = out_se,
        Rcpp::Named("main_measure")       = out_meas,
        Rcpp::Named("scout_n")            = scout.n,
        Rcpp::Named("scout_power")        = scout.power,
        Rcpp::Named("scout_measure")      = scout.meas,
        Rcpp::Named("measure")            = measure,
        Rcpp::Named("measure_value")      = measure_value,
        Rcpp::Named("target_pow")         = target_pow,
        Rcpp::Named("iterations")         = iter_counter,
        Rcpp::Named("identified")         = identified,
        Rcpp::Named("converged")          = converged,
        Rcpp::Named("duration")           = duration);
}