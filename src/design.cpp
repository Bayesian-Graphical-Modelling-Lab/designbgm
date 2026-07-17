#include <string>
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
                        const arma::uword &n_tol = 1, const arma::uword &max_n = 10000) {

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
                arma::mat K_h = random_precision_from_prior(prior, K, nu, G, Kchol);
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
                arma::mat K_h = random_precision_from_prior(prior, K, nu, G, Kchol);
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
                arma::mat K_h = random_precision_from_prior(prior, K, nu, G, Kchol);
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

        // --- phase 2: bisection within bracket ---
        double prob_global_at_n_star = 0.0;  
        arma::vec prob_pw_at_n_star(n_pars, arma::fill::zeros); 

        // (1) Find global bisection 
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

        // (2) Find parameterwise bisection (based on min probability criterion, this ensures all parameters meet the target probability and we find n* as max n* across off-diag parameters)
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

        out["converged_global"] = converged_global;
        out["converged_pw"]     = converged_pw;
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
                K_h = random_precision_from_prior(prior, K, nu, G, Kchol); // random draw from the prior distribution
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
                                    const arma::uword &n_tol = 1, const arma::uword &max_n = 10000) {

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
                arma::mat K_h     = random_precision_from_prior(prior, K, nu, G_temp, Kchol);
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
        
        // --- hypothesis 0: bisect power_h0 ---
        double cur_power = 0.0, 
                cur_error = 1.0, 
                pow0_at_n_star = 0.0, 
                pow1_at_n_star = 0.0, 
                fnr_at_n_star = 1.0, 
                fpr_at_n_star = 1.0;
        arma::uword nl_pow0      = n_min; // because as the prior becomes more informative, the n* decreases and can also be lower than nu
        arma::uword nr_pow0      = warm_start_n(0.01, p, 0.05, pow0); // define upper bracket as the n* for a partial correlation as low as 0.01

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

        bool converged_h0 = (pow0_at_n_star >= pow0) && (nr_pow0 <= max_n);
        if(converged_h0) {
            out["n_star_power_h0"] = nr_pow0;
            out["fpr_at_n_star_power_h0"] = fpr_at_n_star;
            out["power_h0_at_n_star"] = pow0_at_n_star;
        }
        out["converged_h0"] = converged_h0;
        

        // --- hypothesis 1: bisect power_h1 ---
        cur_power = 0.0; cur_error = 1.0;
        arma::uword nl_pow1      = n_min; // because as the prior becomes more informative, the n* decreases and can also be lower than nu
        arma::uword nr_pow1      = warm_start_n(rho_ml, p, 0.05, pow1) * 10; // define upper bracket as the n* for the selected partial correlation and inflate it by 10 

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

        // check convergence and store results 
        bool converged_h1 = (pow1_at_n_star >= pow1) && (nr_pow1 <= max_n);
        if(converged_h1) {
            out["n_star_power_h1"] = nr_pow1;
            out["fnr_at_n_star_power_h1"] = fnr_at_n_star;
            out["power_h1_at_n_star"] = pow1_at_n_star;
        }
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
                arma::mat K_h    = random_precision_from_prior(prior, K, nu, G_temp, Kchol); //rgwishart(1, K, nu, G_temp, 1e-08, 500).slice(0);
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
                                    const arma::uword &max_n = 10000) {

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
                arma::mat K_h     = rgwishart(1, K, nu, G_temp, "direct", 1e-08, 500, 500, R_NilValue).slice(0);  
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
                    arma::mat K_h    = rgwishart(1, K, nu, G_temp, "direct", 1e-08, 500, 500, R_NilValue).slice(0); 
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
                    arma::mat K_h    = rgwishart(1, K, nu, G_temp, "direct", 1e-08, 500, 500, R_NilValue).slice(0);  
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
