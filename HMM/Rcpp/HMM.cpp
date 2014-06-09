//#include <Rcpp.h>
#include <RcppArmadillo.h>
#include <string>     // std::string, std::to_string


// in order to use includes like this it's best to build a package!
//#include "/home/adrian/EnergyAnalytics/HMM/Rcpp/HMM.h"		// class definition

using namespace Rcpp;
using namespace arma;

// Algebra & math is done in Armadillo because of more robust support for this functionality.
// Rcpp is used to interface easily with R.
// Some math functionality is based on "Rcpp sugar"

// [[Rcpp::depends("RcppArmadillo")]]

/* 
	Response distribution density
*/

// C++/Armadillo version
arma::mat p_response(double x_obs, const arma::colvec& z, const arma::mat& beta) {

  arma::mat beta_mu  = beta.cols(0, beta.n_cols-2);
  arma::colvec sigma = beta.col(beta.n_cols-1);
  
  arma::colvec be;
  arma::mat P = zeros<arma::mat>(beta_mu.n_rows, beta_mu.n_rows);
  double mu, p, sd;
  
  for (int i = 0; i<beta_mu.n_rows; i++) {
    sd = sigma(i);
    be = arma::trans(beta_mu.row(i));
    mu = arma::dot(be, z);
    P(i,i) = R::pnorm(x_obs, mu, sd, 1, 0);
  }
  
	return P;
}

// wrapper to export to R
// [[Rcpp::export]]
Rcpp::NumericMatrix p_response_R(double x_obs, 
                                 Rcpp::NumericVector& z_, 
                                 Rcpp::NumericMatrix& beta_) {
                                   
  arma::colvec z = as<arma::colvec>(z_);
  arma::mat beta = as<arma::mat>(beta_);
  arma::mat ret_ = p_response(x_obs, z, beta);
  return(Rcpp::wrap(ret_));
}       

/* 
  Derivatives of response probabilities to response parameters
*/

arma::colvec compute_derivatives_response(arma::colvec &parvec) {
    
  
}


/* 
	Transition distribution density
*/

arma::mat p_transition(const arma::colvec& y, const arma::cube& gamma) {
  
  int K = gamma.n_slices;
  arma::mat G = zeros<arma::mat>(K,K);
  arma::mat gamma_k;
  arma::colvec P_k;
  for (int k = 0; k < K; k++ ) {
    gamma_k = gamma.slice(k);
    P_k = exp(gamma_k * y);
    P_k = P_k / sum(P_k);
    G.row(k) = P_k.t();
  }    
	return G;
}

// wrapper to export to R
// [[Rcpp::export]]
Rcpp::NumericMatrix p_transition_R(Rcpp::NumericVector& y_, Rcpp::List& gamma_) {
                                   
  arma::colvec y = as<arma::colvec>(y_);
  int K          = gamma_.size();
  int p          = (as<arma::mat>(gamma_[0])).n_cols;  
  arma::cube gamma = zeros<arma::cube>(K,p,K);
  for (int k = 0; k<K; k++) {
    arma::mat tmp = gamma_[k];
    gamma.slice(k) = tmp;
  }
  arma::mat P = p_transition(y, gamma);
  return(Rcpp::wrap(P));
}       


/* 
  Transform model parameters to incorporate constraints. (response)
*/

// natural parameters to working parameters
arma::colvec response_nat2work(arma::mat& response) {
  
  // response parameters
  arma::mat beta          = response;
  beta.col(0)             = log(beta.col(0));             // intercept > 0  
  beta.col(beta.n_cols-1) = log(beta.col(beta.n_cols-1)); // sigma^2 > 0
  arma::colvec beta_vec   = arma::reshape(beta, beta.n_rows * beta.n_cols, 1);
  
  return(beta_vec);
}

// wrapper to export to R
// [[Rcpp::export]]
Rcpp::NumericVector response_nat2work_R(Rcpp::NumericMatrix& response_) {
  arma::mat response = as<arma::mat>(response_);
  arma::colvec ret   = response_nat2work(response);
  return(Rcpp::wrap(ret));
}


// working parameters to natural parameters
arma::mat response_work2nat(arma::colvec& parvect, int K) {
  
  // response parameters
  int p = (int)(parvect.n_elem / K);
  arma::mat beta = arma::reshape(parvect, K, p);
  beta.col(0) = exp(beta.col(0));
  beta.col(beta.n_cols-1) = exp(beta.col(beta.n_cols-1));
  
  return(beta);
}

// wrapper to export to R
// [[Rcpp::export]]
Rcpp::NumericMatrix response_work2nat_R(Rcpp::NumericVector& parvec_, int K) {
  arma::colvec parvec   = as<arma::colvec>(parvec_);
  arma::mat ret = response_work2nat(parvec, K);
  return(Rcpp::wrap(ret));
}

/* 
  Transform model parameters to incorporate constraints. (transition)
*/

// natural parameters to working parameters
arma::colvec transition_nat2work(arma::cube& gamma) {
  
  // working parameters have 0 diagonals to remove redundancy due to constraints
  int K = gamma.n_slices;
  int p = (gamma.slice(0)).n_cols;  
  arma::cube wp = zeros<arma::cube>(K-1, p, K);
  for (int k = 0; k < wp.n_slices; k++) {
    arma::mat tmp = gamma.slice(k);
    tmp.shed_rows(k,k);
    wp.slice(k) = tmp;
  }
  arma::colvec ret = arma::reshape(wp, wp.n_slices * wp.n_rows * wp.n_cols, 1, 1);
  return(ret);
}

// wrapper to export to R
// [[Rcpp::export]]
Rcpp::NumericVector transition_nat2work_R(Rcpp::List& gamma_) {
  int K = gamma_.size();
  int p = (as<arma::mat>(gamma_[0])).n_cols;
  arma::cube gamma = zeros<arma::cube>(K, p, K);
  for (int k=0; k<K; k++) {
    arma::mat tmp = gamma_[k];
    gamma.slice(k) = tmp;
  }
  arma::colvec ret = transition_nat2work(gamma);
  return(Rcpp::wrap(ret));
}

// working parameters to natural parameters
// add in rows of 0 for reference parameters (diagonals)
arma::cube transition_work2nat(arma::colvec& parvec, int K, int p) {
  arma::cube gamma(K,p,K);
  int idx = 0;
  for (int k = 0; k < K; k++) {
    arma::colvec sel = parvec(span(idx, idx + (K-1) * p - 1));
    arma::mat tmp = arma::reshape(sel, K-1, p);
    arma::rowvec v = zeros<arma::rowvec>(p);
    tmp.insert_rows(k, v);
    gamma.slice(k) = tmp;
    idx = idx + (K-1) * p;
  }
  return(gamma);
}

// wrapper to export to R
// [[Rcpp::export]]
Rcpp::List transition_work2nat_R(Rcpp::NumericVector& parvec_, int K, int p) {
  arma::colvec parvec = as<arma::colvec>(parvec_);
  arma::cube gamma = transition_work2nat(parvec, K, p);
  Rcpp::List gamma_;
  for (int k=0; k<K; k++) {
    std::stringstream ss; ss << k;
    gamma_[ss.str()] = gamma.slice(k);
  }
  return(Rcpp::wrap(gamma_));
}

// concatenate working parameters into one vector
arma::colvec nat2work(arma::mat& beta, arma::cube& gamma) {
  arma::colvec wp_beta  = response_nat2work(beta);
  arma::colvec wp_gamma = transition_nat2work(gamma);
  arma::colvec wp       = arma::join_rows(wp_beta, wp_gamma);  
  cout<<wp;
  return(wp);  
}

// compute sum of slices of cubes
arma::mat sum_cube(arma::cube& B) {
  arma::mat res = zeros<arma::mat>(B.n_rows, B.n_cols);
  for (int k=0; k < B.n_slices; k++) {
    res = res + B.slice(k);
  }
  return(res);
}

/* 
  Compute log-likelihood function.
	Implements [MacDonald and Zucchinni, 2011], pg. 46.
  Compute gradient and Hessian of forward probabilities.
  Implements [Turner, 2008].
*/

arma::colvec compute_log_likelihood(int K,
                                    arma::colvec& parvec, // all model parameters in column vector format
                                    arma::colvec& x_obs,  // observed data
                                    arma::mat& z_resp,    // response covariates
                                    arma::mat& y_tran){	 // transition covariates
  
  int pr = z_resp.n_cols+1;
  int pt = y_tran.n_cols;
  int T  = x_obs.n_elem;
  int noParams = K*pr + K*(K-1)*pt;
  
  // transform from working parameters to natural parameters
  arma::colvec wp_beta  = parvec(span(0, K*pr-1));
  arma::colvec wp_gamma = parvec(span(K*pr, parvec.n_elem-1));
  arma::mat beta        = response_work2nat(wp_beta, K);
  arma::cube gamma      = transition_work2nat(wp_gamma, K, pt);

  // compute derivatives of response and transition parameters
//  arma::colvec dF = compute_derivatives_response(parvec);
//  arma::colvec dP = compute_derivatives_transition(parvec);

	// compute log-likelihood of the data under the model
	// phi <-- alpha/sum(alpha)
  // use arma objects since Rcpp NumericVector/Matrix do not yet fully implement algebraic operations
  arma::mat P, F;
  arma::mat a = zeros<arma::mat>(noParams, K);
  arma::cube B = zeros<arma::cube>(noParams, noParams, K);
  arma::mat F_NA; F_NA = arma::eye(K, K);
  arma::mat H = zeros<arma::mat>(noParams, noParams);  
  arma::colvec G = zeros<arma::mat>(noParams);
  arma::colvec v;
  arma::colvec phi = arma::ones(K, 1); phi = phi / sum(phi);
  
	double ll = 0, u;
	for (int t = 0; t < T; t++) {

    // treat missing value in response
    if (R_IsNA(x_obs(t))) {
      F = F_NA;
    } else {
      F  = p_response(x_obs(t), (z_resp.row(t)).t(), beta);    
    }
    
		// recursions for log-likelihood    
		P  = p_transition((y_tran.row(t)).t(), gamma);    
		v  = (phi.t() * P * F).t();		// alpha_prime(t)
		u  = sum(v);			  // c(t)
		ll = ll + log(u);		
		phi= v / u;				  // alpha_star(t)

		// recursions for gradient of log-likelihood
    
   // a = compute_a(a, parvecF, )
		// B = 
    
    
	}
  
  G = (1 / u) * sum(a, 1);
  H = (1 / u) * sum_cube(B) + (1 / u)*(1 / u) * sum(a, 1) * (sum(a, 1)).t();
  
  // vectorize gradient and hessian
  arma::colvec H_vec = arma::reshape(H, noParams^2, 1);
  arma::colvec l_vec(1); l_vec(0) = ll;
  
  arma::colvec tmp = arma::join_cols(l_vec, G);
  arma::colvec ret = arma::join_cols(tmp, H_vec);

  return(ret);
};

// [[Rcpp::export]]
Rcpp::List compute_log_likelihood_R(Rcpp::NumericVector& parvec_, 
                                             int K,
                                             Rcpp::NumericVector& x_obs_, 
                                             Rcpp::NumericMatrix& z_resp_, 
                                             Rcpp::NumericMatrix& y_tran_){	
  
  // convert to Arma objects for easier math handling
  arma::colvec parvec= as<arma::colvec>(parvec_);
  arma::colvec x_obs = as<arma::colvec>(x_obs_);
  arma::mat z_resp   = as<arma::mat>(z_resp_);
  arma::mat y_tran   = as<arma::mat>(y_tran_);

  int pr = z_resp.n_cols+1;
  int pt = y_tran.n_cols;
  int noParams = K*pr + K*(K-1)*pt;
  
  arma::colvec ret = compute_log_likelihood(K, parvec, x_obs, z_resp, y_tran);
  
  double ll = ret(0);
  arma::colvec G = ret(span(1, noParams));
  arma::colvec H_vec = ret(span(noParams+1, ret.n_elem-1));
  arma::mat H = arma::reshape(H_vec, noParams, noParams);
  
	return(
		Rcpp::List::create(
			Rcpp::Named("LL")= ll,
			Rcpp::Named("G") = G,
			Rcpp::Named("H") = H)
	);
};

//  // access model parameters
//	int T = x_obs.n_elem;
//	int K = model["nStates"];
//  arma::colvec phi             = as<arma::colvec>(model["initDist"]);
//	Rcpp::List response  		     = model["response"];
//  Rcpp::List transition		     = model["transition"];
//  arma::mat beta               = as<arma::mat>(response["coef"]);
//  arma::mat beta_sd            = as<arma::mat>(response["sder"]);
//
