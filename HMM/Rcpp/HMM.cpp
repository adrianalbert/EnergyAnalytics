//#include <Rcpp.h>
#include <RcppArmadillo.h>

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
arma::mat p_response(double x_obs, arma::colvec z, arma::mat beta) {

  arma::mat beta_mu  = beta.cols(0, beta.n_cols-2);
  arma::colvec sigma = beta.col(beta.n_cols-1);
  
  arma::colvec be;
  arma::mat P(beta_mu.n_rows, beta_mu.n_rows);
  double mu, p, sd;
  
  for (int i = 0; i<beta_mu.n_rows; i++) {
    sd = sigma(i);
    be = arma::trans(beta_mu.row(i));
    mu = arma::dot(be, z);
    P(i,i) = R::dnorm(x_obs, mu, sd, 0);
  }
  
	return P;
}

// wrapper to export to R
// [[Rcpp::export]]
Rcpp::NumericMatrix p_response_R(double x_obs, 
                                 Rcpp::NumericVector z_, 
                                 Rcpp::NumericMatrix beta_) {
                                   
  arma::colvec z = as<arma::colvec>(z_);
  arma::mat beta = as<arma::mat>(beta_);
  arma::mat ret_ = p_response(x_obs, z, beta);
  return(Rcpp::wrap(ret_));
}       

/* 
	Transition distribution density
*/

arma::mat p_transition(arma::colvec y, arma::cube theta) {
  
  int K = theta.n_slices;
  arma::mat theta_k, G = zeros<arma::mat>(K,K);
  arma::colvec P_k;
  for (int k = 0; k < K; k++ ) {
    theta_k = theta.slice(k);
    P_k = exp(theta_k * y);
    P_k = P_k / sum(P_k);
    G.row(k) = P_k;
  }    
	return G;
}

// wrapper to export to R
// [[Rcpp::export]]
Rcpp::NumericMatrix p_transition_R(Rcpp::NumericVector y_, Rcpp::List theta_) {
                                   
  arma::colvec y = as<arma::colvec>(y_);
  int K          = theta_.size();
  int p          = (as<arma::mat>(theta_[0])).n_cols;  
  cout<<p;
  arma::cube theta = zeros<arma::cube>(K,p,K);
  cout<<theta;
  for (int i = 0; i<p; i++) {
    arma::mat tmp = theta_[i];
    cout<<tmp;
    theta.slice(i) = tmp;
  }
  arma::mat P = p_transition(y, theta);
  return(Rcpp::wrap(P));
}       


/* 
  Transform model parameters to incorporate constraints.
*/
//arma::colvec transform_parameters(Rcpp::List model) {
//  
//  // access model parameters
//  int K                 = model["nStates"];
//  arma::colvec phi      = as<arma::colvec>(model["initDist"]);
//	Rcpp::List response  	= model["response"];
//  Rcpp::List transition	= model["transition"];
//  arma::mat beta        = as<arma::mat>(response["coef"]);
//
//  // initial distribution
//
//  // response parameters
//  beta.col(0)             = log(beta.col(0));             // intercept > 0
//  beta.col(beta.n_cols-1) = log(beta.col(beta.n_cols-1)); // sigma^2 > 0
//  
//  // transition parameters
//  
//}

/* 
	Compute log-likelihood function
	Implements [MacDonald and Zucchinni, 2011], pg. 46
*/

//Rcpp::List compute_log_likelihood(int K,
//                                  arma::colvec parvec, 
//                                  arma::colvec x_obs, 
//                                  arma::mat y_resp, 
//                                  arma::mat z_tran){	
//  
//  // transform from working parameters to natural parameters
//  
//  
//	// compute log-likelihood of the data under the model
//	// phi <-- alpha/sum(alpha)
//  // use arma objects since Rcpp NumericVector/Matrix do not yet fully implement algebraic operations
//  arma::mat P, F, B, H;  
//  arma::colvec a, G, v;
//
//	double ll = 0, u;
//	for (int t = 0; t < T; t++) {
//
//    // treat missing value in response
//    
//		// log-likelihood
//    F  = p_response(x_obs(t), (y_resp.row(t)).t(), beta);
//		P  = p_transition(z_tran.row(t), theta);
//		v  = phi * P * F;		// alpha_prime(t)
//		u  = sum(v);			  // c(t)
//		ll = ll + log(u);		
//		phi= v / u;				  // alpha_star(t)
//
//		// recursions for gradient of log-likelihood
//		// a = 
//		// B = 
//	}
//
//	return(
//		Rcpp::List::create(
//			Rcpp::Named("LL")= ll,
//			Rcpp::Named("G") = G,
//			Rcpp::Named("H") = H)
//	);
//};
//
/////* 
//  Compute log-likelihood function
//	Implements [MacDonald and Zucchinni, 2011], pg. 46
//*/
//
//// [[Rcpp::export]]
//Rcpp::List compute_log_likelihood_R(Rcpp::List model, 
//                                    Rcpp::NumericVector x_obs_, 
//                                    Rcpp::NumericMatrix y_resp_, 
//                                    Rcpp::NumericMatrix z_tran_){	
//  
//  // convert to Arma objects for easier math handling
//  arma::colvec x_obs = as<arma::colvec>(x_obs_);
//  arma::mat y_resp   = as<arma::mat>(y_resp_);
//  arma::mat z_tran   = as<arma::mat>(z_tran_);
//
//	// access model parameters
//	int T = x_obs.n_elem;
//	int K = model["nStates"];
//  arma::colvec phi             = as<arma::colvec>(model["initDist"]);
//	Rcpp::List response  		     = model["response"];
//  Rcpp::List transition		     = model["transition"];
//  arma::mat beta               = as<arma::mat>(response["coef"]);
//  arma::mat beta_sd            = as<arma::mat>(response["sder"]);
//
//	// compute log-likelihood of the data under the model
//	// phi <-- alpha/sum(alpha)
//  // use arma objects since Rcpp NumericVector/Matrix do not yet fully implement algebraic operations
//  arma::mat P, F, B, H;  
//  arma::colvec a, G, v;
//  Rcpp::List theta;
//
//	double ll = 0, u;
//	for (int t = 0; t < T; t++) {
//
//    // treat missing value in response
//    
//		// log-likelihood
//    F  = p_response(x_obs(t), (y_resp.row(t)).t(), beta);
//		P  = p_transition(z_tran.row(t), theta);
//		v  = phi * P * F;		// alpha_prime(t)
//		u  = sum(v);			  // c(t)
//		ll = ll + log(u);		
//		phi= v / u;				  // alpha_star(t)
//
//		// recursions for gradient of log-likelihood
//		// a = 
//		// B = 
//	}
//
//	return(
//		Rcpp::List::create(
//			Rcpp::Named("LL")= ll,
//			Rcpp::Named("G") = G,
//			Rcpp::Named("H") = H)
//	);
//};
//
