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
	Build an object to store a HMM model.
*/

//RcppExport SEXP initialize_model(int modelSize, arma::colvec vars_resp, arma::colvec vars_tran) {
//  
//	return(
//		Rcpp::List::create(Rcpp::Named("coef")= 0,
//		Rcpp::Named("se") = 0,
//		Rcpp::Named("df") = 0));
//};

/* 
	Response distribution density
*/

// C++/Armadillo version
arma::mat p_response(double x_obs, arma::colvec z, arma::mat beta) {

  arma::mat beta_mu = beta.cols(0, beta.n_cols-2);
  arma::colvec sigma  = beta.col(beta.n_cols-1);
  
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
Rcpp::NumericMatrix r_p_response(double x_obs, 
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

arma::mat p_transition(arma::colvec y, SEXP theta_) {
  
  Rcpp::List theta(theta_);
	arma::mat G;
  
	return G;
}

/* 
	Compute log-likelihood function
	Implements [MacDonald and Zucchinni, 2011], pg. 46
*/

Rcpp::List compute_log_likelihood(SEXP model_, 
                                  arma::colvec x_obs, 
                                  arma::mat y_resp, 
                                  arma::mat z_tran){	

  Rcpp::List model(model_);

	// access model parameters
	int T = x_obs.n_elem;
	int K = model["nStates"];
	Rcpp::NumericVector initDist = model["initDist"];
	Rcpp::List response  		     = model["response"];
	Rcpp::DataFrame respCoef 	   = response["coef"];
	Rcpp::DataFrame respSder 	   = response["sder"];
	Rcpp::List transition		     = model["transition"];

	// compute log-likelihood of the data under the model
	// phi <-- alpha/sum(alpha)
  // use arma objects since Rcpp NumericVector/Matrix do not yet fully implement algebraic operations
  arma::mat P, F, B, H, beta;  
  arma::colvec phi = as<arma::colvec>(initDist), a, G, v;
  Rcpp::List theta;

	double ll = 0, u;
	for (int t = 0; t < T; t++) {

		// log-likelihood
		P  = p_transition(z_tran.row(t), theta);
		F  = p_response(x_obs(t), y_resp.row(t), beta);
		v  = phi * P * F;		// alpha_prime(t)
		u  = sum(v);			// c(t)
		ll = ll + log(u);		
		phi= v / u;				// alpha_star(t)

		// recursions for gradient of log-likelihood
		// a = 
		// B = 
	}

	return(
		Rcpp::List::create(
			Rcpp::Named("LL")= ll,
			Rcpp::Named("G") = G,
			Rcpp::Named("H") = H)
	);
};


/* 
	Estimate HMM by direct likelihood maximization.
*/

//// [[Rcpp::export]]
//RcppExport SEXP estimate_model(SEXP model_, SEXP x_obs_, SEXP y_resp_, SEXP z_tran_){
//
//
//
//	return(
//		Rcpp::List::create(
//			Rcpp::Named("coef")= 0,
//			Rcpp::Named("se") = 0,
//			Rcpp::Named("df") = 0)            
//	);
//};

