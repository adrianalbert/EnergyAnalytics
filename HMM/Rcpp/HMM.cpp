#include <Rcpp.h>
#include <armadillo>
//#include "HMM.h"		// class definition

using namespace Rcpp;
using namespace arma;

/* 
	Build an object to store a HMM model.
*/

// [[Rcpp::export]]
Rcpp::List initialize_model(int modelSize, 
							Rcpp::NumericVector vars_resp,
							Rcpp::NumericVector vars_tran) {

	return(
		Rcpp::List::create(Rcpp::Named("coef")= 0,
		Rcpp::Named("se") = 0,
		Rcpp::Named("df") = 0));
};

/* 
	Response distribution density
*/

Rcpp::NumericMatrix p_response(
	double x_obs, 
	const Rcpp::NumericVector &beta, 
	const Rcpp::NumericVector &z,
	double sigma) {

	double mu = beta * z;
	double p = dnorm(x_obs, mu, sigma);

	return arma::diagvec(p);
}

/* 
	Transition distribution density
*/

Rcpp::NumericMatrix p_transition(
	const Rcpp::List par_trans, 
	const Rcpp::NumericVector &y) {

	Rcpp::NumericMatrix G;
	int K = par_trans.size();


	for (int i = 0; i < K; i++){
		G[i] = 
	}

	return G;
}

/* 
	Compute log-likelihood function
*/

// [[Rcpp::export]]
double compute_log_likelihood(
	Rcpp::List &model, 
	Rcpp::NumericVector &x_obs, 
	Rcpp::DataFrame &y_resp, 
	Rcpp::DataFrame &z_tran){	


	// access model parameters
	int T = x_obs.size()
	int K = model["nStates"];
	Rcpp::List initDist 		= model["initDist"];
	Rcpp::List response  		= model["response"];
	Rcpp::DataFrame respCoef 	= response["coef"]
	Rcpp::DataFrame respSder 	= response["sder"]
	Rcpp::List transition		= model["transition"];

	// compute log-likelihood of the data under the model
	Rcpp::NumericVector phi = initDist;
	double ll = 0;
	for (t = 0; t < T; t++) {
		G  = p_transition(x_obs[t], z_tran[t])
		P  = p_response(x_obs[t], y_resp[t])
		v  = phi * G * P;
		u  = sum(v);
		ll = ll + log(u);
		phi= v / u;
	}

	return(
		Rcpp::List::create(
			Rcpp::Named("coef")= 0,
			Rcpp::Named("se") = 0,
			Rcpp::Named("df") = 0)
	);
};


/* 
	Compute log-likelihood gradient
*/

/* 
	Compute log-likelihood Hessian
*/

/* 
	Estimate HMM by direct likelihood maximization.
*/

// [[Rcpp::export]]
Rcpp::List estimate_model(
	Rcpp::List model, 
	Rcpp::NumericVector x_obs, 
	Rcpp::DataFrame y_resp, 
	Rcpp::DataFrame z_tran){



	return(
		Rcpp::List::create(
			Rcpp::Named("coef")= 0,
			Rcpp::Named("se") = 0,
			Rcpp::Named("df") = 0)
	);
};