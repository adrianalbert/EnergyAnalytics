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
	Compute forward probabilities
*/


/* 
	Estimate HMM.
*/

// [[Rcpp::export]]
Rcpp::List estimateHMM(
	Rcpp::List model, 
	Rcpp::NumericVector x_obs, 
	Rcpp::DataFrame y_resp, 
	Rcpp::DataFrame z_tran){

	return(
		Rcpp::List::create(Rcpp::Named("coef")= 0,
		Rcpp::Named("se") = 0,
		Rcpp::Named("df") = 0));
};