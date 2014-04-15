#include <Rcpp.h>
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace RcppArmadillo;

/* 
	Inputs
*/	
typedef struct{
	Rcpp::NumericMatrix z_resp;
	Rcpp::NumericMatrix z_tran;
	Rcpp::NumericVector x_obs;	
} HMMInputs;


/* 
	HMM class definition
*/	
class HMM{
private:
	int nStates;
	
	// internal methods
	void compute_backward_probs();	
	void compute_forward_probs();

public:
	// constructor
	void HMM();

	// setter
	int set_model_size();

	// getter
	int get_model_size();

	// method to estimate model
	void estimate();
	
	// method to decode observations given model
	void decode();
	
	// method to access model parameters
	Rcpp::List get_transition_parameters();
	Rcpp::List get_response_parameters();
};