//#include <Rcpp.h>
#include <RcppArmadillo.h>
#include <string>     // std::string, std::to_string
  
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends("RcppArmadillo")]]

/* 
  Compute objective 
*/
  
double compute_objective(arma::mat& Abar,
                         arma::mat& W,
                         arma::mat& U,
                         arma::colvec& g,
                         arma::colvec& q) {

  int N = Abar.n_rows; int tau = Abar.n_cols;
  arma::colvec D = zeros<arma::colvec>(tau);

  for (int t=0; t < tau; t++) {
    D[t] = sum((Abar.col(t).t() * Abar.col(t)) % (U.col(t).t() * U.col(t)));    
    D[t] = D[t] + sum(W.col(t) % (U.col(t) % U.col(t)));
    D[t] = D[t] - 2 * g[t] * sum(Abar.col(t) % U.col(t));
    D[t] = D[t] + g[t]*2;    
  }  
  return(sum(D % q));
};

// [[Rcpp::export]]
Rcpp::NumericVector compute_objective_quad(Rcpp::NumericMatrix& Abar_,
                           Rcpp::List& W_,
                           Rcpp::NumericMatrix& U_,
                           Rcpp::NumericVector& g_,
                           Rcpp::NumericVector& q_) {
  
  arma::colvec q = as<arma::colvec>(q_);
  arma::colvec g = as<arma::colvec>(g_);
  arma::mat Abar = as<arma::mat>(Abar_);
  arma::mat U    = as<arma::mat>(U_);
  arma::mat W    = zeros<arma::mat>(Abar.n_rows, Abar.n_cols);
  for (int i=0; i<W_.size(); i++) {
    arma::mat Wt = W_[i];
    for (int t=0; t<Abar.n_cols; t++) {
      W[i,t] = Wt[t,t];
    };
  };
  double res = compute_objective(Abar, W, U, g, q);
  return(Rcpp::wrap(res));
}    

/*
  Algorithm to optimize non-monotone submodular function. 
  Note that the function is hard coded for now. 
*/  

// [[Rcpp::export]]
Rcpp::List optimize_submodular_LS(Rcpp::List& Omega, Rcpp::List& U_list, Rcpp::List& params)  {
  
  // define objects
  arma::colvec q = as<arma::colvec>(params['q']);
  arma::colvec g = as<arma::colvec>(params['g']);
  arma::mat Abar, U;
  arma::mat W    = zeros<arma::mat>(Abar.n_rows, Abar.n_cols);
  int N = Omega.size(); 
  for (int i=0; i<N; i++) {
    arma::mat Wt = W_[i];
    
  // 
}