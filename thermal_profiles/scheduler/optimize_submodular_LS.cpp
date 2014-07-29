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
    arma::mat tmp = (Abar.col(t) * Abar.col(t).t()) % (U.col(t) * U.col(t).t());
    D[t] = arma::accu(tmp);    
    D[t] = D[t] + sum(W.col(t) % (U.col(t) % U.col(t)));
    D[t] = D[t] - 2 * g[t] * sum(Abar.col(t) % U.col(t));
    D[t] = D[t] + g[t]*g[t];    
  }  
  return(sum(g % g % q) - sum(D % q));
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


  cout<<"Initializing objects...\n";

  // define objects and access data
  
  int N = Omega.size(); 
  arma::colvec q = as<arma::colvec>(params["q"]);
  arma::colvec g = as<arma::colvec>(params["g"]);
  double eps     = 0.01;
  int tau = g.size();
    
  // format arrays
  arma::mat W    = zeros<arma::mat>(N, tau), Wi;
  arma::mat Abar = zeros<arma::mat>(N, tau);
  Rcpp::List Oi; arma::vec Ai;
  for (int i=0; i<N; i++) {
    Oi = Omega[i]; Wi = as<arma::mat>(Oi["w"]); Ai = as<arma::vec>(Oi["a"]);
    arma::vec Wid = Wi.diag();
    W.row(i) = Wid.t();
    Abar.row(i) = Ai.t();
  }

  arma::mat UL  = U_list["UL"];  // list of all schedules
  arma::uvec UA = U_list["UA"];  // indicates which user each schedule belongs to
  UA -= 1; //UA - ones<arma::colvec>(UA.size());
  int Nu = UL.n_rows;

  // initialize sets as 0/1 arrays
  cout<<"Initializing sets..."<<'\n';  
  arma::uvec   A = zeros<arma::uvec>(N);
  arma::uvec   U = ones<arma::uvec>(N) * Nu;
  arma::colvec d = zeros<arma::colvec>(Nu);
  for (int v = 0; v<Nu; v++) {
    int e = UA[v];
    arma::uvec e_vec, v_vec; e_vec<<e; v_vec<<v;
    arma::mat Acur = Abar.rows(e_vec), Wcur = W.rows(e_vec), Ucur = UL.rows(v_vec);
    d[v] = compute_objective(Acur, Wcur, Ucur, g, q);
  }
  arma::uvec tmp = find(d == max(d)); int m = tmp[0];
  A[UA[m]] = 1;
  U[UA[m]] = m;
  
  // main loop
  arma::uvec A1, U1;
  int iter = 0, n = 1;
  bool terminate = 0;
  while(!terminate && n<N) {
    
    iter = iter + 1;
    terminate = 1;    
    tmp = find(A > 0); n = tmp.size();
    arma::uvec e_vec = find(A > 0); 
    arma::uvec v_vec = U.elem(find(U < Nu));
    arma::mat Acur = Abar.rows(e_vec), 
              Wcur = W.rows(e_vec), 
              Ucur = UL.rows(v_vec);
    double obj  = compute_objective(Acur, Wcur, Ucur, g, q);
    double f = 1.0 + eps/(double)pow(n,2);
    Rcout<<"Obj = "<<obj<<"; f="<<f<<"; eps="<<eps<<'\n';
    
    Rcout<<"Step forward: Iteration "<<iter<<"; |A| = "<<n<<'\n';
    // find appropriate element
    bool backstep = 1;
    arma::uvec vvec = zeros<arma::uvec>(N); arma::vec ovec = -ones<arma::vec>(N);
    for (int v = 0; v < Nu; v++) 
      if (A[UA[v]] == 0) {      // if element v is not in A yet
        int e = UA[v]; 
        A1 = A; A1[e] = 1;
        U1 = U; U1[e] = v;
        arma::uvec e_vec1= A1.elem(find(A1 > 0)), v_vec1 = U1.elem(find(U1 < Nu));
        arma::mat Acur1= Abar.rows(e_vec1), Wcur1 = W.rows(e_vec1), Ucur1 = UL.rows(v_vec1);        
        double obj1 = compute_objective(Acur1, Wcur1, Ucur1, g, q);
        if (obj1 > f*obj) {
          vvec[e]=v; ovec[e] =obj1 - obj;
        }        
      }
    arma::uvec evec = find(ovec>0);
    if (evec.size()>0){ // take the maximum
      Rcout<<"|evec|="<<evec.size()<<'\n';
      arma::uvec m = find(ovec == max(ovec));
      int e = as_scalar(m[0]);
      int v = vvec[e];
      A[e] = 1; U[e] = v;
      backstep = 0;
      terminate = 0;
    }
    
    // take step backward
    if (!backstep) continue;
    terminate = 1;
    Rcout<<"Step backward: Iteration "<<iter<<'\n';
        
    ovec = -ones<arma::vec>(N);
    for (int e = 0; e < N; e++) 
      if (A[e] == 1) {      // if element v is in A
        A1 = A; A1[e] = 0;
        U1 = U; U1[e] = Nu;
        arma::uvec e_vec1= A1.elem(find(A1 > 0)), v_vec1 = U1.elem(find(U1 < Nu));
        arma::mat Acur1= Abar.rows(e_vec1), Wcur1 = W.rows(e_vec1), Ucur1 = UL.rows(v_vec1);
        double obj1 = compute_objective(Acur1, Wcur1, Ucur1, g, q);
        if (obj1 > f*obj) {
          ovec[e] = obj - obj1;
        }        
      };                   
    evec = find(ovec>0);
    if (evec.size()>0){ // take the maximum
      Rcout<<"|evec|="<<evec.size()<<'\n';
      arma::uvec m = find(ovec == max(ovec));
      int e = as_scalar(m[0]);
      A[e] = 0; U[e] = Nu;
      terminate = 0;
    }        
  };
  
  Rcout<<"Done"<<'\n';
  
  // return value
  return(Rcpp::List::create(Rcpp::Named("A")=A , Rcpp::Named("U")=U));

}