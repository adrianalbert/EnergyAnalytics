//#include <Rcpp.h>
#include <RcppArmadillo.h>
#include <string>     // std::string, std::to_string
#include <time.h>
  
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

/*
  Given objective at previous iteration, compute objective at current iteration.
*/

double compute_objective_add(double obj0, 
                              arma::mat& A,
                              arma::mat& W,
                              arma::mat& U,                       
                              arma::rowvec& a, 
                              arma::rowvec& w,
                              arma::rowvec& u,
                              arma::colvec& g,
                              arma::colvec& q) {
                                
  double obj = obj0;
  int N = A.n_rows, tau = A.n_cols;
  arma::colvec D = zeros<arma::colvec>(tau);
  
  for (int t=0; t<tau; t++) {
    D[t] = 2 * sum(U.col(t) % A.col(t)) * u[t] * a[t] + u[t]*u[t]*(w[t] + a[t]*a[t]);
    D[t] = D[t] - 2 * g[t] * u[t] * a[t];
  }
  obj = obj - sum(D % q);
  return(obj);                                
}

/*
  Given objective at previous iteration, compute objective at current iteration.
*/

double compute_objective_del(double obj0, 
                            arma::mat& A,
                            arma::mat& W,
                            arma::mat& U,                       
                            int e, 
                            arma::colvec& g,
                            arma::colvec& q) {
                                
  double obj = obj0;
  arma::rowvec a = A.row(e), w = W.row(e), u = U.row(e);
  arma::mat A1 = A, W1 = W, U1 = U;
  A1.shed_row(e); W1.shed_row(e); U1.shed_row(e);
  int N = A1.n_rows, tau = A1.n_cols;
  arma::colvec D = zeros<arma::colvec>(tau);
  
  for (int t=0; t<tau; t++) {
    D[t] = 2 * sum(U1.col(t) % A1.col(t)) * u[t] * a[t] + u[t]*u[t]*(w[t] + a[t]*a[t]);
    D[t] = D[t] - 2 * g[t] * u[t] * a[t];
  }
  obj = obj + sum(D % q);
  return(obj);                                
}

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

  // variables for timings
  clock_t t1,t2;

  // define objects and access data
  
  int N = Omega.size(); 
  arma::colvec q = as<arma::colvec>(params["q"]);
  arma::colvec g = as<arma::colvec>(params["g"]);
  bool verbose   = as<bool>(params["verbose"]);
  double eps     = as<double>(params["eps"]);
  int NOBJ       = as<int>(params["NOBJ"]);
 // double eps     = 0.005, EPS = 1e-5;
  //int NOBJ       = 10;
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
  if (verbose) Rcout<<"Initializing sets..."<<'\n';  
  arma::uvec   A = zeros<arma::uvec>(N);
  arma::uvec   S = zeros<arma::uvec>(N);
  arma::uvec   U = ones<arma::uvec>(N) * Nu;
  arma::colvec O = -ones<arma::colvec>(N);
  arma::colvec d = zeros<arma::colvec>(Nu);
  for (int v = 0; v<Nu; v++) {
    int e = UA[v];
    arma::uvec e_vec, v_vec; e_vec<<e; v_vec<<v;
    arma::mat Acur = Abar.rows(e_vec), Wcur = W.rows(e_vec), Ucur = UL.rows(v_vec);
    d[v] = compute_objective(Acur, Wcur, Ucur, g, q);
  }
  uword m, e; double maxval = d.max(m); e = UA[m];
  A[e] = 1;
  U[e] = m;
  S[0] = e;
  if (verbose) Rcout<<"Initial element: "<<UA[m]<<"; Initial objective = "<<maxval<<endl;
  
  // main loop
  arma::uvec A1, U1;
  int iter = 0, n = 1;
  bool terminate = 0;
  while(!terminate && n<N) {
    
    t1   = clock();    
    iter = iter + 1;
    terminate = 1;    
    arma::uvec e_vec = find(A > 0); n = e_vec.size();
    arma::uvec v_vec = U.elem(find(U < Nu));
    arma::mat Acur = Abar.rows(e_vec), 
              Wcur = W.rows(e_vec), 
              Ucur = UL.rows(v_vec);
    double obj  = compute_objective(Acur, Wcur, Ucur, g, q);
    double f = 1.0 + eps/(double)pow(n,2);
    O[n-1] = obj;
    
    if (verbose) Rcout<<"Step forward: Iteration "<<iter<<"; |A| = "<<n<<"; Objective = "<<obj<<'\n';
    // find appropriate element
    bool backstep = 1;
    arma::uvec vvec = zeros<arma::uvec>(N); arma::vec ovec = -ones<arma::vec>(N);
    for (int v = 0; v < Nu; v++) 
      if (A[UA[v]] == 0) {      // if element v is not in A yet
        int e = UA[v]; 
        A1 = A; A1[e] = 1;
        U1 = U; U1[e] = v;
        arma::rowvec a = Abar.row(e);
        arma::rowvec u = UL.row(v);
        arma::rowvec w = W.row(e);
        double obj1 = compute_objective_add(obj, Acur, Wcur, Ucur, a, w, u, g, q);        
        if (obj1 > f*obj) {
          vvec[e]=v; ovec[e] =obj1; 
          arma::uvec tmp = find(ovec > 0);
          if (tmp.size() > NOBJ) break;
        }        
      }
    arma::uvec evec = find(ovec>0);
    if (evec.size()>0){ // take the maximum
      uword e; maxval = ovec.max(e);
      int v = vvec[e];
      A[e] = 1; U[e] = v; S[n] = e;
      backstep = 0;
      terminate = 0;
    }
    
    // output timing
    t2 = clock();
    float diff ((float)t2-(float)t1); diff /= CLOCKS_PER_SEC;
    if (verbose) Rcout<<"  --> Iteration took: "<<diff<<" seconds"<<endl;
    
    // take step backward
    if (!backstep) continue;
    terminate = 1;
    if (verbose) Rcout<<"Step backward: Iteration "<<iter<<"; |A| = "<<n<<"; Objective = "<<obj<<'\n';
        
    ovec = -ones<arma::vec>(N); int i = 0;    
    for (int e = 0; e < N; e++) 
      if (A[e] == 1) {      // if element v is in A
        double obj1 = compute_objective_del(obj, Acur, Wcur, Ucur, i, g, q);        
        if (obj1 > f*obj) {
          ovec[e] = obj1;
          arma::uvec tmp = find(ovec > 0);
          i += 1;
          if (tmp.size() > NOBJ) break;
        }        
      };                   
    evec = find(ovec>0);
    if (evec.size()>0){ // take the maximum
      uword e; maxval = ovec.max(e);
      A[e] = 0; U[e] = Nu;
      terminate = 0;
    }        
  };
  
  // format returned objects
  arma::uvec tmp = find(A > 0); n = tmp.size();
  arma::uvec e_vec = find(A > 0); 
  arma::uvec v_vec = U.elem(find(U < Nu));
  arma::mat Acur = Abar.rows(e_vec), 
            Wcur = W.rows(e_vec), 
            Ucur = UL.rows(v_vec);
  double obj  = compute_objective(Acur, Wcur, Ucur, g, q);
  obj = sum(g % g % q) - obj;
  arma::uvec S1   = S.elem(find(S>0)) + 1;
  arma::colvec O1 = O.elem(find(O>0));
  O1 = sum(g % g % q) - O1;

  if (verbose) Rcout<<"Done. Final objective = "<<obj<<'\n';
  // return value
  return(Rcpp::List::create(Rcpp::Named("A")=A , 
                            Rcpp::Named("U")=U , 
                            Rcpp::Named("obj")=obj,
                            Rcpp::Named("objvec")=O1,
                            Rcpp::Named("S")=S1));

}