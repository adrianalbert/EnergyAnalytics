# --------------------------------------
# d_err.r
# 
# C-Rcpp implementation of kError algorithm distance computations.
#
# Adrian Albert
# February 2014
# --------------------------------------

library('Rcpp')
library('inline')
library('MASS')
library('parallel')
require(inline)

# ------------------------------------------
# C implementation to speed up calculation
# ------------------------------------------

init = function() {
  
  src_obj = '   
  using namespace Rcpp;
  // input
  
  Rcpp::NumericVector cg(g);
  Rcpp::NumericVector cq(q);
  Rcpp::NumericMatrix cU(U);
  Rcpp::NumericMatrix cA(A);
  Rcpp::List cO(O);
  Rcpp::NumericMatrix cW(W);
  
  int N   = cA.nrow();
  int tau = cA.ncol();
  float obj = 0;
  float v = 0;  

  for (int t = 0; t<tau; t++){
    Rcpp::NumericMatrix Ot = cO[t];
    for (int i = 0; i<N; i++){
      for (int j = 0; j<N; j++) {
        if ( i == j ) {
          v = cW(i,t);
        } else {
          v = Ot(i,j);
          //v = 0;
        }
        obj = obj + cq(t) * cU(i,t) * cU(j,t) * ( v + cA(i,t)*cA(j,t) );
      };
      obj = obj - 2 * cq(t)*cg(t)*cU(i,t)*cA(i,t);
    };
    obj = obj + cq(t) * cg(t) * cg(t);
  };
  
  return wrap(obj);'  

  src_grad = '   
  using namespace Rcpp;
  // input  
  Rcpp::NumericVector cg(g);
  Rcpp::NumericVector cq(q);
  Rcpp::NumericMatrix cU(U);
  Rcpp::NumericMatrix cA(A);
  Rcpp::List cO(O);
  Rcpp::NumericMatrix cW(W);
  
  int N   = cA.nrow();
  int tau = cA.ncol();

  // output
  Rcpp::NumericMatrix grad(N, tau);

  float v = 0;
  for (int t = 0; t<tau; t++){
    Rcpp::NumericMatrix Ot = cO[t];
    for (int i = 0; i<N; i++){
      grad(i,t) = 0;
      for (int j = 0; j<N; j++) {
        if ( i == j ) {
          v = cW(i,t) + cA(i,t) * cA(i,t);
        } else {
          v = Ot(i,j);
          //v = 0;
        }
        grad(i,t) = grad(i,t) + cq(t) * cU(j,t) * v;
      };
      grad(i,t) = grad(i,t) - 2*cq(t)*cg(t)*cA(i,t);
    };
  };
  
  return wrap(grad);'  
  
  src_obj_reg = '   
  using namespace Rcpp;
  // input
  
  Rcpp::NumericVector cg(g);
  Rcpp::NumericVector cq(q);
  Rcpp::NumericMatrix cU(U);
  Rcpp::NumericMatrix cA(A);
  Rcpp::List cO(O);
  Rcpp::NumericMatrix cW(W);
  
  int N   = cA.nrow();
  int tau = cA.ncol();
  float obj = 0;
  float v = 0;  
  
  for (int t = 0; t<tau; t++){
    Rcpp::NumericMatrix Ot = cO[t];
    for (int i = 0; i<N; i++){
      for (int j = 0; j<N; j++) {
        if ( i == j ) {
          v = cW(i,t) + cA(i,t) * cA(i,t);
        } else {
          v = Ot(i,j);
          //v = 0;
        }
        obj = obj + cq(t) * cU(i,t) * cU(j,t) * ( v + cA(i,t)*cA(j,t) );
      };
      obj = obj - 2 * cq(t)*cg(t)*cU(i,t)*cA(i,t);
      obj = obj + fabs(cU(i,t));
    };
    obj = obj + cq(t) * cg(t) * cg(t);
  };

  return wrap(obj);'  
  
  
  fun <- cxxfunction(signature(U = "numeric", g = "numeric",  q = "numeric", A = "numeric", 
                               O = "list", 
                               W = "numeric"),
                     body = src_obj, 
                     plugin = "Rcpp")  
  fun_reg <- cxxfunction(signature(U = "numeric", g = "numeric",  q = "numeric", A = "numeric", 
                               O = "list", 
                               W = "numeric"),
                     body = src_obj_reg, 
                     plugin = "Rcpp")  
  grad <- cxxfunction(signature(U = "numeric", g = "numeric",  q = "numeric", A = "numeric", 
                               O = "list", 
                               W = "numeric"),
                     body = src_grad, 
                     plugin = "Rcpp")  
  return ( list(fun = fun, fun_reg = fun_reg, grad = grad) )

}

# _______________________
# Initialize C functions

cat('Initializing C functions...\n')
func   = init()
f_obj  = func$fun
f_grad = func$grad
f_reg  = func$fun_reg

