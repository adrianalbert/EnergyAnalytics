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

# ------------------------------------------
# Pure R implementations for testing
# ------------------------------------------

# ___________________________________
# Distance function (chi-square)
d_err <- function(x, sx, y, sy) {
  p    = length(x)
  d    = sum((x - y)^2 /sx)
  return(d)
}

# _____________________________________________________________
# Distance function (Mahalanobis) for non-diagonal covariance
d_err_cov <- function(x, S, mu) {
  S1   = solve(S)
  d    = sqrt(t(x - mu) %*% S1 %*% (x - mu))
  return(d)
}

# ________________________
# Compute distance matrix 
d_err_mat <- function(X, SX, M) {
  idx = expand.grid(x=1:nrow(X), m=1:nrow(M))
  D   = lapply(1:nrow(idx), function(i){
    x = X[idx[i,'x'],]
    S = SX[[idx[i,'x']]]
    mu= M[idx[i,'m'],]
    d = d_err_cov(x, S, mu)
    return(d)
  })
  D = matrix(unlist(D), nrow = nrow(X), byrow=F)
}

# ------------------------------------------
# C implementation to speed up calculation
# ------------------------------------------

init_d_err = function() {
  
  # compute KS statistic given a, b
  src_d_err_c <- '
  using namespace Rcpp;
  Rcpp::NumericVector cx(x);
  Rcpp::NumericVector cy(y);
  Rcpp::NumericVector csx(sx);
  Rcpp::NumericVector csy(sy);
  
  Rcpp::NumericVector d = (cx - cy)*(cx - cy) / csx;
  double s = std::accumulate(d.begin(), d.end(), 0.0);
  
  return wrap(s);'
  
  fun <- cxxfunction(signature(x = "numeric",  sx = "numeric", 
                               y = "numeric",  sy= "numeric"),
                     body = src_d_err_c, 
                     plugin = "Rcpp")  
  return (fun)  
}


init_d_err_mat = function() {
  
  # compute KS statistic given a, b
  src_d_err_c <- '
  using namespace Rcpp;
  double d_err_c(NumericVector cx, NumericVector csx, NumericVector cy, NumericVector csy) {
  
    using namespace Rcpp;    
    NumericVector d = (cx - cy)*(cx - cy) / csx;
    double s = std::accumulate(d.begin(), d.end(), 0.0);
    
    return s;
  }'
  
  src_d_mat = '   
  using namespace Rcpp;
  // input
  
  Rcpp::NumericMatrix cX(X);
  Rcpp::NumericMatrix cY(Y);
  Rcpp::NumericMatrix cSX(SX);
  Rcpp::NumericMatrix cSY(SY);
  
  int N1 = cX.nrow();
  int N2 = cY.nrow();
  NumericMatrix D(N1, N2);

  for (int x = 0; x<N1; x++){
    for (int y = 0; y<N2; y++) {
      D(x,y) = d_err_c(cX(x,_), cSX(x,_), cY(y,_), cSY(y,_));
    };
  };
  
  return wrap(D);'  
  fun <- cxxfunction(signature(X = "numeric",  SX = "numeric", 
                               Y = "numeric",  SY = "numeric"),
                     body = src_d_mat, includes = src_d_err_c, 
                     plugin = "Rcpp")  
  return (fun)

}

# _______________________
# Initialize C functions

# cat('Initializing C functions...\n')
# d_err_c     = init_d_err()
# d_err_mat_c = init_d_err_mat()

