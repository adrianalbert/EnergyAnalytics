# --------------------------------------
# d_err.r
# 
# C-Rcpp implementation of kError algorithm distance computations.
#
# Adrian Albert
# September 2013
# --------------------------------------

library('Rcpp')
library('inline')

# ------------------------------------------
# Pure R implementations for testing
# ------------------------------------------

# ___________________________________
# Distance function (chi-square)
d_err <- function(x, sx, y, sy) {
  p    = length(x)
  # stat = sum((x - y)^2 /(sx + sy))
  # d    = pchisq(stat, p)
  d    = sum((x - y)^2 /sx)
  return(d)
}

d_err_mat <- function(X, SX, Y, SY) {
  D = matrix(nrow=nrow(X), ncol=nrow(Y))
  for (x in 1:nrow(X))
    for (y in 1:nrow(Y))
      D[x,y] = d_err(X[x,], SX[x,], Y[y,], SY[y,])
  return(D)
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

cat('Initializing C functions...\n')
#d_err_c     = init_d_err()
d_err_mat_c = init_d_err_mat()

