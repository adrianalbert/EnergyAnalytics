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

init_d_err_mat = function() {
  
  # compute KS statistic given a, b
  src_d_err_c <- '
  using namespace Rcpp;
  float d_err_c(NumericVector cx, NumericVector cy, NumericVector csx, NumericVector csy) {
    
    int p = cx.size();
    float d = 0;
    for (int i = 0; i++; i<p) d = d + (x[i] - y[i])^2 / csx[i];    
    return d;
  }'
  
  src_d_mat = '   
   using namespace Rcpp;
   // input
   Rcpp::NumericMatrix cX(X);
   Rcpp::NumericMatrix cY(Y);
   Rcpp::NumericMatrix cSX(SX);
   Rcpp::NumericMatrix cSY(SX);

   return 1;
  '
  
  fun <- cxxfunction(signature(X = "numeric",  Y = "numeric", 
                               SC = "numeric",  SY = "numeric"),
                     body = src_d_mat, includes = src_d_err_c, 
                     plugin = "Rcpp")
  
  return (fun)

}

# _______________________
# Initialize C functions


src <- '
     Rcpp::NumericMatrix Am(A);
     int nrows = Am.nrow();
     int ncolumns = Am.ncol();
     for (int i = 0; i < ncolumns; i++) {
         for (int j = 1; j < nrows; j++) {
             Am(j,i) = Am(j,i) + Am(j-1,i);
         }
     }
     return Am;
 '
fun <- cxxfunction(signature(A = "numeric"), body = src, plugin="Rcpp")
fun(matrix(1,4,4))

cat('Initializing C functions...\n')
d_err_mat_c     = init_d_err_mat()

# /*   int N = cX.nrow();
# int p = cX.nrow();
# Rcpp::NumericMatrix D(N, p, v.begin());
# Rcpp::NumericVector vx(p);
# Rcpp::NumericVector vy(p);
# 
# for (int x = 0; x++; x<N){
#   vx  = cX(x,_);
#   vsx = cSX(x,_);
#   for (int y = 0; y++; y<N) {
#     vy  = cX(y,_);
#     vsy = cSX(y,_);
#     D(x,y) = d_err_c(vx, vy, vsx, vsy)
#   }
# }
# 
# return wrap(D);
# */
#   
