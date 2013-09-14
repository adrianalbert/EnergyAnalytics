# --------------------------------------
# dks.r
# 
# Computes two-sample KS-distance between two distributions.
#
# Adrian Albert
# February 2013
# --------------------------------------

# ------------------------------------------
# Pure R implementations for testing
# ------------------------------------------

# _______________________________________________________________________
# Compute a shift & scaled version of given CDF
#
# Fy    : input CDF
# x_dom : domain of Fy; Fy is 0 to the left of x_dom and 1 to the right
# x     : coordinate after transformation with scale and shift
# a,b   : scale and shift parameters
cdf.ss = function(Fy, y_dom, x, a, b) {
  
  z = (x - b) / a
  
  # if outside of domain of Fy, it's either 0 or 1
  if (z < min(y_dom)) return(0)
  if (z > max(y_dom)) return(1)
  
  # interpolate value in between
  idx  = which(z >= y_dom)
  i    = idx[length(idx)]
  Fz   = Fy[i] + (Fy[i+1] - Fy[i]) / (y_dom[i+1] - y_dom[i]) * (z - y_dom[i])
  return(Fz)
}

compute.KS = function(theta, Fx, x_dom, Fy, y_dom) {
  
  a = theta[1]
  b = theta[2]
  cdf.diff = sapply(1:length(x_dom), function(i) abs(Fx[i] - cdf.ss(Fy, y_dom, x_dom[i], a, b)))
  idx.max  = which.max(cdf.diff)
  
  res = cdf.diff[idx.max]
  names(res) = 'KS'
  
  return(res)
  
}

# ------------------------------------------
# C implementation to speed up calculation
# ------------------------------------------

init_cdf_ss = function() {

  # shift & scale CDF given a, b
  src <- '
   // input variables
   Rcpp::NumericVector cFy(Fy);
   Rcpp::NumericVector cy_dom(y_dom);
   float cx = Rcpp::as<float>(x);
   float ca = Rcpp::as<float>(a);
   float cb = Rcpp::as<float>(b);
   int n    = cFy.size();
 
   // for output
   Rcpp::NumericVector Fz(1);

   // treat end cases
   float z = (cx - cb) / ca;
   if (z < cy_dom[0]) {
      Fz[0] = 0;
      return Fz;
   }
   if (z > cy_dom[n-1]) {
      Fz[0] = 1;
      return Fz;
   }
    
   // interpolate values in between
   int i = 0;
   while (z > cy_dom[i]) i += 1;
   i = i-1;
   Fz[0] = cFy[i] + (cFy[i+1] - cFy[i]) / (cy_dom[i+1] - cy_dom[i]) * (z - cy_dom[i]);

   return Fz;'
  
  fun <- cxxfunction(signature(Fy = "numeric", y_dom = "numeric", x = "numeric", a = "numeric", b = "numeric"),
                     src, plugin = "Rcpp")
  
  return (fun)
}

init_compute_KS = function() {
  
  # compute KS statistic given a, b
  src_cdf_ss <- '
  using namespace Rcpp;
  float cdf_ss(NumericVector cFy, NumericVector cy_dom, float cx, float ca, float cb) {
    
    int n = cFy.size();
    float Fz;
    // treat end cases
    float z = (cx - cb) / ca;
    if (z < cy_dom[0]) {
      Fz = 0;
      return Fz;
    }
    if (z > cy_dom[n-1]) {
      Fz = 1;
      return Fz;
    }
    
    // interpolate values in between
    int i = 0;
    while (z > cy_dom[i]) i += 1;
    i = i-1;
    Fz = cFy[i] + (cFy[i+1] - cFy[i]) / (cy_dom[i+1] - cy_dom[i]) * (z - cy_dom[i]);
    
    return Fz;
  }'
  
  src_KS = '   
   using namespace Rcpp;
   // input
   NumericVector ctheta(theta);
   NumericVector cFx(Fx);
   NumericVector cx_dom(x_dom);
   NumericVector cFy(Fy);
   NumericVector cy_dom(y_dom);
   int n = cFy.size();
   float a = ctheta[0] , b = ctheta[1];
 
   // compute KS
   float ks_max = 0, cur;
   for (int i=0; i<n; i++) {
     cur = fabs(cFx[i] - cdf_ss(cFy, cy_dom, cx_dom[i], a, b));
     if (cur > ks_max) ks_max = cur;
   }

   return wrap(ks_max);
  '
  
  fun <- cxxfunction(signature(theta = "numeric", 
                               Fx = "numeric",  x_dom = "numeric", 
                               Fy = "numeric",  y_dom = "numeric"),
                     src_KS, includes = src_cdf_ss, 
                     plugin = "Rcpp")
  
  return (fun)

}

# _______________________________________________________________________
# Compute optimum shift and scale parameters

compute.dKS = function(Fx, x_dom, Fy, y_dom, lower = c(0.5, -2), upper = c(2, 2)) {
  
  f_wrap = function(theta) {
    compute.KS.C(theta, Fx, x_dom, Fy, y_dom)
  }
  
  # perform optimization
  theta0 = c(a = 1, b = 0)
  res = optim(theta0, f_wrap, 
              lower=lower, 
              upper=upper, 
              method="L-BFGS-B", 
              control=list(maxit=100))  
  theta = res$par
  dKS   = res$value
  
  return(c(theta[1], theta[2], dKS = dKS))
}

# _______________________
# Initialize C functions

library('Rcpp')
library('inline')

cat('Initializing C functions...\n')
cdf.ss.C     = init_cdf_ss()
compute.KS.C = init_compute_KS()
