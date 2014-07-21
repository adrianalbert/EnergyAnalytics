
library('Rcpp')
library('RcppArmadillo')
library('inline')

ratio_normals_CDF_src <- '
NumericVector y_(y);
double muz = as<double>(mu1), 
       sdz = as<double>(sd1), 
       muw = as<double>(mu2), 
       sdw = as<double>(sd2), 
       rho_= as<double>(rho);
int nElem = y_.size();
double r, s, a, b;
NumericVector P(nElem), yt;
/* get r,s,a,b so that  r*z/w-s=(a+x)/(b+y) */
/* z/w is distributed as s+[(a+x)/(b+y)]/r */

 s = rho_ * sdz / sdw;
 r = sdw / (sdz * sqrt(1 - rho_ * rho_));
 b = muw / sdw;
 a = (muz / sdz - rho_ * muw / sdw) / sqrt(1-rho_ * rho_);
 if (a<0) {a=-a; r=-r;}
 yt = y_ * r - s;
 
 for (int i = 0; i < nElem; i++) P[i] = F_ratio(a, b, yt[i]);

return wrap(P);'

c_functions_src = paste(readLines('zoverw_select.cpp'), collapse = '\n')

ratio_normals_CDF_C <- cxxfunction(signature(y ="numeric", 
                                           mu1 = "double", sd1 = "double", 
                                           mu2 = "double", sd2 = "double",
                                           rho = "double"),
                                 body=ratio_normals_CDF_src,
                                 include=c_functions_src,
                                 plugin = 'RcppArmadillo')

# Hackish simulation version
ratio_normals_CDF = function(y, muz, sdz, muw, sdw, rho = 0) {
  
  ratio_ecdf = ecdf(rnorm(10000, muz, sdz) / rnorm(10000, muw, sdw))
  return(ratio_ecdf(y))
  
}

ratio_normals_abs_CDF = function(y, muz, sdz, muw, sdw, rho = 0) {
  
  ratio_ecdf = ecdf(abs(rnorm(10000, muz, sdz) / rnorm(10000, muw, sdw)))
  return(ratio_ecdf(y))
  
}

# # the single value case
# ratio_normals_CDF(0.01, mu1 = 0.02, sd1 = 1e-4, mu2 = 0.04, sd2 = 3e-4, rho = 0)
# 
# # the vectorised case
# mu1 = 0.0400490748; sd1 = 0.0007254173 
# mu2 = 0.02008236;   sd2 = 0.00121371 
# mud = mu1 - mu2; sdd = sqrt(sd1^2 + sd2^2)
# 
# muz = mud; sdz = sdd;
# muw = mu1; sdw = sd1;
# rho_= 0
# 
# s = rho_ * sdz / sdw;
# r = sdw / (sdz * sqrt(1 - rho_ * rho_));
# b = muw / sdw;
# a = (muz / sdz - rho_ * muw / sdw) / sqrt(1-rho_ * rho_);
# 
# y   = (1:50)/50
# yt = y * r - s;
# 
# P = ratio_normals_CDF(y, mu1 = muz, sd1 = sdz, mu2 = muw, sd2 = sdw, rho = 0)
# ratio_ecdf_o = ecdf(rnorm(10000, muz, sdz) / rnorm(10000, muw, sdw))
# ratio_ecdf_t = ecdf(rnorm(10000, a, 1) / rnorm(10000, b, 1) / r + s)
# ratio_ecdf_n = ecdf(rnorm(10000, a, 1) / rnorm(10000, b, 1))
# plot(y,ratio_ecdf_o(y))
# points(y,ratio_ecdf_t(y), col = 'red')
# points(y,ratio_ecdf_n(yt), col = 'blue')
