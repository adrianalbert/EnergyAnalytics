
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

/* muz=30.5;muw=32;sigz=5;sigw=4;rho=.8;*/

 s = rho_ * sdz / sdw;
 r = sdw / (sdz * sqrt(1 - rho_ * rho_));
 b = muw / sdw;
 a = (muz / sdz - rho_ * muw / sdw) / sqrt(1-rho_ * rho_);
 if (a<0) {a=-a; r=-r;}
 yt = y_ * r - s;
 for (int i = 0; i < nElem; i++) P[i] = F(a, b, yt[i]);

 return wrap(P);'

c_functions_src = paste(readLines('zoverw_select.cpp'), collapse = '\n')

ratio_normals_CDF <- cxxfunction(signature(y ="numeric", 
                                           mu1 = "double", sd1 = "double", 
                                           mu2 = "double", sd2 = "double",
                                           rho = "double"),
                                 body=ratio_normals_CDF_src,
                                 include=c_functions_src,
                                 plugin = 'RcppArmadillo')

# # the single value case
# ratio_normals_CDF(0.1, mu1 = 30.5, sd1 = 5, mu2 = 32, sd2 = 4, rho = .8)
# 
# # the vectorised case
# y   = (1:100)/100
# ratio_normals_CDF(y, mu1 = 2, sd1 = 1e-1, mu2 = 3, sd2 = 4e-1, rho = 0)
