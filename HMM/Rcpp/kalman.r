
library('Rcpp')
library('RcppArmadillo')
library('inline')

kalmanSrc <- '
mat Z = as<mat>(ZS);
Kalman K;
mat Y = K.estimate(Z);
return wrap(Y);'

kalmanClass = paste(readLines('kalman.cpp'), collapse = '\n')

KalmanCpp <- cxxfunction(signature(ZS="numeric"),
                         body=kalmanSrc,
                         include=kalmanClass,
                         plugin="RcppArmadillo")