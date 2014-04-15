# RcppHMM.r
# 
# Tests Rcpp implementation of HMM estimation via direct likelihood maximization
# 
# Author:
# Adrian Albert (adalbert@stanford.edu)
#
# Last modified: 
# April 2014

# -----------------------
# Initializations
# -----------------------

rm(list=ls())

# ________________________
# Libraries

library('Rcpp')
library('RcppArmadillo')
# library('benchmark')

# ________________________
# Load C++ functions

sourceCpp('HMM.cpp')

# -----------------------
# Test on data
# -----------------------
