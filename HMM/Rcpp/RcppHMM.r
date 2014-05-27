# RcppHMM.r
# 
# Tests Rcpp implementation of HMM estimation via direct likelihood maximization
# 
# Author:
# Adrian Albert (adalbert@stanford.edu)
#
# Last modified: 
# May 2014

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

# ________________________
# Set up dummy model

# model size and initial distribution
nStates  = 3
initDist = runif(nStates); initDist = initDist / sum(initDist)

# response coefficients
coefs    = cbind(V0 = runif(nStates), V1 = runif(nStates), Sigma = runif(nStates))
sders    = cbind(V0 = runif(nStates), V1 = runif(nStates), Sigma = runif(nStates))
# response = list(coef = coefs, sder = sders)
response = coefs

# transition coefficients
transition = list()
for (i in 1:nStates) {
  coefs    = cbind(W0 = runif(nStates), W1 = runif(nStates), W2 = runif(nStates), W3 = runif(nStates))
  # sders    = cbind(W0 = runif(nStates), W1 = runif(nStates), W3 = runif(nStates))
  transition[[i]] = coefs # list(coef = coefs, sder = sders)
}

# put together model as list
model = list(nStates = nStates,
             initDist= initDist,
             response= response,
             transition = transition)

# generate some observational data
TT    = 100
x_obs = runif(TT)
z_resp= matrix(runif(TT*2), nrow = TT); colnames(z_resp) = c('V0', 'V1')
y_tran= matrix(runif(TT*4), nrow = TT); colnames(y_tran) = c('W0', 'W1', 'W2', 'W3')

# -----------------------
# Test on data
# -----------------------

# response probability matrix
p_response_R(x_obs[1], z_resp[1,], response)

# transition probability matrix
p_transition_R(y_tran[1,], transition)

# res = compute_log_likelihood(model, x_obs, y_resp, z_tran)
