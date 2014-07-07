# RcppHMM.r
# 
# Tests Rcpp implementation of HMM estimation via direct likelihood maximization
# 
# Author:
# Adrian Albert (adalbert@stanford.edu)
#
# Last modified: 
# June 2014

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
p        = 4
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
  coefs[i,] = 0
  # sders    = cbind(W0 = runif(nStates), W1 = runif(nStates), W3 = runif(nStates))
  transition[[i]] = coefs # list(coef = coefs, sder = sders)
}

# put together model as list
model = list(nStates = nStates,
             initDist= initDist,
             response= response,
             transition = transition)

# generate observational data, some of which is missing
TT    = 100
x_obs = runif(TT); x_obs[sample.int(TT, 10)] = NA
z_resp= matrix(runif(TT*2), nrow = TT); colnames(z_resp) = c('V0', 'V1')
y_tran= matrix(runif(TT*4), nrow = TT); colnames(y_tran) = c('W0', 'W1', 'W2', 'W3')

# -----------------------
# Test on data
# -----------------------

# response probability matrix
FF = p_response_R(x_obs[1], z_resp[1,], response)
FF

# transition probability matrix
p_transition_R(y_tran[1,], transition)

# transform natural parameters to working parameters: response
rwp = response_nat2work_R(response)

# transform natural parameters to working parameters: response
rnp = response_work2nat_R(rwp, nrow(response))

# transform natural parameters to working parameters: transition
twp = transition_nat2work_R(transition)

# transform working parameters to natural parameters: transition
tnp = transition_work2nat_R(twp, nStates, p)

# put together parameters into working vector
parvec = c(rwp, twp)

# compute log-likelihood
res = compute_log_likelihood_R(parvec, nStates, x_obs, z_resp, y_tran)
