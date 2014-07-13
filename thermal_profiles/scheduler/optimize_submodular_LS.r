library('Rcpp')
library('RcppArmadillo')

rm(list = ls())

setwd('~/EnergyAnalytics/thermal_profiles/scheduler/')

# compile functions
sourceCpp('optimize_submodular_LS.cpp')

# ___________________________________________
# Test objective function 

# generate fake data
N = 100
g = 20*runif(24)
q = runif(24)
A = matrix(runif(N*24), ncol = 24)
U = matrix(runif(N*24), ncol = 24)
W = lapply(1:N, function(i) matrix(runif(24*24), ncol = 24))

# test objective
res = compute_objective_quad(A, W, U, g, q)

# ___________________________
# Test optimization algorithm

# generate fake data
N = 1000;   
zeta = sample(2:6, N, replace = T)
usr_names = paste('Name', 1:N)
Omega = lapply(usr_names, function(s) list(a = runif(24), w = matrix(runif(24*24), ncol=24)))
names(Omega) = usr_names
U_list= lapply(zeta, function(z) {
  eta   = sample(1:23, z, replace = T)
  gamma = sample(1:5, z, replace = T)
  U     = lapply(1:z, function(s) {
    u = rep(0, 24)
    u[eta[s]:min(eta[s] + gamma[s], 24)] = 1
    return(u)
  })
  return(U)
})
names(U_list) = usr_names
UL = do.call('rbind', sapply(U_list, function(l) do.call('rbind', l)))
UA = rep(1:N, sapply(U_list, length))

# test optimization
params = list(g = g, q = q, eps = 0.01)
res = optimize_submodular_LS(Omega, list(UL = UL, UA = UA), params)
