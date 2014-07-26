library('Rcpp')
library('RcppArmadillo')

rm(list = ls())

# compile functions
sourceCpp('optimize_submodular_LS.cpp')

# generate fake data
g = runif(24)
q = runif(24)
N = 100;   
zeta = sample(2:6, N, replace = T)
usr_names = paste('Name', 1:N)
Omega = lapply(usr_names, function(s) runif(24))
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

A = matrix(runif(N*24), ncol = 24)
U = matrix(runif(N*24), ncol = 24)
W = lapply(1:N, function(i) matrix(runif(24*24), ncol = 24))
res = compute_objective_quad(A, W, U, g, q)
  
# # define submodular function
# f = function(A, U) {
#   if (length(A) == 0) return(sum(g^2))
#   D = sapply(1:length(A), function(i) A[[i]] * U[[i]])
#   C = sum(g^2) - sum((colSums(D) - g)^2)
#   return(C)
# }
#   
# # optimize submodular function
# res = optimize_submodular(f, Omega, U_list)
