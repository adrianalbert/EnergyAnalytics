# #########################################################################
# kError.r
# --------
#
# kError K-Means variant for clustering vector data with associated errors. 
# This file implements the independent errors case (no covariance across observations).
# 
# Inputs:
# -------
# X:    N x p matrix of N data points in R^p
# S:    N x p matrix of associated measurement errors 
# K:    number of clusters sought / initial assignment of clusters
# iter: maximum number of iterations
# eps:  accuracy threshold for convergence
# 
# Outputs: list with following fields:
# -------
# ASS:  N x 1 cluster assignment vector
# CX:   K x p matrix of centers
# CS:   K x p matrix of center errors
# OBJ:  clustering objective
# 
# Adrian Albert
# Last modified: February 2014.
# #########################################################################

source('d_err.r', chdir = T)

# ___________________________________
# Compute center of cluster
computeCenter <- function(X, S) {
  
  # compute cluster covariance matrix
  S1 = lapply(S, function(s) solve(s))
  SS = S1[[1]]
  if (length(S1)>1) for (i in 2:length(S1)) SS = SS + S1[[i]]
  CS = solve(SS)
  
  # compute cluster mean
  if (length(S1) == 1) X = t(as.matrix(X))
  Sx = lapply(1:length(S1), function(i) S1[[i]] %*% X[i,])
  xx = Sx[[1]]
  if (length(Sx)>1) for (i in 2:length(Sx)) xx = xx + Sx[[i]]
  CX = CS %*% xx
  CX = as.numeric(CX)
  
  return(list(CX = CX, CS = CS))
}

kError <- function(X, S, K, iter = 10, verbose = T) {
  
  # process inputs
  N = nrow(X)
  p = ncol(X)
  if (length(K) > 1) {
    ASS = K 
  } else {
    ASS = sample(1:K, N, replace = TRUE)
  }
  
  # main interation loop
  i = 0; done = 0; obj = 0; obj_vec = c(); 
  while (i < iter & !done ){    
    i       = i + 1
    ASS_old = ASS
    obj_old = obj; 
        
    # M-step
    CX = matrix(nrow=K, ncol=p)
    CS = list()
    for (k in 1:K) {
      idx = which(ASS == k)
      x_k = X[idx,]
      s_k = S[idx]
      res = computeCenter(x_k, s_k)
      CX[k,] = res$CX
      CS[[1+length(CS)]] = res$CS
    }
    
    # E-step  
    D   = d_err_mat(X, S, CX)
    ASS = apply(D, 1, which.min)
    obj = sum(sapply(1:N, function(i) D[i,ASS[i]])); 
    obj_vec = c(obj_vec, obj); 
    
    # check convergence
    per_changed = (sum(ASS != ASS_old) / length(ASS));
    per_obj_chg = abs(obj - obj_old) / obj_old
    done = (per_changed <= 0.01) | (per_obj_chg < 0.005)
  
    if (verbose) print(paste("K =", k, "; Iteration", i, ': Objective =', obj, '; %. changed =', per_changed, '; %obj changed =', per_obj_chg))    
  }
  
  return(list(objective = obj_vec, assignment = ASS, centers = CX, errors = CS))
  
}