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
# Last modified: September 2013.
# #########################################################################

# ___________________________________
# Distance function (chi-square)
d_err <- function(x, sx, y, sy) {
  p    = length(x)
  # stat = sum((x - y)^2 /(sx + sy))
  # d    = pchisq(stat, p)
  d    = sum((x - y)^2 /sx)
  return(d)
}

d_err_mat <- function(X, SX, Y, SY) {
  D = matrix(nrow=nrow(X), ncol=nrow(Y))
  for (x in 1:nrow(X))
    for (y in 1:nrow(Y))
      D[x,y] = d_err(X[x,], SX[x,], Y[y,], SY[y,])
  return(D)
}

# ___________________________________
# Compute center of cluster
computeCenter <- function(X, S) {
  CS = 1 / colSums(1/S) 
  CX = CS * colSums(X / S)
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
  i = 0; done = 0;
  while (i < iter & !done ){    
    i       = i + 1
    ASS_old = ASS
    if (verbose) cat(paste("Iteration", i, '\n'))
    
    # M-step
    CX = matrix(nrow=K, ncol=p)
    CS = matrix(nrow=K, ncol=p)
    for (k in 1:K) {
      idx = which(ASS == k)
      x_k = X[idx,]
      s_k = S[idx,]
      res = computeCenter(x_k, s_k)
      CX[k,] = res$CX
      CS[k,] = res$CS
    }
    # E-step
    D   = d_err_mat(X, S, CX, CS)
    ASS = apply(D, 1, which.min) 
    
    # check convergence
    done = sum((ASS - ASS_old)^2) == 0;
  }
  
  return(list(assignment = ASS, centers = CX, errors = CS))
  
}