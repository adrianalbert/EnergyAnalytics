# compareMCs.r
#
# Compares two MCs using hypothesis testing. 
# 
# Adrian Albert
# Last modified: April 2013.
# -----------------------------------------------------------------------

# ____________________________________
# Pre-defined Markov Chain structures

define_MC = function(K, constrMC = 'telescope') {
  # identity model
  if (constrMC == 'identity') {
    A = diag(1,K)
  }
  # telescope model
  if (constrMC == 'complete') {
    A = matrix(1, nrow = K, ncol = K)
  }
  # telescope model
  if (constrMC == 'telescope') {
    A = matrix(0, nrow = K, ncol = K)
    for (i in 1:K) {
      if (i<K) A[i,i+1] = 1
      A[i,i]   = 1
      A[i,1]   = 1
    }            
  }
  # forward-jump model
  if (constrMC == 'fwdjump') {
    A = matrix(0, nrow = K, ncol = K)
    for (i in 1:K) {
      A[i,i:K] = 1                    
    }    
    A[K,1] = 1
  }
  # outlier state model
  if (constrMC == 'outlier') {
    A = matrix(0, nrow = K, ncol = K)
    for (i in 1:K) {
      if (i<K) A[i,i+1] = 1
      A[i,i]   = 1
      A[i,1]   = 1
      A[i,K]   = 1
      A[K,i]   = 1
    }            
  }
  if (constrMC == 'null') A = matrix(0, ncol = K, nrow = K) else 
    A = A / rowSums(A)
  return(A)
}

# ______________________________________________________________
# Likelihood Ratio test for MCs
# A:  reference model
# B:  data-derived model
# TT: lenght of sequence
# returns p-value of test

test_MC_LR = function(A, B, TT, Nx = NULL) {
  # approximate observed frequency of states if not given
  if (is.null(Nx)) {
    Nx = as.real(eigen(t(B))$vectors[,1])
    Nx = TT * Nx / sum(Nx)
    N  = B * Nx
  }
  
  # compute stationary distribution of reference model A
  pi.A = as.real(eigen(t(A))$vectors[,1])
  pi.A = pi.A / sum(pi.A)
  
  # compute log-likelihood ratio test statistic
  ll   = -2 * sum( (N - TT * pi.A * A) * log(A) )
  
  # perform hypothesis test
  K    = nrow(A)
  df   = K * (K - 1)
  pval = pchisq(ll, df)
  
  return(pval)
}

# ______________________________________________________________
# Chi-squared test for MCs
# A: reference model
# B: data-derived model
# 
# returns p-value of test

test_MC_Chi = function(A, B, TT, Nx = NULL, two.sided = T) {
  # approximate observed frequency of states if not given
  if (is.null(Nx)) {
    Nx = as.real(eigen(t(B))$vectors[,1])
    Nx = TT * Nx / sum(Nx)
    N  = B * Nx
  }
    
  # replace 0 values with very small constants
  eps  = 1e-6
  idx  = which(A<eps)
  if (length(idx)>0) A[idx] = eps
  
  # compute test statistic
  ll   = sum(Nx * (B - A)^2 / A)
  
  # perform hypothesis test
  K    = nrow(A)
  df   = K * (K - 1)
  pval = pchisq(ll, df, lower.tail = two.sided)

  return(pval)
}

# ______________________________________________________________
# Compute most likely type of model for B from list models

library('gtools')
assignMCtypes = function(B, TT, Nx = NULL, models, noPerm = 100) {
  
  # sample some permutations of the nodes in B
  K   = nrow(B)
  res = lapply(1:noPerm, function(j) {
    perm = sample(K, K)    
    pvals= rep(NA,length(models))
    names(pvals) = names(models)
    for (m in names(models)) pvals[m] = test_MC_Chi(models[[m]], B[perm,perm], TT, Nx)
    return(t(pvals))
  })
  
  return(res)  
}

