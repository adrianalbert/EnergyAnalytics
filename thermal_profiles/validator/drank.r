# drank.r
#
# Compute metrics of agreement between rankings using the drank metric.
# Implements "On Rank Correlation and the Distance Between Rankings",
# by Ben Cartarette, 2009.
#
# Adrian Albert
# Last modified: June 2014.
# ---------------------------------------------------------

library('quadprog')

drank = function(y, X) {
  
  # compute mu - vector of means of X
  mu = colMeans(X)
  n  = nrow(X)
  m  = ncol(X)
  
  # sort columns of X according to y
  idx = order(y, decreasing = T)
  X.o = X[,idx]
  
  # compute difference in columns
  X.D = t(diff(t(X.o)))
  muD = diff(mu)
  
  # compute covariance matrix by solving QP
  lambda = 1e-5
  S   = t(X.D - muD) %*% (X.D - muD) / (n-1) + diag(m-1) * lambda
  S1  = solve(S)
  Dmat= n * S1 
  bvec= rep(0, m-1)
  Amat= diag(m-1)
  fit = solve.QP(Dmat, muD, Amat, bvec)
  theta = fit$solution
  d = n * (theta - muD) %*% S1 %*% (theta - muD)
  d = abs(sqrt(d))
  return(list(theta = theta, drank = d))  
}

# estimate p-value of rank distance
drank.pval = function(d, X, B = 100) {
  S = 0
  n = nrow(X)
  for (i in 1:B) {
    idx = sample(n, n, replace = T)
    Xs  = X[idx,]
    mus = colMeans(Xs)
    if (drank(mus, X)$drank >= d) S = S + 1
  }
  p = S / B
  return(p)
}

compute_drank = function(y, X, B = 100) {
  d = drank(y,X)$drank
  p = drank.pval(d, X, B)
  return(c(d = d, p = p))
}