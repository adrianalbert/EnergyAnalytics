# #########################################################################
# test_kError.r
# --------
#
# tests kError implementation. 
# 
# Adrian Albert
# Last modified: February 2014.
# #########################################################################

rm(list = ls())
options(error = recover)
library('ggplot2')
library('reshape')
source('../clustering/kError.r')

# _______________
# Generate data

p  = 24
K  = 3
N  = sample(200:500, K)
mu_orig  = t(sapply(1:K, function(n) runif(p)))
lab_orig = unlist(sapply(1:K, function(k) rep(k, N[k])))
S_orig   = lapply(1:K, function(k) {
  s = 0.1*matrix((-1 + 2*runif(p^2)), ncol = p)
  s = s %*% t(s)
  })
X  = lapply(1:K, function(k) {  
  mvrnorm(N[k], mu_orig[k,], S_orig[[k]])
})
X  = do.call('rbind', X)
S  = lapply(1:K, function(k) {  
  S1 = list()
  for (i in 1:N[k]) {
    S   = S_orig[[k]]
    res = eigen(S)
    l1  = res$values * (0.9 + 0.2 * runif(nrow(S)))
    S1[[length(S1) + 1]]  = res$vectors %*% diag(l1) %*% t(res$vectors)
  }
  return(S1)
})
S = unlist(S, recursive = F)

# __________________________
# Test distance computation

d_err_cov(X[1,], S[[1]], mu_orig[1,])
D = d_err_mat(X, S, mu_orig)

# center computation
res = computeCenter(X, S)

# ________________
# Test clustering 

source('../clustering/kError.r')
res = kError(X, S, K, iter = 10)

# ___________________
# Visualize clusters
if (p == 2) {
  df = as.data.frame(X)
  names(df) = c('x', 'y')
  df$lab = res$assignment
  p = ggplot(df, aes(x, y, color = as.factor(lab)))
  p = p + geom_point()
  print(p)
} else {
  df        = as.data.frame(X)
  df$lab    = as.factor(res$assignment)
  df$Obs    = 1:nrow(df)
  df        = melt(df, id.vars = c('Obs', 'lab'))
  ds        = as.data.frame(t(sapply(S, function(x) diag(x))))  
  ds$Obs    = 1:nrow(ds)
  ds        = melt(ds, id.vars = 'Obs')
  names(ds)[3] = 'se'
  df        = merge(df, ds, by = c('Obs', 'variable'))
  p         = ggplot(df, aes(variable, value, color = lab, group = Obs))
  p         = p + facet_wrap(~lab, ncol = 3)
  p         = p + geom_point() + geom_line()
  p         = p + geom_ribbon(aes(ymin=value-se, ymax=value+se), width=.1, color = 'gray', alpha = 0.1)

  print(p)
}
