# #########################################################################
# test_kError.r
# --------
#
# tests kError implementation. 
# 
# Adrian Albert
# Last modified: September 2013.
# #########################################################################

rm(list = ls())
options(error = recover)
library('ggplot2')
library('reshape')
source('./clustering/kError.r')

# _______________
# Generate data

N  = 100
p  = 20
K  = 3
mu = t(sapply(1:K, function(k) runif(p)))
s  = 0.05*t(sapply(1:(K), function(k) runif(p)))
S  = lapply(1:K, function(k) {
  t(sapply(1:N, function(j) s_cur = s[k,] * runif(p) ))
#   t(sapply(1:N, function(j) s_cur = s[k,] * runif(p) ))
  t(sapply(1:N, function(j) s_cur = s[k,] ))
})
X  = lapply(1:K, function(k) {  
  t(sapply(1:N, function(j) {    
    rnorm(p, mu[k,], S[[k]][j,])
  }))
})
S  = do.call('rbind', S)
X  = do.call('rbind', X)
lab= rep(1:K, each=N)

# ________________
# Test clustering 

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
} 

df        = as.data.frame(X)
df$lab    = res$assignment
df$Obs    = 1:nrow(df)
df        = melt(df, id.vars = c('Obs', 'lab'))
p         = ggplot(df, aes(variable, value, color = as.factor(lab)))
p         = p + geom_point() + geom_line()
print(p)


