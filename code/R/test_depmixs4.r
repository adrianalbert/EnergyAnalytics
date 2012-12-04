# test_depmixs4.r
#
# Tests the depmixs4 package that implements HMM estimation with covariate dependence of the parameters.

rm(list=ls())
library('depmixS4')
data(speed)
 
# simple two-state model
set.seed(1)
mod <- depmix(response = rt ~ 1, data = speed, nstates = 2, trstart = runif(4))
fm  <- fit(mod)

# model with covariates on transition matrix
mod <- depmix(rt ~ 1, data = speed, nstates = 2, family = gaussian(),
              transition = ~ scale(Pacc), instart = runif(2))
fm <- fit(mod, verbose = FALSE)
