# metrics.r
#
# Compute metrics to validate individual model performance. 
#
# Adrian Albert
# Last modified: June 2014.
# ---------------------------------------------------------

library(fpc)
library(fdrtool)
library(VGAM)

# given a data frame containing indications of "states" and (ground truth) dependent and independent variables, 
# compute state-based linear response
compute_ground_truth_response = function(df, dep.var = 'AC', indep.var = 'TemperatureD'){
  # compute ground-truth distributions
  nStates = length(unique(df$state))
  fmla = as.formula(paste(dep.var, indep.var, sep='~'))
  fit  = lapply(1:nStates, function(s) {
    fit = lm(fmla, data = subset(df, state == s))
    sm  = summary(fit)  
    dat = sm$coefficients[,c('Estimate', 'Std. Error')]; dat = as.data.frame(dat); dat$state = s; dat$parameter = rownames(dat)
    return(dat)
  })
  means = sapply(fit, function(l) l$Estimate); rownames(means) = rownames(fit[[1]])
  sderr = sapply(fit, function(l) l[,'Std. Error']); rownames(sderr) = rownames(fit[[1]])    
  
  return(list(means = as.data.frame(means), stderr = as.data.frame(sderr)))
}

# error function
erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1

# CDF of folded normal
F_fnorm <- function(y, mu=0, sigma=1) 1/2*(erf((y+mu)/(sqrt(2)*sigma)) + erf((y-mu)/(sqrt(2)*sigma)))

# expectation of folded normal
E_fnorm <- function(mu=0, sigma=1) sigma * sqrt(2 / pi) * exp(-mu^2/(2*sigma^2)) + mu * (1 - 2*pnorm(-mu/sigma))

# given two Gaussian models for (ground truth and estimated response), compare their relative closeness
compare_magnitudes = function(ground.resp, model.resp, indep.var = 'TemperatureD') {
  
  mu.ground = as.numeric(ground.resp$means[indep.var,]);   mu.model  = as.numeric(model.resp$means[indep.var,])
  s2.ground = as.numeric(ground.resp$stderr[indep.var,]^2); s2.model  = as.numeric(model.resp$stderr[indep.var,]^2)
  nStates   = length(mu.ground)
  
  # compute Bhattcharya distance
  # bd = bhattacharyya.dist(mu.ground, mu.model, diag(s2.ground), diag(s2.model))
  bd = sapply(1:nStates, function(s) bhattacharyya.dist(mu.ground[s], mu.model[s], s2.ground[s], s2.model[s]))
  
  # compute p-value
  mu.diff = mu.ground - mu.model; sd.diff = sqrt(s2.ground + s2.model)
  dx= -0.025 + 0.001*(1:50)
  P = sapply(dx, function(d) F_fnorm(d, mu=mu.diff, sigma=sd.diff))
  P = as.data.frame(t(P)); P$dx = dx;
  E = E_fnorm(mu=mu.diff, sigma=sd.diff)
}

ddf = melt(P, id.vars = 'dx')
ddf = rbind(data.frame(mu = mu.ground, sd = sqrt(s2.ground), type = 'Ground', state=1:length(mu.ground)),
            data.frame(mu = mu.model, sd = sqrt(s2.model), type = 'Model', state=1:length(mu.ground)))

ggplot(ddf, aes(shape = type, color = state)) + stat_function(data = ddf, fun = dnorm, arg = list( mean = mu, sd = sd))
ddf = melt(df[,-1], id.vars = c('state', 'TemperatureF', 'TemperatureD'))

ddf$state = as.factor(ddf$state)
ggplot(ddf, aes(TemperatureF, value, color = variable)) + geom_point(size=2) + facet_grid(state~variable)
