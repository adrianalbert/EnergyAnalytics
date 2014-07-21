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
library(pracma)

# ------------------------------------------
# Statistics functions
# ------------------------------------------

source('~/EnergyAnalytics/thermal_profiles/validator/ratio_normals.r', chdir = T)

# error function
erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1

# CDF of folded normal
F_fnorm <- function(y, mu=0, sigma=1) 1/2*(erf((y+mu)/(sqrt(2)*sigma)) + erf((y-mu)/(sqrt(2)*sigma)))

# expectation of folded normal
E_fnorm <- function(mu=0, sigma=1) sigma * sqrt(2 / pi) * exp(-mu^2/(2*sigma^2)) + mu * (1 - 2*pnorm(-mu/sigma))

# ------------------------------------------
# Metrics extraction functions
# ------------------------------------------

# given a data frame containing indications of "states" and (ground truth) dependent and independent variables, 
# compute state-based linear response
compute_ground_truth_response = function(df, dep.var = 'AC', indep.var = 'TemperatureD'){
  # compute ground-truth distributions
  nStates = length(unique(df$state))
  fmla = as.formula(paste(dep.var, indep.var, sep='~'))
  fit  = lapply(1:nStates, function(s) {
    fit = lm(fmla, data = subset(df, state == s))
    sm  = summary(fit)  
    dat = sm$coefficients[indep.var,c('Estimate', 'Std. Error')]; 
    names(dat) = c("mu", "sd")
    return(dat)
  })
  
  return(fit)
}

# mean and variance of difference of two gaussians
diff_gaussians = function(m1, m2){  
  if (class(m1) == 'numeric' && class(m2) == 'numeric') {
    return(c(mu = m1['mu'] - m2['mu'], sd = sqrt(m1['sd']^2 + m2['sd']^2)))
  }  
  l = list()
  for (i in 1:length(m1)) {
    z = c(m1[[i]]['mu'] - m2[[i]]['mu'], sqrt(m1[[i]]['sd']^2 + m2[[i]]['sd']^2))
    names(z) = c('mu', 'sd')
    l[[length(l) + 1]] = z
  }  
  return(l)
}

# given two Gaussian models for (ground truth and estimated response), compare their relative closeness
# note that we care about both the magnitude and the sign of the relative differences!
compare_gaussian_models = function(ground, estimate, 
                                   support.mag = NULL, support.rel = NULL, support.rel.abs = NULL){
  
  # compute expected relative error 
  expect.rel.err = abs( (ground["mu"] - estimate["mu"]) / ground["mu"] )
  
  # compute distribution of difference
  mu.d = ground["mu"] - estimate["mu"]; sd.d = sqrt(ground["sd"]^2 + estimate["sd"]^2)
  if (is.null(support.mag)) {
    eta  = seq(mu.d - 3 * sd.d, mu.d + 3 * sd.d, length.out = 50)  
  } else {
    eta  = support.mag
  }
  P.mag = pnorm(eta, mean = mu.d, sd = sd.d); 
  
  # compute relative difference
  if (is.null(support.rel)) {
    rel = seq(-1, 1, length.out = 50)
  } else {
    rel  = support.rel
  }
  P.rel = ratio_normals_CDF(rel, mu.d, sd.d, ground["mu"], ground["sd"], 0)
  
  # compute absolute relative difference
  if (is.null(support.rel.abs)) {
    rel.abs = seq(0, 2, length.out = 50)
  } else {
    rel.abs  = support.rel.abs
  }
  P.rel.abs = ratio_normals_abs_CDF(rel.abs, mu.d, sd.d, ground["mu"], ground["sd"], 0)
  
  return(list(P.magnitude = P.mag, support.mag = eta,
              P.relative = P.rel, support.rel = rel,
              P.rel.abs  = P.rel.abs, support.rel.abs = rel.abs,
              expect.rel.err = expect.rel.err))
}

# given two Gaussian state-space response models, compare benchmarks
compute_comparisons_states = function(ground.resp, model.resp) {

  nStates = length(ground.resp) # should be equal to length(model.resp)
  resp.mu.g = sapply(ground.resp, function(l) l["mu"])
  resp.mu.m = sapply(model.resp, function(l) l["mu"])
  resp.sd.g = sapply(ground.resp, function(l) l["sd"])
  resp.sd.m = sapply(model.resp, function(l) l["sd"])
  resp.mu.d = range(resp.mu.g - resp.mu.m)
  resp.sd.d = max(sqrt(resp.sd.g^2 + resp.sd.m^2))
  eta     = seq(resp.mu.d[1] - 3*resp.sd.d, resp.mu.d[2] + 3*resp.sd.d, length.out = 100)
  rel     = seq(-1, 1, length.out = 100)
  rel.abs = seq(0, 1, length.out = 100)
  metrics = lapply(1:nStates, function(i) {
    tmp = compare_gaussian_models(ground.resp[[i]], model.resp[[i]], 
                                  support.mag = eta, 
                                  support.rel = rel, 
                                  support.rel.abs = rel.abs)
  return(list(expect.rel = tmp$expect.rel.err, 
              mag = tmp$P.magnitude, 
              rel = tmp$P.relative,
              rel.abs = tmp$P.rel.abs))
  })
  metr.mag = lapply(metrics, function(l) l$mag); metr.mag = do.call('rbind', metr.mag)
  metr.rel = lapply(metrics, function(l) l$rel); metr.rel = do.call('rbind', metr.rel)
  metr.rab = lapply(metrics, function(l) l$rel.abs); metr.rab = do.call('rbind', metr.rab)
  metr.xpr = lapply(metrics, function(l) l$expect.rel); metr.rxpr = do.call('rbind', metr.xpr)
  
  return(list(P.magnitude = metr.mag, support.mag = eta, 
              P.relative = metr.rel, support.rel = rel,
              P.rel.abs  = metr.rab, support.rab = rel.abs,
              expect.rel = metr.xpr))
}



# # approximate a response curve along a covariate through linear models
binned_response_curve = function(df, bin.out = NULL,
                                 dep.var = 'AC', bin.var = 'TemperatureF', 
                                 indep.var = 'TemperatureD',
                                 nbins = 10){
  
  # temperature bins
  q = quantile(df[,bin.var], probs = seq(0,1,length.out = nbins))
  
  # compute model in each bin
  fmla = as.formula(paste(dep.var, indep.var, sep='~'))
  fit  = lapply(1:(nbins-1), function(b) {
    idx = which(df[,bin.var] >= q[b] & df[,bin.var] <= q[b+1])
    data= df[idx,]
    fit = lm(fmla, data = data)
    sm  = summary(fit)  
    dat = sm$coefficients[indep.var,c('Estimate', 'Std. Error')]; 
    names(dat) = c("mu", "sd")
    return(dat)
  })
    
  # return just intervals and models if no curve is requested
  if (is.null(bin.out)) return(list(intervals = q, models = fit))
  
  # compute curve for output response variable bin.out if requested  
  idx.bin = findInterval(bin.out, q, rightmost.closed = FALSE, all.inside = TRUE)
  
  mu.vec  = sapply(idx.bin, function(b) fit[[b]]['mu']); names(mu.vec) = NULL
  sd.vec  = sapply(idx.bin, function(b) fit[[b]]['sd']); names(sd.vec) = NULL  
  return(list(var = bin.out, mu = mu.vec, sd = sd.vec))
  
}

# function to compute probability distribution of states given model parameters
compute_probability_profile = function(var, tran){
  P = lapply(1:length(var), function(t) {
        prob = lapply(tran, FUN = function(d){
          m = as.matrix(d)
          y = exp(m %*% c(1,var[t]))
          p = y / sum(y)
          return(as.numeric(p))
        })
        prob = do.call('rbind', prob)
        v    = t(eigen(t(prob))$vectors)
        pv   = abs(v)[,1]
        pv = pv / sum(pv)
        return(pv)
  })
  P = do.call('rbind', P)
  return(P)
}

# function to compute probability distribution of states for given temperature ranges
compute_probability_profile_empirical = function(df, 
                                       bin.out = NULL,
                                       bin.var = 'TemperatureF', 
                                       nbins = 10){
  # temperature bins
  q = quantile(df[,bin.var], probs = seq(0,1,length.out = nbins))
  
  # compute model in each bin
  df$state = as.factor(df$state)
  prob  = lapply(1:(nbins-1), function(b) {
    idx = which(df[,bin.var] >= q[b] & df[,bin.var] <= q[b+1])
    st  = table(df[idx,'state'])
    ps  = st / sum(st)
    return(ps)
  })
  
  # return just intervals and models if no curve is requested
  if (is.null(bin.out)) return(list(intervals = q, models = prob))
  
  # compute curve for output response variable bin.out if requested  
  idx.bin = findInterval(bin.out, q, rightmost.closed = FALSE, all.inside = TRUE)  
  probvec = lapply(idx.bin, function(b) prob[[b]])
  probvec = do.call('rbind', probvec)
  return(list(var = bin.out, probvec = probvec, prob.state = prob))
  
}


# function to compute response curve based on state probability profile
model_response_curve = function(var, resp, P) {
  resp.df = do.call('rbind', resp)
  mu_eff = P %*% resp.df[,"mu"]
  sd_eff = sqrt(P %*% (resp.df[,"mu"]^2 + resp.df[,"sd"]^2) - mu_eff^2)
  return(list(var = var, mu = mu_eff, sd = sd_eff))
}

# ____________________________________________________________________________
# Remove apparent bias by taking difference between summer and winter models
# this is very hackish and experimental

# compute "de-biased" response
remove_response_bias = function(mod1, mod2) {
  
  nStates1 = length(mod1$baseload)
  nStates2 = length(mod2$baseload)
  
  # compute expected distance metric
  E_b_s = matrix(0, nrow = nStates1, ncol = nStates2)
  for (i in 1:nStates1) {
    for (j in 1:nStates2) {
      d_mu = (mod1$baseload[[i]]['mu'] + mod1$volatility[[i]]['mu']) - 
        (mod2$baseload[[j]]['mu'] + mod2$volatility[[j]]['sd'])
      d_sd = sqrt(mod1$baseload[[i]]['sd']^2 + mod2$baseload[[j]]['sd']^2 + 
                  mod1$volatility[[i]]['sd']^2 + mod2$volatility[[j]]['sd']^2)
        
      E_b_s[i,j] = E_fnorm(d_mu, d_sd)
    }
  }  
  # match states in the two models according to expected difference
  idx.b = apply(E_b_s, 2, which.min)
  
  # compute and adjust for bias
  bias_state    = mod2$response[idx.b]
  resp1_unb_st  = diff_gaussians(mod1$response, bias_state)
  
  # compute intersection interval for temperature profiles of two models
  r1 = range(mod1$data$TemperatureF); r2 = range(mod2$data$TemperatureF)
  temp_range = trunc(c(max(min(r1), min(r2)), min(max(r1), max(r2))))
  temp_grid  = temp_range[1]:temp_range[2]
  
  # match states in the two models according to relative ranking at different temperatures
#  P1    = compute_probability_profile(temp_grid, mod1$tran)
  P1    = compute_probability_profile_empirical(mod1$data, 
                                              bin.out = temp_grid,
                                              bin.var = 'TemperatureF', 
                                              nbins = 5)
  P1 = P1$probvec
  rank1 = t(apply(P1, 1, order))
#   P2    = compute_probability_profile(temp_grid, mod2$tran)
  P2    = compute_probability_profile_empirical(mod2$data, 
                                              bin.out = temp_grid,
                                              bin.var = 'TemperatureF', 
                                              nbins = 5)
  P2 = P2$probvec
  rank2 = t(apply(P2, 1, order))
  resp1_temp = lapply(1:length(temp_grid), function(t) {
    resp_st_T = diff_gaussians(mod1$response[rank1[t,]], mod2$response[rank2[t,]])
    resp_st_T = resp_st_T[order(rank1[t,])]
    resp_st_T = do.call('rbind', resp_st_T)
    mu_eff = sum(P1[t,] * resp_st_T[,"mu"])
    sd_eff = sqrt(sum(P1[t,] * (resp_st_T[,"mu"]^2 + resp_st_T[,"sd"]^2)) - mu_eff^2)
    c(mu = mu_eff, sd = sd_eff)
  })
  
  a = list(var = temp_grid,
           mu  = sapply(resp1_temp, function(l) l['mu']),
           sd  = sapply(resp1_temp, function(l) l['sd']))
  
  return(list(per_state = resp1_unb_st, per_temperature = a))
}



