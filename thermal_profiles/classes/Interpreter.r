# Interpreter.r
#
# Interpret output from HMM decoder.
# 
# Adrian Albert
# Last modified: October 2013.
# -----------------------------------------------------------------------

library('methods')
library('timeDate')
library('zoo')
require('lmtest')
library('segmented')
library('dummies')
library('depmixS4')
library('R.utils')
library('MASS')
library('lubridate')

source('code/utils/sequence_entropy.r')
source('code/clustering/compareMCs.r')
source('code/viterbi_states.R')
source('code/utils/plot_utils.r')
source('code/utils/acf_ggplot.r')

# clean-up previous definitions of methods for class Interpreter
removeClass('Interpreter')

options(error = recover)

# ________________________
# Class definition
setClass(
  Class = "Interpreter",
  representation = representation(
    UID        = "character",         # unique person ID
    HMM.STATS  = "list"
    )
)

# _______________________________________
# Constructor method for class Interpreter

setMethod(f = "initialize", 
          signature = "Interpreter",
          definition = function(.Object, raw_data, UID, ZIP5, verbose = T) {
            
            .Object@UID     = as.character(UID)
            
            return(.Object)
          })

compute_characteristic_duration = function(P, method = 'decoding'){
  K    = ncol(P)
  P    = as.matrix(P)
  P0   = P
  if (method == 'decoding') {
    done = F
    tau  = rep(1,K)
    while (!done) {
      idx = which(diag(P) > 0.5)
      done = length(idx) == 0
      if (!done) {
        tau[idx] = tau[idx] + 1 
        P = P %*% P0
      }
    }
  } else tau = round(1 / (1 - diag(P)))
  return(tau)
}
            
# _____________________________________________________
# Compute interpretable metrics on the HMM fit

setGeneric(
  name = "computeMetricsHMM",
  def = function(.Object){standardGeneric("computeMetricsHMM")}
)
setMethod('computeMetricsHMM',
          signature  = 'Interpreter',
          definition = function(.Object){
            
            # access model parameters
            fmopt  = .Object@HMM$model
            K_opt  = .Object@HMM$nStates
            stdev  = .Object@HMM$response$stdev
            respm  = .Object@HMM$response$means
            trans  = .Object@HMM$model@transition
            covar  = rownames(respm)[2]
            breakpoints = .Object@breakpoint$psi
              
            # compute transition matrix at baseline (no temperature) values
            P = lapply(1:K_opt, function(j) {
              tmp = fmopt@transition[[j]]@parameters$coefficients[1,]
              tmp = exp(tmp) / sum(exp(tmp))
              idx.0 = which(tmp<1e-4)
              if (length(idx.0)>0) tmp[idx.0] = 0
              tmp = tmp / sum(tmp)
              if (class(tmp) != 'numeric') return(data.frame(tmp)) else return(data.frame(t(tmp)))
            })            
            # format transition probabilities in matrix form
            P      = do.call(rbind, lapply(P, data.frame, stringsAsFactors=FALSE))
            .Object@HMM[['state.duration']] = compute_characteristic_duration(P, method = 'geometric')
            .Object@HMM[['P.trans']]        = P
            
            # compute "stationary" probabilities given dependent variable (temperature)
            # x_var  = .Object@HMM$model@transition[[1]]@x
            x_var  = cbind(rep(1,121), seq(0, 120, by=1))
            colnames(x_var) = colnames(.Object@HMM$model@transition[[1]]@x)
            dep    = lapply(.Object@HMM$model@transition, function(l) {
              y = x_var %*% l@parameters$coefficients 
              y = exp(y) / rowSums(exp(y))
              return(y)
            })
            dep = do.call('cbind', dep)
            dis = lapply(1:nrow(dep), function(i) {
              M = matrix(dep[i,], ncol = sqrt(ncol(dep)), byrow=T)
              p = abs(as.double(eigen(t(M))$vectors[,1]))
              p = p / sum(p)
              return(p)              
            })
            dis = do.call('rbind', dis)
            distr.covar = data.frame(x_var[,2], dis)
            names(distr.covar)[1] = covar
                      
            # compute characteristic time as function of temperature
            x_var  = cbind(rep(1,121), seq(0, 120, by=1))
            colnames(x_var) = colnames(.Object@HMM$model@transition[[1]]@x)
            tau = lapply(1:nrow(x_var), function(i) {
              t = c()
              for (j in 1:.Object@HMM$nStates) {
                l = .Object@HMM$model@transition[[j]]
                y = x_var[i,] %*% l@parameters$coefficients
                y = exp(y) / sum(exp(y))
                t[j] = 1 / (1 - y[j])
              }
              return(t)
            })     
            tau = do.call('rbind', tau)
                        
            # compute state distribution per time-of-day, day-of-week, summer/winter
            df = data.frame(DateTime = as.character(.Object@data.train$timestamps),
                            state = .Object@HMM$states[,1])
            df$DoW  = wday(df$DateTime, label = T)
            df$ToD  = hour(df$DateTime)
            df$Month= month(df$DateTime)
            df$Season= 'Summer'            
            df$Season[which(df$Month %in% c(11,12,1,2,3,4))] = 'Winter'
            df$Weekday = isHoliday(timeDate(as.character(.Object@data.train$timestamps)))
            df$Weekday = factor(df$Weekday, labels = c('Weekday','Weekend/Holiday'))
            tab.dow = table(df$DoW, df$state)
            tab.dow = tab.dow / rowSums(tab.dow)
            tab.tod = table(df$ToD, df$state)
            tab.tod = tab.tod / rowSums(tab.tod)
            tab.wkd = table(df$Weekday, df$state)
            tab.wkd = tab.wkd / rowSums(tab.wkd)
            df.win  = subset(df, Season =='Winter')
            tab.win = table(df.win$ToD, df.win$state)
            tab.win = tab.win / rowSums(tab.win)
            df.sum  = subset(df, Season =='Summer')
            tab.sum = table(df.sum$ToD, df.sum$state)
            tab.sum = tab.sum / rowSums(tab.sum)
            breakdown = list(DoW = tab.dow, ToD = tab.tod, Weekday = tab.wkd, 
                             Winter = tab.win, Summer = tab.sum)
            
            # store metrics
            .Object@HMM$statMetrics             = list()
            .Object@HMM$statMetrics$distr.covar = distr.covar
            .Object@HMM$statMetrics$breakdown   = breakdown
            .Object@HMM$statMetrics$time.temp   = tau
            return(.Object)
          }
)
            

# _____________________________________________________
# Computes covariate contribution to predicted HMM fit

setGeneric(
  name = "computeContributionsHMM",
  def = function(.Object, verbose = T, covar = 'TemperatureF')
    {standardGeneric("computeContributionsHMM")}
)
setMethod('computeContributionsHMM',
          signature  = 'Interpreter',
          definition = function(.Object, verbose=T, covar = 'TemperatureF') {
	    if (verbose) cat(paste('***** Computing HMM components for Interpreter', .Object@UID, '*****\n'))
      
      # first compute estimates due to simple breakpoint model            
	    breaks     = .Object@breakpoint$psi
	    x_var      = .Object@data.train[,covar]
	    respmean   = .Object@breakpoint$slope
	    hard.resp  = sapply(1:nrow(.Object@data.train), function(j){	      
	      brks   = c(-Inf, breaks, Inf)
	      i      = findInterval(x_var[j], brks)
	      r      = respmean[i]
	      if (!is.finite(brks[i])) i = i + 1
	      if (r>0) t = brks[i]-1 else t = brks[i+1]-1
        if(!is.finite(t)) t = brks[i]-1
	      resp   = abs((t - x_var[j]) * r)        
	      return(resp)
	    })   

      # compute estimated temperature response based off stationary probability model
	    state.prob = .Object@HMM$statMetrics$distr.covar[,-1]
	    x_var      = .Object@HMM$statMetrics$distr.covar[,1]
	    resp       = .Object@HMM$response$means
	    Tcur       = .Object@data.train[,covar]
	    dt         = c(0, diff(Tcur))
	    soft.resp  = sapply(1:ncol(resp), function(s){
        mod = .Object@HMM$model@response[[s]][[1]]
        y   = predict(mod) - resp['(Intercept)',s]
	      idx = sapply(Tcur, function(t) which.min(abs(t-x_var)))
        y   =  y * state.prob[idx,s]
	    })   
      soft.resp = rowSums(soft.resp)
	    	    
	    # do not account for low temperature contributions
      adj_means = .Object@HMM$response$means
# 	    type = .Object@HMM$statMetrics$StateType      
#       idx.none = grep('none', type)
 	    #if (length(idx.none) > 0) adj_means[covar,idx.none] = 0
      
      # format covariates
	    v = rownames(.Object@HMM$response$means)
      s = .Object@HMM$states[,-1]	    
	    B = as.data.frame(t(adj_means[,.Object@HMM$states[,1]]))
	    B.avg = as.data.frame(as.matrix(s) %*% as.matrix(t(adj_means)))	 
      data.dum = dummy.data.frame(.Object@data.train[,-which(names(.Object@data.train)=='timestamps')])
      data.dum[,'(Intercept)'] = 1
      if (length(v)>1) X = data.dum[,v] * B.avg else {
        X = data.dum[,v] * B.avg
        X = data.frame(X)
        names(X) = v
	    }

      # aggregate contribution by variable
      vars = setdiff(c(.Object@resp.vars, .Object@addl.vars), '(Intercept)')
      vars.levels   = lapply(vars, function(v) levels(.Object@data.train[,v])[-1])
      comp = lapply(1:length(vars.levels), function(l) {
        cur.vars = paste(vars[l], vars.levels[[l]], sep='')
        if (length(cur.vars)>1) res = rowSums(X[,cur.vars]) else res = X[,cur.vars]
        return(res)
      })
      comp = as.data.frame(do.call('cbind', comp))
      names(comp) = vars
    	    
	    comp$state= .Object@HMM$states[,1]
      if ('(Intercept)' %in% .Object@resp.vars) comp$Intercept = X[,'(Intercept)']
      names(comp)[which(names(comp) == 'Intercept')] = 'Activity'
      
      # compute seasonal effects
	    df.season = data.frame(date = .Object@data.train$timestamps,
                             Soft.Estimate = soft.resp,
                             Activity      = comp$Activity,
                             Hard.Estimate = hard.resp)
      months    = month(df.season$date)
      df.season$Season = 'Summer'
	    df.season$Day    = wday(df.season$date, label=T)
	    df.season$Hour   = hour(df.season$date)
	    df.season$Season[months %in% c(1:4,10:12)] = 'Winter'
      agg              = aggregate(data = df.season, 
                                   cbind(Soft.Estimate, Hard.Estimate, Activity)~Season + Day + Hour,
                                   FUN = mean)
# 	    agg.summer       = aggregate(data = subset(df.season, Season == 'Summer'), 
# 	                                 cbind(Soft.Estimate, Hard.Estimate, Activity)~Day + Hour,
# 	                                 FUN = mean)
# 	    agg.winter       = aggregate(data = subset(df.season, Season == 'Winter'), 
# 	                                 cbind(Soft.Estimate, Hard.Estimate, Activity)~Day + Hour,
# 	                                 FUN = mean)
#       agg.diff         = agg.summer[,-c(1,2)] - agg.winter[,-c(1,2)]
#       agg.diff         = cbind(agg.summer[,c(1,2)], agg.diff)
	    
      # store for later	    
	    .Object@HMM$components           = list()
	    .Object@HMM$components$stat      = aggregate(data = comp, . ~ state, FUN = mean)
	    .Object@HMM$components$soft.resp = soft.resp
	    .Object@HMM$components$hard.resp = hard.resp
	    .Object@HMM$components$agg.tod   = agg
	    .Object@HMM$components$ts        = comp	    
	    .Object@HMM$components$ts$fit    = .Object@HMM$fit
	    .Object@HMM$components$ts$kWh    = .Object@data.train$kWh	
		
	    return(.Object)
})

