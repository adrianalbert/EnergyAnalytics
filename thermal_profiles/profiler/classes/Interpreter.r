# Interpreter.r
#
# Interpret output from HMM decoder.
# 
# Adrian Albert
# Last modified: October 2013.
# -----------------------------------------------------------------------

library('timeDate')
library('dummies')
library('R.utils')
library('lubridate')

# clean-up previous definitions of methods for class Interpreter
removeClass('Interpreter')

options(error = recover)

# ________________________
# Class definition

setClass(
  Class = "Interpreter",
  representation  = representation(
    UID           = "character",         # unique person ID
    benchmarks    = "list",
    temporalStats = "list",
    contributions = "list"
    )
)

# _____________________________________________________
# Compute interpretable metrics on the HMM fit

computeBenchmarks = function(decoder){
            
  # access model parameters
  K_opt  = decoder@HMM$nStates
  stdev  = decoder@HMM$response$stdev
  respm  = decoder@HMM$response$means
  trans  = decoder@HMM$transition
  covar  = rownames(respm)[2]
  
  # compute "stationary" probabilities given dependent variable (temperature)
  # x_var  = decoder@HMM$model@transition[[1]]@x
  x_var  = cbind(rep(1,121), seq(0, 120, by=1))
  idx.rm = c(length(colnames(trans))-1, length(colnames(trans)))
  colnames(x_var) = colnames(trans)[-idx.rm]
  dep    = lapply(1:K_opt, function(s) {
    c = as.matrix(subset(trans, From == s)[,-idx.rm])
    y = x_var %*% t(c)
    y = exp(y) / rowSums(exp(y))
    return(y)
  })
  prob= dep

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
  colnames(x_var) = colnames(trans)[-idx.rm]
  tau = lapply(1:nrow(x_var), function(i) {
    t = c()
    for (j in 1:K_opt) {
      c = as.matrix(subset(trans, From == j)[,-idx.rm])
      y = x_var[i,] %*% t(c)
      y = exp(y) / sum(exp(y))
      t[j] = 1 / (1 - y[j])
    }
    return(t)
  })     
  tau = do.call('rbind', tau)                      
              
  return(list(probProfile = prob, steadyDistr = dis, duration = tau))
}            

# _____________________________________
# Temporal statistics (time of day etc)

computeTemporalStats = function(decoder){
  # compute state distribution per time-of-day, day-of-week, summer/winter
  df = data.frame(DateTime = as.character(decoder@data.train$timestamps),
                  state = decoder@HMM$states[,1])
  df$DoW  = wday(df$DateTime, label = T)
  df$ToD  = hour(df$DateTime)
  df$Month= month(df$DateTime)
  df$Season= 'Summer'            
  df$Season[which(df$Month %in% c(11,12,1,2,3,4))] = 'Winter'
  df$Weekday = isHoliday(timeDate(as.character(decoder@data.train$timestamps)))
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
  
  return(breakdown)
}

# _____________________________________________________
# Computes covariate contribution to predicted HMM fit

computeContributions = function(decoder, verbose=T) {
  	          
      # format covariates
	    v = rownames(decoder@HMM$response$means)
	    adj_means = decoder@HMM$response$means
      s = decoder@HMM$states[,-1]	    
	    B = as.data.frame(t(adj_means[,decoder@HMM$states[,1]]))
	    # B.avg = as.data.frame(as.matrix(s) %*% as.matrix(t(adj_means)))	 
      data.dum = dummy.data.frame(decoder@data.train[,-which(names(decoder@data.train)=='timestamps')])
      data.dum[,'(Intercept)'] = 1
      if (length(v)>1) X = data.dum[,v] * B else {
        X = data.dum[,v] * B
        X = data.frame(X)
        names(X) = v
	    }
    	    
      # compute seasonal effects
	    df.season = as.data.frame(X)
      if ('(Intercept)' %in% names(df.season)) names(df.season) = c('Baseload', names(df.season)[-1])
	    fmla             = as.formula(paste('cbind(', paste(names(df.season), collapse = ','), ') ~ Season + Day + Hour'))
	    fmla.seas        = as.formula(paste('cbind(', paste(names(df.season), collapse = ','), ') ~ Day + Hour'))
	    
      date             = decoder@data.train$timestamps
      months           = month(date)
      df.season$Season = 'Summer'
	    df.season$Day    = wday(date, label=T)
	    df.season$Hour   = hour(date)
	    df.season$Season[months %in% c(1:4,10:12)] = 'Winter'
      
      agg              = aggregate(data = df.season, fmla, FUN = mean)
	    agg.summer       = aggregate(data = subset(df.season, Season == 'Summer'), fmla.seas, FUN = mean)
	    agg.winter       = aggregate(data = subset(df.season, Season == 'Winter'), fmla.seas, FUN = mean)
      agg.diff         = agg.summer[,-c(1,2)] - agg.winter[,-c(1,2)]
      agg.diff         = cbind(agg.summer[,c(1,2)], agg.diff)	    
		
	    return(list(ts = X, agg = agg))
}


# _______________________________________
# Constructor method for class Interpreter

setMethod(f = "initialize", 
          signature = "Interpreter",
          definition = function(.Object, decoder, verbose = T) {
            
            if (verbose) {
              cat(paste('*** Initializing Interpreter (', .Object@UID, ') ***\n', sep=''))
              t0 = proc.time()
            }
            
            .Object@UID     = as.character(decoder@UID)
            
            # compute stationary metrics
            .Object@benchmarks     = computeBenchmarks(decoder)
            
            if (verbose){
              dt = proc.time() - t0
              print(dt)
            }            

            # compute temporal breakdown of states
            .Object@temporalStats  = computeTemporalStats(decoder)
            
            if (verbose){
              dt = proc.time() - t0
              print(dt)
            }            
            
            # compute temporal breakdown of states
            .Object@contributions  = computeContributions(decoder)
            
            if (verbose){
              dt = proc.time() - t0
              print(dt)
            }            
            
            return(.Object)
          })

# ______________________________________
# Method to dump Interpreter model data

setGeneric(
  name = "dumpInterpretedData",
  def = function(.Object, verbose = T, path = NULL){standardGeneric("dumpInterpretedData")}
)
setMethod('dumpInterpretedData',
          signature  = 'Interpreter',
          definition = function(.Object, verbose = T, path = NULL){
            if (verbose) 
              cat(paste('*** Dumping interpreted data (', .Object@UID, ')***\n', sep=''))
            
            data               = list()
            data$UID           = .Object@UID
            data$benchmarks    = .Object@benchmarks
            data$temporalStats = .Object@temporalStats
            data$contributions = .Object@contributions
            
            if (is.null(path)) return(data) else {
              save(list = c('data'), file = paste(path, .Object@UID, '.RData', sep=''))
              return(NULL)
            }
          })


# _______________________________________
# Extract features for classification

setGeneric(
  name = "extractUserFeatures",
  def = function(.Object, verbose = T)
  {standardGeneric("extractUserFeatures")}
)
setMethod('extractUserFeatures',
          signature  = 'Interpreter',
          definition = function(.Object, verbose=T) {
            
            # TODO: extract & format features for classifying ppl
            
          })