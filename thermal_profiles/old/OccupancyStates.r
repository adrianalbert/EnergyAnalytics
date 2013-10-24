# OccupancyStates.r
#
# Extracts occupancy states from kWh time series data.
# - fitting OLS/FGLS models
# - OLS residual analysis
# - ACF analysis
# - HMM analysis w/ covariates
# - plotting: OLS, structure, HMM structure
# - formatted output: print, data frame for further analysis
# 
# Adrian Albert
# Last modified: May 2013.
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

# clean-up previous definitions of methods for class OccupancyStates
removeClass('OccupancyStates')

options(error = recover)

# ________________________
# Class definition
setClass(
  Class = "OccupancyStates",
  representation = representation(
    UID        = "character",         # unique person ID
    ZIPCODE    = "character",         # current premise zipcode
    timestamps = "character",         # vector of timestamps
    weather    = "data.frame",        # weather covariates
    kWh        = "numeric",           # hourly kWh 
    resp.vars  = 'character',         # covariates for analysis (response)
    tran.vars  = 'character',         # covariates for analysis (transition)
    addl.vars  = 'character',         # covariates for analysis (transition)
    OLS        = "list",              # summary of OLS model
    HMM        = "list",              # summary of HMM model
    data.train = 'data.frame',
    data.test  = 'data.frame',
    breakpoint = 'list'    
    )
)

# _______________________________________
# Constructor method for class OccupancyStates

setMethod(f = "initialize", 
          signature = "OccupancyStates",
          definition = function(.Object, raw_data, UID, ZIP5, verbose = T) {
            
            .Object@UID     = as.character(UID)
            .Object@ZIPCODE = as.character(ZIP5)
            
            if (verbose) 
              cat(paste('*** Initializing OccupancyStates (', .Object@UID, ',', .Object@ZIPCODE,') ***\n', sep=''))
                        
            # ______________________________________________
            # Format kWh data to timeseries format

            cons_names  = paste('hkw', 1:24, sep='')
            if (length(intersect(cons_names,names(raw_data)))>0) {
              # columns referring to kWh
              days        = raw_data$date              
              hours       = paste(rep(formatC(0:23, flag=0, width=2), length(days)), ':00:00', sep='')
              days_times  = rep(days, each = 24)
              days_times  = paste(days_times, hours, sep=' ')
              kWh = as.vector(t(raw_data[,cons_names]))
            } else {
              kWh = raw_data$kWh
              days_times = raw_data$date
            }
            
            # remove possible duplicates
            idx_dup = which(duplicated(days_times))
            if (length(idx_dup)>0) {
              raw_data = raw_data[-idx_dup,]
              days_times=days_times[-idx_dup]
            }                        
            .Object@timestamps  = days_times
            
      	    # is everthing else NA?
      	    if (length(na.omit(kWh)) == 0) {
      	      stop('Error: kWh data is all NAs!')
      	    }

            # if data is all zeros
            rc = range(na.omit(kWh))
            if (rc[2] - rc[1] == 0) {
              stop('Error: data supplied to constructor has no variation!')
            }

            .Object@kWh = kWh
                        
            # initialize default weather object
            .Object@weather     = data.frame()
            .Object@OLS         = list()
            .Object@HMM         = list()
            
            return(.Object)
          })


# ____________________________________________________________
# Method to add in exogenous covariates (weather, billing)

setGeneric(
  name = "addWeather",
  def = function(.Object, wthr_times, wthr_data, verbose = T){standardGeneric("addWeather")}
)
setMethod('addWeather',
          signature  = 'OccupancyStates',
          definition = function(.Object, wthr_times, wthr_data, verbose=T){
            if (verbose) 
              cat(paste('*** Adding weather data for OccupancyStates', .Object@UID,' ***\n',sep=''))
            
            wthr_names = names(wthr_data)
            
            # remove duplicate timestamps in weather data, if any
            idx_dup     = which(duplicated(wthr_times))
            if (length(idx_dup) > 0) {
              wthr_times= wthr_times[-idx_dup]
              wthr_data = wthr_data[-idx_dup,]
            }
            
            # remove NAs, if any
            idx.na = which(!complete.cases(wthr_data)) 
            if (length(idx.na)>0){
              wthr_data = wthr_data[-idx.na,]
              wthr_times= wthr_times[-idx.na]
            }
            
            # is there weather data at times that match kWh data?            
            idx_ok      = which(wthr_times %in% .Object@timestamps)            
            if (length(idx_ok) > 0) {
              wthr_data = wthr_data[idx_ok,]
              wthr_times= wthr_times[idx_ok]
            } else {
              stop(paste('!!! No weather data available for user', .Object@UID))
            }
            idx_ok_cons  = which(.Object@timestamps %in% wthr_times)            
                        
            if (class(wthr_data) == 'numeric') {
              wthr_data = data.frame(wthr_data)
              names(wthr_data) = wthr_names
            }
            # update object
            .Object@weather     = wthr_data
            .Object@timestamps  = as.character(wthr_times)
            .Object@kWh = .Object@kWh[idx_ok_cons]

            return(.Object)
          }
)

# ____________________________________________________________
# Print method for class OccupancyStates: display useful stats.

setMethod('show',
          signature  = 'OccupancyStates',
          definition = function(object){
            cat('*** OccupancyStates Object ***\n')
            
            # basic info
            cat(sprintf('UID:         %s\n', object@UID))            
            cat(sprintf('Zipcode:     %s\n', object@ZIPCODE))
            cat(sprintf('Time range:  %s - %s\n', object@timestamps[1], object@timestamps[length(object@timestamps)]))   
            cat(sprintf('Data amount: %s observations\n', length(object@timestamps)))   
            
            # weather info
            cat(paste('Weather vars:  ', paste(names(object@weather), collapse=', '), '\n'))            
            
            # info on model fit 
            if (length(object@OLS) > 0) {
              cat(sprintf('OLS R2 =\t%f; GR2 =\t%f\n', object@OLS$R2, object@OLS$GR2))
            }        
            if (length(object@HMM) > 0) {
      	      cat(sprintf('HMM States:\t%d\n', object@HMM$nStates))
      	      cat(sprintf('HMM MAPE =\t%f; R2 =\t%f; GR2 =\t%f\n', object@HMM$MAPE, object@HMM$R2, object@HMM$GR2))
      	      cat('HMM contributions:\n')
      	      if (!is.null(object@HMM$components$stat))
                print(as.data.frame(round(as.matrix(object@HMM$components$stat), digits=3)))
            }
            
            cat('*** END: OccupancyStates Object ***\n')
          })

# ____________________________________________________
# Auxiliary function for preparing data for modelling

setGeneric(
  name = "prepareData",
  def = function(.Object, train.frac = 1, resp.vars = c(), tran.vars = c(), addl.vars = c()){standardGeneric("prepareData")}
)
setMethod('prepareData',
          signature  = 'OccupancyStates',
          definition = function(.Object, train.frac = 1, resp.vars = c(), tran.vars = c(), addl.vars = c()) {
            
            # add kWh and weather, if any
            data          = data.frame(kWh = .Object@kWh) 
            N             = nrow(data)
            if (length(.Object@weather)>0)
              data        = cbind(data, .Object@weather)
            
            # if no FE defined, create some for ToD, DoW
            wthr_times     = as.POSIXlt(.Object@timestamps)              
            data$DayOfWeek = as.factor(substr(weekdays(wthr_times), 1, 3))
            data$Month     = as.factor(wthr_times$mon)
            data$HourOfDay = as.factor(as.numeric(sapply( sapply( sapply(strsplit(as.character(wthr_times),' '), '[', 2), strsplit, ':' ), '[',1)))
            data$Weekend   = as.factor(as.numeric(isHoliday( timeDate(wthr_times) )))                          
            
            # add in temperature lags 
            data$TemperatureF.1 = c(data$TemperatureF[1], data$TemperatureF[1:(nrow(data)-1)])
            data$TemperatureF.2 = c(data$TemperatureF[1], data$TemperatureF[2], data$TemperatureF[1:(nrow(data)-2)])
            
            # retain only interesting data
            data = subset(data, select = c('kWh', setdiff(resp.vars, '(Intercept)'), setdiff(tran.vars, '(Intercept)'), addl.vars))
            
            # set response and transition covariates
            if (!is.null(resp.vars)) .Object@resp.vars = resp.vars
            if (!is.null(tran.vars)) .Object@tran.vars = tran.vars
            if (!is.null(addl.vars)) .Object@addl.vars = addl.vars
            
            # add in timestamps back
            data$timestamps = .Object@timestamps

            # split up into test and train
            T_train   = trunc(nrow(data) * train.frac)
            T_train   = trunc(T_train / 24) * 24
            idx_train = 1:T_train
            idx_test  = setdiff(1:nrow(data), idx_train)                                
            .Object@data.train = data[idx_train,]
            if (length(idx_test) > 0) .Object@data.test = data[idx_test,] else .Object@data.test = data.frame()                                      
            
            rm(list=c('data'))
            gc()
            return(.Object)
})

# _________________________________
# Method to perform OLS analysis

setGeneric(
  name = "fitOLS",
  def = function(.Object, verbose = T, stats = T){standardGeneric("fitOLS")}
)
setMethod('fitOLS',
          signature  = 'OccupancyStates',
          definition = function(.Object, verbose=T, stats = T){
            
            if (verbose)
              cat(paste('*** OLS analysis for OccupancyStates-SPID (', .Object@UID, ')***\n', sep=''))
            
            # ______________________________________________________________
            # Set-up OLS problem (regress kWh on weather + tod/dow)
            
            # get data
            data.dum = na.omit(.Object@data.train)
            
            # set-up model covariates
            fmla_null    = as.formula("kWh ~ 0")
            if ('(Intercept)' %in% .Object@resp.vars) fmla_full = 'kWh ~ ' else fmla_full = 'kWh ~ -1 + '
            fmla_full    = as.formula(paste(fmla_full, 
                                            paste(setdiff(.Object@resp.vars, '(Intercept)'),collapse='+')))
            # perform regression            
            OLS.null     = lm( fmla_null, data = data.dum )
            OLS.step     = lm( fmla_full, data = data.dum )
            
            # ______________________________________________________________
            # Investigate structure of residual
                      
            # which coefficients are statistically significant at least at 0.1 level?
            summ.fit = summary(OLS.step)
            idx.sigf = which(summ.fit$coefficients[,4] <= 0.1 & abs(summ.fit$coefficients[,1])>1e-6 )
            coef.sig = summ.fit$coefficients[idx.sigf,]
            
	          if (stats) {
              # test for heteroskedasticity using Breusch-Pagan test
              bp.test        = bptest(OLS.step, data = data.dum)
              chi2.95        = qchisq(0.95,bp.test$parameter+1)
              heterosc.BP.95 = bp.test$statistic > chi2.95
      
              # test for serial correlation using Durbin-Watson test
              dw.test    = dwtest(OLS.step, data = data.dum)
              sercorr.DW = dw.test$alternative == "true autocorrelation is greater than 0"
       	    }
            
            # save model result
            .Object@OLS = list()
            .Object@OLS[['fit.summary']] = summ.fit
            .Object@OLS[['fitted.vals']] = predict(OLS.step)
            if (stats) .Object@OLS[['heterosked']]  = heterosc.BP.95
            if (stats) .Object@OLS[['serial.corr']] = sercorr.DW
            .Object@OLS[['coef.signif']] = coef.sig
            .Object@OLS[['all.covars']]  = .Object@resp.vars
            .Object@OLS[['R2']]          = .Object@OLS$fit.summary$adj.r.squared
            .Object@OLS[['LL0']]         = as.numeric(logLik(OLS.null))
            .Object@OLS[['LLf']]         = as.numeric(logLik(OLS.step))
            .Object@OLS[['GR2']]         = GR2(.Object@OLS[['LL0']], as.numeric(logLik(OLS.step)), nrow(.Object@data.train))
            
            rm(list=c('OLS.step', 'data.dum'))
            return(.Object)
          }
)

# ______________________________________________________
# Compute initial temperature breakpoint estimate

setGeneric(
  name = "computeBreakpoint",
  def = function(.Object, verbose = T, break.var = c(TemperatureF=45, TemperatureF=55)) {standardGeneric("computeBreakpoint")}
)
setMethod('computeBreakpoint',
          signature  = 'OccupancyStates',
          definition = function(.Object, verbose=T, break.var = c(TemperatureF=45, TemperatureF=55)) {
            
            if (verbose)
              cat(paste('*** Breakpoint analysis for OccupancyStates-SPID (', .Object@UID, ')***\n', sep=''))

            # compute daily totals
            df.hour = subset(.Object@data.train, select = c('timestamps', 'kWh', names(break.var)))
            df.day  = aggregate(data = df.hour, cbind(kWh, TemperatureF) ~ as.Date(timestamps), FUN = mean)
            df.day$kWh = df.day$kWh

            # first compute linear model 
            fmla_lm = as.formula(paste('kWh ~ ', names(break.var)))
            fit     = lm(fmla_lm, data = df.hour)
            
            # then compute segmented model
            fmla    = as.formula(paste('~',names(break.var)))
            fit.seg = try(segmented(fit, seg.Z = fmla, psi = break.var,
                                    control = seg.control(n.boot = 100, random=T)))
            if (class(fit.seg)[1] == 'try-error') {
              bp     = c(min(df.day[,names(break.var)[1]]), max(df.day[,names(break.var)[1]]))
              slopes = rep(0,3)
            } else {
              bp     = fit.seg$psi[,2]
              slopes = slope(fit.seg)[[1]][,1]
            }
            names(bp) = paste('T', 1:length(bp), sep='.')            
            
            # create additional data columns for breakpoint: train data
#             var.str = paste(names(break.var),c('C','H'), sep='.')
#             diff_1  = .Object@data.train[,names(break.var)[1]] - fit.seg$psi[1,2]
#             diff_2  = .Object@data.train[,names(break.var)[1]] - fit.seg$psi[2,2]
#             .Object@data.train[,var.str[1]] = sapply(diff_1, function(x) min(x,0))
#             .Object@data.train[,var.str[2]] = sapply(diff_2, function(x) max(x,0))
#             var.str = paste(names(break.var),c('C','H'), sep='.')
#             brk.val = .Object@data.train[,names(break.var)[1]]
#             .Object@data.train[,var.str[1]] = as.numeric(brk.val < fit.seg$psi[1,2])
#             .Object@data.train[,var.str[2]] = as.numeric(brk.val >= fit.seg$psi[2,2])
            
            # create additional data columns for breakpoint: test data
#             diff_1  = .Object@data.test[,names(break.var)[1]] - fit.seg$psi[1,2]
#             diff_2  = .Object@data.test[,names(break.var)[1]] - fit.seg$psi[2,2]
#             .Object@data.test[,var.str[1]] = sapply(diff_1, function(x) as.numeric(x>0))
#             .Object@data.test[,var.str[2]] = sapply(diff_2, function(x) as.numeric(x<=0))
#             var.str = paste(names(break.var),c('C','H'), sep='.')
#             brk.val = .Object@data.test[,names(break.var)[1]]
#             .Object@data.test[,var.str[1]] = as.numeric(brk.val < fit.seg$psi[1,2])
#             .Object@data.test[,var.str[2]] = as.numeric(brk.val >= fit.seg$psi[2,2])
            
            # store stats
            # .Object@tran.vars       = c(setdiff(.Object@tran.vars, names(break.var)), var.str)
            # .Object@resp.vars       = c(setdiff(.Object@tran.vars, names(break.var)), var.str)
            .Object@breakpoint        = list(psi = bp, slope = slopes)
                        
            return(.Object)
          }
)

# ______________________________________
# Define HMM model for depmixS4 package

defineHMM = function(data.train, K = 3, type = 'default', 
                     response.vars = NULL, transitn.vars = NULL, 
                     respstart = NULL, trstart = NULL) {
  
  intercept_resp = '(Intercept)' %in% response.vars
  intercept_tran = '(Intercept)' %in% transitn.vars
  response.vars  = setdiff(response.vars, '(Intercept)')
  transitn.vars  = setdiff(transitn.vars, '(Intercept)')
  
  # define response model
  if (intercept_resp) fmla_response = 'kWh ~ 1' else fmla_response = 'kWh ~ -1'
  if (length(response.vars)>0) fmla_response = paste(fmla_response, paste(response.vars, collapse = '+'), sep='+')
  fmla_response = as.formula(fmla_response)            
  
  # define transition model
  if (intercept_tran) fmla_transitn = 'kWh ~ 1' else fmla_transitn = 'kWh ~ -1'
  if (length(transitn.vars)>0) fmla_transitn = paste(fmla_transitn, paste(transitn.vars, collapse = '+'), sep='+')
  fmla_transitn = as.formula(fmla_transitn)            
  
  # set up model for depmixS4
  if (type == 'default') {
    mod <- depmix(response  = fmla_response, 
                  transition= fmla_transitn,
                  data      = data.train, 
                  nstates   = K, 
                  respstart = respstart,
                  trstart   = trstart,
                  family    = gaussian(),
                  instart   = runif(K))
  } else {
    mod = makeModelSelect(data.train, fmla_response, fmla_transitn, K = K)    
  }
  return(mod)
}

# ______________________________________
# Fit HMM using depmixS4 package.

fitHMM = function(mod, nRestarts = 1, verbose = T, maxit = 100){ 

  if (verbose) cat('---> Learning HMM...\n')
    
  # fit given model            
  ok = FALSE
  it = 0
  cur_tol = 1e-3
  while (!ok & it <= nRestarts) {
    
    out <- capture.output(fm  <- try(fit(mod, verbose = T, useC = T, emcontrol = em.control(maxit = maxit, tol = cur_tol))))                
    nlines = length(out)
    
    if (class(fm) != 'try-error') {
      ok = nlines > 2                                   
    } else ok = FALSE
    if (!ok){
      if (verbose) cat('Bad HMM fit; re-estimating...\n')
      cur_tol = cur_tol / 10
    } else {
      if (verbose) cat(paste('Convergence in', nlines*5, 'iterations.\n'))
    }
    it = it + 1
  }
  return(fm)  
}            

# __________________________________________
# Cross-validation for depmixS4 HMM models.

fitHMM.cv = function(data, K = 3, 
                     response.vars = NULL, transitn.vars = NULL, 
                     respstart = NULL, trstart = NULL, 
                     nRestarts = 1, verbose = T, maxit = 100){ 
  if (verbose)
    cat(paste('---> HMM Cross-Validation K =', K,'\n'))
  
  # fit HMM on 1/2 of data
  data_1  = data[seq(1, nrow(data)-1, by=2),]
  mod     = defineHMM(data_1, K = K, 
                      response.vars = response.vars, 
                      transitn.vars = transitn.vars, 
                      respstart     = respstart,
                      trstart       = trstart)

  fm_half = fitHMM(mod, nRestarts = nRestarts, verbose = verbose, maxit = maxit)
  
  if (class(fm_half) == 'try-error') 
    return(list(model.half = NA, metrics = c(MAPE.cv = NA, R2.cv = NA, BIC = NA, AIC = NA)))
  
  # decode the other half
  data_2  = data[seq(2, nrow(data), by=2),]
  data.dum= dummy.data.frame(data_2[,-which(names(data_2)=='timestamps')])
  probs   = viterbi_states(fm_half, data.dum, data.dum$kWh)
  
  # compute Cross-Validation error metrics
  MAPE    = mean(abs((probs$fit - data_2$kWh) / probs$fit))
  R2      = 1 - sum((probs$fit - data_2$kWh)^2) / sum((data_2$kWh - mean(data_2$kWh))^2)
  
  return(list(model.half = fm_half, metrics = c(MAPE.cv = MAPE, R2.cv = R2, BIC = BIC(fm_half), AIC = AIC(fm_half))))
}            

# _________________________________
# Learn HMM and optimum model size

setGeneric(
  name = "learnOccupancyStates",
  def = function(.Object, verbose = T, Kmin = 3, Kmax = 3, nRestarts = 1, maxit = 100, thresh.R2 = 0.8, thresh.MAPE = 0.15)
    {standardGeneric("learnOccupancyStates")}
)
setMethod('learnOccupancyStates',
          signature  = 'OccupancyStates',
          definition = function(.Object, verbose=T, Kmin = 3, Kmax = 3, 
                                nRestarts = 1, maxit = 100, 
                                thresh.R2 = 0.8, thresh.MAPE = 0.15){
            
            if (verbose)
              cat(paste('*** HMM analysis for OccupancyStates-SPID (', .Object@UID, ')***\n', sep=''))
            
            # _____________________________
            # Set up depmixS4 model object
          
            # prepare dataset
            data.train    = .Object@data.train
            response_vars = .Object@resp.vars
            transitn_vars = .Object@tran.vars
            additinl_vars = .Object@addl.vars
            
            # _______________________________________________
            # choose model size K (number of states)

            if (Kmax > Kmin) {
              result = list()
              k      = Kmin
              done   = F
              while (k <= Kmax & !done) {
                result[[as.character(k)]] = fitHMM.cv(data.train, K = k, 
                                                      response.vars = response_vars, transitn.vars = transitn_vars, 
                                                      respstart = NULL, trstart = NULL, 
                                                      nRestarts = nRestarts, verbose = verbose, maxit = maxit)
                R2.cv   = result[[as.character(k)]]$metrics['R2.cv']
                MAPE.cv = result[[as.character(k)]]$metrics['MAPE.cv']
                if (R2.cv > thresh.R2 | MAPE.cv < thresh.MAPE) done = T else k = k + 1
              }
              metric = t(sapply(result, function(l) l$metrics))
              rownames(metric) = names(result)
              if (done) K_opt = k else {            
                # get minimum BIC from good fits 
                BICs   = sapply(result, function(l) l$metrics[3])
                bad.fit= which(!is.finite(BICs))
                if (length(BICs) == length(bad.fit)) {
                  stop('HMM fitting error!')
                }
                Kvec   = Kmin:Kmax
                idx_opt= which.max(metric[,'R2.cv'])             
                K_opt  = Kvec[idx_opt]
              }
            } else { 
              K_opt = Kmin   
              metric = NULL
              result = NULL
            }
                                  
            # _______________________________________________
            # Fit model to full data
            
            trstart = NULL

            # fit model to full data 
            mod    = defineHMM(data.train, K = K_opt, 
                               response.vars = c(response_vars, additinl_vars), 
                               transitn.vars = transitn_vars, 
                               respstart = NULL, 
                               trstart = trstart)
            fm     = fitHMM(mod, nRestarts = nRestarts, verbose = verbose, maxit = maxit)
            
            # ______________________
            # Save computation
            
            .Object@HMM         = list()
            .Object@HMM$nStates = K_opt
            .Object@HMM$model   = fm
            .Object@HMM$cv.metrics = metric
            
            rm(list = c('result', 'fm'))
            gc()
            return(.Object)
          }
)           
         
# _______________________________________________
# Compute estimates of out-of-sample performance

setGeneric(
  name = "computePredictionAccuracy",
  def = function(.Object, test.periods = 5, verbose = T){standardGeneric("computePredictionAccuracy")}
)
setMethod('computePredictionAccuracy',
          signature  = 'OccupancyStates',
          definition = function(.Object, test.periods = 5, verbose = T){
            
            TT  = nrow(.Object@data.test)
            per = (1:TT) %/% (TT %/% test.periods) + 1
            if (length(which(per == test.periods+1))<=2) per[which(per == test.periods+1)] = test.periods
            test.periods = length(unique(per))
            
            data.dum = dummy.data.frame(.Object@data.test[,-which(names(.Object@data.test)=='timestamps')])
            
            stats = list()
            for (p in 1:test.periods) {

              idx_test   = which(per <= p) 
              T_test     = length(idx_test)
              test_data  = data.dum[idx_test,]                        
              
              if (verbose) cat(paste('Test window:', T_test, '\n'))

              # set up Viterbi decoding on test data
              object     = .Object@HMM$model
              probs      = viterbi_states(object, test_data, test_data$kWh)              
              
              # compute decoding performance statistics
              MAPE       = mean(abs((probs$fit - test_data$kWh) / probs$fit))
              R2         = 1 - sum((probs$fit - test_data$kWh)^2) / sum((test_data$kWh - mean(test_data$kWh))^2)
                            
              # store stats
              stats[[as.character(T_test)]] = c(MAPE, R2)
            }
            len_vec = names(stats)
            stats   = do.call('rbind', stats)
            rownames(stats) = len_vec
            colnames(stats) = c('MAPE', 'R2')
            
            .Object@HMM$PredictStats = stats
            return(.Object)
        }
)

# ________________________________________
# Adjust temperature slopes for the states

setGeneric(
  name = "pruneStates",
  def = function(.Object, verbose = T, thresholds = NULL, covar = 'TemperatureF') {standardGeneric("pruneStates")}
)
setMethod('pruneStates',
          signature  = 'OccupancyStates',
          definition = function(.Object, verbose=T, thresholds = NULL, covar = 'TemperatureF') {
            
            # for each state, re-estimate linear model using robust OLS
            fm_opt = .Object@HMM$model      
            states = .Object@HMM$states[,1]
            resp.l = sapply(unique(states), function(s) {
              fmla = fm_opt@response[[s]][[1]]@formula
              data = .Object@data.train[which(states == s),]
              fitr = rlm(fmla, data)
              summ = summary(fitr)
              t    = summ$coefficients[,3]
              pval = 2*pt(-abs(t),df=nrow(data)-1)
              return(c(coef = coef(fitr), err = summ$coefficients[,2], 
                       pval = pval, sd = sqrt(sum(fitr$residuals^2) / nrow(data))))
            })
            
            # extract adjusted coefficients and errors
            respns   = resp.l[grep('coef', rownames(resp.l)),]
            resper   = resp.l[grep('err', rownames(resp.l)),]
            stdev    = resp.l['sd',]            
            pvals    = resp.l[grep('pval', rownames(resp.l)),]
            rownames(respns) = gsub('coef.','',rownames(respns))
            rownames(resper) = gsub('err.','',rownames(resper))
            rownames(pvals)  = gsub('pval.','',rownames(pvals))
            
            # classify states according to temperature dependence
            stat     = atan(respns[covar,]) * 180 / pi
            meaning  = names(thresholds)
            thresholds = c(-Inf, thresholds, Inf)
            names(thresholds) = c(meaning[1], meaning, meaning[length(meaning)])
            intens   = sapply(stat, function(s) {
              return(names(thresholds)[which(abs(s) <= thresholds)[1]])
            })            
            type   = sapply(sign(stat), function(s) if (s<0) return('heating') else return('cooling'))

            # some states however have insignificant temperature dependence
            idx.c = which(abs(pvals[covar,])>0.15)            
            if (length(idx.c)>0) intens[idx.c] = 'none'
            
            # piece state interpretation together
            idx.0 = which(intens =='none')
            if (length(idx.0)>0) type[idx.0] = ''
            intens = paste(intens, type, sep = '.')
            
            .Object@HMM$response$means = respns
            .Object@HMM$response$error = resper
            .Object@HMM$response$sd    = stdev
            .Object@HMM$statMetrics$StateType   = intens
            
            return(.Object)
            
          }
)

# ________________________________________
# Compute statistics for HMM fit

# function to compute Generalized R2
GR2 = function(LL0, LL, N) {
  gr2 = 1 - exp(2*(LL0 - LL) / N)
}

# function to compute R2
R2 = function(y, yfit) {
  r2 = 1 - sum((y - yfit)^2) / sum((y-mean(y))^2)
}

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
            
setGeneric(
  name = "computeStatsHMM",
  def = function(.Object){standardGeneric("computeStatsHMM")}
)
setMethod('computeStatsHMM',
          signature  = 'OccupancyStates',
          definition = function(.Object){
            fm_opt     = .Object@HMM$model
            K_opt      = .Object@HMM$nStates
            
            # ________________________________
            # Response parameters
            
            resp_vars  = colnames(fm_opt@response[[1]][[1]]@x)
            .Object@HMM[['response']] = list()
            stddev     = sapply(1:K_opt, function(j) fm_opt@response[[j]][[1]]@parameters$sd)
            names(stddev) = 1:K_opt
            .Object@HMM[['response']][['stdev']] = stddev
            params = sapply(1:K_opt, function(j) fm_opt@response[[j]][[1]]@parameters$coefficients)
            if (class(params) == 'numeric') params = data.frame(params) else params = data.frame(t(params))
            rownames(params) = 1:K_opt
            colnames(params) = resp_vars
            .Object@HMM[['response']][['means']]  = as.data.frame(t(params))

            # instead of pruneStates
            .Object@HMM$response$sd    = stddev            
            
            # ______________________________
            # Transition parameters
                        
            # format into matrix for output 
            transitn_vars  = colnames(fm_opt@transition[[1]]@x)          
            params = lapply(1:K_opt, function(j) {
              tmp = fm_opt@transition[[j]]@parameters$coefficients
              rownames(tmp) = transitn_vars              
              if (class(tmp) != 'numeric') return(data.frame(tmp)) else return(data.frame(t(tmp)))
            })                          
            params = do.call(cbind, lapply(params, data.frame, stringsAsFactors=FALSE))            
            if (!is.null(ncol(params))) {
              colnames(params) = as.vector(sapply(1:K_opt, function(x) paste(x, '->',1:K_opt,sep='')))
            } else {
              names(params) = as.vector(sapply(1:K_opt, function(x) paste(x, '->',1:K_opt,sep='')))
            }
            if (!is.null(nrow(params))) rownames(params) = transitn_vars
            .Object@HMM[['transition']] = params
            
            # Viterbi states
            .Object@HMM[['states']] = posterior(fm_opt)
            
            # ________________________________________
            # Perform preliminary analysis of errors
            
            # Predicted state means
            fit = sapply(1:K_opt, function(k) predict(fm_opt@response[[k]][[1]])) 
            .Object@HMM[['fit.max']] = sapply(1:nrow(fit), function(j) fit[j,.Object@HMM[['states']][j,1]])
            .Object@HMM[['fit']]     = sapply(1:nrow(fit), function(j) sum(fit[j,] * .Object@HMM[['states']][j,-1]))
            
            # Standard deviation of most likely state
            hmm.sigma = .Object@HMM$response$stdev[.Object@HMM[['states']][,1]]  
            
            # model fit (penalized likelihood)
            .Object@HMM[['BIC']]    = BIC(fm_opt)
            
            # residuals for each state
            residuals = .Object@data.train$kWh - .Object@HMM$fit                          
            .Object@HMM[['residual']] = residuals
            
            # test normality of residuals in each state
            is_normal = sapply(1:.Object@HMM[['nStates']], function(j) {
              idx = which(.Object@HMM[['states']][,1] == j)
              if (length(na.omit(residuals[idx]))>3 & length(na.omit(residuals[idx]))<5000)
                res = shapiro.test(na.omit(residuals[idx])) else res = list(p.value=NA)
              pval= res$p.value
              return(pval)
            })
            .Object@HMM$sw_test = is_normal
                            
            # in-sample fit performance metrics
            .Object@HMM[['MAPE']]    = mean(abs(residuals / .Object@HMM$fit))
            .Object@HMM[['GR2']]     = GR2(.Object@OLS[['LL0']], as.numeric(logLik(fm_opt)), nrow(.Object@data.train))
            .Object@HMM[['R2']]      = R2(.Object@data.train$kWh, .Object@HMM$fit)
            
            # entropy of train set
            .Object@HMM$entropy      = list()
            entropy_gzip_bytes       = gzip_entropy(.Object@HMM$states[,1], file = paste('tmp',.Object@UID,sep=''))
            .Object@HMM$entropy$gzip = -log2(entropy_gzip_bytes*8)
            .Object@HMM$entropy$uncr = uncorr_entropy(.Object@HMM$states[,1])
            .Object@HMM$entropy$rand = rand_entropy(.Object@HMM$states[,1])
            .Object@HMM$entropy$lziv = lempel_ziv_entropy(c(), .Object@HMM$states[,1])$entropy
            .Object@HMM$predict      = list()      
            N = length(.Object@HMM$states[,1])
            .Object@HMM$predict$gzip = compute_predictability(.Object@HMM$entropy$gzip, N)
            .Object@HMM$predict$uncr = compute_predictability(.Object@HMM$entropy$uncr, N)
            .Object@HMM$predict$rand = compute_predictability(.Object@HMM$entropy$rand, N)
            .Object@HMM$predict$lziv = compute_predictability(.Object@HMM$entropy$lziv, N)
  
            rm(list = c('fm_opt'))
            gc()
            return(.Object)
})

# _____________________________________________________
# Compute interpretable metrics on the HMM fit

setGeneric(
  name = "computeMetricsHMM",
  def = function(.Object){standardGeneric("computeMetricsHMM")}
)
setMethod('computeMetricsHMM',
          signature  = 'OccupancyStates',
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
          signature  = 'OccupancyStates',
          definition = function(.Object, verbose=T, covar = 'TemperatureF') {
	    if (verbose) cat(paste('***** Computing HMM components for OccupancyStates', .Object@UID, '*****\n'))
      
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

# _______________________
# Plot OccupancyStates object

setMethod('plot',
          signature  = 'OccupancyStates',
          definition = function(x, type = 'default', interval = NULL, verbose = T, 
                                covar = NULL, PACF = T, highlight = T, separate = F, ...){
            
            if (verbose)
              cat(paste('*** ', type, ' Plot for OccupancyStates (', x@UID, ') ***\n', sep=''))
            
            timestamps = as.POSIXct(x@data.train$timestamps)
            ols.resids = x@OLS[['fit.summary']]$residuals
            ols.resids = naresid(x@OLS$fit.summary$na.action, ols.resids)
            kWh= x@data.train$kWh

            ols.means  = x@OLS[['fitted.vals']]
            hmm.means  = x@HMM$fit
            states     = x@HMM$states[,1]
            hmm.sigma  = x@HMM$response$stdev[states]
            hmm.residual=x@HMM$residual      
            wthr_data   =x@weather
      	    comps.hmm   =x@HMM$components$ts
            if (!is.null(interval)) {
              idx_ok      = which(x@timestamps >= as.POSIXct(interval[1]) & 
                                  x@timestamps <= as.POSIXct(interval[2]))
              timestamps  = timestamps[idx_ok]
              ols.resids  = ols.resids[idx_ok]
              kWh = kWh[idx_ok]
              ols.means   = ols.means[idx_ok]
              hmm.means   = hmm.means[idx_ok]
              hmm.sigma   = hmm.sigma[idx_ok]
              states      = states[idx_ok]
              hmm.residual= hmm.residual[idx_ok]
              wthr_data   = wthr_data[idx_ok,]
      	      comps.hmm   = comps.hmm[idx_ok,]
            }  
                    
            if (type == 'default') {
              plot(timestamps, kWh, 
                   main = paste('kWh OccupancyStates-SPID (', x@UID,')'),
                   xlab = 'Time', ylab = 'kWh', type = 'l', lwd = 2)
            }
            
            if (type == 'weather') {
              if (is.null(wthr_data) | nrow(wthr_data)==0) {
                cat('No weather data!\n')
              } else {
                wthr_names  = c('TemperatureF', 'DewpointF', 'Pressure', 'WindSpeed', 'Humidity', 
                                'HourlyPrecip', 'SolarRadiation')
                wthr = wthr_data[,names(wthr_data) %in% wthr_names]
                wthr$kWh = kWh
                wthr = zoo(wthr, order.by = timestamps)  
                plot(wthr, 
                     main = paste('kWh OccupancyStates-SPID (', x@UID,')'),
                     xlab = 'Time', type = 'l', lwd = 2)                
              }
            }
                        
            if (type == 'OLS-res') {
              layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))            
              plot(timestamps, ols.resids, 
                   main = 'OLS Residuals', ylab = 'Residuals', xlab = 'Time', lwd = 3, type='l')
              plot(density(na.omit(ols.resids)), 
                   main = 'OLS Residuals: Density', ylab = 'pdf', xlab = 'Residuals', lwd = 3, type='l')
              acf(na.omit(ols.resids), 
                   main = 'OLS Residuals', ylab = 'Residuals', xlab = 'Time', lwd = 3)
            }    
            
            if (type == 'OLS-fit') {
              yrange = range(c(ols.means, kWh))
              plot(timestamps, ols.means, 
                   main = 'OLS Fitted Values', ylab = 'Values', xlab = 'Time', 
                   lwd = 3, type='l', col='red', ylim = yrange)
              points(timestamps, kWh, type='b', pch=20, col = 'black')
              legend('topleft', c('OLS Fit', 'Observations'), col=c('red', 'black'))
            }    
                        
            if (type == 'HMM-ts') {
            		title = paste("Zoom-In: States and Observed Emissions (", x@UID,')',sep='')                       
            		plt = plot_hmm_ts(hmm.means, hmm.sigma, states, timestamps, kWh, 
                        				  y.lab = 'kWh', title = title)
            		  return(plt)
            }
                 
            # plot time series of state parameters
            if (type == 'HMM-coefs-ts') {
              covar_state = as.data.frame(t(as.matrix(x@HMM$response$means[covar,])))
              names(covar_state) = covar
              title = paste("Zoom-In: States Parameters (", x@UID,')',sep='')                       
              plt = plot_hmm_coefs_ts(covar_state, states, kWh, timestamps,  
                                      title = title)
              return(plt)
            }
            
            # plot ToD/Dow breakdown of occupancy states
            if (type == 'HMM-state-breakdown') {
              title = paste(covar, " states breakdown (", x@UID,')',sep='')                       
              plt   = plot_state_breakdown(as.factor(states), timestamps, 
                                           state.type =  NULL, #x@HMM$statMetrics$StateType,
                                           title = title)
              return(plt)
            }
            
            # plot heatmap of occupancy states
            if (type == 'HMM-state-heatmap') {
              if (covar == 'States') myMat = as.factor(states) else 
                if (covar == 'kWh') myMat = -kWh else {
                  myMat = -t(x@HMM$response$means[covar,states])
                  colnames(myMat) = NULL
                }
              title = paste(covar, " heatmap (", x@UID,')',sep='')                       
              plot_state_heatmap(myMat, timestamps, title = title)
            }
            
            # plot heatmap of occupancy states (ggplot)
            if (type == 'HMM-state-heatmap2') {
              if (covar == 'States') myMat = as.factor(states) else 
                if (covar == 'kWh') myMat = -kWh else {
                  myMat = -t(x@HMM$response$means[covar,states])
                  colnames(myMat) = NULL
                }
              title = paste(covar, " heatmap (", x@UID,')',sep='')                       
              plt = plot_state_heatmap2(myMat, timestamps, title = title)
              return(plt)
            }
            
            if (type == 'HMM-MC-cov') {
              title             = paste("State Space Diagram (", x@UID, ')',sep='')
              contrib           = x@HMM$components$stat
              contrib$mu        = as.numeric(x@HMM$response$means[1,])
              contrib$sigma2    = as.numeric(x@HMM$response$stdev)
              plt               = plot_HMM_MC_cov(x@HMM$P.trans, x@HMM$state.duration, contrib, title=title)
              return(plt)
            }
            
            if (type =='HMM-pacf') {
              
                plt = acf_ggplot(na.omit(kWh), na.omit(hmm.means), 
                                 title = 'PACF: Empirical vs Model', PACF = PACF)
                return(plt)
            }
                  
            if (type == 'HMM-res') { 
      		    plt = plot_HMM_res(hmm.residual, states, x@HMM$response$stdev, sw_test = x@HMM$sw_test)
      	  	  plt
             }             
                  
            if (type == 'HMM-MC2') {
              title = paste("HMM MC structure (", x@UID, ')',sep='') 
              # organize MC information
              plt = plot_HMM_MC_NET(as.numeric(x@HMM$response$means[1,]), as.numeric(x@HMM$response$stdev), x@HMM$transition)
      	      plot(plt, main = title) 
            }
      
      	    if (type == 'HMM-contrib-ts'){
          		title = paste("Zoom-In: HMM Covariate Contributions (", x@UID,')',sep='')
              covar_state = as.data.frame(t(as.matrix(x@HMM$response$means[covar,])))
              names(covar_state) = covar              
          		plot_components_ts(comps.hmm, timestamps, 
                                 title = title, states = states, 
                                 covars = covar_state)              
            }
        
      	    if (type == 'HMM-contrib-tot') {
              df = x@HMM$components$stat
              title = paste("HMM Covariate Contributions (", x@UID,')',sep='')
              p     = plot_tornado(df, title = title, nrow = 2)		
          		return(p)
      	    }
            
            if (type == 'HMM-dep-covar') {
              title = paste(covar, paste("dependence (", x@UID,')',sep=''))
              marks = x@breakpoint$psi
              names(marks) = c('C', 'H')
              states = x@HMM$states[,1]
#              states = paste(states, x@HMM$statMetrics$StateType[states],sep=':')
#               states = paste(states, ': slope=', 
#                              round(x@HMM$response$means[covar,states]*100, digits=3),'+/-',
#                              round(x@HMM$response$error[covar,states]*100, digits=3),'(/100)',
#                              sep='')
              states = paste(states, ': slope=', 
                             round(x@HMM$response$means[covar,states]*100, digits=3),'(/100)',sep='')
              p     = plot_dep_covar(x@data.train[,covar], x@data.train$kWh, states, 
                                     title = title, x.lab = covar, y.lab = 'kWh', highlight = highlight,
                                     markers = marks, separate = separate)  	
              return(p)
            }            
          
            if (type == 'HMM-trans-prob') {
              title = paste("Transition probabilities (", x@UID,')',sep='')
              marks = x@breakpoint$psi
              names(marks) = c('C', 'H') 
              
              # compute transition probabilities over a predefined grid
              x_var = x@HMM$model@transition[[1]]@x
              x_var = cbind(rep(1,100), seq(min(x_var[,covar]), max(x_var[,covar]), length.out=100))
              colnames(x_var) = colnames(x@HMM$model@transition[[1]]@x)
              dep = lapply(x@HMM$model@transition, function(l) {
                y = x_var %*% l@parameters$coefficients 
                y = exp(y) / rowSums(exp(y))
                return(y)
              })

              # plot 
              x_var = as.data.frame(x_var[,covar])
              names(x_var) = covar
              p     = plot_tran_covar(x_var, dep, 
                                      title = title, x.lab = covar, y.lab = 'P', markers = marks)  	
              return(p)
            }            
          
            if (type == 'HMM-stationary-prob') {
              title = paste("Stationary transition probabilities (", x@UID,')',sep='')
              marks = x@breakpoint$psi
              names(marks) = c('C', 'H') 
              x_var = data.frame(x@HMM$statMetrics$distr.covar[,covar])
              names(x_var) = covar
              dep   = x@HMM$statMetrics$distr.covar[,-which(names(x@HMM$statMetrics$distr.covar) == covar)]
              dep   = list(as.matrix(dep))
              
              # plot 
              p     = plot_tran_covar(x_var, dep,
                                      title = title, x.lab = covar, y.lab = 'P', markers = marks)    
              return(p)
            }            
          
            if (type == 'HMM-err-horiz') {
              title = paste("Out-of-sample decoding performance (", x@UID,')',sep='')
              perf  = x@HMM$PredictStats
              x_var = rownames(x@HMM$PredictStats)
              perf  = data.frame(perf)
              perf$Horizon = rownames(perf)
              df    = melt(perf, id.vars = 'Horizon')
              df$Horizon = as.numeric(df$Horizon)
              # plot 
              plt = ggplot(df, aes(x = Horizon, y = value))
              plt = plt + geom_point(aes(color = variable, shape = variable), show_guide = FALSE, size=4)
              plt = plt + facet_wrap(~variable, ncol = 2, scales = 'free')
              
              # formatting
              plt = plt + 
                theme_bw() +
                theme(panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.background = element_blank(),
                      strip.text.x     = element_text(size=18),
                      axis.text.y      = element_text(size=15),
                      axis.text.x      = element_text(size=15),
                      axis.title.y     = element_text(size=18),
                      axis.title.x     = element_text(size=18),
                      plot.title       = element_text(size=20),
                      legend.text      = element_text(size=15),
                      legend.title     = element_text(size=15),
                      axis.ticks       = element_blank() ) + 
                ylab('Performance') + xlab('Time Horizon [hours]') + 
                theme(plot.title=element_text(family="Times", face="bold", size=20)) + 
                ggtitle( title )
              return(plt)
            }            
            
            if (type == 'HMM-time-temp') {
              title = paste("State duration (", x@UID,')',sep='')
              marks = x@breakpoint$psi
              names(marks) = c('C', 'H') 
              tstart= round(marks[1])
              tstop = round(marks[2])              
              dep   = as.data.frame(x@HMM$statMetrics$time.temp[tstart:tstop,])
              names(dep) = paste('State', 1:ncol(dep))
              dep[,covar] = tstart:tstop
              df.mlt = melt(dep, id.vars = covar)
              
              # plot 
              plt = ggplot(df.mlt, aes_string(x = covar, y = 'value'))
              plt = plt + geom_point(aes(color = variable, shape = variable), size = 4)
              
              # formatting
              plt = plt + 
                theme_bw() +
                theme(panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.background = element_blank(),
                      strip.text.x     = element_text(size=18),
                      axis.text.y      = element_text(size=15),
                      axis.text.x      = element_text(size=15),
                      axis.title.y     = element_text(size=18),
                      axis.title.x     = element_text(size=18),
                      plot.title       = element_text(size=20),
                      legend.text      = element_text(size=15),
                      legend.title     = element_text(size=15),
                      axis.ticks       = element_blank() ) + 
                ylab('Duration') + xlab(covar) + 
                theme(plot.title=element_text(family="Times", face="bold", size=20)) + 
                ggtitle( title )
              return(plt)
            }              
            
            if (type == 'HMM-aggregate-season') {
              title = paste("Seasonal analysis (", x@UID,')',sep='')              
              df.mlt = melt(x@HMM$components$agg.tod, id.vars = c('Day', 'Hour', 'Season'))
              
              # plot 
              plt = ggplot(df.mlt, aes(x = Hour, y = value, color = Season, shape = Day))
              plt = plt + geom_point(size = 4) + geom_line()
              plt = plt + facet_wrap(~variable, scales = 'free')
              
              # formatting
              plt = plt + 
                theme_bw() +
                theme(panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.background = element_blank(),
                      strip.text.x     = element_text(size=18),
                      axis.text.y      = element_text(size=15),
                      axis.text.x      = element_text(size=15),
                      axis.title.y     = element_text(size=18),
                      axis.title.x     = element_text(size=18),
                      plot.title       = element_text(size=20),
                      legend.text      = element_text(size=15),
                      legend.title     = element_text(size=15),
                      axis.ticks       = element_blank() ) + 
                ylab('kWh') + xlab(covar) + 
                theme(plot.title=element_text(family="Times", face="bold", size=20)) + 
                ggtitle( title )
              return(plt)
            }                                  
          })

# _________________________________________
# Method to dump OccupancyStates model data to file
setGeneric(
  name = "dumpComputationToFile",
  def = function(.Object, verbose = T, path = NULL){standardGeneric("dumpComputationToFile")}
)
setMethod('dumpComputationToFile',
          signature  = 'OccupancyStates',
          definition = function(.Object, verbose = T, path = NULL){
          if (verbose) 
            cat(paste('*** Dumping model data to file (', .Object@UID, ')***\n', sep=''))
                    
          # ______________________________
          # Format analysis results: OLS
          
          df     = as.data.frame(t(coef(.Object@OLS$fit.summary)[,1]))
          df$R2  = .Object@OLS$fit.summary$adj.r.squared
          df$GR2 = .Object@OLS$GR2
          df$UID = .Object@UID
          df$ZIPCODE= .Object@ZIPCODE
      	  df$dwTest = .Object@OLS$dw.test
      	  df$bpTest = .Object@OLS$bp.test 
          df.ols = df
          
          # ______________________________
          # Format analysis results: HMM
          
          # stats
          df        = as.data.frame(.Object@HMM$entropy)
          names(df) = paste('entropy',names(df), sep='.')          
          df$UID    = .Object@UID
          df$ZIPCODE= .Object@ZIPCODE
          df$GR2    = .Object@HMM$GR2
          df$R2     = .Object@HMM$R2
          df$MAPE   = .Object@HMM$MAPE
          df$MAPE.cv= .Object@HMM$cv.metrics[as.character(.Object@HMM$nStates),'MAPE.cv']
          df$R2.cv  = .Object@HMM$cv.metrics[as.character(.Object@HMM$nStates),'R2.cv']
          tmp       = as.data.frame(.Object@HMM$predict)
          names(tmp)= paste('predict',names(tmp), sep='.')
          df        = cbind(df, tmp)
          df        = cbind(df, data.frame(t(.Object@breakpoint$psi)))
          df.stat   = df
          
          # response
          
          df        = as.data.frame(t(.Object@HMM$response$means))
          df$Sigma  = .Object@HMM$response$stdev
          df$UID    = .Object@UID
          df$ZIPCODE= .Object@ZIPCODE
          df$State  = 1:.Object@HMM$nStates
          df.resp   = df

          # benchmarks: temperature-dependent state distribution
          x_var     = .Object@HMM$statMetrics$distr.covar[,1]
          df        = as.data.frame(t(.Object@HMM$statMetrics$distr.covar[,-1]))
          names(df) = x_var
          df.distr  = cbind(UID = .Object@UID, ZIPCODE = .Object@ZIPCODE, State = 1:nrow(df), variable = 'Distribution', df)

          # benchmarks: temperature-dependent state duration
          df        = as.data.frame(t(.Object@HMM$statMetrics$time.temp))
          names(df) = x_var
          df.time   = cbind(UID = .Object@UID, ZIPCODE = .Object@ZIPCODE, State = 1:nrow(df), variable = 'Duration', df)
          
          df.tran   = rbind(df.distr, df.time)
          rownames(df.tran) = NULL
          
          # benchmarks: state breakdown by time-of-day
          df = data.frame(matrix(.Object@HMM$statMetrics$breakdown$Summer, ncol = 24, byrow=T))
          names(df) = 0:23
          df.summer = cbind(UID = .Object@UID, ZIPCODE = .Object@ZIPCODE, State = 1:nrow(df), Season = 'Summer', df)
          
          df = data.frame(matrix(.Object@HMM$statMetrics$breakdown$Winter, ncol = 24, byrow=T))
          names(df) = 0:23
          df.winter = cbind(UID = .Object@UID, ZIPCODE = .Object@ZIPCODE, State = 1:nrow(df), Season = 'Winter', df)
                  
          df.season = rbind(df.summer, df.winter)
          
      	  # ________________________
      	  # Component contributions
  	 
       	  df.comp  = .Object@HMM$components$agg.tod
      	  df.comp  = cbind(UID = .Object@UID, ZIPCODE = .Object@ZIPCODE, df.comp)

          # ______________________________
          # Dump to file
            
      	  if (!is.null(path)) {
      		  write.csv(df.ols,    file = paste(path,.Object@UID, '_ols.csv', sep=''), quote = F, row.names = F)
      		  write.csv(df.resp,   file = paste(path,.Object@UID, '_hmm_resp.csv', sep=''), quote = F, row.names = F)
      		  write.csv(df.tran,   file = paste(path,.Object@UID, '_hmm_tran.csv', sep=''), quote = F, row.names = F)
      		  write.csv(df.season, file = paste(path,.Object@UID, '_hmm_seas.csv', sep=''), quote = F, row.names = F)
      		  write.csv(df.comp,   file = paste(path,.Object@UID, '_hmm_comp.csv', sep=''), quote = F, row.names = F)
      	  }
          
          return(list(OLS = df.ols, 
                      HMM.response = df.resp, HMM.transition = df.tran, HMM.stats = df.stat,
                      HMM.season = df.season, HMM.comp = df.comp))
})

# _________________________________________
# Method to dump OccupancyStates model data to file
setGeneric(
  name = "dumpComputation",
  def = function(.Object, verbose = T, path = NULL){standardGeneric("dumpComputation")}
)
setMethod('dumpComputation',
          signature  = 'OccupancyStates',
          definition = function(.Object, verbose = T, path = NULL){
            if (verbose) 
              cat(paste('*** Dumping model data to RData file (', .Object@UID, ')***\n', sep=''))
            
            # trim down object for saving to disc
            tmp         = list()
            tmp$UID     = .Object@UID
            tmp$ZIP     = .Object@ZIPCODE
            tmp$OLS     = .Object@OLS[c('heterosked', 'serial.corr', 'R2', 'GR2')]
            tmp$breakpoint = .Object@breakpoint
            tmp$HMM     = .Object@HMM[c('nStates', 'cv.metrics', 'PredictStats',
                                        'response', 'BIC', 'sw_test', 'MAPE',
                                        'GR2', 'R2', 'state.duration', 'time.temp',
                                        'P.trans', 'statMetrics', 'fit')]
            tmp$HMM$states = data.frame(states = .Object@HMM$states[,1],
                                       DateTime = .Object@data.train$timestamps)
            tmp$HMM$entropy = unlist(.Object@HMM$entropy)
            tmp$HMM$predict = unlist(.Object@HMM$predict)
            tmp$HMM$comps   = .Object@HMM$components
                      
            # save to disc if requested
            if (!is.null(path)) 
              save(file = paste(path, .Object@UID, '.RData', sep=''), list = 'tmp')
            
            return(tmp)
          })

