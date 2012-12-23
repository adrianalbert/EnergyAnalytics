# Person.r
#
# The main purposes of  this class is to encapsulate pre-processing and model fitting functions where
# the unit of analysis is a (person,premise) tuple. It should have slots for 
# - reading in unformatted data and formatting it accordingly
# - fitting OLS/FGLS models
# - OLS residual analysis
# - ACF analysis
# - HMM analysis w/ covariates
# - plotting: OLS, structure, HMM structure
# - formatted output: print, data frame for further analysis
# - save RData file
# 
# Adrian Albert
# Last modified: December 2012.
# -----------------------------------------------------------------------

# Person class: data & analysis for a unique tuple (PER_ID, SP_ID)
# -----------------------------------------------------------------

library('methods')
library('zoo')
library('lubridate')
require('lmtest')

# clean-up previous definitions of methods for class Person
removeClass('Person')

options(error = recover)

# ________________________
# Class definition

setClass(
  Class = "Person",
  representation = representation(
    UID        = "numeric",           # unique person ID
    ZIPCODE    = "numeric",           # current premise zipcode
    N_OBS      = "numeric",           # number of data points
    N_DAYS     = "numeric",           # number of days (= N_OBS * 24)
    timestamps = "character",         # vector of timestamps
    weather    = "data.frame",        # weather covariates
    consumption= "numeric",           # hourly consumption 
    OLS        = "list",              # summary of OLS model
    HMM        = "list",               # summary of HMM model
    data.tmp   = 'data.frame'
    )
)

# all covariates
hourly_vars   = paste('HourOfDay', 0:23, sep='.')
monthly_vars  = paste('Month', 1:12, sep = '.')
weekly_vars   = paste('DayOfWeek', c('Mon', 'Tue', 'Wed', 'Thu', 'Fri', 'Sat', 'Sun'), sep='.')
wthr_vars_all = c('TemperatureF', 'WindSpeed', 'Humidity', 'HourlyPrecip', 'SolarRadiation')
trend_vars    = c('Trend.8', 'Trend.24', 'Trend.4320')
holiday_vars  = c('IsHoliday.0', 'IsHoliday.1')
all_covars    = c(hourly_vars, monthly_vars, weekly_vars, wthr_vars_all, trend_vars, holiday_vars)

# _______________________________________
# Constructor method for class Person

setMethod(f = "initialize", 
          signature = "Person",
          definition = function(.Object, raw_data, UID, ZIP5, res = 60, verbose = T, log = F, long = T) {
            
            .Object@UID     = UID
            .Object@ZIPCODE = ZIP5
            
            if (verbose) 
              cat(paste('*** Initializing Person (', .Object@UID, ',', .Object@ZIPCODE,') ***\n', sep=''))
                        
            # ______________________________________________
            # Format consumption data to timeseries format

            cons_names  = paste('hkw', 1:(24*60/res), sep='')
            if (long) {                                            # long format (Nx1)                          
              # columns referring to consumption
              days        = raw_data$date              
              hours       = paste(rep(formatC(0:23, flag=0, width=2), length(days)), ':00:00', sep='')
              days_times  = rep(days, each = 24)
              days_times  = paste(days_times, hours, sep=' ')
              consumption = as.vector(t(raw_data[,cons_names]))
            } else {
              consumption = raw_data[,cons_names]         # day-profile format (Nx24)
              days_times  = raw_data$date 
            }
            
            # remove possible duplicates
            idx_dup = which(duplicated(days_times))
            if (length(idx_dup)>0) {
              raw_data = raw_data[-idx_dup,]
              days_times=days_times[-idx_dup]
            }                        
            .Object@timestamps  = days_times
            
	    # is everthing else NA?
	    if (length(na.omit(consumption)) == 0) {
	      stop('Error: consumption data is all NAs!')
	    }

            # if data is all zeros
            rc = range(na.omit(consumption))
            if (rc[2] - rc[1] == 0) {
              stop('Error: data supplied to constructor has no variation!')
            }

            # is the log of consumption desired?
            if (log) {
              tmp = log(consumption)
              tmp[which(!is.finite(tmp))] = 0
              .Object@consumption = tmp
            } else
              .Object@consumption = consumption
                        
            # initialize default weather object
            .Object@weather     = data.frame()
          
            # _________________________________________
            # Compute summary information of interest
            
            .Object@N_OBS = length(.Object@consumption) 
            .Object@N_DAYS  = .Object@N_OBS / 24
            
            return(.Object)
          })


# ____________________________________________________________
# Method to add in exogenous covariates (weather, billing)

setGeneric(
  name = "addWeather",
  def = function(.Object, wthr_data, verbose = T){standardGeneric("addWeather")}
)
setMethod('addWeather',
          signature  = 'Person',
          definition = function(.Object, wthr_data, verbose=T){
            if (verbose) 
              cat(paste('*** Adding weather data for Person', .Object@UID,' ***\n',sep=''))
            
            # select covariates of interest
            wthr_times  = wthr_data$date
            wthr_names  = c('TemperatureF', 'DewpointF', 'Pressure', 'WindSpeed', 'Humidity', 
                            'HourlyPrecip', 'SolarRadiation')
            time_names  = c('DayOfWeek', 'HourOfDay', 'Month', 'IsHoliday')
            wthr_data   = wthr_data[,intersect(names(wthr_data), c(wthr_names, time_names))]            
              
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
            
            # is there weather data at times that match consumption data?            
            idx_ok      = which(wthr_times %in% .Object@timestamps)            
            if (length(idx_ok) > 0) {
              wthr_data = wthr_data[idx_ok,]
              wthr_times= wthr_times[idx_ok]
            } else {
              stop(paste('!!! No weather data available for user', .Object@UID))
            }
            idx_ok_cons  = which(.Object@timestamps %in% wthr_times)            
                                                
            # update object
            .Object@weather     = wthr_data
            .Object@timestamps  = wthr_times
            .Object@consumption = .Object@consumption[idx_ok_cons]
            .Object@N_OBS       = length(.Object@consumption) 
            .Object@N_DAYS      = round(.Object@N_OBS / 24)

            rm(list = c('wthr_data', 'wthr_names', 'wthr_times'))
            gc()
            return(.Object)
          }
)

# ____________________________________________________________
# Print method for class Person: display useful stats.

setMethod('show',
          signature  = 'Person',
          definition = function(object){
            cat('*** Person Object ***\n')
            
            # basic info
            cat(paste('Person ID', object@UID, '\n'))            
            cat(sprintf('# Days = \t%d; \n# Observations =%d\n', object@N_DAYS, object@N_OBS))
            cat(sprintf('Zipcode = \t%d\n', object@ZIPCODE))
            cat(sprintf('# Time range: \t%s - %s\n', object@timestamps[1], object@timestamps[object@N_OBS]))   
            
            # info on data quality
            wthr_names  = c('TemperatureF', 'DewpointF', 'Pressure', 'WindSpeed', 'Humidity', 
                            'HourlyPrecip', 'SolarRadiation')
            wthr_names  = intersect(wthr_names, names(object@weather))
            wthr_na_col = sapply(wthr_names, function(j){
              length(which(is.na(object@weather[,j])))
            })
            kwh_na      = length(which(is.na(object@consumption)))
            
            cat(paste(wthr_names, '\t:', wthr_na_col/nrow(object@weather)*100, '% NAs\n', sep=''))
            cat(paste('kWh\t\t:', kwh_na/nrow(object@weather)*100, '% NAs\n', sep=''))
            
            # info on model fit 
            if (!is.null(object@OLS)) {
              cat(sprintf('OLS MARE =\t%f; OLS pval =\t%f\n', object@OLS$MARE, mean(object@OLS$p.values)))
	      cat('OLS contributions:\n')
	      cat(paste(names(object@OLS$components$stat), ':\t', round(object@OLS$components$stat, digits=3), '\n')) 
            }        
            if (!is.null(object@HMM)) {
              cat(sprintf('HMM MARE =\t%f; HMM pval =\t%f\n', object@HMM$MARE, mean(object@HMM$p.values)))	
	      cat(sprintf('HMM States:\t%d\n', object@HMM$nStates))
	      cat('HMM contributions:\n')
	      print(round(object@HMM$components$stat, digits=3))
            }
            
            cat('*** END: Person Object ***\n')
          })

# ____________________________________________________
# Auxiliary function for preparing data for modelling

setGeneric(
  name = "prepareData",
  def = function(.Object, trends = NULL){standardGeneric("prepareData")}
)
setMethod('prepareData',
          signature  = 'Person',
          definition = function(.Object, trends = NULL) {
            
            # add consumption and weather, if any
            data          = data.frame(consumption = .Object@consumption)                      
            time_vars_all = c('HourOfDay', 'DayOfWeek', 'Month', 'IsHoliday')
            wthr_vars_all = c('TemperatureF', 'WindSpeed', 'Humidity', 'HourlyPrecip', 'SolarRadiation')
            wthr_vars     = intersect(c(wthr_vars_all, time_vars_all), names(.Object@weather))            
            if (length(wthr_vars)>0) data = cbind(data, .Object@weather[,wthr_vars])
            data$time     = 1:nrow(data)
            
            # create dummies
            data.dum      = dummy.data.frame(data=data, sep='.', all=T)
            rm.col        = sapply(time_vars_all, function(c){
              idx = which(regexpr(c, names(data.dum))>0)[1]
            })
            data.dum      = data.dum[,-rm.col]

            # add seasonal trends
            if (!is.null(trends)) {
              trend_vars = paste('Trend',trends,sep='.')            
              for (DT in trends) {
                data.dum[,paste('Trend',DT,sep='.')] = sin(2 * pi * data.dum$time / DT)
              }
            }
            
            .Object@data.tmp = data.dum
            return(.Object)
})

# _________________________________
# Method to perform OLS analysis

library('timeDate')
library('dummies')
setGeneric(
  name = "fitOLS",
  def = function(.Object, verbose = T, stats = T){standardGeneric("fitOLS")}
)
setMethod('fitOLS',
          signature  = 'Person',
          definition = function(.Object, verbose=T, stats = T){
            
            if (verbose)
              cat(paste('*** OLS analysis for Person-SPID (', .Object@UID, ')***\n', sep=''))
            
            # ______________________________________________________________
            # Set-up OLS problem (regress consumption on weather + tod/dow)
            
            # get data
            data.dum = na.omit(.Object@data.tmp)
            
            # set-up model covariates
            predictors = setdiff(names(data.dum), c('consumption', 'time'))
            fmla_null = as.formula("consumption ~ 1")
            fmla_full = as.formula(paste('consumption ~ ', paste(predictors,collapse='+')))
            
            if (verbose) print(fmla_full)
            
            # perform regression            
#             OLS.null = lm( fmla_null, data = data.dum )
            OLS.step     = lm( fmla_full, data = data.dum )
            fmla_full_ok = as.formula(paste('consumption ~ ', paste(names(na.omit(coef(OLS.step)[-1])),collapse='+')))
            OLS.step     = lm( fmla_full_ok, data = data.dum )
#             OLS.step = step(OLS.null, direction = 'forward', trace = 0, 
#                               scope = list(upper = fmla_full, lower = fmla_null))
            
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
            .Object@OLS[['all.covars']]  = predictors
            rel_mare = abs(residuals(OLS.step) / data.dum$consumption)
            rel_mare[!is.finite(rel_mare)] = 0
            .Object@OLS[['MARE']]        = mean(rel_mare)
            .Object@OLS[['p.values']]    = pnorm(residuals(OLS.step))
                      
            rm(list=c('OLS.step', 'data.dum'))
            return(.Object)
          }
)

# ______________________________________________________
# Computes covariate constribution to predicted OLS fit
setGeneric(
  name = "computeContributionsOLS",
  def = function(.Object, verbose = T) {standardGeneric("computeContributionsOLS")}
)
setMethod('computeContributionsOLS',
          signature  = 'Person',
          definition = function(.Object, verbose=T) {
	    if (verbose) cat(paste('***** Computing OLS components for Person', .Object@UID, '*****\n'))
	    # compute contributions by covariate
	    b  = .Object@OLS$fit.summary$coefficients[,1]
	    X  = as.data.frame(t(b[-1] * t(.Object@data.tmp[,names(b)[-1]])))	    
	    
	    # aggregate DoW, ToD, Month
	    idx.tod  = names(X)[which(names(X) %in% hourly_vars)]
	    idx.dow  = names(X)[which(names(X) %in% weekly_vars )]
	    idx.mnth = names(X)[which(names(X) %in% monthly_vars)]
	    idx.trnd = names(X)[which(names(X) %in% trend_vars)]
	    idx.wthr = names(X)[which(names(X) %in% wthr_vars_all)]
	    if (length(idx.tod)>1) X.tod    = rowSums(X[,idx.tod]) else X.tod = X[,idx.tod]
	    if (length(idx.dow)>1) X.dow    = rowSums(X[,idx.dow]) else X.dow = X[,idx.dow]
	    if (length(idx.mnth)>1) X.mnth  = rowSums(X[,idx.mnth]) else X.mnth = X[,idx.mnth]
	    if (length(idx.mnth)>1) X.trnd  = rowSums(X[,idx.trnd]) else X.trnd = X[,idx.trnd]
	    Xagg     = X[,idx.wthr]
	    if (length(idx.tod)>0) Xagg  = cbind(Xagg, data.frame(ToD = X.tod))
	    if (length(idx.dow)>0) Xagg  = cbind(Xagg, data.frame(DoW = X.dow))
	    if (length(idx.mnth)>0) Xagg = cbind(Xagg, data.frame(Month = X.mnth))
	    if (length(idx.trnd)>0) Xagg = cbind(Xagg, data.frame(Trend = X.trnd))
	    Xagg$mean= b[1]
	    
	    tot     = rowSums(abs(Xagg))
	    Xr      = as.data.frame(as.matrix(Xagg) / tot*100)
	    # store for later	    
	    .Object@OLS$components      = list()
	    .Object@OLS$components$stat = colMeans(Xr)
	    .Object@OLS$components$ts   = Xagg	
	    .Object@OLS$components$ts$fit             = .Object@OLS$fitted.vals
	    .Object@OLS$components$ts$consumption     = .Object@data.tmp$consumption	
	    rm(list=c('X', 'b', 'Xr', 'Xagg'))
	    gc()
	    return(.Object)
})

# _________________________________
# Method to perform HMM analysis339.94 

library('depmixS4')
library('R.utils')
setGeneric(
  name = "fitHMM",
  def = function(.Object, verbose = T, Kmin = 3, Kmax = 3, constrMC = NULL, 
                 response_vars = NULL, transitn_vars = NULL, ols_vars = T)
    {standardGeneric("fitHMM")}
)
setMethod('fitHMM',
          signature  = 'Person',
          definition = function(.Object, verbose=T, Kmin = 3, Kmax = 3, constrMC = NULL,
                                response_vars = NULL, transitn_vars = NULL, ols_vars = T){
            
            if (verbose)
              cat(paste('*** HMM analysis for Person-SPID (', .Object@UID, ')***\n', sep=''))
            
            # _____________________________
            # Set up depmixS4 model object
          
            # prepare dataset
            data.dum      = .Object@data.tmp
            existing_vars = names(data.dum)
                        
            # see which variables have been included that are reflected in the data as well
            response_vars = intersect(response_vars, existing_vars)            
            transitn_vars = intersect(transitn_vars, existing_vars)
            
            # do we include just variables that are significant in regression?
            if (ols_vars) {
              ols_vars_names= rownames(.Object@OLS[['fit.summary']]$coefficients)
#               ols_vars_names= rownames(.Object@OLS[['coef.signif']])
              if (length(ols_vars_names)>0) 
                response_vars = intersect(response_vars, ols_vars_names)             
            }

            # define response model
            fmla_response = 'consumption ~ 1'
            if (length(response_vars)>0) 
              fmla_response = paste('consumption ~ ',paste(response_vars,collapse='+'))
            fmla_response   = as.formula(fmla_response)
            
            if (verbose) print(fmla_response)
            
            # define transition model
            fmla_transitn = '~ 1'
            if (length(transitn_vars)>0) 
              fmla_transitn = paste('~ ',paste(transitn_vars,collapse='+'))
            fmla_transitn   = as.formula(fmla_transitn)
            
            # _______________________________________
            # choose model size K (number of states)
            
            # set up and fit model
            compute_K = function(K) {
              
              # initialize state parameters with OLS estimates?
              respstart = NULL
              if (ols_vars) {   
                vars_ok = response_vars
                if ('(Intercept)' %in% ols_vars_names) {
                  vars_ok = c('(Intercept)',vars_ok)
                  ols_vars_coef = .Object@OLS[['fit.summary']]$coefficients[vars_ok]
                  #ols_vars_coef = .Object@OLS[['coef.signif']][vars_ok,1]
                } else {
#                   ols_vars_coef = .Object@OLS[['coef.signif']][vars_ok,1]  
                  ols_vars_coef = .Object@OLS[['fit.summary']]$coefficients[vars_ok]
                  ols_vars_coef = c(0,ols_vars_coef)
                }
                respstart     = rep(c(ols_vars_coef,0), K)             
#		respstart     = runif(length(respstart )) * respstart
              }
              
              # place structure on transition matrix?
              if (!is.null(constrMC)) {
                # define constraints on transition matrix: telescope model
                if (constrMC == 'telescope') {
                  A = matrix(0, nrow = K, ncol = K)
                  for (i in 1:K) {
                    if (i<K) A[i,i+1] = 1
                    A[i,i]   = 1
                    A[i,1]   = 1
                  }            
                }
                # define forward-jump model
                if (constrMC == 'fwdjump') {
                  A = matrix(0, nrow = K, ncol = K)
                  for (i in 1:K) {
                    A[i,i:K] = 1                    
                  }    
                    A[K,1] = 1
                }
                A = A / rowSums(A)
                A = as.numeric(A)    
                trstart = A
              } else trstart = NULL
                
              # unconstrained model
              
              ok = FALSE
              it = 0
	      cur_tol = 5e-4
              while (!ok & it <= 1) {
                mod <- depmix(response  = fmla_response, 
                              transition= fmla_transitn,
                              data      = data.dum, 
                              nstates   = K, 
                              respstart = respstart,
                              trstart   = trstart,
                              family    = gaussian(),
                              instart   = runif(K))
                out <- capture.output(fm  <- try(fit(mod, verbose = T, useC = T, emcontrol = em.control(maxit = 50, tol = cur_tol))))
                nlines = length(out)
                if (class(fm) != 'try-error') {
                  ok = nlines > 2                                   
                } else ok = FALSE
                if (!ok){
                  BIC_cur = Inf                  
                  if (verbose) cat('Bad HMM fit; re-estimating...\n')
		  cur_tol = cur_tol / 10
                } else {
                  if (verbose) cat(paste('Convergence in', nlines*5, 'iterations.\n'))
                  BIC_cur = BIC(fm)
                }
                it = it + 1
              }
              return(list(model = fm, BIC = BIC_cur, nlines = nlines))
            }
                     
            result = lapply(Kmin:Kmax, compute_K)
            BICs   = sapply(result, function(l)l[['BIC']])
            nlines = sapply(result, function(l)l[['nlines']])
            idx_opt= which.min(BICs)            
            K_opt  = (Kmin:Kmax)[idx_opt]
            fm_opt = result[[idx_opt]][['model']]
            
            if (!is.finite(BICs[idx_opt])) stop('HMM fitting error!')
            
            # ______________________
            # Save computation
            
            .Object@HMM = list()
            
            .Object@HMM$converg = nlines[idx_opt]
            
            # response parameters
            .Object@HMM[['response']] = list()
            stddev = sapply(1:K_opt, function(j) fm_opt@response[[j]][[1]]@parameters$sd)
            names(stddev) = 1:K_opt
            .Object@HMM[['response']][['stdev']] = stddev
            params = sapply(1:K_opt, function(j) fm_opt@response[[j]][[1]]@parameters$coefficients)
            if (class(params) == 'numeric') params = data.frame(params) else params = data.frame(t(params))
            rownames(params) = 1:K_opt
            colnames(params) = c('Intercept', response_vars)
            .Object@HMM[['response']][['means']]  = as.data.frame(t(params))
            
            # transition parameters
            .Object@HMM[['transition']] = list()
            params = lapply(1:K_opt, function(j) {
              tmp = fm_opt@transition[[j]]@parameters$coefficients
              if (class(tmp) != 'numeric') return(data.frame(tmp)) else return(data.frame(t(tmp)))
            })
            params = do.call(cbind, lapply(params, data.frame, stringsAsFactors=FALSE))
            if (!is.null(ncol(params))) {
              colnames(params) = as.vector(sapply(1:K_opt, function(x) paste(x, '->',1:K_opt,sep='')))
            } else {
              names(params) = as.vector(sapply(1:K_opt, function(x) paste(x, '->',1:K_opt,sep='')))
            }
            if (!is.null(nrow(params))) rownames(params) = c('Intercept', transitn_vars)
            .Object@HMM[['transition']] = params
            
            # Viterbi states
            .Object@HMM[['states']] = posterior(fm_opt)
            
            # other parameters
            .Object@HMM[['nStates']]    = fm_opt@nstates
            
            # ________________________________________
            # Perform preliminary analysis of errors
            
            # Predicted state means
            means = sapply(1:K_opt, function(k) predict(fm_opt@response[[k]][[1]])) 
            .Object@HMM[['means']] = sapply(1:nrow(means), function(j) means[j,.Object@HMM[['states']][j,1]])

            # Confidence p-values
            hmm.sigma = .Object@HMM$response$stdev[.Object@HMM[['states']][,1]]  
            z.scores  = abs(.Object@consumption - .Object@HMM$means) / sqrt(hmm.sigma)
            .Object@HMM[['p.values']]   = pnorm(z.scores)

            # model fit (penalized likelihood)
            .Object@HMM[['BIC']]    = BIC(fm_opt)
            
            # residuals for each state
            residuals = .Object@consumption - .Object@HMM[['means']]                          
            .Object@HMM[['residual']] = residuals
            
            # test normality of residuals in each state
            is_normal = sapply(1:.Object@HMM[['nStates']], function(j) {
              idx = which(.Object@HMM[['states']][,1] == j)
              res = shapiro.test(residuals[idx])
              pval= res$p.value
              return(pval)
            })
	    #df = data.frame(p.values = .Object@HMM[['p.values']], state = .Object@HMM[['states']][,1])
	    #zz = aggregate(p.values ~ state, data = df, FUN = mean, na.rm = T)
            .Object@HMM$sw_test = is_normal#zz$p.values
                            
            # accuracy of prediction
            .Object@HMM[['MARE']] = mean(abs(residuals / .Object@HMM[['means']]))
            
            rm(list = c('result', 'data.dum', 'fm_opt'))
            gc()
            return(.Object)
          }
)           

# _____________________________________________________
# Computes covariate contribution to predicted HMM fit

setGeneric(
  name = "computeContributionsHMM",
  def = function(.Object, verbose = T)
    {standardGeneric("computeContributionsHMM")}
)
setMethod('computeContributionsHMM',
          signature  = 'Person',
          definition = function(.Object, verbose=T) {
	    if (verbose) cat(paste('***** Computing HMM components for Person', .Object@UID, '*****\n'))
	    # compute contributions by covariate
	    v = rownames(.Object@HMM$response$means)[-1]
	    B = as.data.frame(t(.Object@HMM$response$means[,.Object@HMM$states[,1]]))
	    X = .Object@data.tmp[,v] * B[,-1]
	    
	    # aggregate DoW, ToD, Month
	    idx.tod  = names(X)[which(names(X) %in% hourly_vars)]
	    idx.dow  = names(X)[which(names(X) %in% weekly_vars )]
	    idx.mnth = names(X)[which(names(X) %in% monthly_vars)]
	    idx.trnd = names(X)[which(names(X) %in% trend_vars)]
	    idx.wthr = names(X)[which(names(X) %in% wthr_vars_all)]
	    if (length(idx.tod)>1) X.tod    = rowSums(X[,idx.tod]) else X.tod = X[,idx.tod]
	    if (length(idx.dow)>1) X.dow    = rowSums(X[,idx.dow]) else X.dow = X[,idx.dow]
	    if (length(idx.mnth)>1) X.mnth  = rowSums(X[,idx.mnth]) else X.mnth = X[,idx.mnth]
	    if (length(idx.trnd)>1) X.trnd  = rowSums(X[,idx.trnd]) else X.trnd = X[,idx.trnd]
	    Xagg     = X[,idx.wthr]
	    if (length(idx.tod)>0) Xagg  = cbind(Xagg, data.frame(ToD = X.tod))
	    if (length(idx.dow)>0) Xagg  = cbind(Xagg, data.frame(DoW = X.dow))
	    if (length(idx.mnth)>0) Xagg = cbind(Xagg, data.frame(Month = X.mnth))
	    if (length(idx.trnd)>0) Xagg = cbind(Xagg, data.frame(Trend = X.trnd))
	    Xagg$mean= B[,1]

	    tot     = rowSums(abs(Xagg))
	    Xr      = as.data.frame(as.matrix(Xagg) / tot * 100)
	    Xr$state= .Object@HMM$states[,1]
	    # store for later	    
	    .Object@HMM$components      = list()
	    .Object@HMM$components$stat = aggregate(data = Xr, . ~ state, FUN = mean, na.rm = T)
	    .Object@HMM$components$ts   = Xagg		    
	    .Object@HMM$components$ts$fit             = .Object@HMM$means
	    .Object@HMM$components$ts$consumption     = .Object@data.tmp$consumption	
		
	    return(.Object)

})

# _______________________
# Plot Person object

setMethod('plot',
          signature  = 'Person',
          definition = function(x, type = 'default', interval = NULL, verbose = T, ...){
            
            if (verbose)
              cat(paste('*** ', type, ' Plot for Person (', x@UID, ') ***\n', sep=''))
            
            timestamps = as.POSIXct(x@timestamps)
            ols.resids = x@OLS[['fit.summary']]$residuals
            ols.resids = naresid(x@OLS$fit.summary$na.action, ols.resids)
            consumption= x@consumption

            ols.means  = x@OLS[['fitted.vals']]
            hmm.means  = x@HMM$means
            states     = x@HMM$states[,1]
            hmm.sigma  = x@HMM$response$stdev[states]
            hmm.residual=x@HMM$residual      
            wthr_data   =x@weather
            p.vals      =x@HMM$p.values
	    comps.ols   =x@OLS$components$ts
	    comps.hmm   =x@HMM$components$ts
            if (!is.null(interval)) {
              idx_ok      = which(x@timestamps >= as.POSIXct(interval[1]) & 
                                  x@timestamps <= as.POSIXct(interval[2]))
              timestamps  = timestamps[idx_ok]
              ols.resids  = ols.resids[idx_ok]
              consumption = consumption[idx_ok]
              ols.means   = ols.means[idx_ok]
              hmm.means   = hmm.means[idx_ok]
              hmm.sigma   = hmm.sigma[idx_ok]
              states      = states[idx_ok]
              hmm.residual= hmm.residual[idx_ok]
              wthr_data   = wthr_data[idx_ok,]
	      comps.ols   = comps.ols[idx_ok,]
	      comps.hmm   = comps.hmm[idx_ok,]
	      p.vals      = p.vals[idx_ok]
            }  

            if (type == 'default') {
print(length(consumption))
              plot(timestamps, consumption, 
                   main = paste('Consumption Person-SPID (', x@UID,')'),
                   xlab = 'Time', ylab = 'kWh', type = 'l', lwd = 2)
            }
            
            if (type == 'weather') {
              if (is.null(wthr_data) | nrow(wthr_data)==0) {
                cat('No weather data!\n')
              } else {
                wthr_names  = c('TemperatureF', 'DewpointF', 'Pressure', 'WindSpeed', 'Humidity', 
                                'HourlyPrecip', 'SolarRadiation')
                wthr = wthr_data[,names(wthr_data) %in% wthr_names]
                wthr$kWh = consumption
                wthr = zoo(wthr, order.by = timestamps)  
                plot(wthr, 
                     main = paste('Consumption Person-SPID (', x@UID,')'),
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
              yrange = range(c(ols.means, consumption))
              plot(timestamps, ols.means, 
                   main = 'OLS Fitted Values', ylab = 'Values', xlab = 'Time', 
                   lwd = 3, type='l', col='red', ylim = yrange)
              points(timestamps, consumption, type='b', pch=20, col = 'black')
              legend('topleft', c('OLS Fit', 'Observations'), col=c('red', 'black'))
            }    
            
	    if (type == 'OLS-contrib-ts') {
		title = paste("Zoom-In: OLS Covariate Contributions (", x@UID,')',sep='')
		p     = plot_components_ts(comps.ols, timestamps, title = title)
		
		return(p)
	    }

	    if (type == 'OLS-contrib-tot') {
		df = as.data.frame(t(x@OLS$components$stat))
		df$state = 1
		title = paste("OLS Covariate Contributions (", x@UID,')',sep='')
		p     = plot_tornado(df, title = title)
		
		return(p)
	    }

	    if (type == 'OLS-contrib-ts') {
		title = paste("Zoom-In: OLS Covariate Contributions (", x@UID,')',sep='')
		p     = plot_components_ts(comps.ols, timestamps, title = title)
		
		return(p)
	    }

            if (type == 'HMM-ts') {
		title = paste("Zoom-In: States and Observed Emissions (", x@UID,')',sep='')   
              
		plt = plot_hmm_ts(hmm.means, hmm.sigma, states, timestamps, consumption, 
				  y.lab = 'kWh', title = title)
                return(plt)
            }
            
            
            # plot state space network for a given HSMM using ggplot
            if (type == 'HMM-MC') {
		title = paste("State Space Diagram (", x@UID, ')',sep='')
		plt   = plot_HMM_MC(as.numeric(x@HMM$response$means[1,]),
				    as.numeric(x@HMM$response$stdev), x@HMM$transition, title=title)
                return(plt)
            }
            
            # plot confidence interval
            if (type == 'HMM-ci') { 
		title = paste("Avg. Prediction Confidence (", x@UID, ')',sep='') 
		plt   = plot_HMM_CI(p.vals, states, timestamps, title = title)
                return(plt)
            }
            
            if (type =='HMM-acf') {
              plt = acf_ggplot(na.omit(consumption), 
                               na.omit(hmm.means), 
                               title = 'ACF: Empirical vs Model')
              return(plt)
            }
            
            if (type == 'HMM-res') { 
		plt = plot_HMM_res(hmm.residual, states, x@HMM$response$stdev, sw_test = x@HMM$sw_test)
		plt
	    }             
            
            if (type == 'HMM-MC2') {
              title = paste("Avg. Prediction Confidence (", x@UID, ')',sep='') 
              # organize MC information
              plt = plot_HMM_MC_NET(as.numeric(x@HMM$response$means[1,]), as.numeric(x@HMM$response$stdev), x@HMM$transition)
	      plot(plt, main = title) 
            }

	    if (type == 'HMM-contrib-ts'){
		title = paste("Zoom-In: HMM Covariate Contributions (", x@UID,')',sep='')
		p     = plot_components_ts(comps.hmm, timestamps, title = title)
		return(p)
            }

	    if (type == 'HMM-contrib-tot') {
		df = x@HMM$components$stat
		title = paste("HMM Covariate Contributions (", x@UID,')',sep='')
		p     = plot_tornado(df, title = title, nrow = 2)
		
		return(p)
	    }
})


# _________________________________________
# Method to dump Person model data to file
setGeneric(
  name = "dumpComputationToFile",
  def = function(.Object, verbose = T, path = './', dump = TRUE){standardGeneric("dumpComputationToFile")}
)
setMethod('dumpComputationToFile',
          signature  = 'Person',
          definition = function(.Object, verbose = T, path = './', dump = TRUE){
          if (verbose) 
            cat(paste('*** Dumping model data to file (', .Object@UID, ')***\n', sep=''))

          # define all variables in the model          
          hourly_vars   = paste('HourOfDay', 0:23, sep='.')
          monthly_vars  = paste('Month', 1:12, sep = '.')
          weekly_vars   = paste('DayOfWeek', c('Mon', 'Tue', 'Wed', 'Thu', 'Fri', 'Sat', 'Sun'), sep='.')
          wthr_vars_all = c('TemperatureF', 'WindSpeed', 'Humidity', 'HourlyPrecip', 'SolarRadiation')
          trend_vars    = c('Trend.8', 'Trend.24')
          holiday_vars  = c('IsHoliday.0', 'IsHoliday.1')
          all_covars    = c(hourly_vars, monthly_vars, weekly_vars, wthr_vars_all, trend_vars, holiday_vars)
          
          # ______________________________
          # Format analysis results: OLS
          
          df     = as.data.frame(t(.Object@OLS$coef.signif[,1]))
          df[,setdiff(all_covars, names(df))] = NA
          df     = df[,all_covars]
          df$R2  = .Object@OLS$fit.summary$adj.r.squared
          df$MARE= .Object@OLS$MARE
          df$pval= mean(.Object@OLS$p.values)
          df$UID = .Object@UID
          df$ZIPCODE= .Object@ZIPCODE
	  df$dwTest = .Object@OLS$dw.test
	  df$bpTest = .Object@OLS$bp.test 
          df.ols = df
          
          # ______________________________
          # Format analysis results: HMM
          
          # response
          df        = as.data.frame(t(.Object@HMM$response$means))
          df[,setdiff(all_covars, names(df))] = NA
          df$Sigma  = .Object@HMM$response$stdev
          df$pval   = mean(.Object@HMM$p.values)
          df$MARE   = .Object@HMM$MARE
          df$UID    = .Object@UID
          df$ZIPCODE= .Object@ZIPCODE
          df$State  = 1:.Object@HMM$nStates
          df.resp = df

          # transition
          df = as.data.frame(t(.Object@HMM$transition))
          df$Transit   = rownames(df)
          rownames(df) = NULL
          df$UID       = .Object@UID
          df$ZIPCODE   = .Object@ZIPCODE
          df.tran      = df
          
	  # ________________________
	  # Component contributions
	 
    	  ols.comp         = as.data.frame(t(.Object@OLS$components$stat))
	  ols.comp$UID     = .Object@UID
	  ols.comp$ZIPCODE = .Object@ZIPCODE

     	  hmm.comp         = .Object@HMM$components$stat
	  hmm.comp$UID     = .Object@UID
	  hmm.comp$ZIPCODE = .Object@ZIPCODE

          # ______________________________
          # Dump to file
          
	  if (dump) {
		  write.csv(df.ols, file = paste(path,.Object@UID, '_ols.csv', sep=''), quote = F, row.names = F)
		  write.csv(df.resp, file = paste(path,.Object@UID, '_hmm_resp.csv', sep=''), quote = F, row.names = F)
		  write.csv(df.tran, file = paste(path,.Object@UID, '_hmm_tran.csv', sep=''), quote = F, row.names = F)
	  }
          return(list(OLS = df.ols, HMM.response = df.resp, HMM.transition = df.tran, OLS.comp = ols.comp, HMM.comp = hmm.comp))
})


 # _______________________
 # Test constructor


if (0 == 1) {
 source('~/EnergyAnalytics/code/R/utils/sql_utils.r')
 source('~/EnergyAnalytics/code/R/utils/plot_utils.r')
 source('~/EnergyAnalytics/code/R/utils/timing.r')
 
 raw_data  = run.query("select * from pge_res_final3_unique LIMIT 0,1000", user = 'adalbert', password = 'adrian')
 raw_data  = subset(raw_data, UID == raw_data$UID[1])
 test      = new(Class='Person', raw_data, unique(raw_data$UID), unique(raw_data$ZIP5), log = F)

 query     = paste("SELECT * FROM ZIP_", unique(raw_data$ZIP5), sep='')
 wthr_data = run.query(query, db = 'PGE_WEATHER', user = 'adalbert', password = 'adrian')
 test      = addWeather(test, wthr_data)
 test      = prepareData(test, trends = c(8,24, 24*30*6))


 # _______________________
 # Test modeling methods
 
 test          = fitOLS(test)
 test          = fitHMM(test, constrMC = NULL, Kmin = 6, Kmax = 6, ols_vars = T,
                        response_vars = c(hourly_vars, wthr_vars_all, trend_vars))
 
 # ________________________________________
 # Test formatting and aggregation methods

test	       = computeContributionsOLS(test, verbose = T)
test	       = computeContributionsHMM(test, verbose = T)

 
# # _______________________
# # Test plotting methods
# 
# png(filename = 'kwh.png', width = 1200, height = 800)
# plot(test, interval = c('2010-12-01', '2010-12-10'))
# dev.off()
# 
# png(filename = 'weather.png', width = 1200, height = 800)
# plot(test, type = 'weather', interval = c('2010-12-01', '2010-12-10'))
# dev.off()
# 
# png(filename = 'OLS.png', width = 1200, height = 800)
# plot(test, type = 'OLS-fit', interval = c('2010-12-01', '2010-12-10'))
# dev.off()

png(filename = 'OLS-contrib-ts.png', width = 1200, height = 800)
print(plot(test, type = 'OLS-contrib-ts'), interval = c('2010-12-01', '2010-12-10'))
dev.off()

png(filename = 'OLS-contrib-tot.png', width = 600, height = 600)
print(plot(test, type = 'OLS-contrib-tot'), interval = c('2010-12-01', '2010-12-10'))
dev.off()

png(filename = 'HMM_ts.png', width = 1200, height = 800)
print(plot(test, type = 'HMM-ts', interval = c('2010-12-01', '2010-12-10')))
dev.off()
 
png(filename = 'HMM_MC.png', width = 1200, height = 800)
print(plot(test, type = 'HMM-MC', interval = c('2010-12-01', '2010-12-10')))
dev.off()

png(filename = 'HMM_ci.png', width = 1200, height = 800)
print(plot(test, type = 'HMM-ci', interval = c('2010-12-01', '2010-12-10')))
dev.off()

png(filename = 'HMM_MC2.png', width = 1200, height = 800)
plot(test, type = 'HMM-MC2', interval = c('2010-12-01', '2010-12-10'))
dev.off()
 
png(filename = 'HMM_res.png', width = 1200, height = 800)
print(plot(test, type = 'HMM-res', interval = c('2010-12-01', '2010-12-10')))
dev.off()

png(filename = 'HMM-contrib-ts.png', width = 1200, height = 800)
print(plot(test, type = 'HMM-contrib-ts'), interval = c('2010-12-01', '2010-12-10'))
dev.off()

png(filename = 'HMM-contrib-tot.png', width = 1200, height = 600)
print(plot(test, type = 'HMM-contrib-tot', nrow = 2), interval = c('2010-12-01', '2010-12-10'))
dev.off()

}
