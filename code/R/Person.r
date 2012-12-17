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
library('data.table')

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
            }        
            if (!is.null(object@HMM)) {
              cat(sprintf('HMM MARE =\t%f; HMM pval =\t%f\n', object@HMM$MARE, mean(object@HMM$p.values)))
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
  def = function(.Object, verbose = T){standardGeneric("fitOLS")}
)
setMethod('fitOLS',
          signature  = 'Person',
          definition = function(.Object, verbose=T){
            
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
                      
            # test for heteroskedasticity using Breusch-Pagan test
            require('lmtest')
            bp.test        = bptest(OLS.step, data = data.dum)
            chi2.95        = qchisq(0.95,bp.test$parameter+1)
            heterosc.BP.95 = bp.test$statistic > chi2.95
      
            # test for serial correlation using Durbin-Watson test
            dw.test    = dwtest(OLS.step, data = data.dum)
            sercorr.DW = dw.test$alternative == "true autocorrelation is greater than 0"
            
            # which coefficients are statistically significant at least at 0.1 level?
            summ.fit = summary(OLS.step)
            idx.sigf = which(summ.fit$coefficients[,4] <= 0.1 & abs(summ.fit$coefficients[,1])>1e-6 )
            coef.sig = summ.fit$coefficients[idx.sigf,]
            
            # save model result
            .Object@OLS = list()
            .Object@OLS[['fit.summary']] = summ.fit
            .Object@OLS[['fitted.vals']] = predict(OLS.step)
            .Object@OLS[['heterosked']]  = heterosc.BP.95
            .Object@OLS[['serial.corr']] = sercorr.DW
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

# _________________________________
# Method to perform HMM analysis

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
              while (!ok & it <= 1) {
                mod <- depmix(response  = fmla_response, 
                              transition= fmla_transitn,
                              data      = data.dum, 
                              nstates   = K, 
                              respstart = respstart,
                              trstart   = trstart,
                              family    = gaussian(),
                              instart   = runif(K))
                out <- capture.output(fm  <- try(fit(mod, verbose = T, useC = T, emcontrol = em.control(maxit = 50, tol = 5e-3))))
                nlines = length(out)
                if (class(fm) != 'try-error') {
                  ok = nlines > 2                                   
                } else ok = FALSE
                if (!ok){
                  BIC_cur = Inf                  
                  if (verbose) cat('Bad HMM fit; re-estimating...\n')
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
            .Object@HMM$sw_test = is_normal
                            
            # accuracy of prediction
            .Object@HMM[['MARE']] = mean(abs(residuals / .Object@HMM[['means']]))
            
            rm(list = c('result', 'data.dum', 'fm_opt'))
            gc()
            return(.Object)
          }
)
            

# _______________________
# Plot Person object

library(igraph)
library('useful')
library('reshape')
library('RColorBrewer')

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
            
            if (type == 'default') {
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
            
            if (type == 'HMM-ts') {    
              
              # construct plotting data frames
              df.hmm = data.frame(Mean   = hmm.means, 
                                  Sigma  = hmm.sigma, 
                                  State  = states,                                  
                                  Time   = format(timestamps, "%a,%m/%d %H:00"))              
              df.obs = data.frame(Observed = consumption, 
                                  Time   = format(timestamps, "%a,%m/%d %H:00"),
                                  Sigma    = rep(0, length(consumption)))
              
              df.mlt.obs = melt(df.obs, id.vars = c('Time', 'Sigma'))
              df.mlt.hmm = melt(df.hmm, id.vars = c('Time', 'Sigma', 'State'), measure.vars = c('Mean'))
              df.mlt.hmm$variable = paste(df.mlt.hmm$variable, df.mlt.hmm$State)
              df.mlt.hmm$State = NULL
              
              # plot current zoom-in
              pd <- position_dodge(.1)
              plt = ggplot(df.mlt.hmm, aes(x = as.numeric(Time), y = value)) + 
                geom_errorbar(aes(ymin = value-Sigma, ymax=value+Sigma), 
                              colour="black", width=.1, position=pd, alpha=0.2) +
                geom_line(position=pd) + geom_point(position=pd, aes(colour=variable)) + 
                geom_line(col='black', data = df.mlt.obs, alpha = 0.4) + 
                geom_point(col='black', data = df.mlt.obs, alpha = 0.4) + 
                scale_x_continuous(breaks = seq(1,length(df.mlt.hmm$Time), length.out=10), 
                                   labels = format(df.mlt.hmm$Time[seq(1,length(df.mlt.hmm$Time), length.out=10)], 
                                                   format = "%a,%m/%d %H:00")) + 
                theme_bw() +
                theme(panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(),
                       panel.background = element_blank(),
                       axis.title.x = element_blank(),
                       panel.background = element_rect(fill = "transparent",colour = NA),
                       axis.ticks = element_blank() ) + 
                ylab("kWh") + 
                theme(plot.title=element_text(family="Times", face="bold", size=20)) + 
                ggtitle( paste("Zoom-In: States and Observed Emissions (", x@UID,')',sep=''))
              
              return(plt)
            }
            
            
            # plot state space network for a given HSMM using ggplot
            if (type == 'HMM-MC') {
              df = data.frame(Mean   = as.numeric(x@HMM$response$means[1,]),
                              Sigma  = as.numeric(x@HMM$response$stdev),
                              State  = 1:x@HMM$nStates)
              state.grid = expand.grid(State.i = as.factor(1:1:x@HMM$nStates), 
                                       State.j = as.factor(1:1:x@HMM$nStates))
              trans = sapply(1:nrow(state.grid), function(s) {
                i = state.grid[s,1]
                j = state.grid[s,2]
                tr= paste(i,'->', j, sep='')
                ret = x@HMM$transition[1,tr]
              })
                            
              df.P = state.grid
              df.P$P = trans
              df.P$Mean.i  = as.numeric(as.character(df.P$State.i))
              df.P$Mean.j  = as.numeric(as.character(df.P$State.j))
              df.P$Sigma.i = as.numeric(x@HMM$response$stdev[as.character(df.P$State.i)])
              df.P$Sigma.j = as.numeric(x@HMM$response$stdev[as.character(df.P$State.j)])
              
              plt = ggplot(df) + 
                geom_point(aes(x = State, y = Sigma, color = State, size = Mean)) + 
                #scale_size(to = c(2, 12)) + 
                geom_segment(aes(x = Mean.i, xend = Mean.j, 
                                 y = Sigma.i, yend = Sigma.j, 
                                 size=P, alpha=P), data=df.P) +
                theme_bw() +
                theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     panel.background = element_blank(),
              #        axis.title.y = element_blank(),
              #        axis.title.x = element_blank(),
                     legend.position = "none",
                     panel.background = element_blank(),
                     axis.ticks = element_blank()) + 
                ylab('Sigma') + xlab('Mean') + 
                ggtitle(paste("State Space Diagram (", x@UID, ')',sep=''))
              return(plt)
            }
            
            # plot confidence interval
            if (type == 'HMM-ci') {  
              df = data.frame( Prob = p.vals, State = as.factor(states), 
                               Weekday = weekdays(timestamps), Hour = hour(timestamps) )
              df.st = aggregate(Prob ~ State, FUN = mean, data = df)
              names(df.st) = c('State')
              df.mu = aggregate(Prob ~ Weekday + Hour, FUN = mean, data = df)
              names(df.mu) = c('Weekday', 'Hour', 'Mean')
              df.sd = aggregate(Prob ~ Weekday + Hour, FUN = sd, data = df)
              names(df.sd) = c('Weekday', 'Hour', 'Sigma')
              df.sd.mu = merge(df.sd, df.mu)
              
              pd <- position_dodge(.1)
              plt = ggplot(df.sd.mu, aes(x = Hour, y = Mean, color = Weekday)) + 
                geom_errorbar(aes(ymin = Mean-Sigma, ymax=Mean+Sigma), 
                              colour="black", width=.1, position=pd, alpha=0.2) +
                geom_line(position=pd) + geom_point(position=pd) + 
                theme_bw() +
                theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     panel.background = element_blank(),
                     panel.background = element_rect(fill = "transparent",colour = NA),
                     axis.ticks = element_blank()) + 
                     ylab('Confidence') + xlab('Hour') + 
                ggtitle( paste("Avg. Prediction Confidence (", x@UID, ')',sep=''))
              
              return(plt)
            }
            
            if (type =='HMM-acf') {
              plt = acf_ggplot(na.omit(consumption), 
                               na.omit(hmm.means), 
                               title = 'ACF: Empirical vs Model')
              return(plt)
            }
            
            if (type == 'HMM-res') {              
              df = data.frame(State = factor(states), Residual = hmm.residual)
              plt = ggplot(df) + 
                geom_density(aes(x = Residual, colour = State), size=2) + 
                theme(legend.position = c(0.9,0.7)) + 
                theme_bw() + facet_wrap( ~ State, nrow=1, scales = 'free') + 
                theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     panel.background = element_blank(),
                     panel.background = element_rect(fill = "transparent",colour = NA),
                     axis.ticks = element_blank()) + 
                     ylab('pdf') + xlab('Residual') + 
                ggtitle("HMM Residuals")
              return(plt)
            }
            
            if (type == 'HMM-qq') {
              df      = data.frame(State = factor(states), Residual = hmm.residual)
              pvector = x@HMM$sw_test
              levels(df$State) = paste(levels(df$State), round(pvector, digits=5), sep=':')
              for (j in 1:x@HMM[['nStates']]) {
                idx = which(df$Residual == j)
                df$Residual[idx] = df$Residual[idx] / x@HMM$response$stdev[j]^2
              }
              plt = ggplot(df, aes(sample = Residual)) + facet_wrap(~State) +  
                stat_qq(geom = "point", size = 2, position = "identity") + 
                        #dparams = x@HMM$response$stdev) + 
                theme(legend.position = c(0.9,0.7)) + 
                theme_bw() + facet_wrap( ~ State, nrow=1) + 
                theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     panel.background = element_blank(),
                     panel.background = element_rect(fill = "transparent",colour = NA),
                     axis.ticks = element_blank()) + 
                     ylab('Sample Quantiles') + xlab('Theoretical Quantiles')
              return(plt)
            }
            
            if (type == 'HMM-MC2') {
              
              # organize MC information
              df = data.frame(Mean   = as.numeric(x@HMM$response$means[1,]),
                              Sigma  = as.numeric(x@HMM$response$stdev),
                              State  = 1:x@HMM$nStates)              
              state.grid = expand.grid(State.i = 1:x@HMM$nStates, State.j = 1:x@HMM$nStates)
              trans = sapply(1:nrow(state.grid), function(s) {
                i = state.grid[s,1]
                j = state.grid[s,2]
                tr= paste(i,'->', j, sep='')
                ret = x@HMM$transition[1,tr]
              })
                        
              # form igraph object
              g = graph.empty()              
              g = add.vertices(g, nrow(df), State=as.character(df$State), Mean =df$Mean, Sigma = df$Sigma)              
              g = add.edges(g, t(as.matrix(state.grid)), P = trans)                            
              
              # add graph plotting attributes
              scale <- function(v, a, b) {
                v <- v-min(v) ; v <- v/max(v) ; v <- v * (b-a) ; v+a
              }              
              pallete    = colorRampPalette(brewer.pal(9,"Blues"))(100)
              V(g)$color = 'grey'              
              V(g)$size  <- scale(V(g)$Mean, 10, 20)
              E(g)$color <- pallete[scale(E(g)$P, 10,100)]
              
              # plot
              plot(g, main = 'MC Structure')              
            }
})


# _________________________________________
# Method to dump Person model data to file
setGeneric(
  name = "dumpComputationToFile",
  def = function(.Object, verbose = T, path = './'){standardGeneric("dumpComputationToFile")}
)
setMethod('dumpComputationToFile',
          signature  = 'Person',
          definition = function(.Object, verbose = T, path = './'){
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
          df$Transit= rownames(df)
          rownames(df) = NULL
          df$UID    = .Object@UID
          df$ZIPCODE= .Object@ZIPCODE
          df.tran = df
          
          # ______________________________
          # Dump to file
          
          write.csv(df.ols, file = paste(path,.Object@UID, '_ols.csv', sep=''), quote = F, row.names = F)
          write.csv(df.resp, file = paste(path,.Object@UID, '_hmm_resp.csv', sep=''), quote = F, row.names = F)
          write.csv(df.tran, file = paste(path,.Object@UID, '_hmm_tran.csv', sep=''), quote = F, row.names = F)

          return(list(OLS = df.ols, HMM.response = df.resp, HMM.transition = df.tran))
})


# # _______________________
# # Test constructor
# 
# source('~/Dropbox/ControlPatterns/code/R/utils/sql_utils.r')
# source('~/Dropbox/ControlPatterns/code/R/utils/timing.r')
# raw_data  = run.query("select * from pge_res_final3_unique LIMIT 0,1000")
# raw_data  = subset(raw_data, UID == raw_data$UID[1])
# test      = new(Class='Person', raw_data, unique(raw_data$UID), unique(raw_data$ZIP5))
# 
# query     = paste("SELECT * FROM ZIP_", unique(raw_data$ZIP5), sep='')
# wthr_data = run.query(query, db = 'PGE_WEATHER')
# test      = addWeather(test, wthr_data)
# test      = prepareData(test, trends = c(8,24, 24*30*6))
# 
# # _______________________
# # Test modeling methods
# 
# hourly_vars   = paste('HourOfDay', 0:23, sep='.')
# monthly_vars  = paste('Month', 1:12, sep = '.')
# weekly_vars   = paste('DayOfWeek', c('Mon', 'Tue', 'Wed', 'Thu', 'Fri', 'Sat', 'Sun'), sep='.')
# wthr_vars_all = c('TemperatureF', 'WindSpeed', 'Humidity', 'HourlyPrecip', 'SolarRadiation')
# trend_vars    = c('Trend.8', 'Trend.24', 'Trend.4320')
# holiday_vars  = c('IsHoliday.0', 'IsHoliday.1')
# all_covars    = c(hourly_vars, monthly_vars, weekly_vars, wthr_vars_all, trend_vars, holiday_vars)
# 
# test          = fitOLS(test)
# test          = fitHMM(test, constrMC = NULL, Kmin = 2, Kmax = 6, ols_vars = T,
#                        response_vars = c(hourly_vars, wthr_vars_all, trend_vars))
# 
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
# 
# png(filename = 'HMM_ts.png', width = 1200, height = 800)
# print(plot(test, type = 'HMM-ts', interval = c('2010-12-01', '2010-12-10')))
# dev.off()
# 
# png(filename = 'HMM_MC.png', width = 1200, height = 800)
# print(plot(test, type = 'HMM-MC', interval = c('2010-12-01', '2010-12-10')))
# dev.off()
# 
# png(filename = 'HMM_ci.png', width = 1200, height = 800)
# print(plot(test, type = 'HMM-ci', interval = c('2010-12-01', '2010-12-10')))
# dev.off()
# 
# png(filename = 'HMM_MC2.png', width = 1200, height = 800)
# plot(test, type = 'HMM-MC2', interval = c('2010-12-01', '2010-12-10'))
# dev.off()
# 
# png(filename = 'HMM_qq.png', width = 1200, height = 800)
# print(plot(test, type = 'HMM-qq', interval = c('2010-12-01', '2010-12-10')))
# dev.off()
# 
# png(filename = 'HMM_res.png', width = 1200, height = 800)
# print(plot(test, type = 'HMM-res', interval = c('2010-12-01', '2010-12-10')))
# dev.off()
