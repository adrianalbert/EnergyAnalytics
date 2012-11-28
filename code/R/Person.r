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
# Last modified: November 2012.
# -----------------------------------------------------------------------

# Person class: data & analysis for a unique tuple (PER_ID, SP_ID)
# -----------------------------------------------------------------

library('methods')
library('zoo')
library('lubridate')
removeClass('Person')

# ________________________
# Class definition

setClass(
  Class = "Person",
  representation = representation(
    PER_ID     = "numeric",           # unique person ID
    SP_ID      = "numeric",           # unique service point ID
    ZIPCODE    = "numeric",           # current premise zipcode
    N_OBS      = "numeric",           # number of data points
    N_DAYS     = "numeric",           # number of days (= N_OBS * 24)
    attr_const = "data.frame",        # Attributes of this tuple that are constant in time
    attr_vary  = "data.frame",        # Attributes of this tuple that vary in time
    timestamps = "POSIXct",           # vector of timestamps
    weather    = "zoo",               # weather covariates
    consumption= "zoo",               # hourly consumption 
    OLS        = "list",              # summary of OLS model
    HMM        = "list"               # summary of HMM model
    )
)

# _______________________________________
# Constructor method for class Person

setMethod(f = "initialize", 
          signature = "Person",
          definition = function(.Object, raw_data, res = 60, verbose = T, log = F) {
            
            .Object@PER_ID = unique(raw_data$PER_ID)
            .Object@SP_ID  = unique(raw_data$SP_ID)
            .Object@ZIPCODE= unique(raw_data$ZIP5)
            
            if (verbose) 
              cat(paste('*** Initializing Person-SP ( ', .Object@PER_ID, ',', .Object@SP_ID,') ***\n', sep=''))
                        
            # ______________________________________________
            # Clean up data: remove some obvious errors
            
            # remove very unfrequent Proxy IDs
            tab_psa_id = table(raw_data$PSA_ID) / nrow(raw_data)
            infreq_ids = which(tab_psa_id < 0.05)
            if (length(infreq_ids) > 0) {
              idx_rm   = which(raw_data$PSA_ID %in% names(infreq_ids))
              raw_data = raw_data[-idx_rm,]
            }

            # ______________________________________________
            # Format consumption data to timeseries format

            # columns referring to consumption
            days        = as.POSIXct(raw_data$date,tz="PST8PDT", '%Y-%m-%d')
            days        = seq(min(days), max(days), by = 3600*24)
            # create a row of hourly values for each day
            # sapply returns an array of numeric epoch seconds (for origin '1970-01-01')
            dateMat     = sapply(days,FUN=function(x) x + (0:(24-1) * 3600))            
            # flatten into a vector and re-convert into date objects
            .Object@timestamps  = as.POSIXct(as.vector(dateMat),origin='1970-01-01')
            idx_order   = order(.Object@timestamps)
            .Object@timestamps  = .Object@timestamps[idx_order]
            cons_names  = paste('hkw', 1:(24*60/res), sep='')
            kwhMat      = raw_data[idx_order,cons_names]
            consumption = as.vector(t(kwhMat))
            consumption = consumption[idx_order]
            
            # snap to common time vector
            orig        = .Object@timestamps[1]
            dest        = .Object@timestamps[length(.Object@timestamps)]
            time_vec    = seq(orig, dest, by=60*60)  
            if (log) {
              tmp = log(consumption)
              tmp[which(!is.finite(tmp))] = 0
              .Object@consumption = zoo(tmp, order.by = time_vec)
            } else
              .Object@consumption = zoo(consumption, order.by = time_vec)
            
            # initialize default weather object
            .Object@weather     = zoo(NULL)

            # _______________________________________________
            # Format attributes data: varying vs not varying
            
            # columns referring to person attributes (should not change in time)
            attr_const_names = c('CRSCODE', 'ENDUSE', 'CARE', 'STATUS', 'COMMTYPE', 
                                'FERA', 'CLIMSMRT', 'ACCTTYPE', 'DRPROG', 'CEEPROG', 
                                'max_total_duration',                   
                                'renter', 'sfo', 'owned_premises', 'PRIORITY', 'ua_uc_ru', 'median_income',
                                'median_income_quantiles', 'ro', 'vro', 
                                'good_sample_period', 'vro_total_good_sample_period', 'num_prems_for_vro',
                                'SP_ID', 'SA_TYPE', 'RSCHED', 'PSA_ID', 'SERCITY', 'GCOUNTY', 'DIVOFF', 
                                'CLIMATE', 'PREMTYPE', 'SM_SPST', 'NETMETER', 'MSTRMTR',
                                'ZIP5', 'AREA', 'WTHRSTN', 'CECCLMZN', 'DTONPREM_DATE', 'total_duration', 
                                'num_good_people_at_prem')
            # get unique values for attributes          
            # should be only one value per person!
            attr_levels = sapply(attr_const_names, function(c) { length(unique(raw_data[,c])) })
            mult_levels = which(attr_levels > 1)
            attr_vary_names = c()
            if (length(mult_levels)>0) {
              attr_vary_names  = c(attr_vary_names, attr_const_names[mult_levels])
              attr_const_names = attr_const_names[-mult_levels]
            }            
            # set person-specific attributes
            .Object@attr_const = as.data.frame(t(sapply(raw_data[,attr_const_names], unique)))
            # set premise-specific attributes
            .Object@attr_vary  = raw_data[,attr_vary_names]
            
            # _____________________________
            # Check for object consistency 
            
            #validObject(.Object)
            
            # _________________________________________
            # Compute summary information of interest
            
            .Object = computeStats(.Object)
            
            return(.Object)
          })

# __________________________________
# Validator method for class Person

validityPerson = function(object) {
  error = (length(object@timestamps) != length(object@consumption))
  if (error) cat('Timestamps and consumption have different lengths!\n')
  return (!error)
}
setValidity("Person", validityPerson)

# __________________________________
# Compute summary stats for object

setGeneric(
  name = "computeStats",
  def = function(.Object, verbose=T){standardGeneric("computeStats")}
)
setMethod('computeStats',
          signature  = 'Person',
          definition = function(.Object, verbose=T){
            if (verbose) 
              cat(paste('*** Computing stats for person (', .Object@PER_ID, ',', .Object@SP_ID, ') ***\n',sep=''))
            .Object@N_OBS   = length(.Object@consumption) 
            .Object@N_DAYS  = .Object@N_OBS / 24
            return(.Object)
          }
)

# ____________________________________________________________
# Method to add in exogenous covariates (weather, billing)

setGeneric(
  name = "addWeather",
  def = function(.Object, wthr_data, verbose = T, imputate = T){standardGeneric("addWeather")}
)
setMethod('addWeather',
          signature  = 'Person',
          definition = function(.Object, wthr_data, verbose=T, imputate = T){
            if (verbose) 
              cat(paste('*** Adding weather data for Person-SPID (', .Object@PER_ID, ',', .Object@SP_ID, ') ***\n',sep=''))
            
            # select covariates of interest
            wthr_time0  = as.POSIXct(wthr_data$date)
            wthr_names  = c('TemperatureF', 'DewpointF', 'Pressure', 'WindSpeed', 'Humidity', 
                            'HourlyPrecip', 'SolarRadiation')
            wthr_data   = wthr_data[,wthr_names]            
              
            # remove duplicate timestamps in weather data, if any
            idx_dup     = which(duplicated(wthr_time0))
            if (length(idx_dup) > 0) {
              wthr_time0= wthr_time0[-idx_dup]
              wthr_data = wthr_data[-idx_dup,]
            }
            
            # is there weather data at times that match consumption data?            
            idx_ok      = which(wthr_time0 >= min(.Object@timestamps) & 
                                wthr_time0 <= max(.Object@timestamps))
            if (length(idx_ok) == 0) return(.Object) else {
              wthr_data = zoo(wthr_data[idx_ok,], order.by = wthr_time0[idx_ok])
              wthr_time0= wthr_time0[idx_ok]
            }
                                    
            # create cleaned weather data object
            wthr_time   = as.POSIXct(seq(min(.Object@timestamps), max(.Object@timestamps), by=3600))
            wthr_data   = merge.zoo(wthr_data, zoo(order.by = wthr_time), all = c(F,T))
            names(wthr_data) = wthr_names
            
            # match dates for weather and consumption (assume both are ordered by date/time)
            .Object@weather = merge.zoo(wthr_data, zoo(0, order.by = as.POSIXct(.Object@timestamps)), all = c(F,T))
            
            .Object@weather = .Object@weather[,-ncol(.Object@weather)]
            
            # imputate missing values
            if (imputate == T) .Object = imputateMissingValues(.Object)
            
            rm(list = c('wthr_data', 'wthr_names', 'wthr_time'))
            return(.Object)
          }
)

# ____________________________________________________________
# Print method for class Person: display useful stats.

setMethod('show',
          signature  = 'Person',
          definition = function(object){
            cat('*** Person Object ***\n')
            cat(paste('Person ID', object@PER_ID, '; Service Point ID', object@SP_ID, '\n'))            
            cat(sprintf('# Days = %d; # Observations = %d\n', object@N_DAYS, object@N_OBS))
            cat(sprintf('Zipcode = %d\n', object@ZIPCODE))
            cat(sprintf('# Time range: %s - %s\n', object@timestamps[1], object@timestamps[object@N_OBS]))            
            cat('*** END: Person Object ***\n')
          })

# ___________________________________________
# Imputate missing observations using EM/SVD

library('Amelia')
library('imputation')
setGeneric(
  name = "imputateMissingValues",
  def = function(.Object,verbose = T){standardGeneric("imputateMissingValues")}
)
setMethod('imputateMissingValues',
          signature  = 'Person',
          definition = function(.Object, verbose=T){
            if (verbose) 
              cat(paste('*** Imputating missing values for Person-SPID (', .Object@PER_ID,',', .Object@SP_ID, ')***\n', sep=''))
            
              X.imp    = coredata(.Object@weather)
              sd_col   = apply(X.imp, 2, function(x) sd(x,na.rm=T))
              rm.vars  = which(sd_col == 0 | is.na(sd_col))
              if (length(rm.vars) == ncol(X.imp)) {
                print('All non-NA covariate values are constant - cannot impute!')
                .Object@weather = zoo(NULL)
              } else {
                if (length(rm.vars)>0) X.imp = X.imp[,-rm.vars]
          #     k        = cv.SVDImpute(as.matrix(X.imp))$k
          #     X.ok.imp = SVDImpute(as.matrix(X.imp), k)    
          #     X.ok.imp = as.data.frame(X.ok.imp$x)
                bds        = matrix(nrow=ncol(X.imp), ncol=3)
                bds[,1]    = 1:ncol(X.imp)
                bds[,2:ncol(bds)] = t(sapply(1:ncol(X.imp), function(i) range(X.imp[,i], na.rm=T)))    
                res   = try(amelia(m=1, x = X.imp, bounds = bds)$imputations[[1]])
                if (class(res) == 'try-error') X.ok.imp = X.imp else X.ok.imp = res
                .Object@weather = zoo(X.ok.imp, order.by = index(.Object@weather))
              }
            return(.Object)
          }
)

# ____________________________________________________
# Auxiliary function for preparing data for modelling

prepareData = function(.Object, trends = NULL) {
  
  # at a minimum, data contains consumption data
  data             = data.frame(consumption = coredata(.Object@consumption))
  time_vec         = .Object@timestamps         
  
  # add time Fixed Effects
  data$Day.Of.Week = as.factor(dayOfWeek(timeDate(time_vec)))
  data$Month       = as.factor(month(timeDate(time_vec)))
  data$Hour.Of.Day = as.factor(as.numeric(sapply( sapply( sapply(strsplit(as.character(time_vec),' '), '[', 2), strsplit, ':' ), '[',1)))
  data$Is.Holiday  = as.factor(as.numeric(isHoliday( timeDate(time_vec) )))
  data$time        = as.numeric(time_vec - time_vec[1]) / 3600
  data.dum         = dummy.data.frame(data=data, sep='.', all=T)
  hourly_vars      = names(data.dum)[which(regexpr('Hour.Of.Day', names(data.dum))>0)]
  monthly_vars     = names(data.dum)[which(regexpr('Month', names(data.dum))>0)]
  weekly_vars      = names(data.dum)[which(regexpr('Day.Of.Week', names(data.dum))>0)]
  data.dum[,hourly_vars[1]] = NULL
  data.dum$Is.Holiday.0     = NULL
  data.dum[,weekly_vars[1]] = NULL
  data.dum[,monthly_vars[1]]= NULL

  # add weather covariates
  wthr_vars_all    = c('TemperatureF', 'WindSpeed', 'Humidity', 'HourlyPrecip', 'SolarRadiation')
  wthr_vars        = intersect(wthr_vars_all, names(.Object@weather))            
  if (length(wthr_vars)>0) data.dum = cbind(data.dum, coredata(.Object@weather[,wthr_vars]))
  
  # add seasonal trends
  if (!is.null(trends)) {
    for (DT in trends) {
      data.dum[,paste('Trend',DT,sep='.')] = sin(2 * pi * data.dum$time / DT)
    }
  }
  
  return( list(data.dum = data.dum, 
               hourly_vars = hourly_vars[-1], 
               monthly_vars = monthly_vars[-1],
               weekly_vars = weekly_vars[-1], 
               holiday_vars = 'Is.Holiday.1',
               trends = trends, 
               weather_vars = wthr_vars))
}
# _________________________________
# Method to perform OLS analysis

library('timeDate')
library('dummies')
setGeneric(
  name = "fitOLS",
  def = function(.Object, verbose = T, covariates = NULL, trends = NULL){standardGeneric("fitOLS")}
)
setMethod('fitOLS',
          signature  = 'Person',
          definition = function(.Object, verbose=T, covariates = NULL, trends = NULL){
            
            if (verbose)
              cat(paste('*** OLS analysis for Person-SPID (', .Object@PER_ID,',', .Object@SP_ID, ')***\n', sep=''))
            
            # ______________________________________________________________
            # Set-up OLS problem (regress consumption on weather + tod/dow)
            
            # indicators for tod/dow
            res       = prepareData(.Object, trends = trends)
            data.dum  = res$data.dum
            
            # set-up model covariates
            if (!is.null(trends)) trend_vars = paste('Trend',trends,sep='.')
            if (is.null(covariates)) {
              predictors = setdiff(names(data.dum), c('consumption', 'time'))
            } else {
              predictors = c(trends, covariates)
              predictors = intersect(predictors, names(data.dum))
            }
            fmla_null = as.formula("consumption ~ 1")
            fmla_full = as.formula(paste('consumption ~ ',
                                         paste(predictors,collapse='+')))
            
            # perform regression            
            OLS.null = lm( fmla_null, data = na.omit(data.dum) )
            OLS.step = step(OLS.null, direction = 'forward', trace = 0, 
                              scope = list(upper = fmla_full, lower = fmla_null))
            fmla_fin = as.formula(paste( 'consumption ~ ', 
                                  paste(names(coef(OLS.step))[-1], collapse = '+')))
            OLS.step = lm( fmla_fin, data = na.omit(data.dum), na.action = na.exclude )
            
            # ______________________________________________________________
            # Investigate structure of residual
                      
            # test for heteroskedasticity using Breusch-Pagan test
            require('lmtest')
            bp.test        = bptest(OLS.step, data = na.omit(data.dum))
            chi2.95        = qchisq(0.95,bp.test$parameter+1)
            heterosc.BP.95 = bp.test$statistic > chi2.95
      
            # test for serial correlation using Durbin-Watson test
            dw.test    = dwtest(OLS.step, data = na.omit(data.dum))
            sercorr.DW = dw.test$alternative == "true autocorrelation is greater than 0"
            
            # which coefficients are statistically significant at least at 0.1 level?
            summ.fit = summary(OLS.step)
            idx.sigf = which(summ.fit$coefficients[,4] <= 0.1)
            coef.sig = summ.fit$coefficients[idx.sigf,]
            
            # save model result
            .Object@OLS = list()
            .Object@OLS[['fit.summary']] = summ.fit
            .Object@OLS[['fitted.vals']] = predict(OLS.step)
            .Object@OLS[['heterosked']]  = heterosc.BP.95
            .Object@OLS[['serial.corr']] = sercorr.DW
            .Object@OLS[['coef.signif']] = coef.sig
            .Object@OLS[['all.covars']]  = predictors
                      
            rm(list=c('OLS.null', 'OLS.step', 'res', 'data.dum'))
            return(.Object)
          }
)

# _________________________________
# Method to perform HMM analysis

library('depmixS4')
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
              cat(paste('*** HMM analysis for Person-SPID (', .Object@PER_ID,',', .Object@SP_ID, ')***\n', sep=''))
            
            # _____________________________
            # Set up depmixS4 model object
          
            included_vars= c(response_vars, transitn_vars)
          
            # have any trends been included in model specification?
            trend_vars = c(included_vars[which(regexpr('Trend', included_vars)>0)])
            if (length(trend_vars)>0) {
              trends = as.numeric(sapply(1:length(trend_vars), 
                                         function(j) strsplit(trend_vars[j], '\\.')[[1]][2]))
            } else {
              trends     = c(10,24)
              trend_vars = paste('Trend', trends, sep='.')
            }
            
            # prepare dataset
            res          = prepareData(.Object, trends = trends)
            data.dum     = res$data.dum
            hourly_vars  = res$hourly_vars
            monthly_vars = res$monthly_vars
            weekly_vars  = res$weekly_vars
            holiday_vars = res$holiday_vars
            wthr_vars    = res$weather_vars
                        
            # see which variables have been included that are reflected in the data as well
            if (length(included_vars)>0) {
              hourly_vars   = intersect(included_vars, hourly_vars)
              monthly_vars  = intersect(included_vars, monthly_vars)
              weekly_vars   = intersect(included_vars, weekly_vars)
              holiday_vars  = intersect(included_vars, holiday_vars)
              wthr_vars     = intersect(included_vars, wthr_vars)
            }
            existing_vars = c(hourly_vars, monthly_vars, weekly_vars, holiday_vars, wthr_vars, trend_vars)
            response_vars = intersect(response_vars, existing_vars)            
            transitn_vars = intersect(transitn_vars, existing_vars)
            
            # do we include just variables that are significant in regression?
            if (ols_vars) {
              ols_vars_names= rownames(.Object@OLS[['coef.signif']])
              if (length(ols_vars_names)>0) 
                response_vars = intersect(response_vars, ols_vars_names)             
            }

            # define response model
            fmla_response = 'consumption ~ 1'
            if (length(response_vars)>0) 
              fmla_response = paste('consumption ~ ',paste(response_vars,collapse='+'))
            fmla_response   = as.formula(fmla_response)
            print(fmla_response)
            
            # define transition model
            fmla_transitn = '~ 1'
            if (length(transitn_vars)>0) 
              fmla_transitn = paste('~ ',paste(transitn_vars,collapse='+'))
            fmla_transitn   = as.formula(fmla_transitn)
            
            # _______________________________________
            # choose model size K (number of states)
            
            # set up and fit model
            minBIC  = Inf
            vec_BIC = c()
            K   = Kmin
            K_opt = Kmin
#             for (K in Kmin:Kmax) {
              
              # initialize state parameters with OLS estimates?
              respstart = NULL
              if (ols_vars) {   
                vars_ok = response_vars
                if ('(Intercept)' %in% ols_vars_names) {
                  vars_ok = c('(Intercept)',vars_ok)
                  ols_vars_coef = .Object@OLS[['coef.signif']][vars_ok,1]
                } else {
                  ols_vars_coef = .Object@OLS[['coef.signif']][vars_ok,1]  
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
              set.seed(1)
              mod <- depmix(response  = fmla_response, 
                            transition= fmla_transitn,
                            data      = na.exclude(data.dum), 
                            nstates   = K, 
                            respstart = respstart,
                            trstart   = trstart,
                            family    = gaussian(),
                            instart   = runif(K))
              fm  <- fit(mod, verbose = verbose, useC = T, 
                         emcontrol = em.control(maxit = 50, tol = 1e-3))
              
              # retain measures of fit
              vec_BIC[as.character(K)] = BIC(fm)
#               if (BIC(fm) < minBIC) {
#                 fm_opt = fm
#                 minBIC = BIC(fm)
#                 K_opt  = K
#               }
#             }

            K_opt  = K
            fm_opt = fm
            
            # ______________________
            # Save computation
            
            .Object@HMM = list()
            
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
            
            # Predicted state means
            means = sapply(1:K_opt, function(k) predict(fm_opt@response[[k]][[1]])) 
            .Object@HMM[['means']] = sapply(1:nrow(means), function(j) means[j,.Object@HMM[['states']][j,1]])

            # other parameters
            .Object@HMM[['nStates']]    = fm_opt@nstates
            .Object@HMM[['BIC']]        = vec_BIC
            
            # ________________________________________
            # Perform preliminary analysis of errors
            
            # which observations were missing (NAs) in the original data?
            idx_na = attributes(na.exclude(data.dum))$na.action
            idx_ok = (1:length(.Object@consumption))[-as.numeric(idx_na)]
                        
            # Confidence p-values
            hmm.sigma = .Object@HMM$response$stdev[.Object@HMM[['states']][,1]]  
            if (length(idx_na)>0) {
              z.scores  = abs(coredata(.Object@consumption)[-idx_na] - .Object@HMM$means) / sqrt(hmm.sigma)
            } else 
              z.scores  = abs(coredata(.Object@consumption) - .Object@HMM$means) / sqrt(hmm.sigma)
            .Object@HMM[['p.values']]   = pnorm(z.scores)
            
            # position of deleted observations
            .Object@HMM[['NA']] = as.numeric(idx_na)
            
            # residuals for each state
            if (length(idx_na) > 0) {
              residuals = coredata(.Object@consumption)[-idx_na] - .Object@HMM[['means']]
            } else 
              residuals = coredata(.Object@consumption) - .Object@HMM[['means']]                          
            .Object@HMM[['residual']] = residuals
            
            # test normality of residuals in each state
            is_normal = sapply(1:.Object@HMM[['nStates']], function(j) {
              idx = which(.Object@HMM[['states']][,1] == j)
              res = shapiro.test(residuals[idx])
              pval= res$p.value
              return(pval)
            })
            .Object@HMM$sw_test = is_normal
                               
            return(.Object)
          }
)
            

# _______________________
# Plot Person object

library('linkcomm')
library('useful')
library('reshape')

setMethod('plot',
          signature  = 'Person',
          definition = function(x, type = 'default', interval = NULL, verbose = T, ...){
            
            if (verbose)
              cat(paste('*** ', type, ' Plot for Person-SPID (', x@PER_ID,',', x@SP_ID, ') ***\n', sep=''))
            
            timestamps = x@timestamps
            ols.resids = x@OLS[['fit.summary']]$residuals
            ols.resids = naresid(x@OLS$fit.summary$na.action, ols.resids)
            consumption= x@consumption
            idx.na     = as.numeric(x@HMM[['NA']])
            if (length(idx.na)>0) {
              consumption = consumption[-idx.na]
              timestamps  = timestamps[-idx.na]
            }
            ols.means  = x@OLS[['fitted.vals']]
            hmm.means  = x@HMM$means
            states     = x@HMM$states[,1]
            hmm.sigma  = x@HMM$response$stdev[states]
            hmm.residual=x@HMM$residual            
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
              plot(timestamps, ols.means, 
                   main = 'OLS Fitted Values', ylab = 'Values', xlab = 'Time', 
                   lwd = 3, type='l', col='red')
              points(timestamps, consumption, type='b', pch=20, col = 'black')
              legend('topleft', c('OLS Fit', 'Observations'), col=c('red', 'black'))
            }    
            
            if (type == 'default') {
              plot(consumption, 
                   main = paste('Consumption Person-SPID (', x@PER_ID,',', x@PER_ID, ')'),
                   xlab = 'Time', ylab = 'kWh', type = 'l', lwd = 2)
            }
            
            if (type == 'weather') {
              tmp = consumption
              if (length(x@weather)>0) tmp = merge(tmp, x@weather)
              plot(tmp, 
                   main = paste('Consumption Person-SPID (', x@PER_ID,',', x@PER_ID, ')'),
                   xlab = 'Time', ylab = 'kWh', type = 'l', lwd = 2)
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
                ggtitle( paste("Zoom-In: States and Observed Emissions (", x@PER_ID, ',', x@SP_ID,')',sep=''))
              
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
                ggtitle(paste("State Space Diagram (", x@PER_ID, ',', x@SP_ID, ')',sep=''))
              return(plt)
            }
            
            # plot confidence interval
            if (type == 'HMM-ci') {
              # compute z-scores of observations
              z.scores = abs(consumption - hmm.means) / sqrt(hmm.sigma)
              p.vals   = pnorm(z.scores)
  
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
                ggtitle( paste("Avg. Prediction Confidence (", x@PER_ID, ',', x@SP_ID, ')',sep=''))
              
              return(plt)
            }
            
            if (type =='HMM-acf') {
              plt = acf_ggplot(coredata(consumption), 
                               hmm.means, 
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
            }
})

# _______________________
# Test constructor
# 
# rm(list=ls())
# source('~/Dropbox/ControlPatterns/code/R/utils/sql_utils.r')
# source('~/Dropbox/ControlPatterns/code/R/utils/timing.r')
# raw_data  = run.query("select * from pge_res_final3_unique WHERE PER_ID = 8420562867 AND SP_ID = 5534067894 ORDER BY date")
# test      = new(Class='Person', raw_data)
# zips      = unique(raw_data$ZIP5)
# query     = paste("SELECT * FROM weather_60 WHERE zip5 IN (", 
#                   paste(zips,collapse=','), ')')
# wthr_data = run.query(query, db = 'PGE_WEATHER')
# test      = addWeather(test, wthr_data, imputate = T)
# 
# # _______________________
# # Test modeling methods
# 
# test      = fitOLS(test)
# test      = fitHMM(test, constrMC = F, Kmin = 2, Kmax = 5)
# 
# # _______________________
# # Test plotting methods
# 
# plot(test, type = 'OLS-fit', interval = c('2010-12-01', '2010-12-10'))
# plot(test, type = 'HMM-ts', interval = c('2010-12-01', '2010-12-20'))
# plot(test, type = 'HMM-MC', interval = c('2010-12-01', '2010-12-20'))
# plot(test, type = 'HMM-ci', interval = c('2010-12-01', '2010-12-10'))
