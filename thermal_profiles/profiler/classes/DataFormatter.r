# DataFormatter.r
#
# Creates a object that sets up the data needed for analysis.
# 
# Adrian Albert
# Last modified: October 2013.
# -----------------------------------------------------------------------

library('methods')
library('timeDate')
library('zoo')
library('R.utils')
library('lubridate')

# clean-up previous definitions of methods for class DataFormatter
removeClass('DataFormatter')

options(error = recover)

# ________________________
# Class definition
setClass(
  Class = "DataFormatter",
  representation = representation(
    UID        = "character",         # unique time series ID
    timestamps = "character",         # vector of timestamps
    covariates = "data.frame",        # covariates
    obs        = "numeric"           # observations 
    )
)

# _______________________________________
# Constructor method for class DataFormatter

setMethod(f = "initialize", 
          signature = "DataFormatter",
          definition = function(.Object, raw_data, UID, verbose = T) {
            
            .Object@UID     = as.character(UID)
                        
            if (verbose) {
              cat(paste('*** Initializing DataFormatter (', .Object@UID,') ***\n', sep=''))
              t0 = proc.time()
            }
                        
            # ______________________________________________
            # Format observation data to timeseries format

            if (ncol(raw_data)>2) {  # is it wide format? (sequence of days)
              # columns referring to kWh
              days        = raw_data$date       
              cons_names  = paste('obs', 1:24, sep='')
              names(raw_data) = c('date', cons_names)              
              hours       = paste(rep(formatC(0:23, flag=0, width=2), length(days)), ':00:00', sep='')
              days_times  = rep(days, each = 24)
              days_times  = paste(days_times, hours, sep=' ')
              obs         = as.vector(t(raw_data[,cons_names]))
            } else {                 # if not, it is long format
              obs = raw_data$obs
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
      	    if (length(na.omit(obs)) == 0) {
      	      stop('Error: observations data is all NAs!')
      	    }

            # if data is all zeros
            rc = range(na.omit(obs))
            if (rc[2] - rc[1] == 0) {
              stop('Error: data supplied to constructor has no variation!')
            }

            .Object@obs = obs
                        
            # initialize default weather object
            .Object@covariates  = data.frame()
            
            if (verbose) {
              dt = proc.time() - t0;
              print(dt)
            }
            
            return(.Object)
          })


# ____________________________________________________________
# Method to add in exogenous covariates (weather, billing)

setGeneric(
  name = "addCovariates",
  def = function(.Object, covar_times, covar_data, verbose = T){standardGeneric("addCovariates")}
)
setMethod('addCovariates',
          signature  = 'DataFormatter',
          definition = function(.Object, covar_times, covar_data, verbose=T){
            if (verbose) {
              cat(paste('*** Adding covariates data for DataFormatter', .Object@UID,' ***\n',sep=''))
              t0 = proc.time()
            }
            
            covar_names = names(covar_data)
            
            # remove duplicate timestamps in weather data, if any
            idx_dup     = which(duplicated(covar_times))
            if (length(idx_dup) > 0) {
              covar_times= covar_times[-idx_dup]
              covar_data = covar_data[-idx_dup,]
            }
            
            # remove NAs, if any
            idx.na = which(!complete.cases(covar_data)) 
            if (length(idx.na)>0){
              covar_data = covar_data[-idx.na,]
              covar_times= covar_times[-idx.na]
            }
            
            # is there weather data at times that match kWh data?            
            idx_ok      = which(covar_times %in% .Object@timestamps)            
            if (length(idx_ok) > 0) {
              covar_data = covar_data[idx_ok,]
              covar_times= covar_times[idx_ok]
            } else {
              stop(paste('!!! No covariate data available for user', .Object@UID))
            }
            idx_ok_cons  = which(.Object@timestamps %in% covar_times)            
                        
            if (class(covar_data) == 'numeric') {
              covar_data = data.frame(covar_data)
              names(covar_data) = covar_names
            }
            # update object
            .Object@covariates  = covar_data
            .Object@timestamps  = as.character(covar_times)
            .Object@obs         = .Object@obs[idx_ok_cons]

            if (verbose) {
              dt = proc.time() - t0;
              print(dt)
            }
            
            return(.Object)
          }
)


# ____________________________________________________
# Method to extract data to pass on to other modules

setGeneric(
  name = "extractFormattedData",
  def = function(.Object){standardGeneric("extractFormattedData")}
)
setMethod('extractFormattedData',
          signature  = 'DataFormatter',
          definition = function(.Object) {
            
          # add kWh and weather, if any
          data          = data.frame(obs = .Object@obs) 
          if (length(.Object@covariates)>0)
            data        = cbind(data, .Object@covariates)
          
          covr_times     = as.POSIXlt(.Object@timestamps)              
          data$DayOfWeek = as.factor(substr(weekdays(covr_times), 1, 3))
          data$Month     = as.factor(covr_times$mon)
          data$HourOfDay = as.factor(as.numeric(sapply( sapply( sapply(strsplit(as.character(covr_times),' '), '[', 2), strsplit, ':' ), '[',1)))
          data$Weekend   = as.factor(as.numeric(isHoliday( timeDate(covr_times) )))                          
                                            
          
          return(list(UID        = .Object@UID,
                      data       = data,
                      timestamps = .Object@timestamps))     
          })      

# ____________________________________________________________
# Print method for class DataFormatter: display useful stats.

setMethod('show',
          signature  = 'DataFormatter',
          definition = function(object){
            cat('*** DataFormatter Object ***\n')
            
            # basic info
            cat(sprintf('UID:         %s\n', object@UID))            
            cat(sprintf('Time range:  %s - %s\n', object@timestamps[1], object@timestamps[length(object@timestamps)]))   
            cat(sprintf('Data amount: %s observations\n', length(object@timestamps)))   
            
            # weather info
            cat(paste('Covariates:  ', paste(names(object@covariates), collapse=', '), '\n'))            
            cat('*** END: DataFormatter Object ***\n')
          })

