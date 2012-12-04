# Person.r
#
# The main purposes of  this class is to encapsulate pre-processing and model fitting functions where
# the unit of analysis is a person. It should have slots for 
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

# Person class: data & analysis for a unique PER_ID
# -------------------------------------------------

library('methods')
removeClass('Person')

# ________________________
# Class definition

setClass(
  Class = "Person",
  representation = representation(
    PER_ID     = "numeric",           # unique person ID
    N_OBS      = "numeric",           # number of data points
    N_DAYS     = "numeric",           # number of days (=N_OBS * 24)
    PREMISE    = "numeric",           # vector of meters/premises associated with this person
    ZIPCODE    = "character",         # vector of zipcodes associated with this person
    TIMES_MOVED= "numeric",           # how many times they moved  
    pers_attr  = "data.frame",        # Person characteristics (fixed over time)
    timestamps = "POSIXlt",           # vector of timestamps
    prem_attr  = "data.frame",        # Premise characteristics (can change in time)  
    weather    = "data.frame",        # weather covariates
    consumption= "numeric"         # hourly consumption 
    )
)

# _______________________________________
# Constructor method for class Person

setMethod(f = "initialize", 
          signature = "Person",
          definition = function(.Object, raw_data, res = 60, verbose = T) {
            
            .Object@PER_ID = unique(raw_data$PER_ID)
            
            if (verbose) cat(paste('*** Initializing Person', .Object@PER_ID, '***\n'))
                        
            # ______________________________________________
            # Format consumption data to timeseries format

            # columns referring to consumption
            days        = as.POSIXct(raw_data$date,tz="PST8PDT", '%Y-%m-%d')
            # create a row of hourly values for each day
            # sapply returns an array of numeric epoch seconds (for origin '1970-01-01')
            dateMat     = sapply(days,FUN=function(x) x + (0:(24-1) * 3600))            
            # flatten into a vector and re-convert into date objects
            .Object@timestamps  = as.POSIXlt(as.vector(dateMat),origin='1970-01-01')
            idx_order   = order(.Object@timestamps)
            .Object@timestamps  = .Object@timestamps[idx_order]
            cons_names  = paste('hkw', 1:(24*60/res), sep='')
            kwhMat      = raw_data[idx_order,cons_names]
            .Object@consumption = as.vector(t(kwhMat))
            .Object@consumption = .Object@consumption[idx_order]
            
            # _______________________
            # Format attributes data
            
            # columns referring to person attributes (fixed in time)
            pers_attr_names = c('CRSCODE', 'ENDUSE', 'CARE', 'STATUS', 'COMMTYPE', 
                                'FERA', 'CLIMSMRT', 'ACCTTYPE', 'DRPROG', 'CEEPROG', 
                                'max_total_duration',                   
                                'renter', 'sfo', 'owned_premises', 'PRIORITY', 'ua_uc_ru', 'median_income',
                                'median_income_quantiles', 'ro', 'vro', 
                                'good_sample_period', 'vro_total_good_sample_period', 'num_prems_for_vro') 
            # columns referring to premise attributes (can change with time)
            prem_attr_names = c('SP_ID', 'SA_TYPE', 'RSCHED', 'PSA_ID', 'SERCITY', 'GCOUNTY', 'DIVOFF', 
                                'CLIMATE', 'PREMTYPE', 'SM_SPST', 'NETMETER', 'MSTRMTR',
                                'FIRSTDTE_DATE', 'LASTDTE_DATE', 'earliest_start', 'latest_end', 
                                'SA_START_DATE', 'SA_END_DATE', 'DTONPREM_DATE', 'SA_DURATION',                 
                                'SM_START_DATE', 'SM_END_DATE', 'SM_DURATION',
                                'ZIP5', 'AREA', 'WTHRSTN', 'CECCLMZN', 'DTONPREM_DATE', 'total_duration', 
                                'num_good_people_at_prem')
            # get unique values for attributes          
            # should be only one value per person!
            attr_levels = sapply(pers_attr_names, function(c) { length(unique(raw_data[,c])) })
            mult_levels = which(attr_levels > 1)
            if (length(mult_levels)>0) {
              prem_attr_names = c(prem_attr_names, pers_attr_names[mult_levels])
              pers_attr_names = pers_attr_names[-mult_levels]
            }            
            # set person-specific attributes
            .Object@pers_attr = as.data.frame(t(sapply(raw_data[,pers_attr_names], unique)))
            # set premise-specific attributes
            .Object@prem_attr = raw_data[,prem_attr_names]
            
            # _____________________________
            # Check for object consistency 
            
            validObject(.Object)
            
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
            if (verbose) cat(paste('*** Computing stats for person', .Object@PER_ID, '***\n'))
            .Object@PREMISE = unique(.Object@prem_attr[,'SP_ID'])
            .Object@ZIPCODE = as.character(unique(.Object@prem_attr[,'ZIP5']))
            .Object@TIMES_MOVED = length(.Object@PREMISE) 
            .Object@N_OBS   = length(.Object@consumption) 
            .Object@N_DAYS  = .Object@N_OBS / 24
            return(.Object)
          }
)

# ____________________________________________________________
# Method to add in exogenous covariates (weather, billing)

setGeneric(
  name = "addWeather",
  def = function(.Object, wthr_data, verbose = T){standardGeneric("addWeather")}
)
setMethod('addWeather',
          signature  = 'Person',
          definition = function(.Object, wthr_data, verbose=T){
            if (verbose) cat(paste('*** Adding weather data to person', .Object@PER_ID, '***\n'))
            # match dates for weather and consumption (assume both are ordered by date/time)
            start_time= .Object@timestamps[1]
            stop_time = .Object@timestamps[.Object@N_OBS] 
            wthr_data$date = as.POSIXlt(wthr_data$date)
            idx_ok    = which(wthr_data$date >= start_time &
                              wthr_data$date <= stop_time &
                              wthr_data$zip5 %in% .Object@ZIPCODE)
            if (length(idx_ok) == 0) 
              stop(paste('No weather data for person', .Object@PER_ID))
            wthr_vars = c('TemperatureF', 'DewpointF', 'Pressure', 'WindSpeed', 'Humidity', 
                          'Clouds', 'HourlyPrecip', 'SolarRadiation')
            .Object@weather = wthr_data[idx_ok, wthr_vars]
            
            return(.Object)
          }
)

# ____________________________________________________________
# Print method for class Person: display useful stats.

setMethod('show',
          signature  = 'Person',
          definition = function(object){
            cat('*** Person Object ***\n')
            cat(paste('ID', object@PER_ID, '\n'))            
            cat(sprintf('# Days = %d; # Observations = %d\n', object@N_DAYS, object@N_OBS))
            cat(sprintf('# Premises changes = %d; # Zipcodes changed = %d\n', 
                        object@TIMES_MOVED, length(object@ZIPCODE)))
            cat(sprintf('# Time range: %s - %s\n', object@timestamps[1], object@timestamps[object@N_OBS]))            
            cat('*** END: Person Object ***\n')
          })

# _________________________________
# Method to perform OLS analysis

setGeneric(
  name = "analyzeOLS",
  def = function(.Object,verbose = T){standardGeneric("analyzeOLS")}
)
setMethod('analyzeOLS',
          signature  = 'Person',
          definition = function(.Object, verbose=T){
            if (verbose) cat(paste('*** Performing OLS analysis for Person', .Object@PER_ID, '***\n'))
            
            # set-up OLS problem (regress consumption on weather + tod/dow)
            
            # investigate 
            
            return(.Object)
          }
)



# _______________________
# Test constructor

# rm(list=ls())
# source('~/Dropbox/ControlPatterns/code/R/utils/sql_utils.r')
# source('~/Dropbox/ControlPatterns/code/R/utils/timing.r')
# raw_data  = run.query("select * from pge_res_final3_unique WHERE PER_ID = 8420562867")
# large_owners  = run.query("select * from pge_res_final3_unique WHERE vro = 1")
# rich_owners   = run.query("select * from pge_res_final3_unique WHERE ro = 1")
# test      = new(Class='Person', raw_data)
# zips      = unique(raw_data$ZIP5)
# query     = paste("SELECT * FROM weather_60 WHERE zip5 IN (", 
#                   paste(zips,collapse=','), ')')
# wthr_data = run.query(query, db = 'PGE_WEATHER')
# test      = addWeather(test, wthr_data)


