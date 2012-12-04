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

# class definition
setClass(
  Class = "Person",
  representation = representation(
    PER_ID     = "numeric",           # unique person ID
    SP_ID      = "numeric",           # vector of meters/premises associated with this person
    ZIPCODE    = "numeric",           # vector of zipcodes associated with this person
    TIMES_MOVED= "numeric",           # how many times they moved    
    pers_attr  = "data.frame",        # Person characteristics (fixed over time)
    timestamps = "character",         # vector of timestamps
    prem_attr  = "data.frame",        # Premise characteristics (can change in time)  
    weather    = "data.frame",        # weather covariates
    consumption= "data.frame"         # hourly consumption 
    )
)

# define validation function for class Person


# define method to populate attributes from R data frame
setMethod(f = "initialize", 
          signature = "Person",
          definition = function(.Object, raw_data, res = 60) {
            
            .Object@PER_ID = unique(raw_data$PER_ID)
            
            # _______________________
            # Format attributes data
            
            # columns referring to person attributes
            pers_attr_names = c('CRSCODE', 'ENDUSE', 'CARE', 'STATUS', 'COMMTYPE', 
                                'FERA', 'CLIMSMRT', 'ACCTTYPE', 'DRPROG', 'CEEPROG', 
                                'max_total_duration',                   
                                'renter', 'sfo', 'owned_premises', 'PRIORITY', 'ua_uc_ru', 'median_income',
                                'median_income_quantiles', 'ro', 'vro', 
                                'good_sample_period', 'vro_total_good_sample_period', 'num_prems_for_vro') 
            # columns referring to premise attributes
            prem_attr_names = c('SA_TYPE', 'RSCHED', 'PSA_ID', 'SERCITY', 'GCOUNTY', 'DIVOFF', 
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
            
            # ______________________________________________
            # Format consumption data to timeseries format

            # columns referring to consumption
            cons_names  = paste('hkw', 1:(24*60/res), sep='')
            kwhMat      = raw_data[,cons_names]
            consumption = as.vector(t(kwhMat))
            days        = as.POSIXct(raw_data$date,tz="PST8PDT", '%Y-%m-%d')
            # create a row of hourly values for each day
            # sapply returns an array of numeric epoch seconds (for origin '1970-01-01')
            dateMat     = sapply(days,FUN=function(x) x + (0:(24-1) * 3600))            
            # flatten into a vector and re-convert into date objects
            timestamps  = as.POSIXlt(as.vector(dateMat),origin='1970-01-01')
            
            # _____________________________
            # Check for object consistency 
            
            # _________________________________________
            # Compute summary information of interest
            
            return(.Object)
          })

test = new(Class='Person', raw_data)

# Method to add in exogenous covariates (weather, billing)
setGeneric(
  name = "add.weather",
  def = function(object){standardGeneric("add.weather")}
)
setMethod('add.weather',
          signature  = 'Person',
          definition = function(.Object, wthr_data, verbose=T){
            if (verbose) cat(paste('*** Adding weather data to person', .Object@PER_ID, '***'))
            
          }
)

ResDataClass = function(sp_id,zip=NULL,w=NULL,db='pge_res'){
  query = paste(
    'SELECT 
      zip5,DATE,
    hkw1, hkw2, hkw3, hkw4, hkw5, hkw6, hkw7, hkw8, hkw9, hkw10,hkw11,hkw12,
    hkw13,hkw14,hkw15,hkw16,hkw17,hkw18,hkw19,hkw20,hkw21,hkw22,hkw23,hkw24 
    FROM',conf.meterTable(zip),'WHERE sp_id =',sp_id,'ORDER BY DATE')
  raw = run.query(query,conf.meterDB())
  
  if(length(raw)==0) stop(paste('No data found for sp_id',sp_id))
  
  zipcode = raw[1,1]
  kwMat = raw[,3:26]
  # reshape the kW readings into a vector matching the dates
  kw    = as.vector(t(kwMat))
  
  days = as.POSIXct(raw[,2],tz="PST8PDT", '%Y-%m-%d')
  # create a row of hourly values for each day
  daySteps = 24
  dtDay = daySteps/24 * 60 * 60 # in seconds
  # sapply returns an array of numeric epoch seconds (for origin '1970-01-01')
  dateMat = sapply(days,FUN=function(x) x + (0:(daySteps-1) * dtDay))
  
  # flatten into a vector and re-convert into date objects
  dates = as.POSIXlt(as.vector(dateMat),origin='1970-01-01')
  
  if (is.null(w)) w = WeatherClass(zipcode) # todo: pass in the dates to interpolate ,dates)
  
  # TODO: clear out obviously bad readings
  #keepers   = which(kw > 0)
  
  obj = list (
    dates = dates,
    kw  = kw,
    kwMat = kwMat,
    days = days,
    zipcode = zipcode,
    weather = w,
    tout = w$resample(dates),
    get = function(x) obj[[x]],
    # Not sure why <<- is used here
    # <<- searches parent environments before assignment
    # http://stat.ethz.ch/R-manual/R-patched/library/base/html/assignOps.html
    set = function(x, value) obj[[x]] <<- value,
    props = list()
  )
  
  obj$norm = function(data) {
    # divide data by the 97th %ile
    return (data/quantile(data,0.97,na.rm=TRUE))
  }
  
  obj$add = function(name, value) {
    obj[[name]] = value
  }
  
  # note how list manipulation requires the use of assign
  # not sure why values can't be set in place, but it 
  # appears to have to do with variable scoping
  obj$addProp = function(name, value) {
    p <- obj$props
    p[[name]] <- value
    assign('props', p, envir=obj)
  }
  
  obj$matchDates = function(newDates) {
    a = approx(obj$dates, obj$tout, newDates, method="linear" )[[2]]
    return(a)
  }
  
  #obj <- list2env(obj)
  class(obj) = "ResDataClass"
  return(obj)
}

# Weather class: data & preprocessing for a zipcode
# -------------------------------------------------

