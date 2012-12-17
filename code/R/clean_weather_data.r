# clean_weather_data.r
#
# Perform imputation on weather data
#
# Adrian Albert
#
# Last modified: December 2012.

rm(list=ls())

options(error = recover)

setwd('~/EnergyAnalytics/code/R')
source('~/EnergyAnalytics/code/R/utils/timing.r')
source('~/EnergyAnalytics/code/R/utils/sql_utils.r')

info_file   = '~/EnergyAnalytics/data/per_sp_id_info.RData'
prof_file   = '~/EnergyAnalytics/Rprof_weather.out'
dump_path   = '~/EnergyAnalytics/data/pge_data/weather/'

# _____________________________
# Preliminary setup

# disconnect all MySQL connections
library('RMySQL')
all_cons <- dbListConnections(MySQL())
for(con in all_cons) dbDisconnect(con)

# get unique zipcodes
zips = read.csv('~/EnergyAnalytics/data/metadata/zipcode_res.csv')
zips$ZIP_ID = 1:nrow(zips)
write.csv(zips, file = '~/EnergyAnalytics/data/metadata/zipcode_res.csv', row.names = F)

# get data from database...
query         = paste("SELECT * FROM weather_60") 
wthr_data_all = run.query(query, db = 'PGE_WEATHER')

# ___________________________________________
# Function to imputate missing observations using EM/SVD

library('Amelia')
library('imputation')
library('zoo')
imputateMissingValues = function(X.imp, verbose=T){
  sd_col   = apply(X.imp, 2, function(x) sd(x,na.rm=T))
  rm.vars  = which(sd_col == 0 | is.na(sd_col))
  if (length(rm.vars) == ncol(X.imp)) {
    print('All non-NA covariate values are constant - cannot impute!')
    X.ok.imp = X.imp
  } else {
    if (length(rm.vars)>0) X.imp = X.imp[,-rm.vars]
    bds        = matrix(nrow=ncol(X.imp), ncol=3)
    bds[,1]    = 1:ncol(X.imp)
    bds[,2:ncol(bds)] = t(sapply(1:ncol(X.imp), function(i) range(X.imp[,i], na.rm=T)))    
    res   = try(amelia(m=1, x = X.imp, bounds = bds)$imputations[[1]])
    if (class(res) == 'try-error') X.ok.imp = X.imp else X.ok.imp = res
  }
  return(X.ok.imp)
}

# __________________________________________________
# Function to clean up data for a given zipcode

library('timeDate')
library('lubridate')
clean_weather = function(zip, path = dump_path) {
  
  wthr_data = subset(wthr_data_all, zip5 == zip)
  
  # select covariates of interest
  wthr_times  = as.character(wthr_data$date)
  wthr_posix  = as.POSIXct(wthr_times)
  wthr_names  = c('TemperatureF', 'DewpointF', 'Pressure', 'WindSpeed', 'Humidity', 
                  'HourlyPrecip', 'SolarRadiation')
  wthr_data   = wthr_data[,wthr_names]            
    
  # remove duplicate timestamps in weather data, if any
  idx_dup     = which(duplicated(wthr_times))
  if (length(idx_dup) > 0) {
    wthr_times = wthr_times[-idx_dup]
    wthr_data  = wthr_data[-idx_dup,]
  }
  
  # temporal indicators
  Day.Of.Week = as.factor(dayOfWeek(timeDate(wthr_times)))
  Month       = as.factor(month(timeDate(wthr_times)))
  Hour.Of.Day = as.factor(as.numeric(sapply( sapply( sapply(strsplit(as.character(wthr_times),' '), '[', 2), strsplit, ':' ), '[',1)))
  Is.Holiday  = as.factor(as.numeric(isHoliday( timeDate(wthr_times) )))
    
  # imputate missing values  
  wthr_clean = imputateMissingValues(wthr_data)
  wthr_clean = cbind(data.frame(date = wthr_times, Month = Month, DayOfWeek = Day.Of.Week, 
                                HourOfDay = Hour.Of.Day, IsHoliday = Is.Holiday), 
                     wthr_clean)

  # open up connection to database
  db.con <- dbConnect(dbDriver("MySQL"), 
                      host = "sssl-cyclops.stanford.edu",
                      user = "adrian", password = "gambit", 
                      dbname = 'PGE_WEATHER')

  # write cleand-up data back to database
  dbWriteTable(db.con, paste("ZIP", zip, sep='_'), wthr_clean, overwrite = T)
  
  # disconnect from database
  dbDisconnect(db.con)

  # dump data in a file
  write.csv(file = paste(path, paste(zip, 'csv', sep='.'), sep=''), wthr_clean, row.names = F)  
  
  return(NULL)
}

# __________________________________________________
# Run through all zipcodes

wrapper = function(zip_id){
  zip = zips$ZIP5[zip_id]
  cat(paste('Processing Zipcode', zip, '(', zip_id, '/', nrow(zips), ')', '\n'))
  
  res = try(clean_weather(zip))
  if (class(res) == 'try-error') {
    cat(paste('Error processing zipcode', zip, '\n'))
  }
}

# 
# for (zip_id in 1:nrow(zips)) {
#   wrapper(zip_id)  
# }

# parallel execution using parallel package
library(parallel)

cat(paste('Parallel computation on', detectCores(), 'cores\n'))

ptm <- proc.time()
Rprof(filename = paste(prof_file, sep=''), 
      interval = 0.02, memory.profiling = T)
result = mclapply(1:nrow(zips), FUN = wrapper, mc.silent = F, 
                   mc.preschedule = FALSE, mc.cores = 8)
Rprof(NULL)
dt = proc.time() - ptm

profiled = summaryRprof(filename=prof_file, memory = 'none')
