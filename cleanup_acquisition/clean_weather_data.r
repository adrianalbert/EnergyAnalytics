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

prof_file   = '~/EnergyAnalytics/Rprof_weather.out'
dump_path   = '~/EnergyAnalytics/data/pge_data/weather/'
plots_path  = '~/Dropbox/ControlPatterns/plots/weather_clean/'

# create directories specified in profiles if they don't already exist
if (!file.exists(plots_path)) dir.create(plots_path, recursive = TRUE)
if (!file.exists(dump_path)) dir.create(dump_path, recursive = TRUE)

# _____________________________
# Preliminary setup

# get unique zipcodes
zips = read.csv('~/EnergyAnalytics/data/metadata/zipcode_res.csv')
zips$ZIP_ID = 1:nrow(zips)
write.csv(zips, file = '~/EnergyAnalytics/data/metadata/zipcode_res.csv', row.names = F)

# disconnect all MySQL connections
library('RMySQL')
all_cons <- dbListConnections(MySQL())
for(con in all_cons) dbDisconnect(con)

# delete existing tables of processed weather data
db.con <- dbConnect(dbDriver("MySQL"), 
                    host = "sssl-cyclops.stanford.edu",
                    user = "adrian", password = "xmenxmen", 
                    dbname = 'PGE_WEATHER')
for (zip in zips$ZIP5) dbRemoveTable(db.con, paste('ZIP', zip, sep='_'))
dbDisconnect(db.con)

# get data from database...
query         = paste("SELECT * FROM weather_60") 
wthr_data_all = run.query(query, db = 'PGE_WEATHER')

# _______________________________________________________
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

# ________________________________________________________
# Function to perform analysis and do plots of imputation

perform_analysis = function(orig, clean, timestamps) {
    
  # plots of densities for each covariate (before/after)
  covars = names(clean)
  orig   = orig[,covars]
  orig$Type  = 'orig'
  clean$Type = 'clean'
  df = rbind(melt(orig, id.vars = c('Type')), melt(clean, id.vars = c('Type')))
  plt1 = ggplot(df, aes(x = value, color = Type)) + 
    facet_wrap(~variable, nrow=2, scales = 'free') + 
    geom_density(size=2) + 
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.text.x      = element_text(size=18),                  
          axis.text.y      = element_text(size=18),  
          strip.text.x     = element_text(size=18),                                   
          legend.text      = element_text(size=16),                         
          plot.title       = element_text(size=18),          
          axis.title.x = element_blank(),
          panel.background = element_blank(),
          axis.ticks = element_blank()) + 
    ggtitle('Density: Raw and Cleaned')
  
  # table of # of points imputed per covariate
  nr.NA = sapply(covars, function(v) length(which(is.na(orig[,v]))))
  levels(df$variable) = paste(levels(df$variable), nr.NA[levels(df$variable)])
  
  # time series plot: original and estimate
  
  orig$date = timestamps
  clean$date= timestamps
  df = rbind(melt(orig, id.vars = c('Type', 'date')), melt(clean, id.vars = c('Type', 'date')))
  plt2 = ggplot(df, aes(x = date, y = value, color = Type)) + 
    facet_wrap(~variable, ncol = 1, scales = 'free') + 
    geom_line() + 
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.text.x      = element_text(size=18),                  
          axis.text.y      = element_text(size=18),                         
          legend.text      = element_text(size=16),                         
          strip.text.x     = element_text(size=18),                         
          plot.title       = element_text(size=18),
          axis.title.x     = element_blank(),
          panel.background = element_blank(),
          axis.ticks = element_blank()) + 
   ggtitle(paste('Weather: Actual and Imputed'))
  
  return(list(nr.NA/nrow(orig), plt1, plt2))
}

# __________________________________________________
# Function to clean up data for a given zipcode

source('~/EnergyAnalytics/code/R/utils/gbmImpute.r')
library('imputation')
library('dummies')
library('timeDate')
library('lubridate')
clean_weather = function(zip) {
  
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
  
  # clean up weather covariates (remove unreasonable values)
  ranges = list(TemperatureF = c(2, 125), 
                Pressure     = c(25, 35),
                Humidity     = c(5, 100))
  for (r in intersect(names(ranges),names(wthr_data))){
    idx = which(wthr_data[,r]<ranges[[r]][1] | wthr_data[,r]>ranges[[r]][2]) 
    if (length(idx)>0) wthr_data[idx,r] = NA
  }

  # remove covariates that have less than 1/2 of the values
  na.number = sapply(names(wthr_data), function(v) {
    return(length(which(is.na(wthr_data[,v]))))
  })
  idx = which(na.number < nrow(wthr_data) / 2)                   
  if (length(idx) > 0) wthr_data = wthr_data[,-idx]              
  
  # temporal indicators
  Day.Of.Week = as.factor(dayOfWeek(timeDate(wthr_times)))
  Month       = as.factor(month(timeDate(wthr_times)))
  Hour.Of.Day = as.factor(as.numeric(sapply( sapply( sapply(strsplit(as.character(wthr_times),' '), '[', 2), strsplit, ':' ), '[',1)))
  Is.Holiday  = as.factor(as.numeric(isHoliday( timeDate(wthr_times) )))
    
  # imputate missing values (at least one covariate value in an observation tuple needs to be defined)
  wthr_clean = imputateMissingValues(wthr_data)
  wthr_covar = names(wthr_clean)
  wthr_clean = cbind(data.frame(Month = Month, HourOfDay = Hour.Of.Day), wthr_clean)

  # use gbm to estimate gaps (no covariate value available for a given observation)
  X.dum = dummy.data.frame(wthr_clean)
  X.imp = gbmImpute(X.dum, max.iters = 2, cv.fold = 4, verbose = F)$x 
  wthr_clean = X.imp[,wthr_covar]
  
  # plots 
  res = perform_analysis(wthr_data, wthr_clean, wthr_posix)
  png(paste(plots_path, paste(zip, 'density.png', sep='_'), sep=''), width=1200, height = 800)
  print(res[[2]])
  dev.off()
  png(paste(plots_path, paste('zip', 'ts.png', sep='_'), sep=''), width=2000, height = 1200)
  print(res[[3]])
  dev.off()
  
  # add in dummies
  wthr_clean = cbind(data.frame(date = wthr_times, Month = Month, HourOfDay = Hour.Of.Day, 
                                IsHoliday = Is.Holiday, DayOfWeek = Day.Of.Week), 
                     wthr_clean)
  
  # open up connection to database
  db.con <- dbConnect(dbDriver("MySQL"), 
                      host = "sssl-cyclops.stanford.edu",
                      user = "adrian", password = "xmenxmen", 
                      dbname = 'PGE_WEATHER')

  # write cleand-up data back to database
  dbWriteTable(db.con, paste("ZIP", zip, sep='_'), wthr_clean, overwrite = T)
  
  # disconnect from database
  dbDisconnect(db.con)

  # dump data in a file
  write.csv(file = paste(dump_path, paste(zip, 'csv', sep='.'), sep=''), wthr_clean, row.names = F)  
  
  df = res[[1]][wthr_names]  
  return(df)
}

# __________________________________________________
# Run through all zipcodes

wrapper = function(zip){
  zip = zips$ZIP5[zip_id]
  cat(paste('Processing Zipcode', zip, '(', zip_id, '/', nrow(zips), ')', '\n'))
  
  res = try(clean_weather(zip))
  if (class(res) == 'try-error') {
    cat(paste('Error processing zipcode', zip, '\n'))
  }
  
  return(res)
}

stats = data.frame()
for (zip_id in 21:30) {
  res = wrapper(zip_id)  
  if (nrow(stats)>0) stats = res else stats = rbind(stats, res)
}

# parallel execution using parallel package
library(parallel)

cat(paste('Parallel computation on', detectCores(), 'cores\n'))

ptm <- proc.time()
Rprof(filename = paste(prof_file, sep=''), 
      interval = 0.02, memory.profiling = T)
result = mclapply(30:60, FUN = wrapper, mc.silent = F, 
                   mc.preschedule = FALSE, mc.cores = 8)
Rprof(NULL)
dt = proc.time() - ptm

profiled = summaryRprof(filename=prof_file, memory = 'none')
