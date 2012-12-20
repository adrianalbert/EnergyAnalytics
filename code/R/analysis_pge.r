#!/usr/bin/Rscript

# analysis_pge.r
# 
# Typical OLS analysis for PGE data:
# - regression residuals: argue non-Gaussian
# - ACF: argue non-iid
# - different models for different time-of-day/holiday/week-weekend: argue different states
# - two users: argue different responses
#
# HMM analysis (w/ covariates):
# - 
# Adrian Albert
# Last modified: December 2012.

# -----------------------------------------------
# Define constants, load libraries and functions
# -----------------------------------------------

rm(list=ls())

options(error = recover)

setwd('~/EnergyAnalytics/code/R')
source('~/EnergyAnalytics/code/R/utils/timing.r')
source('~/EnergyAnalytics/code/R/utils/sql_utils.r')
source('~/EnergyAnalytics/code/R/utils/acf_ggplot.r')
source('~/EnergyAnalytics/code/R/Person.r')
source('~/EnergyAnalytics/code/R/personAnalysis.r')

# default options for this run
NOBS_THRESH = 180   # discard users with less than a 1/2 year at a given premise
N_CORES     = 1     # number of cores to run process on 
DATA_ACCESS = 'sql' # flag to indicate type of media used for data
PLOTS       = FALSE
VERBOSE     = TRUE
USER        = 'adalbert'
PASSWORD    = 'adrian'
ZIP_START   = 2
ZIP_STOP    = 2
PROFILE     = F

plots_path  = '~/Dropbox/ControlPatterns/plots/'
save_path   = '~/EnergyAnalytics/fits/'
kwh_path    = ''
wthr_path   = ''
prof_file   = '~/EnergyAnalytics/Rprof.out'
zips_file   = '~/EnergyAnalytics/data/metadata/zipcode_res.csv'

# if in batch mode, parse command-line arguments
library(utils)
myargs = commandArgs(trailingOnly = TRUE)
cat(paste('Arguments:', paste(myargs,sep='; '), '\n'))

if (length(myargs)==1) {
  source(myargs)
}

if (length(myargs)>=2) {
  ZIP_START = as.numeric(myargs[1])
  ZIP_STOP = as.numeric(myargs[2])
  if (length(myargs)>=3) nProc  = as.numeric(myargs[3]) else nProc = detectCores()
  if (length(myargs)==5) {
    user     = myargs[4]
    password = myargs[5]
    cat(paste('Credentials:', user, password), '\n')
  }
} 

# ------------------------------------------------
# Load identification information on users
# ------------------------------------------------

# disconnect all MySQL connections
library('RMySQL')
#all_cons <- dbListConnections(MySQL())
#for(con in all_cons) dbDisconnect(con)

# get zipcodes
zips = read.csv(zips_file)
zips = zips$ZIP5

# ------------------------------------------
# Wrapper to send computation to each core
# ------------------------------------------

analysis_wrapper <- function(iZip, zips, type = 'sql') {
  
  sink(paste(save_path, "log_", iZip,".txt", sep=''), append=F, type = c('output', 'message'))
  
  cat(paste('--------- Processing Zipcode', iZip, '/', length(zips), '---------\n'))
  
  if (type == 'sql') {
    # retrieve current zip consumption data from DB
    query     = paste("SELECT * FROM pge_res_60_", zips[iZip], ' ORDER BY date', sep='')
    raw_data  = run.query(query, db = 'PGE_SAM', user = USER, password = PASSWORD)  
    
    # retrieve current zip weather data from DB
    query     = paste("SELECT * FROM ZIP_", zips[iZip], ' ORDER BY date', sep='')
    cur_wthr  = run.query(query, db = 'PGE_WEATHER', user = USER, password = PASSWORD)  
    
    # check data size in memory
    cat(paste('Weather data size for zip', zips[iZip], '=', object.size(cur_wthr) / 1024^2, '\n'))
    cat(paste('Data size for zip', zips[iZip], '=', object.size(raw_data)/1024^2, '\n'))

  } else {
    
    # TODO: read data from files!
  }
  
  # perform calculations on each user in the dataset
  response.df = data.frame()
  transitn.df = data.frame()
  ols.coefs   = data.frame()
  iUser       = 0
  all_users   = unique(raw_data$UID)
  for (usr in all_users) {     
    iUser = iUser + 1

    # _______________________________________
    # Perform analysis flow for current user
    
    cur_data  = subset(raw_data, UID == usr)
    if (nrow(cur_data) < NOBS_THRESH) {
      cat(paste('---> Too few datapoints! (', nrow(cur_data), ')\n', sep=''))
      next
    }
    
    # perform analysis flow
    ptm <- proc.time()
    res = try(personAnalysis(cur_data, cur_wthr, usr, zips[iZip],
                             plots = PLOTS, verbose = VERBOSE,
                             plots_path = plots_path, fits_path = save_path,
                             transitn.df = transitn.df, response.df = response.df, 
                             ols.coefs = ols.coefs))
    if (class(res) != 'try-error') {
      transitn.df = res$transition
      response.df = res$response
      ols.coefs   = res$ols.coefs
    } else {
      cat(paste('Estimation error at person', usr, '\n'))
    }
    dt = proc.time() - ptm
    
    cat(paste('User', iUser, '/', length(all_users), '( Zip', iZip,'): dt =', dt[3],'\n'))
    
    # remove this in final run
    # if (iUser > 10) break
  }
    
  # save temporary data
  write.csv(response.df, file = paste(save_path, 'response_zip_', iZip, '.csv', sep=''))
  write.csv(transitn.df, file = paste(save_path, 'transitn_zip_', iZip, '.csv', sep=''))
  write.csv(ols.coefs,   file = paste(save_path, 'olscoefs_zip_', iZip, '.csv', sep=''))  
    
  # reset output connection
  sink(NULL)
  
  rm(list = c('cur_data', 'cur_wthr'))
  gc()
  return(list(response = response.df, transition = transitn.df, ols.coefs = ols.coefs))
}

# ------------------------------------------------
# Run analysis script
# ------------------------------------------------

library(parallel)

cat('****** Running analysis on PGE data ******\n')
cat(paste('Computation on zips', ZIP_START, '-', ZIP_STOP, '\n'))

if (PROFILE) {
  ptm <- proc.time()
  Rprof(filename = paste(prof_file, sep=''), interval = 0.02, memory.profiling = T)
  result = analysis_wrapper(2, zips, user = USER, password = PASSWORD)
  Rprof(NULL)
  dt = proc.time() - ptm
  profiled = summaryRprof(filename=prof_file, memory = 'none')
} else {
  # parallel execution using parallel package
  ptm <- proc.time()
  result = mclapply(ZIP_START:ZIP_STOP, FUN = analysis_wrapper, zips, user = USER, password = PASSWORD, 
	            mc.silent = F, mc.preschedule = FALSE, mc.cores = N_PROC)
  dt = proc.time() - ptm
}

