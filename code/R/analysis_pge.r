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

N_USERS     = 100   # send batches of 50 users to each core
NOBS_THRESH = 180   # discard users with less than a 1/2 year at a given premise

plots_path  = '~/Dropbox/ControlPatterns/plots/'
save_path   = '~/EnergyAnalytics/fits/'
info_file   = '~/EnergyAnalytics/data/per_sp_id_info.RData'
prof_file   = '~/EnergyAnalytics/Rprof.out'

# ------------------------------------------------
# Load identification information on users
# ------------------------------------------------

# disconnect all MySQL connections
library('RMySQL')
all_cons <- dbListConnections(MySQL())
for(con in all_cons) dbDisconnect(con)

# get zipcodes
zips = read.csv('~/EnergyAnalytics/data/metadata/zipcode_res.csv')
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
    raw_data  = run.query(query, db = 'PGE_SAM')  
    
    # retrieve current zip weather data from DB
    query     = paste("SELECT * FROM ZIP_", zips[iZip], ' ORDER BY date', sep='')
    cur_wthr  = run.query(query, db = 'PGE_WEATHER')  
    
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
                             plots = T, verbose = T,
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
# Some preliminary analysis
# ------------------------------------------------

library(utils)
# # 
# ptm <- proc.time()
# Rprof(filename = paste(prof_file, sep=''), interval = 0.02, memory.profiling = T)
# result = analysis_wrapper(20)
# Rprof(NULL)
# # Stop the clock
# dt = proc.time() - ptm
# 
# profiled = summaryRprof(filename=prof_file, memory = 'none')

library(parallel)

# if in batch mode, parse arguments
args = commandArgs(TRUE)
if (length(args)>0) {
  ch_min = args[1]
  ch_max = args[2]
  if (length(args)>=3) nProc  = args[3] else nProc = detectCores()
} else {
  ch_min = 1
  ch_max = 20
  nProc  = detectCores()
}
cat(paste('Computation on zips', ch_min, '-', ch_max, '\n'))

# parallel execution using parallel package
ptm <- proc.time()
Rprof(filename = paste(prof_file, sep=''), 
      interval = 0.02, memory.profiling = T)
result = mclapply(ch_min:ch_max, FUN = analysis_wrapper, zips, mc.silent = F, 
                   mc.preschedule = FALSE, mc.cores = nProc)
Rprof(NULL)
dt = proc.time() - ptm

profiled = summaryRprof(filename=prof_file, memory = 'none')
