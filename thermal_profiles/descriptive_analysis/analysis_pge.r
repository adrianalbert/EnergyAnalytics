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

# load classes, utils
setwd('~/EnergyAnalytics/code/R')
source('~/EnergyAnalytics/code/R/utils/timing.r')
source('~/EnergyAnalytics/code/R/utils/sql_utils.r')
source('~/EnergyAnalytics/code/R/utils/acf_ggplot.r')
source('~/EnergyAnalytics/code/R/Person.r')
source('~/EnergyAnalytics/code/R/personAnalysis.r')

# _____________________________
# Default options for this run

NOBS_THRESH = 90   # discard users with less than a 1/4 year at a given premise
N_PROC      = 1     # number of cores to run process on 
DATA_ACCESS = 'sql' # flag to indicate type of media used for data
PLOTS       = FALSE
VERBOSE     = TRUE
USER        = 'adrian'
PASSWORD    = 'xmenxmen'
ZIP_START   = 1
ZIP_STOP    = 100
PROFILE     = TRUE
LOG_OUTPUT  = FALSE
SAVE_TO_SQL = TRUE
OVERWRITE   = TRUE

# directories
plots_path  = '~/Dropbox/ControlPatterns/plots2/'
save_path   = '~/EnergyAnalytics/fits2/'
logs_path   = '~/EnergyAnalytics/logs2/'
kwh_path    = ''
wthr_path   = ''
prof_file   = '~/EnergyAnalytics/Rprof.out'
zips_file   = '~/EnergyAnalytics/data/metadata/zipcode_res.csv'

# _______________________________________________
# If in batch mode, parse command-line arguments

library(utils)
myargs = commandArgs(trailingOnly = TRUE)
cat(paste('Arguments:', paste(myargs,sep='; '), '\n'))

# Was a profile file specified?
if (length(myargs)==1) {
  cat(paste('Loading profile from file...', myargs, '\n'))
  source(myargs)
}

# Are parameters given in command line?
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

# create directories specified in profiles if they don't already exist
if (!file.exists(plots_path)) dir.create(plots_path, recursive = TRUE)
if (!file.exists(save_path)) dir.create(save_path, recursive = TRUE)
if (!file.exists(logs_path)) dir.create(logs_path, recursive = TRUE)

# _________________________________________
# Load identification information on users

# disconnect all MySQL connections
library('RMySQL')
#all_cons <- dbListConnections(MySQL())
#for(con in all_cons) dbDisconnect(con)

# get zipcodes
zips = read.csv(zips_file)
zips = zips$ZIP5

# re-analyze zipcodes already processed? 
db.con <- dbConnect(dbDriver("MySQL"), 
                    host = "sssl-cyclops.stanford.edu",
                    user = USER, password = PASSWORD, 
                    dbname = 'PGE_FITS')
tables = dbListTables(db.con)
if (!OVERWRITE) {
   # open up connection to database
   if (length(tables)>0) {
     tmp = strsplit(tables, '_')
     tmp = unique(sapply(1:length(tmp), function(j) tmp[[j]][1]))
     processed.zips = substring(tmp, 4, nchar(tmp))	
   } else processed.zips = c()
   zips.todo = setdiff(1:length(zips), as.numeric(processed.zips))
   cat(sprintf('Already processed: %d zips\n', length(processed.zips)))
   cat(sprintf('Left to process:   %d zips\n', length(zips.todo)))
} else { 
  cat('Deleting tables from database...\n')
  zips.todo = 1:length(zips)
  for (tab in tables) dbRemoveTable(db.con, tab)
}
dbDisconnect(db.con)

# ------------------------------------------
# Wrapper to send computation to each core
# ------------------------------------------

analysis_wrapper <- function(iZip, zips, type = 'sql', SAVE_TO_SQL = FALSE) {
   # check whether zipcode is already processed
   db.con <- dbConnect(dbDriver("MySQL"), 
                        host = "sssl-cyclops.stanford.edu",
                        user = USER, password = PASSWORD, 
                        dbname = 'PGE_FITS')
   tables = dbListTables(db.con)
   dbDisconnect(db.con)
   if (length(tables)>0) {
     tmp = strsplit(tables, '_')
     tmp = unique(sapply(1:length(tmp), function(j) tmp[[j]][1]))
     processed.zips = substring(tmp, 4, nchar(tmp))  
   } else processed.zips = c()   

   # if already processed, terminate
   if (iZip %in% processed.zips) {
	  stop(sprintf('Already processed zip %d\n', iZip))	
   }
    
  # log output?
  if (LOG_OUTPUT) sink(paste(logs_path, "log_", iZip,".txt", sep=''), append=F, type = c('output', 'message'))

  cat(paste('--------- Processing Zipcode', iZip, '/', length(zips.todo), '---------\n'))
  
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
  
  # assign each (PER_ID, PREM_ID) tuple a unique UID 
   
  # perform calculations on each user in the dataset
  hmmstats.df = data.frame()
  response.df = data.frame()
  transitn.df = data.frame()
  ols.coefs   = data.frame()
  ols.comps   = data.frame()
  hmm.comps   = data.frame()
  iUser       = 0
  all_users   = unique(raw_data$UID)
  for (usr in all_users) {     
    iUser = iUser + 1
    cat(paste('Analyzing user', iUser, '/', length(all_users), '( Zip', iZip,')\n'))

    # _______________________________________
    # Perform analysis flow for current user
    
    cur_data  = subset(raw_data, UID == usr)
    if (nrow(cur_data) < NOBS_THRESH) {
      cat(paste('---> Too few datapoints at user ', usr, '! (', nrow(cur_data), ')\n', sep=''))
      next
    }
    
    # perform analysis flow
    ptm <- proc.time()
    res = try(personAnalysis(cur_data, cur_wthr, usr, zips[iZip], NOBS_THRESH = NOBS_THRESH,
                             plots = PLOTS, verbose = VERBOSE,
                             plots_path = plots_path, fits_path = save_path,
                             transitn.df = transitn.df, response.df = response.df, 
                             hmmstats.df = hmmstats.df, 
                             ols.coefs = ols.coefs, ols.comps = ols.comps, hmm.comps = hmm.comps))
    if (class(res) != 'try-error') {
      transitn.df = res$transition
      response.df = res$response
      hmmstats.df = res$hmmstats
      ols.coefs   = res$ols.coefs
      ols.comps   = res$ols.comps
      hmm.comps   = res$hmm.comps
    } else {
      cat(paste('Estimation error at person', usr, '\n'))
    }
    dt = proc.time() - ptm
    
    cat(paste('User', iUser, '/', length(all_users), '( Zip', iZip,'): dt =', dt[3],'\n'))
    
    # remove this in final run
    # if (iUser > 10) break
  }
    
  # save computation results to disc
  write.csv(response.df, file = paste(save_path, 'response_zip_', iZip, '.csv', sep=''))
  write.csv(transitn.df, file = paste(save_path, 'transitn_zip_', iZip, '.csv', sep=''))
  write.csv(hmmstats.df, file = paste(save_path, 'hmmstats_zip_', iZip, '.csv', sep=''))
  write.csv(ols.coefs,   file = paste(save_path, 'olscoefs_zip_', iZip, '.csv', sep=''))  
  write.csv(ols.comps,   file = paste(save_path, 'olscomps_zip_', iZip, '.csv', sep=''))  
  write.csv(hmm.comps,   file = paste(save_path, 'hmmcomps_zip_', iZip, '.csv', sep=''))  
     
  # store computation data to sql, if requested
  if (SAVE_TO_SQL) {
    cat(paste('***** Storing computation to database for zip ', iZip, '*****\n'))

    # open up connection to database
    db.con <- dbConnect(dbDriver("MySQL"), 
                        host = "sssl-cyclops.stanford.edu",
                        user = USER, password = PASSWORD, 
                        dbname = 'PGE_FITS')

    # write fits
    dbWriteTable(db.con, paste("ZIP", iZip, '_OLS_COEF', sep=''), ols.coefs, overwrite = T)
    dbWriteTable(db.con, paste("ZIP", iZip, '_HMM_RESP', sep=''), response.df, overwrite = T)
    dbWriteTable(db.con, paste("ZIP", iZip, '_HMM_TRAN', sep=''), transitn.df, overwrite = T)  
    dbWriteTable(db.con, paste("ZIP", iZip, '_HMM_STAT', sep=''), hmmstats.df, overwrite = T)  
    # write component stats
    dbWriteTable(db.con, paste("ZIP", iZip, '_OLS_COMP', sep=''), ols.comps, overwrite = T)
    dbWriteTable(db.con, paste("ZIP", iZip, '_HMM_COMP', sep=''), hmm.comps, overwrite = T)
    # disconnect from database
    dbDisconnect(db.con)
  }

  # reset output connection
  if (LOG_OUTPUT) sink(NULL)
  
  rm(list = c('cur_data', 'cur_wthr'))
  gc()
  return(paste('Done zip', iZip,'\n'))
}

# ------------------------------------------------
# Run analysis script
# ------------------------------------------------

library(parallel)

cat('****** Running analysis on PGE data ******\n')
cat('Computation on zips\n')
print(zips.todo[ZIP_START:ZIP_STOP])

if (PROFILE) Rprof(filename = paste(prof_file, sep=''), interval = 0.02, memory.profiling = T)

# for (zip in zips.todo[ZIP_START:ZIP_STOP]) {
#   analysis_wrapper(zip, zips)
# }

# parallel execution using parallel package
ptm <- proc.time()
result = mclapply(zips.todo[ZIP_START:ZIP_STOP], FUN = analysis_wrapper, zips, SAVE_TO_SQL = SAVE_TO_SQL, 
             mc.silent = F, mc.preschedule = FALSE, mc.cores = N_PROC)
dt = proc.time() - ptm
print(dt)

if (PROFILE) { 
  Rprof(NULL)
  profiled = summaryRprof(filename=prof_file, memory = 'none')
}

