#!/usr/bin/Rscript

# analysis_pge_selection.r
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
# Last modified: May 2013.

# -----------------------------------------------
# Define constants, load libraries and functions
# -----------------------------------------------

rm(list=ls())

options(error = recover)

# load classes, utils
setwd('~/Dropbox/OccupancyStates/')
source('code/utils/timing.r')
source('code/utils/acf_ggplot.r')
source('code/OccupancyStates.r')
source('code/occupancyAnalysis.r')

# _____________________________
# Default options for this run

NOBS_THRESH = 100    # at least 100 days for analysis
N_USR_CHUNK = 20      # number of users in a chunk
N_PROC      = 6     # number of cores to run process on 
PLOTS       = TRUE
VERBOSE     = TRUE
PROFILE     = FALSE
LOG_OUTPUT  = T

# directories
plots_path  = 'plots/plots_sel/'
save_path   = 'fits/fits_sel/'
logs_path   = 'logs/logs_sel/'
prof_file   = 'Rprof.out'

# create directories specified in profiles if they don't already exist
if (!file.exists(plots_path)) dir.create(plots_path, recursive = TRUE)
if (!file.exists(save_path))  dir.create(save_path, recursive = TRUE)
if (!file.exists(logs_path))  dir.create(logs_path, recursive = TRUE)

# _________________________________________
# Load selection data (Bakersfield zipcode)

load('data/selection_consumption.RData')

tmp = data.frame()
for (zip in names(weather.ok)) {
  tmp2 = weather.ok[[zip]]
  tmp2$ZIP5 = zip
  if (nrow(tmp)==0) tmp = tmp2 else tmp = rbind(tmp, tmp2)
}
weather.ok = tmp

# create chunks for parallel processing
usr_chunks = subset(good.overlap.info, select = c('UID', 'ZIP5', 'PER_ID', 'SP_ID'))
usr_chunks = usr_chunks[with(usr_chunks, order(ZIP5, UID)), ]
usr_chunks$CHUNK = (1:nrow(usr_chunks)) %/% N_USR_CHUNK + 1

# ------------------------------------------
# Wrapper to send computation to each core
# ------------------------------------------

analysis_wrapper <- function(chunkID) {
    
  # create folder for each chunk
  
  save_path_cur  = paste(save_path, chunkID, '/', sep='')
  if (!file.exists(save_path_cur)) dir.create(save_path_cur, recursive = TRUE)  
  
  # log output?
  if (LOG_OUTPUT) sink(paste(logs_path, "log_", chunkID,".txt", sep=''), append=F, type = c('output', 'message'))

  cat(paste('--------- Processing chunk', chunkID, '---------\n'))
    
  # get current data 
  this_chunk  = subset(usr_chunks, CHUNK == chunkID)
  chunk_zips  = unique(this_chunk$ZIP5)
  chunk_uids  = unique(this_chunk$UID)
  chunk_wthr  = subset(weather.ok, ZIP5 %in% chunk_zips)
  names(chunk_wthr)[1] = 'date'
  chunk_data  = subset(consumption.ok, UID %in% chunk_uids)
  
  # perform calculations on each user in the dataset
  hmmstats.df = list()
  response.df = list()
  transitn.df = list()
  ols.coefs   = list()
  hmm.comps   = list()
  hmmseasn.df = list()
  iUser       = 0
  all_users   = unique(chunk_data$UID)
  
  for (usr in all_users) {     
    iUser = iUser + 1
    cat(paste('Analyzing user', iUser, '/', length(all_users), '\n'))

    # select data for current user
    cur_data  = subset(chunk_data, UID == usr)
    zip       = cur_data$ZIP5[1]
    cur_wthr  = subset(chunk_wthr, ZIP5 == zip)
    
    # _______________________________________
    # Perform analysis flow for current user
    
    if (nrow(cur_data) < NOBS_THRESH) {
      cat(paste('---> Too few datapoints at user ', usr, '! (', nrow(cur_data), ')\n', sep=''))
      next
    }
    
    # plot 10% of the time
    plots_path_cur = paste(plots_path, chunkID, '/', sep='')
    if (!PLOTS) plots_path_cur = NULL else
        if(runif(1) >= 0.3) plots_path_cur = NULL    

    # perform analysis flow
    ptm <- proc.time()
    
    res = try(occupancyAnalysis(cur_data, cur_wthr, usr, zip, 
                                Kmin = 4, Kmax = 6, 
                                NOBS_THRESH = NOBS_THRESH,
                                verbose = VERBOSE,
                                thresh.R2 = 0.85, thresh.MAPE = 0.12,
                                resp.vars = c('(Intercept)', 'TemperatureF'),
                                tran.vars = c('(Intercept)', 'TemperatureF'),
                                addl.vars = c(),
                                plots_path = plots_path_cur, dump_path = NULL))  
    
    if (class(res) != 'try-error') {
      transitn.df[[iUser]] = res[[3]]$HMM.transition
      response.df[[iUser]] = res[[3]]$HMM.response
      hmmstats.df[[iUser]] = res[[3]]$HMM.stats
      ols.coefs[[iUser]]   = res[[3]]$OLS
      hmm.comps[[iUser]]   = res[[3]]$HMM.comp
      hmmseasn.df[[iUser]] = res[[3]]$HMM.season
    } else {
      cat(paste('Estimation error at person', usr, '\n'))
    }
    dt = proc.time() - ptm
    
    cat(paste('User', iUser, '/', length(all_users), ': dt =', dt[3],'\n'))  
    
    rm(list = c('res'))
  }
    
  print(paste('Length = ',length(transitn.df)))
  
  transitn.df = do.call('rbind', transitn.df)
  response.df = do.call('rbind', response.df)
  hmmstats.df = do.call('rbind', hmmstats.df)
  ols.coefs   = do.call('rbind', ols.coefs)
  hmm.comps   = do.call('rbind', hmm.comps)
  hmmseasn.df = do.call('rbind', hmmseasn.df)
  
  # save computation results to disc
  write.csv(response.df, file = paste(save_path_cur, 'response_ID_', chunkID, '.csv', sep=''))
  write.csv(transitn.df, file = paste(save_path_cur, 'transitn_ID_', chunkID, '.csv', sep=''))
  write.csv(hmmstats.df, file = paste(save_path_cur, 'hmmstats_ID_', chunkID, '.csv', sep=''))
  write.csv(ols.coefs,   file = paste(save_path_cur, 'olscoefs_ID_', chunkID, '.csv', sep=''))  
  write.csv(hmm.comps,   file = paste(save_path_cur, 'hmmcomps_ID_', chunkID, '.csv', sep=''))  
  write.csv(hmmseasn.df, file = paste(save_path_cur, 'hmmseasn_ID_', chunkID, '.csv', sep=''))  
  
  # reset output connection
  if (LOG_OUTPUT) sink(NULL)
  
  rm(list = c('chunk_data', 'chunk_wthr'))
  gc()
  return(paste('Done chunk', chunkID,'\n'))
}

# ------------------------------------------------
# Run analysis script
# ------------------------------------------------

library(parallel)
library(utils)

cat('****** Running analysis on PGE data ******\n')

if (PROFILE) Rprof(filename = paste(prof_file, sep=''), interval = 0.02, memory.profiling = F)

# options(error = recover)
# for (cur_chunk_ID in unique(usr_chunks$CHUNK)[1:1]) { 
#   analysis_wrapper(cur_chunk_ID)
# }

# parallel execution using parallel package
ptm <- proc.time()
result = mclapply(unique(usr_chunks$CHUNK), FUN = analysis_wrapper, 
             mc.silent = F, mc.preschedule = TRUE, mc.cores = N_PROC)
dt = proc.time() - ptm
print(dt)

if (PROFILE) { 
  Rprof(NULL)
  profiled = summaryRprof(filename=prof_file, memory = 'none')
}

