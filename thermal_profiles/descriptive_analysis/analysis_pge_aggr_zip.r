#!/usr/bin/Rscript

# analysis_pge_aggr_zip.r
# 
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
N_PROC      = 10     # number of cores to run process on 
PLOTS       = TRUE
VERBOSE     = TRUE
PROFILE     = FALSE
LOG_OUTPUT  = T

# directories
plots_path  = 'plots/aggregate/zips/'
save_path   = 'fits/aggregate/'
logs_path   = 'logs/aggregate/'
prof_file   = 'Rprof.out'

# create directories specified in profiles if they don't already exist
if (!file.exists(plots_path)) dir.create(plots_path, recursive = TRUE)
if (!file.exists(save_path))  dir.create(save_path, recursive = TRUE)
if (!file.exists(logs_path))  dir.create(logs_path, recursive = TRUE)

# _________________________________________
# Load aggregate data per zipcode

load('data/zip_aggregate.RData')

# set up data in the format taken by the analysis wrapper
zips  = names(aggrZip)
noUsr = sapply(aggrZip, function(l) unique(l$Count))
aggrZip = do.call('rbind', aggrZip)
aggrZip = cbind(UID = rep(1:length(zips), each = 365), ZIP5 = rep(zips, each = 365), aggrZip[,-ncol(aggrZip)])
aggrZip$CHUNK = aggrZip$UID %/% N_USR_CHUNK + 1
usr_chunks = aggrZip
weather.ok = weather
names(weather.ok)[1] = 'ZIP5'

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
  chunk_data  = this_chunk[,-nrow(this_chunk)]
  
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

# for (cur_chunk_ID in unique(usr_chunks$CHUNK)[1:2]) { 
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

