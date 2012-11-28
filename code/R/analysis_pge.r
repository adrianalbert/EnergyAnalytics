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
# Last modified: November 2012.

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
NOBS_THRESH = 365   # discard users with less than a year at a given premise

plots_path  = '~/EnergyAnalytics/plots/'
save_path   = '~/EnergyAnalytics/fits/'
info_file   = '~/EnergyAnalytics/data/per_sp_id_info_100.RData'

# ------------------------------------------------
# Load identification information on users
# ------------------------------------------------

# disconnect all MySQL connections
library('RMySQL')
all_cons <- dbListConnections(MySQL())
for(con in all_cons) dbDisconnect(con)

if (!file.exists(info_file)) {
  
  cat('Loading info data from MySQL database...\n')
  
  # SELECT * from pge_res_final3_unique WHERE total_duration >= 365 
  #     INTO OUTFILE 'table.csv'
  #     FIELDS TERMINATED BY ',' OPTIONALLY ENCLOSED BY '"'
  #     LINES TERMINATED BY '\n';
  
  # get a list of unique person ids
  pers_info     = run.query("select PER_ID, SP_ID, date from pge_res_final3_unique", db = 'pge_res')
  person_col    = pers_info$PER_ID
  person_ids    = sort(unique(person_col))
  
  # how many unique (PER_ID,SP_ID) tuples where people have lived for at least 1 year?
  tuples        = paste(pers_info$PER_ID, pers_info$SP_ID, sep=',')
  tuples_unq    = unique(tuples)
  length(tuples_unq)
  
  # divide up data into chunks
  ids_vec    <- seq_along(tuples_unq)
  chunks_vec <- split(tuples_unq, ceiling(ids_vec/N_USERS))
  
  # save to RData file
  save(file = info_file, list = c('chunks_vec', 'person_ids', 'tuples'))
} else {  
  cat(paste('Loading info data from file', info_file, '...\n'))  
  load(info_file)
}

# ------------------------------------------
# Wrapper to send computation to each core
# ------------------------------------------

analysis_wrapper <- function(chunk_id) {
  
  sink(paste(plots_path, "log_", chunk_id,".txt", sep=''), append=TRUE)
  
  cat(paste('--------- Processing Chunk', chunk_id, '---------\n'))
  
  tic()
  per_ids   = sapply(chunks_vec[[chunk_id]], function(x) strsplit(x, ',')[[1]][1])
  sp_ids    = sapply(chunks_vec[[chunk_id]], function(x) strsplit(x, ',')[[1]][2])  
  # retrieve current zip consumption data from DB
  query     = paste("SELECT * FROM pge_res_final3_unique WHERE PER_ID IN (", paste(per_ids, collapse = ','), 
                    ') AND SP_ID IN (', paste(sp_ids, collapse = ','), 
                    ') ORDER BY date' )
  raw_data  = run.query(query, db = 'pge_res')
  zips      = unique(raw_data$ZIP5)
  
  # retrieve corresponding weather data for current zip
  query     = paste("SELECT * FROM weather_60 WHERE zip5 IN (", 
                    paste(zips, collapse = ','), ')')
  wthr_data = run.query(query, db = 'PGE_WEATHER')
  
  time_db_sec  = toc()
  
  all_users = chunks_vec[[chunk_id]]
  # perform calculations on each user in the dataset
  response.df = data.frame()
  transitn.df = data.frame()
  ols.coefs   = data.frame()
  iUser       = 0
  for (usr_info in all_users) {     

    # _________________________________
    # Select data for current person
    
    cur_PER_ID= strsplit(usr_info,',')[[1]][1]
    cur_SP_ID = strsplit(usr_info,',')[[1]][2]    
    cur_data  = subset(raw_data, PER_ID == cur_PER_ID & SP_ID == cur_SP_ID)
    zip       = unique(cur_data$ZIP5)    
    cur_wthr  = subset(wthr_data, zip5 == zip)
    
    # discard users with too few data points
    if (nrow(cur_data) < NOBS_THRESH) {
      cat(paste('--- Too few datapoints! (', nrow(cur_data), ')\n', sep=''))
      next
    }

    # perform analysis flow
    ptm <- proc.time()
    res = try(personAnalysis(cur_data, cur_wthr, plots = F, verbose = F,
                             plots_path = plots_path, 
                             transitn.df = transitn.df, response.df = response.df, 
                             ols.coefs = ols.coefs))
    if (class(res) != 'try-error') {
      transitn.df = res$transition
      response.df = res$response
      ols.coefs   = res$ols.coefs
    } else {
      cat(paste('Estimation error at person', cur_PER_ID, ',', cur_SP_ID, '\n'))
    }
    dt = proc.time() - ptm
    
    iUser = iUser + 1
    cat(paste('User', iUser, '/', length(all_users), '( Chunk', chunk_id,'): dt =', dt[1,],'\n'))
  }
    
  # save temporary data
  write.csv(response.df, file = paste(save_path, 'response_chunk_', chunk_id, '.csv', sep=''))
  write.csv(transitn.df, file = paste(save_path, 'transitn_chunk_', chunk_id, '.csv', sep=''))
  write.csv(ols.coefs,   file = paste(save_path, 'olscoefs_chunk_', chunk_id, '.csv', sep=''))  
    
  # reset output connection
  sink()
  
  rm(list = c('cur_data', 'cur_wthr'))
  return(list(response = response.df, transition = transitn.df, ols.coefs = ols.coefs))
}

# ------------------------------------------------
# Some preliminary analysis
# ------------------------------------------------

library(utils)

# ptm <- proc.time()
# Rprof(filename = paste('~/Dropbox/ControlPatterns/fits/Rprof.out', sep=''), 
#       interval = 0.02, memory.profiling = T)
# result = analysis_wrapper(2)
# Rprof(NULL)
# # Stop the clock
# dt = proc.time() - ptm
# 
# profiled.1 = summaryRprof(filename='~/Dropbox/ControlPatterns/fits/Rprof.out', memory = 'none')

# # ___________________________
# # Profiling using rprof
# 
# library('profr')
# 
# # read profiling results
# profiled.2 = parse_rprof('~/Dropbox/ControlPatterns/fits/Rprof.out', interval=0.02)
# 
# # dump to file
# write.table(profiled.2, file = '~/Dropbox/ControlPatterns/fits/profiling.csv', quote = F, sep = '\t\t')
# 
# # profiling plot
# png('~/Dropbox/ControlPatterns/plots/profiling.png', height=600, width=800)
# ggplot.profr(profiled.2)
# dev.off()

# if in batch mode, parse arguments
args = commandArgs(TRUE)
if (length(args)>0) {
  ch_min = args[1]
  ch_max = args[2]
  if (length(args)>=3) nProc  = args[3] else nProc = 6
} else {
  ch_min = 1
  ch_max = 10#length(chunks_vec)
  nProc  = detectCores()
}

# intialize logs
for (chunk in 1:nProc){
  writeLines(c(""), paste(plots_path, "log_", chunk, ".txt", sep=''))
}

# parallel execution using parallel package
library(parallel)
ptm <- proc.time()
Rprof(filename = paste('~/Dropbox/ControlPatterns/fits/Rprof.out', sep=''), 
      interval = 0.02, memory.profiling = T)
result = mclapply(ch_min:ch_max, FUN = analysis_wrapper, mc.silent = F, 
                   mc.preschedule = FALSE, mc.cores = nProc)
Rprof(NULL)
dt = proc.time() - ptm

profiled.1 = summaryRprof(filename='~/Dropbox/ControlPatterns/fits/Rprof.out', memory = 'none')
