# analysis_OLS_pge.r
# 
# Typical OLS analysis for PGE data:
# - regression residuals: argue non-Gaussian
# - ACF: argue non-iid
# - different models for different time-of-day/holiday/week-weekend: argue different states
# - two users: argue different responses
#
# Adrian Albert
# Last modified: November 2012.

# -----------------------------------------------
# Define constants, load libraries and functions
# -----------------------------------------------

rm(list=ls())

source('~/Dropbox/ControlPatterns/code/R/utils/timing.r')
source('~/Dropbox/ControlPatterns/code/R/utils/sql_utils.r')
source('~/Dropbox/ControlPatterns/code/R/Person.r')

N_USERS = 5000

# ------------------------------------------------
# Perform some initial analysis
# ------------------------------------------------

# get a list of unique person ids
pers_info     = run.query("select PER_ID, SP_ID, date, from pge_res_final3_unique", db = 'pge_res')
pers_info_ok  = run.query("select PER_ID, SP_ID, date, from pge_res_final3_unique WHERE total_duration >= 365", db = 'pge_res')
person_col = pers_info$PER_ID
person_ids = sort(unique(person_col))

# how many unique (PER_ID,SP_ID) tuples?
tuples = paste(pers_info$PER_ID, pers_info$SP_ID, sep=',')
tuples_unique = unique(tuples)
length(tuples_unique)

# ------------------------------------------------
# Distribute computation across nodes 
# ------------------------------------------------
# Break down computation on a per-zipcode basis

# split users for distribution across cores
tmp <- seq_along(person_ids)
person_chunks <- split(person_ids, ceiling(tmp/N_USERS))

# wrapper to send computation to each core
analysis_wrapper <- function(chunk_id) {
  
  tic()
  
  # retrieve current consumption data from DB
  all_users = person_chunks[[chunk_id]]
  query     = paste("SELECT * FROM PGE_SAM WHERE PER_ID IN (", 
                    paste(all_users,collapse=','), ')')
  raw_data  = run.query(query)
  
  # retrieve corresponding weather data for current zips
  zips      = unique(raw_data$ZIP5)
  query     = paste("SELECT * FROM weather_60 WHERE zip5 IN (", 
                    paste(zips,collapse=','), ')')
  wthr_data = run.query(query, db = 'PGE_WEATHER')
  
  time_sec  = toc()
  
  # perform calculations on each user in the dataset
  for (user_id in all_users) {   
    
    # get data for current chunk of people
    cur_data  = subset(raw_data, PER_ID == user_id)
    cur_wthr  = subset(wthr_data, zip5 %in% unique(cur_data$ZIP5))
    user      = new(Class='Person', cur_data)  
    user      = addWeather(user, cur_wthr)
    
    # perform OLS analysis on current user
    
    # perform HMM analysis on current user
    
    # save analysis results
    
    # clear used variables
    rm(c('cur_data', 'cur_wthr', 'user'))
  }
  
  return(TRUE)
}

# parallel execution using multicore
library(multicore)
result = mclapply(1:length(person_chunks), FUN = analysis_wrapper,
                   mc.preschedule = TRUE, mc.cores = 6)


