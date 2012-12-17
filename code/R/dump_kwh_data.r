# dump_kwh_data.r
#
# Dumps PGE consumption data for easier access (per PER_ID/SP_ID tuple).
# Also, create key table for (PER_ID, SP_ID) characteristics. This table contains
# actual identification information; files only contain timestamps & consumption.
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
prof_file   = '~/EnergyAnalytics/Rprof_kwh.out'
dump_path   = '~/EnergyAnalytics/data/pge_data/consumption/'
info_file   = '~/EnergyAnalytics/data/user_info.RData'

N_USERS     = 100   # send batches of 100 users to each core
NOBS_THRESH = 365   # discard users with less than a year at a given premise

# ------------------------------------------------
# Load identification information on users
# ------------------------------------------------

# disconnect all MySQL connections
library('RMySQL')
all_cons <- dbListConnections(MySQL())
for(con in all_cons) dbDisconnect(con)

if (!file.exists(info_file)) {
  
  cat('Loading info data from MySQL database...\n')
  
  # get a list of unique person ids
  pers_ids      = run.query("select PER_ID, SP_ID, UID, SM_DURATION from pge_res_final3_unique", db = 'pge_res')
  tuples        = paste(pers_ids$PER_ID, pers_ids$SP_ID, pers_ids$UID, sep=',')
  tuples_unq    = unique(tuples)
  
  print(length(tuples_unq))

  # get user information table 
  pers_info     = run.query("select * from pge_res_final", db = 'pge_res')
  tmp_tuples    = paste(pers_info$PER_ID, pers_info$SP_ID, pers_info$UID, sep=',')
  idx_ok        = which(tmp_tuples %in% tuples_unq)
  pers_info     = pers_info[idx_ok,]
  
  # divide up data into chunks by zipcode
  zips          = unique(pers_info$ZIP)
  
  # save to RData file
  save(file = info_file, list = c('pers_info', 'zips'))
} else {  
  cat(paste('Loading info data from file', info_file, '...\n'))  
  load(info_file)
}

# ------------------------------------------
# Wrapper to send computation to each core
# ------------------------------------------

wrapper <- function(iZip, path = dump_path) {
  
  cat(paste('--------- Processing Zipcode', iZip, '/', length(zips), '---------\n'))
  
  # retrieve current zip consumption data from DB
  query     = paste("SELECT * FROM pge_res_60_", zips[iZip], ' ORDER BY date', sep='')
  raw_data  = run.query(query, db = 'PGE_SAM')
  
  # only retain those users with more than NTHRESH observations  
  raw_data  = subset(raw_data, SM_DURATION >= NOBS_THRESH)
    
  # dump current zipcode to file
  col_names = c('UID', 'date', paste('hkw', 1:24, sep=''))
  raw_data  = raw_data[,col_names]    
  write.csv(file = paste(path, iZip, '.csv', sep=''), raw_data, row.names = F)    

  rm(list = c('raw_data'))
  gc()
  return(NULL)
}

# ------------------------------------------------
# Dump data to disk
# ------------------------------------------------

# parallel execution using parallel package
library(parallel)
ptm <- proc.time()
Rprof(filename = paste(prof_file, sep=''), interval = 0.02, memory.profiling = T)
result = mclapply(1:length(zips), FUN = wrapper, mc.silent = F, mc.preschedule = FALSE, mc.cores = 6)
Rprof(NULL)
dt = proc.time() - ptm
print(dt)
profiled = summaryRprof(filename=prof_file, memory = 'none')

