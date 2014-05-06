# test_pecan_download.r
#
# Download Pecan St data.
# 
# Adrian Albert
#
# Last modified: May 2014.

# -------------------------------
# Initializations ... 
# -------------------------------

rm(list=ls())

options(error = recover)

library('lubridate')
library('zoo')
library('ggplot2')
library('reshape')
library('dummies')
library('RPostgreSQL')
library('imputation')
library('parallel')

# directory to save files
DUMP_PATH = '~/energy-data/pecan_street/'
USGE_PATH = paste(DUMP_PATH, 'usage-orig', sep = '/')
META_PATH = paste(DUMP_PATH, 'metadata', sep = '/')
dir.create(file.path(USGE_PATH))    
dir.create(file.path(DUMP_PATH, '/metadata'))    

# ------------------------------------------
# Parser logic
# ------------------------------------------

parse_data = function(data, selection = NULL) {
  
  # select columns of interest
  if (!is.null(selection)) data = subset(data, select = selection)  
  
  # only store columns that have actual data
  col.na = which(sapply(data, function(x)all(is.na(x))) | sapply(data, function(x)all(x==0)))
  if (length(col.na)>0) data = data[,-col.na]
  
  # if too little data return null
  if (ncol(data)<=3) return(NULL)
  
  return(data)  
}

# ------------------------------------------
# Connect to database, get UIDs
# ------------------------------------------

# disconnect all old connections, if any
all_cons <- dbListConnections(PostgreSQL())
for(con in all_cons) dbDisconnect(con)

## loads the PostgreSQL driver
drv <- dbDriver("PostgreSQL")
dbListConnections(drv)

## Open a connection
con <- dbConnect(drv,
                 user= "rajagopal", password="smartenergy",
                 host = 'db.wiki-energy.org', port = 5432,
                 dbname="postgres")

# list all tables
dbListTables(con)

# ------------------------------------------
# Save meta-data on users
# ------------------------------------------

for (yr in c('2012', '2013')) {
  # get survey data
  query    = paste("select * from \"PecanStreet_RawData\".participant_survey_", yr, sep = '')
  cur_data = dbGetQuery(con, query)   
  filepath = file.path(DUMP_PATH, paste('metadata/survey_', yr, '.csv', sep=''))
  write.csv(cur_data, file = filepath)  
  
  # get survey data description
  query    = paste("select * from \"PecanStreet_RawData\".participant_survey_", yr, "_field_descriptions", sep = '')
  cur_data = dbGetQuery(con, query)   
  filepath = file.path(DUMP_PATH, paste('metadata/survey_desc_', yr, '.csv', sep=''))
  write.csv(cur_data, file = filepath)    
}


# ------------------------------------------
# Extract user IDs first - so that data can
# be extracted per user (not all in one go)
# ------------------------------------------

uids  = list()
for (yr in 2012:2014) {
  query    = paste("select distinct dataid from \"PecanStreet_RawData\".egauge_minutes_", yr, sep = '')
  cur_data = dbGetQuery(con, query)       
  uids[[as.character(yr)]] = unique(cur_data$dataid)    
}

# ------------------------------------------
# Parse and dump time series data
# ------------------------------------------

# pull raw time series data 
for (yr in names(uids)) {
  cur_uids = uids[[yr]]
  yr       = as.numeric(yr)
  dir.create(file.path(USGE_PATH, yr))    
  res      = lapply(1:length(cur_uids), function(i) {       
    uid = cur_uids[i]
    cat(paste('Fetching user ', i, '/', length(cur_uids), ' (ID=', uid, ')', sep = ''))
    
    # get raw data for current user
    query    = paste("select * from \"PecanStreet_RawData\".egauge_minutes_", yr, " where dataid = ", uid, sep = '')
    raw_data = dbGetQuery(con, query)       
    if (nrow(raw_data) == 0) {
      cat('...no data for this user!\n')
      return(-1)      
    }
    
    # parse data 
    raw_data_ok = parse_data(raw_data) 
    if (is.null(raw_data_ok)) {
      cat('...not enough data!\n')
      return(-1)
    }
    
    # save to CSV file; if file exists, append
    filename = paste(uid, '.csv', sep = '')
    filepath = file.path(USGE_PATH, yr, filename)
    write.csv(raw_data_ok, file = filepath, row.names = F)
    cat('...done\n')
  
    # return columns available for each user
    return(names(data))
  })  
}

# pull curated data
yr = '2013'
dir.create(file.path(USGE_PATH, yr))    
for (group in 1:2) {
  cur_uids = uids[[yr]]
  
  res      = lapply(1:11, function(mo) {       
    cat(paste('Fetching group ', group, ' month ', mo, '/', 11, sep = ''))
    
    query    = paste("select * from \"PecanStreet_CuratedSets\".group", group, "_disaggregated_", yr, "_", sprintf("%02d",10), sep='')
    cur_data = dbGetQuery(con, query)    
    
    # parse data 
    cur_data_ok = parse_data(cur_data) 
    if (is.null(cur_data_ok)) {
      cat('...not enough data!\n')
    }
    
    return(cur_data)
  })  
  res = do.call('rbind', res)
  
  # save to CSV file; if file exists, append
  filename = paste('group', group, '_', yr, '.csv', sep = '')
  filepath = file.path(USGE_PATH, filename)
  write.csv(res, file = filepath, row.names = F)
  cat('...done\n')  
}

# ------------------------------------------
# Disconnect from database
# ------------------------------------------

## Closes the connection
dbDisconnect(con)

## Frees all the resources on the driver
dbUnloadDriver(drv)
