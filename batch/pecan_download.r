# test_pecan_download.r
#
# Download Pecan St data.
# 
# Adrian Albert
#
# Last modified: March 2014.

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

# some interesting components
select.keep  = c('dataid', 'localminute', 'use')
select.AC    = c("air1", "air2", "air3", "airwindowunit1")
select.HV    = c("furnace1", "furnace2", "heater1", "housefan1")
select.light = c("lights_plugs1", "lights_plugs2", "lights_plugs3", "lights_plugs4", "lights_plugs5", "lights_plugs6",
                 "outsidelights_plugs1", "outsidelights_plugs2")

# directory to save files
DUMP_PATH = '~/energy-data/pecan_street/'
dir.create(file.path(DUMP_PATH, '/usage'))    
dir.create(file.path(DUMP_PATH, '/metadata'))    

# ------------------------------------------
# Parser logic
# ------------------------------------------

parse_data = function(data, selection = NULL) {
  
  # select columns of interest
  if (!is.null(selection)) data = subset(data, select = selection)  
  
  # only store columns that have actual data
  col.na = which(sapply(data, function(x)all(is.na(x))))
  if (length(col.na)>0) data = data[,-col.na] else return(NULL)
  
  # if too little data return null
  if (ncol(data)<=3) return(NULL)
  
  # aggregate at 15-minute and hourly levels
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

## get UIDs for users in Pecan data consumption time series
years = c(2012, 2013, 2014)
uids  = list()
for (yr in years) {
  query = paste("select distinct dataid from \"PecanStreet_RawData\".egauge_minutes_", yr, sep = '')
  uids[[as.character(yr)]] <- dbGetQuery(con, query)
}

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
# Parse and dump time series data
# ------------------------------------------

# pull time series data 
for (yr in names(uids)) {
  cur_uids = uids[[yr]]$dataid
  yr       = as.numeric(yr)
  res      = lapply(1:length(cur_uids), function(i) {       
    uid = cur_uids[i]
    cat(paste('Fetching user ', i, '/', length(cur_uids), ' (ID=', uid, ')\n', sep = ''))
    
    # get current data    
    query    = paste("select * from \"PecanStreet_RawData\".egauge_minutes_", yr, " where dataid = ", uid, sep = '')
    cur_data = dbGetQuery(con, query)   
        
    # parse data 
    ok_data = parse_data(cur_data, selection = c(select.keep, select.AC, select.HV, select.light))
    if is.null(ok_data) return(-1)
    
    # save to CSV file; if file exists, append
    filename = paste('usage/', uid, '.csv', sep = '')
    filepath = file.path(DUMP_PATH, filename)
    if (file.exists(filepath)) {
      write.table(ok_data, file = filepath, col.names=FALSE, sep=",", append=TRUE)        
    } else {
      write.csv(ok_data, file = filepath)
    }    
  
    # return columns available for each user
    return(setdiff(names(data), select.keep))
  })  
}

# ------------------------------------------
# Disconnect from database
# ------------------------------------------

## Closes the connection
dbDisconnect(con)

## Frees all the resources on the driver
dbUnloadDriver(drv)
