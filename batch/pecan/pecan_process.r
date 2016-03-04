# pecan_process.r
#
# Parse & integrate Pecan St data.
# Clean: remove unreasonable values, 
#        fill in small gaps, 
#        remove columns with many missing values.
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
library('parallel')

setwd('~/EnergyAnalytics/')
source('./utils/aggregate_data.r')
source('./utils/weather/clean_weather_data.r')

# directory to save files
DATA_PATH = '~/S3L_server/energy-data/pecan_street/'
IN_PATH   = paste(DATA_PATH, 'usage-orig', sep = '')
OUT_PATH  = paste(DATA_PATH, 'usage-processed', sep = '')
METADATA_PATH = '~/S3L_server/energy-data/pecan_street/metadata/'
dir.create(file.path(OUT_PATH))    

# ----------------------------------------------
# Assign each user a unique human-readable name
# ----------------------------------------------

files_all = list.files(path=IN_PATH, full.names = T, recursive = T)
users_all = unique(sapply(files_all, function(s) strsplit(tail(strsplit(s, '/')[[1]], 1), '\\.')[[1]][1]))

# load baby names
# we're going to name each user for later easiness of use
baby_names  = read.csv('~/S3L_server/energy-data/pecan_street/baby-names.csv')
user_names  = data.frame(ID = users_all, name = unique(baby_names$name)[1:length(users_all)])

# save names to file
write.csv(user_names, file = paste(METADATA_PATH, 'user_names_ids.csv', sep = '/'))

# ------------------------------------------
# Logic to process data 
# ------------------------------------------

process_data = function(raw, dateCol = 'date', method = 'IRMI') {
  
  # remove nonsensical values
  df = sapply(raw[,-which(names(raw) == dateCol)], function(x) {
    q   = quantile(x, probs = c(0.005, 0.995), na.rm = T)
    idx = which(x < q[1] | x > q[2])
    if (length(idx)>0) x[idx] = NA
    return(x)
  })
  df = as.data.frame(df)
  dat= as.POSIXct(raw[,dateCol])
  
  # remove columns with all missing data/zeroes
  col.na = which(sapply(df, function(x)all(is.na(x) | x==0)))
  if (length(col.na)>0) df = df[,-col.na] 
  
  # have we removed all but one column?
  if (!is.data.frame(df)) return(NULL)

  # remove columns with many NAs
  nr.na = sapply(df, function(x)length(which(is.na(x))))
  nr.na = which(nr.na >= 0.9*nrow(df))
  if (length(nr.na)>0) df = df[,-nr.na] 
  
  # have we removed all but one column?
  if (!is.data.frame(df)) return(NULL)
  
  # remove complete NAs row-wise
  idx = apply(df, 1, function(x) sum(as.numeric(is.na(x))))
  idx = which(idx == ncol(df))
  if (length(idx)>0) {
    df  = df[-idx, ]
    dat = dat[-idx]
  }
  
  # imputate missing data if still any
  no.nas = apply(df, 2, function(x) length(which(is.na(x))))
  if (sum(no.nas) > 0) {    
    imp = imputateMissingValues(df, method = method) 
  } else imp = df
  
  # add date column back in
  imp = cbind(dat, imp)
  names(imp)[1] = 'date'
  
  return(imp)
}

# ------------------------------------------
# Aggregate data in given directory
# ------------------------------------------

files.proc  = list.files(OUT_PATH, full.names = TRUE, recursive = T)
uids.done   = sapply(files.proc, function(x) {
  tmp = strsplit(x, '/')[[1]]
  uid = tmp[length(tmp)]
  uid = strsplit(uid, '_')[[1]][1]
  yr  = tmp[length(tmp)-1]
  return(paste(yr, uid, sep = '/'))
})
files.input = list.files(IN_PATH, full.names = TRUE, recursive = T)
uids.input  = sapply(files.input, function(x) {
  tmp = strsplit(x, '/')[[1]]
  uid = tmp[length(tmp)]
  uid = strsplit(uid, '\\.')[[1]][1]
  yr  = tmp[length(tmp)-1]
  return(paste(yr, uid, sep = '/'))
})

idx = which(!(uids.input %in% uids.done))

if (length(idx) == 0) stop("All files have already been processed!")

res      = mclapply(1:length(files.input[idx]), 
                    mc.cores = 4, 
                    function(i) {       
  file = files.input[idx[i]]  
  cat(paste('Processing ', file, ' : ', idx[i], '/', length(files.input[idx]), '\n', sep = ''))
  year = strsplit(file, '/')[[1]]
  year = year[length(year)-1]
  
  # get current data  
  data = read.csv(file)   
  uid  = data$dataid[1]
  data = data[,-1]  
  if (uid %in% uids.done) {
    cat('Already processed!\n')
    return(0)
  }
  
  # process minute-level data
  data_01_proc = process_data(data, method = 'NONE', dateCol = 'localminute')  
  
  # is there enough data for this user?
  if (is.null(data_01_proc)) {
    cat('Not enough data for this user!\n')
    return(NA);
  }
  
  # if enough minute-level data, aggregate data to higher resolutions
  data_60 = toHourly(data, dateCol = 'localminute')
  data_15 = toXmin(data, min = 15, dateCol = 'localminute')
  
  # clean data
  data_60_proc = process_data(data_60, method = 'NONE')
  data_15_proc = process_data(data_15, method = 'NONE')
  
  # is there enough data for this user?
  if (is.null(data_60_proc) || is.null(data_15_proc)) {
    cat('Not enough data for this user!\n')
    return(NA);
  }
  
  # save to file
  dir.create(file.path(OUT_PATH, year))      
  write.csv(data_01_proc, file = paste(OUT_PATH, paste(year, '/', uid, '_minute.csv', sep = ''), sep = '/'), row.names = F)
  write.csv(data_15_proc, file = paste(OUT_PATH, paste(year, '/', uid, '_15mins.csv', sep = ''), sep = '/'), row.names = F)
  write.csv(data_60_proc, file = paste(OUT_PATH, paste(year, '/', uid, '_hourly.csv', sep = ''), sep = '/'), row.names = F)
  
  return(0)
})  
