# pecan_process.r
#
# Parse & integrate Pecan St data.
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
library('VIM')

setwd('~/EnergyAnalytics/')
source('./utils/aggregate_data.r')
source('./utils/weather/clean_weather_data.r')

# directory to save files
DATA_PATH = '~/energy-data/pecan_street/'
IN_PATH   = paste(DATA_PATH, 'usage-orig', sep = '')
OUT_PATH  = paste(DATA_PATH, 'usage-processed', sep = '')
dir.create(file.path(OUT_PATH))    

# ------------------------------------------
# Logic to process data 
# ------------------------------------------

process_data = function(raw, thresh = 50, dateCol = 'date') {
  
  # remove nonsensical values
  df = sapply(raw[,-which(names(raw) == dateCol)], function(x) {
    idx = which(x < 0 | x > thresh)
    x[idx] = NA
    return(x)
  })
  df = as.data.frame(df)
  dat= as.POSIXct(raw[,dateCol])
  
  # remove columns with too much missing data
  col.na = which(sapply(df, function(x)all(is.na(x))) | sapply(df, function(x)all(x==0)))
  if (length(col.na)>0) df = df[,-col.na] 
  
  nr.na = sapply(df, function(x)length(which(is.na(x))))
  nr.na = which(nr.na >= nrow(df)/3)
  if (length(nr.na)>0) df = df[,-nr.na] 
  
  if (!is.data.frame(df)) return(NULL)
  
  # remove complete NAs row-wise
  idx = apply(df, 1, function(x) sum(as.numeric(is.na(x))))
  idx = which(idx == ncol(df))
  if (length(idx)>0) {
    df  = df[-idx, ]
    dat = dat[-idx]
  }
  
  # imputate missing data  
  no.nas = apply(df, 2, function(x) length(which(is.na(x))))
  if (sum(no.nas) > 0) {    
#     imp = irmi(df) 
    imp = imputateMissingValues(df) 
  } else imp = df
  imp[,dateCol] = dat
  
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
                    mc.cores = 5, 
                    function(i) {       
  file = files.input[idx[i]]  
  cat(paste('Processing ', file, ' : ', idx[i], '/', length(files.input[idx]), sep = ''))
  year = strsplit(file, '/')[[1]]
  year = year[length(year)-1]
  
  # get current data  
  data = read.csv(file)   
  uid  = data$dataid[1]
  data = data[,-1]  
  if (uid %in% uids.done) {
    cat('Already processed!\n')
    return(0)
  } else cat('\n')
  
  #source('~/EnergyAnalytics/utils/aggregate_data.r')
  # aggregate data
  data_60 = toHourly(data, dateCol = 'localminute')
  data_15 = toXmin(data, dateCol = 'localminute', min = 15)
  
  # clean data
  data_60 = process_data(data_60)
  data_15 = process_data(data_15)
  
  # is there enough data for this user?
  if (is.null(data_60) || is.null(data_15)) return(NA);
  
  # save to file
  dir.create(file.path(OUT_PATH, year))      
  # write.csv(data, file = paste(OUT_PATH, paste(year, '/', uid, '_minute.csv', sep = ''), sep = '/'), row.names = F)
  write.csv(data_15, file = paste(OUT_PATH, paste(year, '/', uid, '_15mins.csv', sep = ''), sep = '/'), row.names = F)
  write.csv(data_60, file = paste(OUT_PATH, paste(year, '/', uid, '_hourly.csv', sep = ''), sep = '/'), row.names = F)
  
  return(0)
})  
