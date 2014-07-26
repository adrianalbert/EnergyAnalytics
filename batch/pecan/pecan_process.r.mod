# pecan_process.r
#
# Parse & integrate Pecan St data.
# Clean: remove unreasonable values, 
#        fill in small gaps, 
#        remove columns with many missing values.
# 
# Adrian Albert
#
# Last modified: June 2014.

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
DATA_PATH = '~/energy-data/pecan_street/'
IN_PATH   = paste(DATA_PATH, 'usage-orig', sep = '')
OUT_PATH  = paste(DATA_PATH, 'usage-processed', sep = '')
dir.create(file.path(OUT_PATH))    

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
  nr.na = which(nr.na >= nrow(df)/2)
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
# Concatenate data across multiple years
# ------------------------------------------

concatenate_data = function(file_list) {
  
  # load data
  data_list = lapply(file_list, function(f) {
    # get current data  
    data = read.csv(f)   
    uid  = data$dataid[1]
    data = data[,-1]  
    return(data)
  })
  
  if (length(data_list) == 1) return(data_list)
  
  # make sure all end-uses are aligned
  all_cols = unique(unlist(lapply(data_list, names)))
  data = data_list[[1]]
  for (i in 2:length(data_list)) {
    tmp = data_list[[i]]
    cur_cols = setdiff(names(data), names(tmp))
    new_cols = setdiff(names(data_list[[i]]), names(data))
    if (length(cur_cols)>0)  {
      tmp[,cur_cols] = NA
    }
    if(length(new_cols)>0) {
      data[,new_cols]= NA;
    }
    data = rbind(data, tmp)
  }
  
  return(data)
}

# ------------------------------------------
# Get files to be processed
# ------------------------------------------

files.proc  = list.files(OUT_PATH, full.names = TRUE, recursive = T)
uids.done   = unique(sapply(files.proc, function(x) {
  tmp = strsplit(x, '/')[[1]]
  uid = tmp[length(tmp)]
  uid = strsplit(uid, '_')[[1]][1]
  return(uid)
}))
files.input = list.files(IN_PATH, full.names = TRUE, recursive = T)
uids.input  = unique(sapply(files.input, function(x) {
  tmp = strsplit(x, '/')[[1]]
  uid = tmp[length(tmp)]
  uid = strsplit(uid, '\\.')[[1]][1]
  return(uid)
}))

idx = which(!(uids.input %in% uids.done))

if (length(idx) == 0) stop("All files have already been processed!")

# ------------------------------------------
# Aggregate data in given directory
# ------------------------------------------

res      = mclapply(1:length(uids.input[idx]), 
                    mc.cores = 6, 
                    function(i) {       
  uid = uids.input[idx[i]]  
  cat(paste('Processing uid ', uid, ' : ', i, '/', length(files.input[idx]), '\n', sep = ''))
  files = files.input[grep(paste('/', uid, '.csv', sep=''), files.input)]
  
  # put together data from multiple years
  data = concatenate_data(files)  
  
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
  write.csv(data_01_proc, file = paste(OUT_PATH, paste(uid, '_minute.csv', sep = ''), sep = '/'), row.names = F)
  write.csv(data_15_proc, file = paste(OUT_PATH, paste(uid, '_15mins.csv', sep = ''), sep = '/'), row.names = F)
  write.csv(data_60_proc, file = paste(OUT_PATH, paste(uid, '_hourly.csv', sep = ''), sep = '/'), row.names = F)
  
  return(0)
})  