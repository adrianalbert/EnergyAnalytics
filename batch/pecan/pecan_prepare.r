# pecan_prepare.r
#
# Prepares pecan data for statistical analysis.
#
# Adrian Albert
# Last modified: June 2014.
# ---------------------------------------------------------

rm(list = ls())
options(error = recover)

library('parallel')

# __________________________________________________
# Initializations...

setwd('~/EnergyAnalytics/batch/')
source('../utils/aggregate_data.r')
source('./pecan/define_categories_pecan.r')

# select one year of data to use
DATA_PATH     = '~/S3L_server/energy-data/pecan_street/usage-processed/'
DATA_OUT_PATH = '~/S3L_server/energy-data/pecan_street/usage-select'
dir.create(DATA_OUT_PATH)

# load weather data
weather.hourly = read.csv('~/S3L_server/energy-data/pecan_street/weather/weather_hourly.csv')
weather.15mins = read.csv('~/S3L_server/energy-data/pecan_street/weather/weather_15mins.csv')
weather.hourly$date = as.POSIXct(as.character(weather.hourly$date))
weather.15mins$date = as.POSIXct(as.character(weather.15mins$date))

# load user names
user_names = read.csv('~/S3L_server/energy-data/pecan_street/metadata/user_names_ids.csv')

# __________________________________________________
# Access usage data

# list all data files by year/uid
files.proc  = list.files(DATA_OUT_PATH, full.names = TRUE, recursive = T)
uids.done   = unique(sapply(files.proc, function(x) {
  tmp = strsplit(x, '/')[[1]]
  uid = tmp[length(tmp)]
  uid = strsplit(uid, '\\.')[[1]][1]
  return(uid)
}))
files.input = list.files(DATA_PATH, full.names = TRUE, recursive = T)
uids.input  = unique(sapply(files.input, function(x) {
  tmp = strsplit(x, '/')[[1]]
  uid = tmp[length(tmp)]
  uid = strsplit(uid, '_')[[1]][1]
  return(uid)
}))

# only process those files not processed before
idx.ok   = which(!(uids.input %in% uids.done))
if (length(idx.ok) == 0) stop("All files have already been processed!")

# select files with appropriate resolution 
files_01 = files.input[grep('minute',files.input)]
files_15 = files.input[grep('15mins',files.input)]
files_60 = files.input[grep('hourly', files.input)]

# ______________________________________________________
# Combine end uses by functionality and add in weather

# function to add 
to_categories = function(homeData, dateCol = 'date') {
  names(homeData) = tolower(names(homeData))
  
  # what kind of data does this user have?
  hv.cols       = intersect(names(homeData), select.HV)
  ac.cols       = intersect(names(homeData), select.AC)
  user.cols     = intersect(names(homeData), select.user)
  alwOn.cols    = intersect(names(homeData), select.alwOn)
  sched.cols    = intersect(names(homeData), select.sched)
  light.cols    = intersect(names(homeData), select.light)
  
  # if no HVAC data, not interesting
  if (length(hv.cols) == 0 && length(ac.cols) == 0) return(NULL)
  if (!(select.total %in% names(homeData))) return(NULL)
  
  # aggregate components
  AC_kwh        = add.columns(homeData, select.AC)
  HV_kwh        = add.columns(homeData, select.HV)
  total_kwh     = add.columns(homeData, select.total)
  user_kwh      = add.columns(homeData, select.user)
  light_kwh     = add.columns(homeData, select.light)
  sched_kwh     = add.columns(homeData, select.sched)
  alwOn_kwh     = add.columns(homeData, select.alwOn)
  nonHVAC_kwh   = total_kwh - AC_kwh - HV_kwh  
  
  # create new dataset
  data_new = subset(homeData, select = c(dateCol, select.total))
  names(data_new)[names(data_new)==select.total[1]] <- "total"
  if (length(ac.cols)>0) data_new = cbind(data_new, C = AC_kwh)
  if (length(hv.cols)>0) data_new = cbind(data_new, H = HV_kwh)
  if (length(user.cols)>0) data_new = cbind(data_new, user = user_kwh)
  if (length(light.cols)>0) data_new = cbind(data_new, lights = light_kwh)
  if (length(alwOn.cols)>0) data_new = cbind(data_new, always_on = alwOn_kwh)
  if (length(sched.cols)>0) data_new = cbind(data_new, scheduled = sched_kwh)
  data_new = cbind(data_new, nonHVAC = nonHVAC_kwh)
  
  return(data_new)
}

# ------------------------------------------
# Concatenate data across multiple years
# ------------------------------------------

concatenate_data = function(data_list, dateCol = 'date') {
  
  if (length(data_list) == 1) return(data_list[[1]])
  
  # make sure all end-uses are aligned
  all_cols = unique(unlist(lapply(data_list, names)))
  i = 1
  while (is.null(data_list[[i]]) & i < length(data_list)) i = i + 1; 
  if (i<=length(data_list)) data = data_list[[i]] else return(NULL)
  for (i in 2:length(data_list)) {
    tmp = data_list[[i]]
    if(is.null(tmp)) next
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
  
  if (!is.null(data)) data[,dateCol] = as.POSIXct(data[,dateCol])
  return(data)
}

# ------------------------------------------
# Prepare data for analysis
# ------------------------------------------

dir.create(file.path(DATA_OUT_PATH, '01min'), recursive = T)
dir.create(file.path(DATA_OUT_PATH, '15min'), recursive = T)
dir.create(file.path(DATA_OUT_PATH, '60min'), recursive = T)

res = mclapply(1:length(uids.input[idx.ok]),
            mc.cores = 6,
            function(i) {

# for (i in 1:length(uids.input[idx.ok])) {
  # current user info
  uid = uids.input[idx.ok[i]]  
  user_name = as.character(user_names$name[which(user_names$ID == uid)])
  cat(paste('Processing uid ', uid, ' : ', i, '/', length(idx.ok), '\n', sep = ''))
              
  # some gymnastics to get appropriate data files for this user  
  #files_01_i = files_01[grep(paste('/', uid, '_minute.csv', sep=''), files_01)]
  #files_15_i = files_15[grep(paste('/', uid, '_15mins.csv', sep=''), files_15)]
  files_60_i = files_60[grep(paste('/', uid, '_hourly.csv', sep=''), files_60)]

  # read in data at different resolutions across years...
  # aggregate individual appliances into categories of interest
  # also flag users who don't have HVAC data for removal
#   data_list_01 = lapply(files_01_i, function(f) {
#     # get current data  
#     data = read.csv(f)   
#     data = to_categories(data)
#     return(data)
#   })
#   data_list_15 = lapply(files_15_i, function(f) {
#     # get current data  
#     data = read.csv(f)   
#     data = to_categories(data)
#     return(data)
#   })
  data_list_60 = lapply(files_60_i, function(f) {
    # get current data  
    data = read.csv(f)   
    data = to_categories(data)
    return(data)
  })
    
  # concatenate data across years
#   data_01_sel = concatenate_data(data_list_01)
#   data_15_sel = concatenate_data(data_list_15)
  data_60_sel = concatenate_data(data_list_60)  
  
  # if the user has no HVAC data, skip
  # if (is.null(data_01_sel) || is.null(data_15_sel) || is.null(data_60_sel)) {
  if (is.null(data_60_sel)) {
    cat('No HVAC data for this user!\n')
    return(NULL)   
  }  
  
  # skip users that don't have enough data, at least a month of hourly data
  if (nrow(data_60_sel) < 30 * 24) {
    cat('Not enough samples (less than one month!) !\n')
    return(NULL)
  }  
  
  # add in weather (temperature) data at the resolutions available
  # data_15_sel = merge(data_15_sel, subset(weather.15mins, select = c('date', 'TemperatureF')), by = 'date')
  data_60_sel = merge(data_60_sel, subset(weather.hourly, select = c('date', 'TemperatureF')), by = 'date')  
  names(data_60_sel)[names(data_60_sel)=="TemperatureF"] <- "Temperature"
  
  # skip users that don't have enough data, at least a month of hourly data
  if (nrow(data_60_sel) < 30 * 24) {
    cat('Not enough samples (less than one month!) !\n')
    return(NULL)
  }  
  
  # save to disc  
#   write.csv(data_01_sel, file = paste(DATA_OUT_PATH, paste('01min/', uid, '.csv', sep = ''), sep = '/'), row.names = F)
#   write.csv(data_15_sel, file = paste(DATA_OUT_PATH, paste('15min/', uid, '.csv', sep = ''), sep = '/'), row.names = F)
  write.csv(data_60_sel, file = paste(DATA_OUT_PATH, paste('60min/', uid, '.csv', sep = ''), sep = '/'), row.names = F)    
  
  return(1)
})
#}