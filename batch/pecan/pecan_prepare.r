# pecan_prepare.r
#
# Prepares pecan data for statistical analysis.
#
# Adrian Albert
# Last modified: May 2014.
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
DATA_PATH_ALL = '~/energy-data/pecan_street/usage-orig/'
DATA_PATH     = '~/energy-data/pecan_street/usage-processed/'
METADATA_PATH = '~/energy-data/pecan_street/metadata/'
DATA_OUT_PATH = '~/energy-data/pecan_street/usage-select'

# load weather data
weather.hourly = read.csv('~/energy-data/pecan_street/weather/weather_hourly.csv')
weather.15mins = read.csv('~/energy-data/pecan_street/weather/weather_15mins.csv')
weather.hourly$date = as.POSIXct(as.character(weather.hourly$date))
weather.15mins$date = as.POSIXct(as.character(weather.15mins$date))

# __________________________________________________
# Assign each user a unique human-readable name

files_all = list.files(path=DATA_PATH_ALL, full.names = T, recursive = T)
users_all = unique(sapply(files_all, function(s) strsplit(tail(strsplit(s, '/')[[1]], 1), '\\.')[[1]][1]))

# load baby names
# we're going to name each user for later easiness of use
baby_names  = read.csv('~/Dropbox/OccupancyStates/data/baby-names.csv')
user_names  = data.frame(ID = users_all, name = unique(baby_names$name)[1:length(users_all)])

# save names to file
write.csv(user_names, file = paste(METADATA_PATH, 'user_names_ids.csv', sep = '/'))

# __________________________________________________
# Access usage data

# list all data files by year/uid
files.input = list.files(path=DATA_PATH, full.names = T, recursive = T)
uids.input  = sapply(files.input, function(x) {
  tmp = strsplit(x, '/')[[1]]
  uid = tmp[length(tmp)]
  uid = strsplit(uid, '_')[[1]][1]
  yr  = tmp[length(tmp)-1]
  return(paste(yr, uid, sep = '/'))
})

files.proc  = list.files(DATA_OUT_PATH, full.names = TRUE, recursive = T)
uids.done   = sapply(files.proc, function(x) {
  tmp = strsplit(x, '/')[[1]]
  uid = tmp[length(tmp)]
  uid = strsplit(uid, '\\.')[[1]][1]
  yr  = tmp[length(tmp)-1]
  return(paste(yr, uid, sep = '/'))
})
uids.done = unique(uids.done)

# only process those files not processed before
idx.ok   = which(!(uids.input %in% uids.done))
files    = files.input[idx.ok]

# select files with appropriate resolution 
files_01 = files[grep('minute',files)]
files_15 = files[grep('15mins',files)]
files_60 = files[grep('hourly', files)]

# select those users for which enough data is available
usersVec = data.frame(UID = sapply(files_15, function(s) strsplit(tail(strsplit(s, '/')[[1]], 1), '_')[[1]][1]),
                      year= sapply(files_15, function(s) tail(strsplit(s, '/')[[1]], 2)[1]))

# ______________________________________________________
# Combine end uses by functionality and add in weather

#!!! TODO !!!
# Get implement "rules" to determine how to assign appliances to different end-uses
# e.g., if it's summer the HV probably is not on etc.

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
  if (length(ac.cols)>0) data_new = cbind(data_new, AC = AC_kwh)
  if (length(hv.cols)>0) data_new = cbind(data_new, HV = HV_kwh)
  if (length(user.cols)>0) data_new = cbind(data_new, user = user_kwh)
  if (length(light.cols)>0) data_new = cbind(data_new, lights = light_kwh)
  if (length(alwOn.cols)>0) data_new = cbind(data_new, always_on = alwOn_kwh)
  if (length(sched.cols)>0) data_new = cbind(data_new, scheduled = sched_kwh)
  data_new = cbind(data_new, nonHVAC = nonHVAC_kwh)
  
  return(data_new)
}

res = mclapply(1:nrow(usersVec),
               mc.cores = 5,
            function(i) {

  # current user info
  user   = user_names[which(user_names$ID == usersVec[i,'UID']),]
  year   = usersVec[i,'year']
  # select those users for which enough data is available
  
  dir.create(file.path(DATA_OUT_PATH, '01min', year), recursive = T)
  dir.create(file.path(DATA_OUT_PATH, '15min', year), recursive = T)
  dir.create(file.path(DATA_OUT_PATH, '60min', year), recursive = T)
                                                  
  # some gymnastics to get appropriate data files for this user  
  files_01_i = files_01[grep(paste(year, '/', user$ID, '_minute.csv', sep=''), files_01)]
  files_15_i = files_15[grep(paste(year, '/', user$ID, '_15mins.csv', sep=''), files_15)]
  files_60_i = files_60[grep(paste(year, '/', user$ID, '_hourly.csv', sep=''), files_60)]

  # read in data at different resolutions...
  cat(paste(user$ID, '...\n'))
  data_01 = read.csv(files_01_i); data_01$date = as.POSIXct(data_01$date); 
  data_15 = read.csv(files_15_i); data_15$date = as.POSIXct(data_15$date);
  data_60 = read.csv(files_60_i); data_60$date = as.POSIXct(data_60$date);
  
  # skip users that don't have enough data, at least a month of hourly data
  if (nrow(data_60) < 30 * 24) {
    cat('Not enough samples (less than one month!) !\n')
    return(NULL)
  }
  
  # aggregate individual appliances into categories of interest
  # also flag users who don't have HVAC data for removal
  data_01_sel = to_categories(data_01)
  data_15_sel = to_categories(data_15)
  data_60_sel = to_categories(data_60)
  
  # if the user has no HVAC data, skip
  if (is.null(data_01_sel) || is.null(data_15_sel) || is.null(data_60_sel)) {
    cat('No HVAC data for this user!\n')
    return(NULL)   
  }
  
  # add in weather (temperature) data at the resolutions available
  data_15_sel = merge(data_15_sel, subset(weather.15mins, select = c('date', 'TemperatureF')), by = 'date')
  data_60_sel = merge(data_60_sel, subset(weather.hourly, select = c('date', 'TemperatureF')), by = 'date')  
  
  # skip users that don't have enough data, at least a month of hourly data
  if (nrow(data_60) < 30 * 24) {
    cat('Not enough samples (less than one month!) !\n')
    return(NULL)
  }  
  
  # save to disc  
  write.csv(data_01_sel, file = paste(DATA_OUT_PATH, paste('01min/', year, '/', user$ID, '.csv', sep = ''), sep = '/'), row.names = F)
  write.csv(data_15_sel, file = paste(DATA_OUT_PATH, paste('15min/', year, '/', user$ID, '.csv', sep = ''), sep = '/'), row.names = F)
  write.csv(data_60_sel, file = paste(DATA_OUT_PATH, paste('60min/', year, '/', user$ID, '.csv', sep = ''), sep = '/'), row.names = F)    
  
  return(1)
})
