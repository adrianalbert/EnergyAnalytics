#!/usr/bin/Rscript
# 
# align_weather_data.r
#
# Align weather data to 15 min/hourly grid. Perform minimal cleaning.
#
# Adrian Albert
#
# Last modified: February 2012.

rm(list=ls())

options(error = recover)

raw_path    = '~/EnergyAnalytics/data/pge_data/weather/raw/'
path_15min  = '~/EnergyAnalytics/data/pge_data/weather/aligned/15_min/'
path_hourly = '~/EnergyAnalytics/data/pge_data/weather/aligned/hourly/'
plots_raw_path  = '~/Dropbox/ControlPatterns/plots/weather/raw/'
plots_path  = '~/Dropbox/ControlPatterns/plots/weather/aligned/'

# create directories specified in profiles if they don't already exist
if (!file.exists(path_15min)) dir.create(path_15min, recursive = TRUE)
if (!file.exists(path_hourly)) dir.create(path_hourly, recursive = TRUE)
if (!file.exists(plots_raw_path)) dir.create(plots_raw_path, recursive = TRUE)
if (!file.exists(plots_path)) dir.create(plots_path, recursive = TRUE)

# get station lat/long
stations = read.csv('~/EnergyAnalytics/data/metadata/station_lat_long_elev.csv')

THRESH_OBS = 30 * 24 * 4 # only consider stations with at least one month of observations

# _______________________________________________________
# Function for preliminary cleaning of the data

# Only very beningn modifications - leave data as intact as possible for subsequent analysis.

prelim_cleaning = function(data) {
  # clean up weather covariates (remove unreasonable values)
  ranges = list(TemperatureF = c(2, 125), 
                Pressure     = c(27, 33),
                Humidity     = c(1, 100), 
                HourlyPrecip = c(0, Inf), 
                WindSpeed    = c(0, 70), 
                DewpointF    = c(0, 125))
  for (r in intersect(names(ranges),names(data))){
    idx = which(data[,r]<ranges[[r]][1] | data[,r]>ranges[[r]][2]) 
    if (length(idx)>0) data[idx,r] = NA
  }
  
  # remove covariates that have less than 1/2 of the values
  na.number = sapply(names(data), function(v) {
    return(length(which(is.na(data[,v]))))
  })
  idx = which(na.number > nrow(data) / 2)                   
  if (length(idx) > 0) data = data[,-idx]              
  
  # remove covariates that have no variation
  sd_col   = apply(data[,-1], 2, function(x) sd(x,na.rm=T))
  rm.vars  = which(sd_col == 0 | is.na(sd_col))
  if (length(rm.vars)>0) data = data[,-(rm.vars+1)]
  
  # remove any duplicate time stamps
  dup = which(duplicated(data$TimeDate))
  if (length(dup)>0) data = data[-dup,]  
  
  return(data)
}

# _______________________________________________________
# Function to align every time stamp to a 15-minute grid
# This is a simple nearest-neighbor interpolation

library(zoo) 
align_timestamps = function(data) {
  
  # snap to 15-minute grid (nearest-neighbor)
  timestamps_raw = as.POSIXct(as.character(data$TimeDate))
  timestamps_60m = as.POSIXct(strftime(timestamps_raw, format="%Y-%m-%d %H:00:00"))
  minutes_vec    = strftime(timestamps_raw, format="%M")
  hours_vec      = strftime(timestamps_raw, format="%Y-%m-%d %H")
  snap.minutes   = sapply(as.numeric(minutes_vec), function(m) {
    min.int = c(0,15,30,45,60)
    # min.str = c('00', '15', '30', '45', '60')
    idx     = which.min(abs(m-min.int))
    return(min.int[idx])
  })
  #   snap.15min.time= as.POSIXct(paste(hours_vec, snap.minutes, '00', sep=':'))
  snap.15min.time= timestamps_60m + snap.minutes * 60
  data$TimeDate  = snap.15min.time
  
  # do 15-min averages
  data_15m       = aggregate(.~TimeDate, data, FUN = mean, na.action = na.pass)
  
#   # define uniform time axis
#   timestamps_15m = as.POSIXct(seq(min(timestamps_60m), max(timestamps_60m)+3600, by = 15 * 60))
#   tmp            = data.frame(TimeDate = timestamps_15m)
#   data_15m       = merge(data_15m, tmp, by = 'TimeDate', all = T)
  
  # aggregate hourly dataset
  data_hour          = data_15m
  data_hour$TimeDate = as.POSIXct(strftime(data_hour$TimeDate, format="%Y-%m-%d %H:00:00"))
  data_hour          = aggregate(.~TimeDate, data_hour, FUN = mean, na.action = na.pass)
    
  return(list(min15 = data_15m, hourly = data_hour))
}

# _______________________________________________________
# Function to compute quality stats for current data

compute_stats = function(data, this_station) {
  
  # how many NA's are there per covariate?
  stats = sapply(names(data)[-1], function(v) length(which(is.na(data[,v]))))

  df = data.frame(Station = this_station, 
                  StartDate = data$TimeDate[1], 
                  EndDate = data$TimeDate[nrow(data)],
                  NoObs   = nrow(data))
  df = cbind(df, as.data.frame(t(stats)))
  
  return(df)
}

# __________________________________________________________________
# Function to align timestamps & minimal cleaning per station

wrapper = function(i) {
  
  st = stations$Station[i]
  cat(paste('Processing Station', st, ':', i, '/', nrow(stations), '\n'))
  
  # load data from appropriate station file
  load_path = paste(raw_path, st, '.csv', sep='')
  if (!file.exists(load_path)) {
    cat(paste('Data for station', st, 'not found!\n'))
    return(NULL)
  }
  data     = read.csv(pipe(paste("cut -f1,2,3,4,7,9,10 -d,", load_path)))
  if (nrow(data) < THRESH_OBS ) {
    cat(paste('--> Not enough data at current station!\n'))
    return(NULL)
  }
  names(data)[1] = 'TimeDate'
  
  # Select covariates of interest & rename with simpler names
  wthr_names         = c('TimeDate', 'TemperatureF', 'DewpointF', 'Pressure', 'WindSpeed', 'Humidity', 
                         'HourlyPrecip')  
  names(wthr_names)  = c('TimeDate', 'TemperatureF', 'DewpointF', 'PressureIn', 'WindSpeedMPH', 'Humidity', 
                         'HourlyPrecipIn')  
  data        = data[,which(names(data) %in% names(wthr_names))]
  names(data) = wthr_names[names(data)]
  
  # perform minimal data cleaning 
  data = prelim_cleaning(data)

  # plot aligned time series 
  tmp = zoo(data[,-1], order.by = data[,1])
  png(paste(plots_raw_path, st,'_raw.png', sep = ''), height = 600, width = 1600)
  plot(tmp, main = paste(st,'(raw)'))
  dev.off()
  
  # align timestamps to a 15-min grid
  res  = align_timestamps(data)
  
  # compute stats about quality
  stats_15 = compute_stats(res$min15, st)
  stats_60 = compute_stats(res$hourly, st)
  stats    = list(min15 = stats_15, min60 = stats_60)
  
  # save data to file
  write.csv(res$min15,  file = paste(path_15min, st, '.csv', sep=''), row.names = F)
  write.csv(res$hourly, file = paste(path_hourly, st, '.csv', sep=''), row.names = F)
  
  # plot aligned time series 
  tmp15 = zoo(res$min15[,-1], order.by = res$min15[,1])
  png(paste(plots_path, st,'_15min.png', sep = ''), height = 600, width = 1600)
  plot(tmp15, main = paste(st,'(15 min)'))
  dev.off()
  tmp60 = zoo(res$hourly[,-1], order.by = res$hourly[,1])
  png(paste(plots_path, st,'_60min.png', sep = ''), height = 600, width = 1600)
  plot(tmp60, main = paste(st,'(15 min)'))
  dev.off()
  
  rm(list = c('data', 'res'))
  return(stats)
}

# _________________________
# Iterate through stations

library('parallel')

# result = list()
# for (i in 1:3) {
#   result[[i]] = wrapper(i)
# }

result = mclapply(1:length(stations$Station), FUN = wrapper, mc.silent = F, 
                  mc.preschedule = FALSE, mc.cores = 8)

# _____________________________
# Save stats to RData file

stats_15   = data.frame()
stats_60   = data.frame()
stat_names = c('Station', 'StartDate', 'EndDate', 'NoObs', 
               'TemperatureF', 'DewpointF', 'Pressure', 'WindSpeed', 'Humidity', 'HourlyPrecip')
for (res in result) {  
  df = res$min15
  df[,setdiff(stat_names, names(df))] = NA 
  if (nrow(stats_15)==0) stats_15 = df else stats_15 = rbind(stats_15, df)
  df = res$min60
  df[,setdiff(stat_names, names(df))] = NA 
  if (nrow(stats_60)==0) stats_60 = df else stats_60 = rbind(stats_60, df)
}

write.csv(stats_15,   
          file = '~/EnergyAnalytics/data/pge_data/weather/station_stats_15min.csv', row.names = F)
write.csv(stats_60,   
          file = '~/EnergyAnalytics/data/pge_data/weather/station_stats_60min.csv', row.names = F)

