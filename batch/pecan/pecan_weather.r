# pecan_weather.r
#
# Obtain, parse, and integrate weather data for Pecan Street users. 
# 
# Adrian Albert
#
# Last modified: May 2014.

# -------------------------------
# Initializations ... 
# -------------------------------

rm(list=ls())

options(error = recover)

library('ggplot2')
library('reshape')
library('VIM')
library('weatherData')

source('~/EnergyAnalytics/utils/weather//wunderground.R')
source('~/EnergyAnalytics/utils/aggregate_data.r')
source('~/EnergyAnalytics/utils/weather/clean_weather_data.r')

# directory to save files
OUT_PATH  = '~/energy-data/pecan_street/weather/'
dir.create(file.path(OUT_PATH))    

# ------------------------------------------
# Get weather data for 2012-2014 for Austin
# ------------------------------------------

station='KTXAUSTI90' # this is Austin, TX
data  = getWUData(station, '2012-01-01', '2014-04-01')

# time average to hourly data
weather_15m = toXmin(data, dateCol = 'Time', min = 15)
weather_60m = toHourly(data, dateCol = 'Time')
  
write.csv(weather_15m, row.names = F, file = paste(OUT_PATH, '/', 'weather_15mins.csv', sep=''))
write.csv(weather_60m, row.names = F, file = paste(OUT_PATH, '/', 'weather_hourly.csv', sep=''))

# ------------------------------------------
# Clean weather data
# ------------------------------------------

# weather_15m = read.csv('~/energy-data/pecan_street/weather/weather_15mins.csv')
# weather_60m = read.csv('~/energy-data/pecan_street/weather/weather_hourly.csv')

wthr_15_clean = clean_weather_data(weather_15m, dateCol = 'date', addToD = F)
wthr_60_clean = clean_weather_data(weather_60m, dateCol = 'date', addToD = F)

write.csv(wthr_15_clean, row.names = F, file = paste(OUT_PATH, '/', 'weather_15mins.csv', sep=''))
write.csv(wthr_60_clean, row.names = F, file = paste(OUT_PATH, '/', 'weather_hourly.csv', sep=''))
