# pecan_weather.r
#
# Obtain, parse, and integrate weather data for Pecan Street users. 
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
library('parallel')
library('VIM')
library('weatherData')

source('~/EnergyAnalytics/utils/weather//wunderground.R')
source('~/EnergyAnalytics/utils/aggregate_data.r')

# directory to save files
DATA_PATH = '~/energy-data/pecan_street/metadata/'
OUT_PATH  = '~/energy-data/pecan_street/weather/'
dir.create(file.path(OUT_PATH))    

# ------------------------------------------
# Function to compute distance on Earth
# ------------------------------------------

# Points are of the form (lon,lat)

geodetic.distance <- function(point1, point2) 
{ 
  R <- 6371 
  p1rad <- point1 * pi/180 
  p2rad <- point2 * pi/180 
  d <- sin(p1rad[2])*sin(p2rad[2])+cos(p1rad[2])*cos(p2rad[2])*cos(abs(p1rad[1]-p2rad[1]))  
  d <- acos(d) 
  R*d 
} 

# ------------------------------------------
# Read in location information on users
# ------------------------------------------

# location data is available only for the 2012 survey!
data_2012 = read.csv(paste(DATA_PATH, 'survey_2012.csv', sep = '/'))
data_2012 = subset(data_2012, select = c('dataid', 'longitude', 'latitude', 'zip_code'))

# get closest stations to each user
airports = mclapply(1:nrow(data_2012), 
                    mc.cores = 5,
                    function(i) {
  lat = data_2012$latitude
  lon = data_2012$longitude
  dis = sapply(1:nrow(USAirportWeatherStations), function(j) {
      lat.1 = USAirportWeatherStations[i,'Lat']
      lon.1 = USAirportWeatherStations[i,'Lon']
      d     = geodetic.distance(c(lon, lat), c(lon.1, lat.1))
  })
  # idx = order(dis)[1:5]
  idx = which.min(dis)
  return(USAirportWeatherStations[idx, 'airportCode'])
})
stations = data.frame(dataid = data_2012$dataid, airport = airports)

# > table(stations$airport)
# 
# PARL 
# 163 
# looks like all available users are closest to the PARL airport station!
# however PARL has no data!
# let's use the KATT airport station which is in Austin, TX

# ------------------------------------------
# Get weather data for 2012-2014 for Austin
# ------------------------------------------

station='KTXAUSTI90' # this is ausitn Tx
data  = getWUData(station, '2012-01-01', '2014-03-01')

# time average to hourly data
weather_15m = toXmin(data, dateCol = 'Time', min = 15)
weather_60m = toHourly(data, dateCol = 'Time')
  
write.csv(weather_15m, row.names = F, file = paste(OUT_PATH, '/', 'weather_15mins.csv', sep=''))
write.csv(weather_60m, row.names = F, file = paste(OUT_PATH, '/', 'weather_hourly.csv', sep=''))
