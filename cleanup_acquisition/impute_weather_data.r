# impute_weather_data.r
#
# Perform imputation on weather data.
#
# Adrian Albert
#
# Last modified: February 2012.

rm(list=ls())

options(error = recover)

setwd('~/EnergyAnalytics/code/R')
source('~/EnergyAnalytics/code/R/utils/timing.r')
source('~/EnergyAnalytics/code/R/utils/sql_utils.r')

input_path   = '~/EnergyAnalytics/data/pge_data/weather/aligned/hourly/'
output_path  = '~/EnergyAnalytics/data/pge_data/weather/imputed/hourly/'
plots_path   = '~/Dropbox/ControlPatterns/plots/weather_imputed/'

# create directories specified in profiles if they don't already exist
if (!file.exists(plots_path)) dir.create(plots_path, recursive = TRUE)
if (!file.exists(output_path)) dir.create(output_path, recursive = TRUE)

# _____________________________
# Preliminary setup

library('zipcode')
data(zipcode)

# get unique zipcodes
zips = read.csv('~/EnergyAnalytics/data/metadata/zipcode_res.csv')
zips$ZIP_ID = 1:nrow(zips)

# get station lat/long
stations = read.csv('~/EnergyAnalytics/data/metadata/station_lat_long_elev.csv')

# _______________________________________
# Analyze data quality across stations

all_files = list.files(path = input_path, pattern = '.csv', all.files = FALSE)
station_quality = data.frame()
for (file in all_files) {
  cur_data     = read.csv(paste(input_path, file, sep=''))
  
  # some general info on the station data
  this_station = strsplit(file, '\\.')[[1]][1]
  df = data.frame(Station = this_station, 
                  StartDate = cur_data$TimeDate[1], 
                  EndDate = cur_data$TimeDate[nrow(cur_data)])
  
  # how many NA's are there per covariate?
  summary(cur_data)
  
}

# _____________________________
# Compute distance matrices

library('sp')

# zip-zip distance
zips.gis            = subset(zipcode, zip %in% zips$ZIP5)
zip.lat.long        = as.matrix(zips.gis[,c('latitude', 'longitude')])
zips.dist           = spDists(zip.lat.long, zip.lat.long)
rownames(zips.dist) = zips.gis$zip
colnames(zips.dist) = zips.gis$zip

# station-station distance
station.lat.long       = as.matrix(stations[,c('Lat', 'Long')])
station.dist           = spDists(station.lat.long, station.lat.long)
rownames(station.dist) = stations$Station
colnames(station.dist) = stations$Station

# zip-station distance
zip.st.dist           = spDists(zip.lat.long, station.lat.long)
rownames(zip.st.dist) = zips.gis$zip
colnames(zip.st.dist) = stations$Station

# _______________________________________
# Impute missing values for a given zip

impute_values = function(zip, time_vec) {
  
  for (var in )
  
}

# _______________________________________
# Iterate over zipcodes

for (zip in zips.gis$zip) {
  
  # construct required time vector for current zipcode
  info     = subset(zips, ZIP5 == zip)
  times_60 = seq(as.POSIXct(as.character(info$MIN.sm_start_date.)), as.POSIXct(as.character(info$MAX.sm_end_date.)), by = 3600) 
  times_15 = seq(as.POSIXct(as.character(info$MIN.sm_start_date.)), as.POSIXct(as.character(info$MAX.sm_end_date.)), by = 15*60) 
  
  # build large matrix of all stations per covariate
  
  
  # find nearest stations to each zipcode in a 3-miles radius
  st_radius  = zip.st.dist[zip, zip.st.dist[zip,] < 1]
  station.ok = colnames(zip.st.dist)[order(st_radius)]
  distances  = zip.st.dist[zip, station.ok]
  elevations = subset(stations, Station %in% station.ok & zip == zip)$Elev
  
  # build regression model for each covariate
  
  
}

