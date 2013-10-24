#!/usr/bin/Rscript

# analysis_simulated.r
#
# HMM analysis of simulated building data.
# 
# - 
# Adrian Albert
# Last modified: May 2013.

# -----------------------------------------------
# Define constants, load libraries and functions
# -----------------------------------------------

# load classes, utils
rm(list = ls())

options(error = recover)

setwd('~/Dropbox/OccupancyStates/')
source('code/OccupancyStates.r')
source('code/occupancyAnalysis.r')

plots_path_root = 'plots/simulation/'

# ------------------------------------------
# Load simulation data
# ------------------------------------------

# ____________________
# Consumption data

csv.files = list.files(path = 'data/simulation', pattern = '*.csv', full.names = T)
data = list()
for (f in csv.files) {
  cat(paste('Loading simulation data from file ', f, '\n'))
  # read from file
  data[[f]] = read.csv(file = f)
  # convert J to kWh
  data[[f]][,-1] = data[[f]][,-1] / (3600 * 1000)
  # process timestamp
  data[[f]]$Date.Time = paste('2010', as.character(data[[f]]$Date.Time), sep = '/')
  data[[f]]$Date.Time = as.POSIXct(data[[f]]$Date.Time)
  names(data[[f]])    = c('Date.Time', 'Gas.HVAC', 'kWh.HVAC', 'kWh', 'Gas', 'Cooling')
  data[[f]]$kWh       = data[[f]]$kWh
  # aggregate to hourly data
  data[[f]]$Date.Time = paste(substr(data[[f]]$Date.Time, 1, 13), '00', sep = ':')
  data[[f]]$Date.Time = as.POSIXct(data[[f]]$Date.Time)
  data[[f]]           = aggregate(data[[f]][,-1], by = list(data[[f]]$Date.Time), FUN=sum, na.rm=TRUE)
  names(data[[f]])[1] = 'date'
  data[[f]]$date = as.character(data[[f]]$date)
}

# __________________
# Weather data

# read in data
weather = read.csv(file = 'data/simulation/CZ01RV2.epw', header = F, skip = 8)
header = read.csv('data/simulation/epw_specification.txt')

# epw files have 35 columns. The first 5 are date and time related.
# The 6th is the crazy metadata, 7th is dry bulb, 9th RH, and the rest are other weather stats
# the year 2010 is hard coded to give the simulation a specific REAL timeframe
dates  = paste('2010', weather[,2], weather[,3], sep = '-')
times  = paste(weather[,4], '00:00', sep = ':')
dates_times = paste(dates, times)
names(weather) = header$EpWColumn
weather     = data.frame(date = as.character(as.POSIXct(dates_times)), TemperatureF = (weather[,7] * 9/5+32), Humidity = weather[,9])

# save data for later processing
save(file = 'data/simulation/simulation_processed.RData', list = c('data', 'weather'))

# ------------------------------------------------
# Run analysis script
# ------------------------------------------------

cat('****** Running analysis on simulated data ******\n')

source('code/OccupancyStates.r')
source('code/occupancyAnalysis.r')

names(data) = c('CZ011-No-Schedule-NoSB-NoGains',
                'CZ012-No-Schedule-2CSB-NoGains',
                'CZ013-No-Schedule-4CSB-NoGains',
                'CZ014-Schedule-7to6Wkdys-NoSB-NoGains',
                'CZ015-Schedule-7to6Wkdys-4CSB-NoGains',
                'CZ016-Schedule-7to6Wkdys-4CSB-Gains',
                'CZ017-No-Schedule-NoSB-Gains')
simulation = list()
for (i in names(data)[-c(4,5)]) {
  
  # create folder for each test case
  plots_path = paste(plots_path_root, i, '/', sep='')
  if (!file.exists(plots_path)) dir.create(plots_path, recursive = TRUE)  
  
  cur_data = data[[i]][, c('date', 'kWh')]
  cur_wthr = weather
  simulation[[i]] = occupancyAnalysis(cur_data, cur_wthr, i, 1, 
                                      Kmin = 3, Kmax = 5, 
                                      verbose = T, 
                                      thresh.R2 = 0.9,
                                      resp.vars = c('(Intercept)', 'TemperatureF'),
                                      tran.vars = c('(Intercept)', 'TemperatureF'),
                                      addl.vars = c(),
                                      plots_path = plots_path, dump_path = plots_path)  
}


# ------------------------------------------------
# Produce plots
# ------------------------------------------------

cat('****** Running analysis on simulated data ******\n')
for (hid in names(data)[-c(4,5)]) {
  # create folder for each test case
  plots_path = paste(plots_path_root, hid, '/', sep='')
  if (!file.exists(plots_path)) dir.create(plots_path, recursive = TRUE)
  
  # format data
  df  = subset(data[[hid]], select = c('kWh', 'kWh.HVAC', 'Cooling'))
  df$TemperatureF = weather$TemperatureF
  df1 = melt(df[,-which(names(data[[hid]]) == 'date')], id.vars = c('TemperatureF'))
  df2 = melt(df, id.vars = c('date'))  
  
  # _____________________________
  # scatterplot with temperature
  
  plt = ggplot(df1, aes(TemperatureF, value, color = variable)) + geom_point(size = 2)
  plt = plt + facet_wrap(~variable, scales = 'free')
  plt = plt + theme_bw() +
    theme(panel.grid.major  = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          strip.text.x     = element_text(size=15),
          axis.title.x     = element_text(size = 18),
          axis.text.y      = element_text(size=15), 
          axis.text.x      = element_text(size=18),
          axis.title.y     = element_text(size=18),
          axis.title.x     = element_text(size=18),
          plot.title       = element_text(size=20),            
          legend.position  = 'none',
          axis.ticks = element_blank() ) + 
    theme(plot.title=element_text(family="Times", face="bold", size=20)) + 
    ggtitle( paste('Temperature Response:', hid )) + ylab('kWh')
  
  png(paste(plots_path,hid,'_temperature_profile.png', sep=''), width = 1600, height = 600, res = 200)
  print(plt)
  dev.off()
  
  # _________________________________
  # Time series view 
  
  plt = ggplot(df2, aes(date, value, color = variable)) + geom_point(size = 1.5) + geom_line()
  plt = plt + facet_wrap(~variable, scales = 'free', ncol = 1)
  plt = plt + theme_bw() +
    theme(panel.grid.major  = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          strip.text.x     = element_text(size=15),
          axis.title.x     = element_text(size = 18),
          axis.text.y      = element_text(size=15), 
          axis.text.x      = element_text(size=18),
          axis.title.y     = element_text(size=18),
          axis.title.x     = element_text(size=18),
          plot.title       = element_text(size=20),            
          legend.position  = 'none',
          axis.ticks = element_blank() ) + 
    theme(plot.title=element_text(family="Times", face="bold", size=20)) + 
    ggtitle( paste('Consumption Readings:', hid) )
  
  png(paste(plots_path,hid,'temperature_time_series.png', sep=''), width = 1800, height = 1400, res = 200)
  print(plt)
  dev.off()
  
}

# -----------------------------------------------
# Perform analysis on individual houses
# -----------------------------------------------

source('code/OccupancyStates.r')
source('code/occupancyAnalysis.r')

simulation = list()
for (i in names(redd)[c(1,2,3,5,6)]) {
  
  # create folder for each test case
  plots_path = paste(plots_path_root, i, '/', sep='')
  if (!file.exists(plots_path)) dir.create(plots_path, recursive = TRUE)
  
  # add mains to obtain 
  cur_data = redd[[i]][, c('date', 'mains')]
  names(cur_data) = c('date', 'kWh')
  cur_data$date = as.character(cur_data$date)
  cur_wthr = redd[[i]][, c('date', 'TemperatureF')]
  cur_wthr$date = as.character(cur_wthr$date)
  simulation[[i]] = occupancyAnalysis(cur_data, cur_wthr, i, 1, 
                                      Kmin = 2, Kmax = 4, 
                                      verbose = T, 
                                      NOBS_THRESH = 0, 
                                      resp.vars = c('(Intercept)', 'TemperatureF'),
                                      tran.vars = c('(Intercept)', 'TemperatureF'),
                                      addl.vars = c(),
                                      plots_path = plots_path, dump_path = plots_path)  
}

