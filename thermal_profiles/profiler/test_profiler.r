# test_profiler.r
#
# Test for Scheduler class.
# 
# Adrian Albert
# Last modified: February 2014.
# -----------------------------------------------------------------------

rm(list = ls())
options(error = recover)
setwd('~/EnergyAnalytics/thermal_profiles/profiler/')

# _________________________
# Initializations....

source('classes/DataFormatter.r')
source('classes/StateDecoder.r')
source('classes/Interpreter.r')
source('classes/Visualizer.r')

library(utils)
library('segmented')

source('stateProcessorWrapper.r')
source('stateVisualizerWrapper.r')

# set plots directory
PLOTS_PATH = '~/Dropbox/OccupancyStates/plots/bakersfield/'
dir.create(file.path(PLOTS_PATH))

DUMP_PATH = '~/Dropbox/OccupancyStates/fits/bakersfield/'
dir.create(file.path(DUMP_PATH))

# load consumption and weather data
load('~/Dropbox/OccupancyStates/data/selection_consumption.RData')

# get test data
# UID       = as.numeric('3284167') # unique(consumption.ok$UID)[1]
UID       = as.numeric('3675267') # unique(consumption.ok$UID)[1]
raw_data  = subset(consumption.ok, UID == as.numeric('3675267'))
ZIP       = as.character(raw_data$ZIP5[1])
raw_data  = raw_data[,-c(1:4)]
wthr_data = weather.ok[[ZIP]]
names(wthr_data)[1] = 'date'
timesteps = wthr_data$date
wthr_data = wthr_data[,c('date', 'TemperatureF')]
wthr_data$TemperatureD = wthr_data$TemperatureF - 65

#wthr_data$date = NULL

# Rprof(filename = 'Rprof.out', interval = 0.02, memory.profiling = F)

# define model learning controls
controls = list(
  Kmin = 4, Kmax = 4, 
  maxit = 50, nRestarts = 5, tol = 1e-4,
  thresh.R2 = 0.85, thresh.MAPE = 0.10,
  test.periods = 12)

# learn model
source('stateProcessorWrapper.r')
res = stateProcessorWrapper(raw_data, wthr_data, 'Bob', 
                            controls = controls,
                            train.frac = 0.95, 
                            verbose = F, 
                            dump_path = DUMP_PATH)

# visualizations
source('stateVisualizerWrapper.r')
vis.interval = 3*24
interval = c(timesteps[1], timesteps[1+vis.interval])   
vis = stateVisualizerWrapper(res$decoder, res$interpreter,
                             interval = interval,
                             plots_path = PLOTS_PATH)

# Rprof(NULL)
# profiled = summaryRprof(filename='Rprof.out', memory = 'none')
