# _________________________
# Test methods for profiler

rm(list = ls())

options(error = recover)

source('classes/DataFormatter.r')
source('classes/StateDecoder.r')
source('classes/Interpreter.r')
source('classes/Visualizer.r')

library(utils)
library('segmented')

setwd('~/Dropbox/EnergyAnalytics/thermal_profiles/profiler/')
source('stateProcessorWrapper.r')
source('stateVisualizerWrapper.r')

PLOTS_PATH = './'

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
#wthr_data$date = NULL

# Rprof(filename = 'Rprof.out', interval = 0.02, memory.profiling = F)

# estimate a "breakpoint" indoors temperature
# fit  = lm('use ~ TemperatureF', data = cbind(obs = rawwthr_data)    
# fmla = as.formula(paste('~', 'TemperatureF'))
# fit.seg <- try(segmented(fit, seg.Z = fmla, psi = 75))
# psi  = fit.seg$psi[1,2]        
wthr_data$TemperatureD = wthr_data$TemperatureF - 75

# define model learning controls
controls = list(
  Kmin = 4, Kmax = 4, 
  maxit = 50, nRestarts = 5, tol = 1e-4,
  thresh.R2 = 0.85, thresh.MAPE = 0.10,
  test.periods = 12)

dir.create(file.path(PLOTS_PATH, UID))
# learn model
res = stateProcessorWrapper(raw_data, wthr_data, 'Bob', 
                            controls = controls,
                            train.frac = 0.95, 
                            verbose = F)

# visualizations
source('stateVisualizerWrapper.r')
vis.interval = 3*24
interval = c(timesteps[1], timesteps[1+vis.interval])   
vis = stateVisualizerWrapper(res$decoder, res$interpreter,
                             interval = interval,
                             plots_path = paste(PLOTS_PATH, '/', UID, '/', sep=''))

# Rprof(NULL)
# profiled = summaryRprof(filename='Rprof.out', memory = 'none')
