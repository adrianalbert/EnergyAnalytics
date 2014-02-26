# run_profiler_bakersfield.r
#
# Runs the Profiler for the Bakersfield data.
# 
# Adrian Albert
# Last modified: February 2014.
# -----------------------------------------------------------------------

rm(list = ls())
options(error = recover)
setwd('~/EnergyAnalytics/thermal_profiles/profiler/')

LOG_OUTPUT = TRUE

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

LOGS_PATH = '~/Dropbox/OccupancyStates/logs/bakersfield/'
dir.create(file.path(LOGS_PATH))

# load consumption and weather data
load('~/Dropbox/OccupancyStates/data/selection_consumption.RData')
consumption.ok$UID = as.character(consumption.ok$UID)
uids = sort(unique(consumption.ok$UID))
will.plot = runif(length(uids))

# Rprof(filename = 'Rprof.out', interval = 0.02, memory.profiling = F)

# define model learning controls
controls = list(
  Kmin = 4, Kmax = 4, 
  maxit = 50, nRestarts = 5, tol = 1e-4,
  thresh.R2 = 0.85, thresh.MAPE = 0.10,
  test.periods = 12)

# wrapper for profiling & visualization functionality
run_wrapper = function(i) {
  
  set.seed(0)
  uid = uids[i]
  cat(paste(i, '...', sep=''))
  if (i %% 20 == 0) cat('\n')
  
  if (LOG_OUTPUT) sink(paste(LOGS_PATH, "log_", i,".txt", sep=''), append=F, type = c('output', 'message'))
  cat(paste('Analyzing user', i, '/', length(uids), '\n'))
  
  # get test data
  raw_data  = subset(consumption.ok, UID == uid)
  ZIP       = as.character(raw_data$ZIP5[1])
  raw_data  = raw_data[,-c(1:4)]
  wthr_data = weather.ok[[ZIP]]
  names(wthr_data)[1] = 'date'
  timesteps = wthr_data$date
  wthr_data = wthr_data[,c('date', 'TemperatureF')]
  wthr_data$TemperatureD = wthr_data$TemperatureF - 65

  # learn model
  res = stateProcessorWrapper(raw_data, wthr_data, uid, 
                              controls = controls,
                              train.frac = 0.95, 
                              verbose = F, 
                              dump_path = DUMP_PATH)
  
  # visualizations
  if (will.plot[i] < 0.1) {
    vis.interval = 3*24
    interval = c(timesteps[1], timesteps[1+vis.interval])   
    vis = stateVisualizerWrapper(res$decoder, res$interpreter,
                                 interval = interval,
                                 plots_path = PLOTS_PATH)    
  }  
  
  # reset output connection
  if (LOG_OUTPUT) sink(NULL)    
}

# run profiler in parallel
ptm <- proc.time()
result = mclapply(1:length(uids), 
                  FUN = run_wrapper, 
                  mc.silent = F, mc.preschedule = TRUE, mc.cores = 5)
dt = proc.time() - ptm
print(dt)
