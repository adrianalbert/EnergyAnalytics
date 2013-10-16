# _______________________________________
# Test methods for class OccupancyStates

rm(list = ls())

options(error = recover)

source('classes/DataFormatter.r')
source('classes/StateDecoder.r')
library(utils)

load('~/Dropbox/OccupancyStates/data/selection_consumption.RData')

# get test data
UID       = as.character(20854) # unique(consumption.ok$UID)[1]
raw_data  = subset(consumption.ok, UID == UID)[,-c(1:4)]
ZIP       = as.character(93301)# as.character(raw_data$ZIP5[1])
wthr_data = weather.ok[[ZIP]]
names(wthr_data)[1] = 'date'
timesteps = wthr_data$date
wthr_data = wthr_data[,-which(names(wthr_data) == 'date')]

Rprof(filename = 'Rprof.out', interval = 0.02, memory.profiling = F)

# construct DataFormatter object
formatter = new(Class='DataFormatter', raw_data, UID)  
formatter = addCovariates(formatter, timesteps, wthr_data)    
good.data = extractFormattedData(formatter)    

# define model learning controls
controls = list(
  Kmin = 2, Kmax = 4, 
  maxit = 100, nRestarts = 5,
  thresh.R2 = 0.85, thresh.MAPE = 0.15,
  test.periods = 12)

# initialize model
source('classes/StateDecoder.r')
decoder   = new(Class='StateDecoder', 
                good.data$data, good.data$timestamps, good.data$UID,
                train.frac = 0.9, 
                tran.vars = c('(Intercept)', 'TemperatureF'), 
                resp.vars = c('(Intercept)', 'TemperatureF'),
                controls = controls)  
# HMM analysis
decoder   = learnStateDecoder(decoder, verbose = T)
show(decoder)

# perform interpretation and feature extraction
interpreter  = new(class = "Interpreter", decoder)
features     = extractClassificationFeatures(interpreter)

Rprof(NULL)
profiled = summaryRprof(filename='Rprof.out', memory = 'none')
