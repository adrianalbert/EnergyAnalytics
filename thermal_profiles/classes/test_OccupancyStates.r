# _______________________________________
# Test methods for class OccupancyStates

rm(list = ls())

options(error = recover)

setwd('~/Dropbox/OccupancyStates/')
source('code/OccupancyStates.r')
source('code/occupancyAnalysis.r')
library(utils)

load('data/selection_consumption.RData')

# get test data
UID       = as.character(20854) # unique(consumption.ok$UID)[1]
raw_data  = subset(consumption.ok, UID == UID)
ZIP       = as.character(93301)# as.character(raw_data$ZIP5[1])
wthr_data = weather.ok[[ZIP]]
names(wthr_data)[1] = 'date'

Rprof(filename = 'Rprof.out', interval = 0.02, memory.profiling = F)
res = occupancyAnalysis(raw_data, wthr_data, UID, ZIP, NOBS_THRESH = 90, 
                        verbose = T, plots_path = './', dump_path = './',
                        resp.vars = c('(Intercept)', 'TemperatureF'),
                        tran.vars = c('(Intercept)', 'TemperatureF'),
                        addl.vars = c(),
                        Kmin = 3, Kmax = 5)
Rprof(NULL)
profiled = summaryRprof(filename='Rprof.out', memory = 'none')
