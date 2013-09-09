# _______________________________________
# Test methods for class OccupancyStates

rm(list = ls())

options(error = recover)

setwd('~/Dropbox/OccupancyStates/')
source('code/OccupancyStates.r')

load('data/selection_consumption.RData')

# get test data
raw_data  = subset(consumption.ok, UID == unique(consumption.ok$UID)[1099])
wthr_data = weather.ok[[as.character(raw_data$ZIP5[1])]]
names(wthr_data)[1] = 'date'

# populate OccupancyStates object
test      = new(Class='OccupancyStates', raw_data, unique(raw_data$UID), unique(raw_data$ZIP5))
test      = addWeather(test, wthr_data$date, wthr_data[,-1])
test      = prepareData(test, train.frac = 0.8, 
                        resp.vars = c('(Intercept)', 'TemperatureF'),
                        tran.vars = c('(Intercept)', 'TemperatureF'),
                        addl.vars = c('HourOfDay'))

# _______________________
# OLS analysis

test      = fitOLS(test)
test      = computeBreakpoint(test)

# _______________________
# HMM analysis

test      = learnOccupancyStates(test, Kmin = 4, Kmax = 6)
test      = computePredictionAccuracy(test, test.periods = 12)
test      = computeStatsHMM(test)
test      = computeContributionsHMM(test, verbose = T)

thresholds = c(none = 0.1, low = 0.3, high = 0.8)
test      = computeMetricsHMM(test, thresholds = thresholds)

source('code/OccupancyStates.r')
# dump coefficients to file
res1      = dumpComputationToFile(test, path = NULL)
res2      = dumpComputation(test, path = NULL)

# _______________________
# Test plotting methods

PLOT = T

source('code/utils/plot_utils.r')

if (PLOT) {
  png(filename = 'kwh.png', width = 1200, height = 800)
  plot(test, interval = c('2010-12-01', '2010-12-10'))
  dev.off()
  
  png(filename = 'weather.png', width = 1200, height = 800)
  plot(test, type = 'weather', interval = c('2010-12-01', '2010-12-10'))
  dev.off()
  
  png(filename = 'OLS.png', width = 1200, height = 800)
  plot(test, type = 'OLS-fit', interval = c('2010-12-01', '2010-12-10'))
  dev.off()
  
  png(filename = 'HMM_ts.png', width = 1400, height = 600, res = 100)
  print(plot(test, type = 'HMM-ts', interval = c('2010-12-01', '2010-12-10')))
  dev.off()
  
  png(filename = 'HMM_MC_cov.png', width = 1200, height = 800, res = 100)
  print(plot(test, type = 'HMM-MC-cov', interval = c('2010-12-01', '2010-12-10')))
  dev.off()
  
  png(filename = 'HMM_res.png', width = 1400, height = 400, res = 100)
  plot(test, type = 'HMM-res')
  dev.off()
  
  png(filename = 'HMM-coefs-ts.png', width = 1400, height = 800, res = 100)
  plot(test, type = 'HMM-coefs-ts', interval = c('2010-12-01', '2010-12-10'),
       covar = c('TemperatureF'))
  dev.off()
  
  png(filename = 'HMM-consumption-heatmap.png', width = 1400, height = 800, res = 100)
  plot(test, type = 'HMM-state-heatmap', covar = 'kWh')
  dev.off()
  
  png(filename = 'HMM-temperature-heatmap.png', width = 1400, height = 800, res = 100)
  plot(test, type = 'HMM-state-heatmap', covar = 'TemperatureF')
  dev.off()
  
  png(filename = 'HMM-state-heatmap.png', width = 1400, height = 800, res = 100)
  plot(test, type = 'HMM-state-heatmap', covar = 'States')
  dev.off()
  
  source('code/OccupancyStates.r')  
  source('code/utils/plot_utils.r')
  png(filename = 'HMM-state-breakdown.png', width = 1400, height = 800, res = 100)
  print(plot(test, type = 'HMM-state-breakdown'))
  dev.off()
  
  png(filename = 'HMM-state-heatmap2.png', width = 2000, height = 1000, res = 200)
  print(plot(test, type = 'HMM-state-heatmap2', covar = 'States'))
  dev.off()
  
  png(filename = 'HMM-contrib-ts.png', width = 1400, height = 800, res = 100)
  plot(test, type = 'HMM-contrib-ts', interval = c('2010-12-01', '2010-12-10'),
       covar = c('TemperatureF'))
  dev.off()
  
  png(filename = 'HMM-contrib-tot.png', width = 800, height = 600, res = 100)
  print(plot(test, type = 'HMM-contrib-tot'))
  dev.off()
  
  png(filename = 'HMM-dep-covar.png', width = 1000, height = 600, res = 100)
  print(plot(test, type = 'HMM-dep-covar', covar = 'TemperatureF'))
  dev.off()
  
  source('code/utils/plot_utils.r')
  png(filename = 'HMM-trans-prob.png', width = 1600, height = 800, res = 100)
  print(plot(test, type = 'HMM-trans-prob', covar = 'TemperatureF'))
  dev.off()
  
  png(filename = 'HMM-dep-covar-sep.png', width = 1000, height = 600, res = 100)
  print(plot(test, type = 'HMM-dep-covar', covar = 'TemperatureF', separate = T))
  dev.off()
  
  png(filename = 'HMM-dep-covar-orig.png', width = 1000, height = 600, res = 100)
  print(plot(test, type = 'HMM-dep-covar', covar = 'TemperatureF', highlight = F))
  dev.off()
  
  png(filename = 'HMM-pacf.png', width = 1000, height = 600, res = 100)
  print(plot(test, type='HMM-pacf', PACF = T))
  dev.off()
}