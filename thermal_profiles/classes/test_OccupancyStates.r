# _______________________________________
# Test methods for class OccupancyStates

rm(list = ls())

options(error = recover)

source('classes/DataFormatter.r')
source('classes/StateDecoder.r')
source('classes/Interpreter.r')
source('classes/Visualizer.r')

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
  Kmin = 4, Kmax = 4, 
  maxit = 100, nRestarts = 5,
  thresh.R2 = 0.85, thresh.MAPE = 0.10,
  test.periods = 12)

# initialize model
source('classes/StateDecoder.r')
decoder   = new(Class='StateDecoder', 
                good.data$data, good.data$timestamps, good.data$UID,
                train.frac = 0.95, 
                tran.vars = c('(Intercept)', 'TemperatureF'), 
                resp.vars = c('(Intercept)', 'TemperatureF'),
                controls = controls)  
# HMM analysis
decoder   = learnStateDecoder(decoder, verbose = T)
data.d    = dumpDecodedData(decoder)
show(decoder)

# perform interpretation and feature extraction
interpreter  = new(Class = "Interpreter", decoder)
data.i       = dumpInterpretedData(interpreter)
features     = extractUserFeatures(interpreter)

# visualize decoded & interpreted data
interval = c(timesteps[1], timesteps[1+24*7])   
visualizer   = new(Class = "Visualizer", decoder, interpreter,
                   interval = interval, covar = 'TemperatureF')

plots_path = './'

# just the data
png(paste(plots_path, visualizer@UID, '_raw.png', sep=''), width=1400, height=600, res = 100)
plot(visualizer)
dev.off()

# HMM fit
p1 = plot(visualizer, type='HMM-ts')
png(paste(plots_path, visualizer@UID, '_HMM_fit.png', sep=''), width=1400, height=600, res = 100)
print(p1)
dev.off()

# heatmap plots
source('classes/Visualizer.r')
png(paste(plots_path, visualizer@UID, '_HMM-consumption-heatmap.png', sep=''), width = 1400, height = 800, res = 100)
plot(visualizer, type = 'HMM-heatmap-obs')
dev.off()


# heatmap plots
png(paste(plots_path, user@UID, '_HMM-consumption-heatmap.png', sep=''), width = 1400, height = 800, res = 100)
plot(user, type = 'HMM-state-heatmap', covar = 'kWh')
dev.off()

png(paste(plots_path, user@UID, '_HMM-temperature-heatmap.png', sep=''), width = 1400, height = 800, res = 100)
plot(user, type = 'HMM-state-heatmap', covar = 'TemperatureF')
dev.off()

png(paste(plots_path, user@UID, '_HMM-state-heatmap.png', sep=''), width = 1400, height = 800, res = 100)
plot(user, type = 'HMM-state-heatmap', covar = 'States')
dev.off()

png(paste(plots_path, user@UID,'_HMM-state-breakdown.png', sep=''), width = 1400, height = 800, res = 100)
print(plot(user, type = 'HMM-state-breakdown'))
dev.off()

png(paste(plots_path, user@UID, '_HMM-state-heatmap2.png', sep=''), width = 2000, height = 1000, res = 200)
print(plot(user, type = 'HMM-state-heatmap2', covar = 'States'))
dev.off()      

# coefficient time series
png(paste(plots_path, user@UID, '_HMM-coefs-ts.png', sep=''), width = 1400, height = 800, res = 100)
print(plot(user, type = 'HMM-coefs-ts', interval = interval, covar = c('TemperatureF')))
dev.off()      

# HMM structure
png(paste(plots_path, user@UID, '_HMM_MC_cov.png', sep=''), width=1200, height=800, res = 100)
print(plot(user, type='HMM-MC-cov'))
dev.off()

# HMM PACF
png(paste(plots_path, user@UID, '_HMM_pacf.png', sep=''), width=1000, height=600, res = 100)
print(plot(user, type='HMM-pacf', PACF=T))
dev.off()

# HMM residuals
png(paste(plots_path, user@UID, '_HMM_res.png', sep=''), width=1200, height=500, res = 100)
plot(user, type = 'HMM-res')
dev.off()

# HMM covariate contributions 
png(paste(plots_path, user@UID, '_HMM_contrib_ts.png', sep=''), width = 1400, height = 800, res = 100)
plot(user, type = 'HMM-contrib-ts', interval=interval, covar = c('TemperatureF'))
dev.off()

# HMM covariate contributions (total)
png(paste(plots_path, user@UID, '_HMM_contrib_tot.png', sep=''), width = 800, height = 400, res = 100)
print(plot(user, type = 'HMM-contrib-tot', nrow=2))
dev.off()      

# temperature dependence
png(paste(plots_path, user@UID, '_HMM-dep-covar.png', sep = ''), width = 1000, height = 600, res = 150)
print(plot(user, type = 'HMM-dep-covar', covar = 'TemperatureF'))
dev.off()      

png(paste(plots_path, user@UID, 'HMM-dep-covar-sep.png', sep=''), width = 1000, height = 600, res = 150)
print(plot(user, type = 'HMM-dep-covar', covar = 'TemperatureF', separate = T))
dev.off()

png(paste(plots_path, user@UID, '_HMM-dep-covar-orig.png', sep = ''), width = 1000, height = 600, res = 150)
print(plot(user, type = 'HMM-dep-covar', covar = 'TemperatureF', highlight = F))
dev.off()

png(paste(plots_path, user@UID,'HMM-trans-prob.png',sep=''), width = 1600, height = 800, res = 100)
print(plot(user, type = 'HMM-trans-prob', covar = 'TemperatureF'))
dev.off()

png(paste(plots_path, user@UID,'HMM-stationary-prob.png',sep=''), width = 1000, height = 400, res = 100)
print(plot(user, type = 'HMM-stationary-prob', covar = 'TemperatureF'))
dev.off()

png(paste(plots_path, user@UID,'HMM-err-horiz.png',sep=''), width = 1000, height = 400, res = 100)
print(plot(user, type = 'HMM-err-horiz'))
dev.off()

png(paste(plots_path, user@UID,'HMM-time-temp.png',sep=''), width = 1000, height = 400, res = 100)
print(plot(user, type = 'HMM-time-temp', covar = 'TemperatureF'))
dev.off()

png(paste(plots_path, user@UID,'HMM-aggregate-season.png',sep=''), width = 1600, height = 600, res = 100)
print(plot(user, type = 'HMM-aggregate-season'))
dev.off()

Rprof(NULL)
profiled = summaryRprof(filename='Rprof.out', memory = 'none')
