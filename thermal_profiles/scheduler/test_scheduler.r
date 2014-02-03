# test_scheduler.r
#
# Test for Scheduler class.
# 
# Adrian Albert
# Last modified: October 2013.
# -----------------------------------------------------------------------
  
rm(list = ls())
options(error = recover)
library('lubridate')

source('classes/DataImporter.r')
source('classes/Scheduler.r')

plots_path = '~/Dropbox/OccupancyStates/plots/scheduling/'

# read in data
importer = new(Class='DataImporter', path = '~/Dropbox/OccupancyStates/fits/fits_sel/')
importer

# load in temperature forecast
weather  = read.csv('~/Dropbox/OccupancyStates/data/weather_stanford_10_29_2013.csv')
weather  = subset(weather, select = c('Time', 'TemperatureF'))
weather$Hour = hour(as.POSIXct(weather$Time))
weather  = aggregate(TemperatureF ~ Hour, data = weather, FUN = mean)
weather$Hour = NULL
weather.low = 25 + weather
weather.high = weather
weather.high$TemperatureF = weather.low$TemperatureF * 1.1

# generation data
generation = data.frame(Generation = 50 * weather$TemperatureF * (1 + 0.1*runif(23)))
# generation = data.frame(Generation = rep(1200,23))

# initialize scheduler
scheduler.l = new(Class = "Scheduler", importer, profile = weather.low, target = generation, budget = 10)
# scheduler.h = new(Class = "Scheduler", importer, profile = weather.high, target = generation)

# solve deterministic proglem
scheduler.l = solveDeterministic(scheduler.l)

source('classes/Scheduler.r')
# plot schedule heatmap
selected.rnd = importer@STATES_PAR$UID[sample(1:nrow(importer@STATES_PAR), 20)]
png(filename = paste(plots_path, 'heatmap_select_20.png', sep=''), height = 800, width = 800, res = 150)
plot(scheduler.l, selected = selected.rnd)
dev.off()

# plot selected schedules (Alice and Bob)
selected = c(Alice = 3284167, Bob = 3675267)
png(filename = paste(plots_path, 'profile_selected_2.png', sep=''), height = 600, width = 1200, res = 150)
plot(scheduler.l, selected = selected, type = 'profile')
dev.off()

# plot selected schedules (randomly selected users)
png(filename = paste(plots_path, 'profile_selected_10.png', sep=''), height = 600, width = 1200, res = 150)
plot(scheduler.l, selected = selected.rnd[1:10], type = 'profile')
dev.off()

# plot density of total daily effort
png(filename = paste(plots_path, 'profile_density.png', sep=''), height = 600, width = 1200, res = 150)
plot(scheduler.l, type = 'density')
dev.off()

# plot density of total daily effort
png(filename = paste(plots_path, 'profile_density_hourly.png', sep=''), height = 1000, width = 2000, res = 150)
plot(scheduler.l, type = 'density-hourly')
dev.off()

# plot density of total daily effort
png(filename = paste(plots_path, 'profile_inputs.png', sep=''), height = 600, width = 2000, res = 150)
plot(scheduler.l, type = 'inputs')
dev.off()

# plot aggregate profiles (Generation, Baseload, Controlled)
png(filename = paste(plots_path, 'aggregate_profile.png', sep=''), height = 800, width = 1600, res = 150)
plot(scheduler.l, type = 'aggregates')
dev.off()

