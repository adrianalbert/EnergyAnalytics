# test_scheduler.r
#
# Test for Scheduler class.
# 
# Adrian Albert
# Last modified: February 2014.
# -----------------------------------------------------------------------
  
rm(list = ls())
# options(warn=2)
options(error = recover)
setwd('~/EnergyAnalytics/thermal_profiles/scheduler/')

# ____________________________________________________
# Initializations....

library('lubridate')

source('classes/DataImporter.r')
source('classes/Scheduler.r')

PLOTS_PATH = '~/Dropbox/OccupancyStates/plots/scheduling/'
DATA_PATH  = '~/Dropbox/OccupancyStates/fits/bakersfield/'

# load in temperature forecast
weather  = read.csv('~/Dropbox/OccupancyStates/data/weather_stanford_10_29_2013.csv')
weather  = subset(weather, select = c('Time', 'TemperatureF'))
weather$Hour = hour(as.POSIXct(weather$Time))
weather  = aggregate(TemperatureF ~ Hour, data = weather, FUN = mean)
weather$Hour = NULL
weather$TemperatureF = weather$TemperatureF + 25

# Goal data
Goal = data.frame(Goal = 50 * weather$TemperatureF * (1 + 0.1*runif(23)))
# Goal = data.frame(Goal = rep(1200,23))

# ____________________________________________________
# Format data for scheduler

# read in data
source('classes/DataImporter.r')
importer = new(Class='DataImporter', path = DATA_PATH)
importer

# produce inputs for scheduler
inputs = computeProfiles(importer, x_vec = weather$TemperatureF)

save.image(file = '~/Dropbox/OccupancyStates/data/bakersfield_processed_profiles.RData')

# # ____________________________________________________
# # Set up scheduling problem
# 
# # initialize scheduler
# setup_info = list(cvx.setup.dir = "/usr/local/MATLAB/cvx/",
#                   budget        = 10,
#                   DT            = 5)
#                   
# scheduler = new(Class = "Scheduler", importer, profile = weather, target = Goal, setup = setup_info)
# 
# # solve deterministic proglem
# scheduler = solveDeterministic(scheduler)
# 
# ____________________________________________________
# Produce plots 

# # plot schedule heatmap
# selected.rnd = importer@STATES_PAR$UID[sample(1:nrow(importer@STATES_PAR), 20)]
# png(filename = paste(PLOTS_PATH, 'heatmap_select_20.png', sep=''), height = 800, width = 800, res = 150)
# plot(scheduler, selected = selected.rnd)
# dev.off()
# 
# # plot selected schedules (Alice and Bob)
# selected = c(Alice = 3284167, Bob = 3675267)
# png(filename = paste(PLOTS_PATH, 'profile_selected_2.png', sep=''), height = 600, width = 1200, res = 150)
# plot(scheduler, selected = selected, type = 'profile')
# dev.off()
# 
# # plot selected schedules (randomly selected users)
# png(filename = paste(PLOTS_PATH, 'profile_selected_10.png', sep=''), height = 600, width = 1200, res = 150)
# plot(scheduler, selected = selected.rnd[1:10], type = 'profile')
# dev.off()
# 
# # plot density of total daily effort
# png(filename = paste(PLOTS_PATH, 'profile_density.png', sep=''), height = 600, width = 1200, res = 150)
# plot(scheduler, type = 'density')
# dev.off()
# 
# # plot density of total daily effort
# png(filename = paste(PLOTS_PATH, 'profile_density_hourly.png', sep=''), height = 1000, width = 2000, res = 150)
# plot(scheduler, type = 'density-hourly')
# dev.off()
# 
# # plot density of total daily effort
# png(filename = paste(PLOTS_PATH, 'profile_inputs.png', sep=''), height = 600, width = 2000, res = 150)
# plot(scheduler, type = 'inputs')
# dev.off()
# 
# # plot aggregate profiles (Goal, Baseload, Controlled)
# png(filename = paste(PLOTS_PATH, 'aggregate_profile.png', sep=''), height = 800, width = 1600, res = 150)
# plot(scheduler, type = 'aggregates')
# dev.off()
# 
