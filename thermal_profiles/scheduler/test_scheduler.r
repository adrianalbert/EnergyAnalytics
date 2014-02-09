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
weather = rbind(weather, weather$TemperatureF[nrow(weather)]*0.95)

# Goal data
Goal = data.frame(Goal = 50 * weather$TemperatureF * (1 + 0.1*runif(23)))
# Goal = data.frame(Goal = rep(1200,23))

# load in baby names
baby_names  = read.csv('~/Dropbox/OccupancyStates/data/baby-names.csv')

# ____________________________________________________
# Format data for scheduler

# read in data
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
# Produce example plots 

# plot performance stats
png(filename = paste(PLOTS_PATH, 'performance_cross_validation.png', sep=''), height = 400, width = 1000, res = 150)
plot(importer)
dev.off()

# sample example profiles and assign names to users
uids_sel        = sample(names(inputs), 3)
names(uids_sel) = sample(baby_names$name, 3)
selection       = inputs[uids_sel]
names(selection)= names(uids_sel)

# plot example profiles
png(filename = paste(PLOTS_PATH, 'example_profiles.png', sep=''), height = 800, width = 1600, res = 180)
plot(importer, type = 'profiles', inputs = selection)
dev.off()

# plot covariance matrices
plts = plot(importer, type = 'covmat', inputs = selection)
png(filename = paste(PLOTS_PATH, 'example_covmat_response.png', sep=''), height = 600, width = 1600, res = 180)
plts[['response']]
dev.off()
png(filename = paste(PLOTS_PATH, 'example_covmat_baseload.png', sep=''), height = 600, width = 1600, res = 180)
plts[['baseload']]
dev.off()

# plot example propagated distributions
png(filename = paste(PLOTS_PATH, 'example_distribution.png', sep=''), height = 600, width = 2000, res = 180)
plot(importer, type = 'distributions', inputs = selection)
dev.off()

# plot example state space charaterization (a, b, sigma)
source('classes/DataImporter.r')
png(filename = paste(PLOTS_PATH, 'example_state_space.png', sep=''), height = 400, width = 1600, res = 180)
plot(importer, type = 'state-space', inputs = uids_sel)
dev.off()

# temperature & goal plot
df = cbind(weather, 
           Goal = 10 * c(rep(3, 6), rep(2, 8), rep(5, 4), rep(4,6)), 
           Reductions = 10 * c(rep(3, 6)*(0.9 + 0.2*runif(6)), rep(2, 8)*(0.9 + 0.2*runif(8)), rep(5, 4)*(0.9 + 0.2*runif(4)), rep(4,6)*(0.9 + 0.2*runif(6))),
           Time = 1:nrow(weather))
names(df) = c('Temperature T(t) [F]', 'Goal g(t) [kWh]', paste('Reductions', expression(Delta), '(t) [kWh]'), 'Time')
df = melt(df, id.vars = 'Time')

# construct plot
plt = ggplot(df, aes(x = Time, y = value, color = variable)) + 
  geom_point(size = 2.5) + geom_line(size=1.5) + 
  facet_wrap(~variable, nrow = 1, scales = 'free')
plt = plt + theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.text.x     = element_text(size=18),
        axis.text.y      = element_text(size=18), 
        axis.text.x      = element_text(size=18),
        axis.title.y     = element_text(size=18),
        axis.title.x     = element_text(size=18),
        plot.title       = element_text(size=20),            
        legend.text      = element_text(size=18),
        legend.title     = element_text(size=18),
        legend.position  = 'none',
        axis.ticks = element_blank()) + 
  ggtitle('Input Profiles') + xlab('Hour of Day')
  
png(filename = paste(PLOTS_PATH, 'inputs_profiles.png', sep=''), height = 450, width = 1800, res = 180)
plt
dev.off()

# ____________________________________________________
# Produce scheduling plots



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
