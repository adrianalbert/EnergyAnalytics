# test_scheduler.r
#
# Test for Scheduler class.
# 
# Adrian Albert
# Last modified: July 2014.
# -----------------------------------------------------------------------
  
rm(list = ls())
options(error = recover)
setwd('~/EnergyAnalytics/thermal_profiles/scheduler/')

# ____________________________________________________
# Initializations....

library('lubridate')
library('parallel')
source('../../clustering/kError.r', chdir = T)

source('classes/DataImporter.r', chdir = T)
source('classes/Scheduler.r', chdir = T)

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
Goal = data.frame(Goal = 5 * c(rep(3, 6), rep(2, 8), rep(5, 4), rep(4,6)))
# Goal = data.frame(Goal = rep(1200,23))

# load in baby names
baby_names  = read.csv('~/Dropbox/OccupancyStates/data/baby-names.csv')

# summer rate for PG&E
## summer rate from Residential Time-of-Use Schedule E-6 and Rate
prate  = 0.28719
pprate = 0.17528
oprate = 0.10074
ratevec= c(rep(oprate, 10), rep(prate, 11), rep(pprate, 3))
ratevec0 = rep(1,24)

# ____________________________________________________
# Import and format data for scheduler

# read in data
importer = new(Class='DataImporter', path = DATA_PATH)
importer

# produce inputs for scheduler
inputs = computeProfiles(importer, x_vec = weather$TemperatureF)

# give the users some more readable names
user_names = unique(as.character(baby_names$name))[1:length(inputs)]
names(user_names) = names(inputs)
names(inputs) = user_names

# sample example profiles and assign names to users
uids_sel        = user_names[c(10, 20, 30)]
selection       = inputs[uids_sel]

# ____________________________________________________
# Set up scheduling problem for a given budget

# initialize scheduler
setup_info = list(DT            = 5,
                  goal          = Goal,
                  tou_rates     = ratevec)

# initialize scheduler
scheduler = new(Class = "Scheduler", inputs, setup = setup_info)

# compute schedules in the most general case (tailored schedule for each user)
options = list(budget = 10, gamma = 0.01)
scheduler = solveSchedules(scheduler, options = options)

# plot example effort profiles
png(filename = paste(PLOTS_PATH, 'effort-profiles.png', sep=''), height = 1000, width = 2000, res = 180)
plot(scheduler, type = 'effort-profiles')
dev.off()
 
# plot matching between aggregate profile and goal
png(filename = paste(PLOTS_PATH, 'goal_match.png', sep=''), height = 600, width = 1000, res = 180)
plot(scheduler, type = 'goal-match')
dev.off()

# ____________________________________________________
# Function to compute cluster-based agreement 

compute_clust_agg = function(Abar, W, U, ASS, g, q) {
  
  U1    = U[ASS,]; W1 = W[ASS]
  Delta.bar = colSums(Abar * U1)
  Delta.Var = matrix(0, nrow = ncol(Abar), ncol = ncol(Abar))
  for (i in 1:length(W)) {
    Delta.Var = Delta.Var + diag(U1[i,]) %*% W[[i]] %*% diag(U1[i,])
  }
  EC = sum((Delta.bar - g)^2 * q)
  
  return(list(mu = Delta.bar, var = Delta.Var, EC = EC))
}

# ____________________________________________________
# Compute cluster-based solution segmentation 

# options 
setup_info = list(DT            = 5,
                  goal          = Goal,
                  tou_rates     = ratevec)

# prepare data 
cur_inputs = inputs
Abar = do.call('rbind', lapply(cur_inputs, function(l) as.numeric(l$a$mu[,1])))
rownames(Abar)= names(cur_inputs)
Abar[which(Abar < 0)] = 0
W = lapply(cur_inputs, function(l) diag(diag(l$a$covmat)))

# group together all analysis logic
perform_analysis_cluster = function(k) {
  
  # cluster profiles
  print(k)
  fit = kError(Abar, W, k, iter = 100)
  
  # scheduling inputs
  inputs.clust = lapply(1:k, function(k) {
    a = list(mu = as.matrix(fit$centers[k,]),
             covmat = fit$errors[[k]])
    return(list(a = a))
  })
  names(inputs.clust) = 1:k
  
  # initialize scheduler
  scheduler.clust = new(Class = "Scheduler", inputs.clust, setup = setup_info)
  
  # compute schedules in the most general case (tailored schedule for each user)
  Nr = table(fit$assignment)
  options = list(budget = 10, gamma = 0.01, Nr = Nr, presaved = FALSE)
  scheduler.clust = solveSchedules(scheduler.clust, options = options)
  
  # compute aggregate by applying to each profile its class-based schedule
  r = compute_clust_agg(Abar, W, scheduler.clust@OUTPUT$U, fit$assignment, setup_info$goal, setup_info$tou_rates)
  
  PLOTS_PATH_CUR = paste(PLOTS_PATH, 'varying_K/K_', k, '/', sep = '')
  dir.create(PLOTS_PATH_CUR, recursive = T)
  
  # plot schedules per classs
  pdf(file = paste(PLOTS_PATH_CUR, 'cluster-effort-profiles.pdf', sep=''), height = 8, width = 16)
  plot(scheduler.clust, type = 'effort-profiles')
  dev.off()
  
  # plot match between goal and aggregate
  pdf(file = paste(PLOTS_PATH_CUR, 'cluster-agg-goal.pdf', sep=''), height = 4, width = 6)
  plot(scheduler.clust, type = 'goal-match', compare = r)
  dev.off()
  
  # plot comparison of # users selected from each class
  pdf(file = paste(PLOTS_PATH_CUR, 'cluster-selected.pdf', sep=''), height = 4, width = 8)
  plot(scheduler.clust, type = 'number-selected')
  dev.off()
  
  return(list(nr        = scheduler.clust@OUTPUT$nr,
              Delta.bar = r$mu,
              Delta.var = r$var,
              EC        = r$EC)) 
}


# ____________________________________________________
# Produce example plots: data inputs

# plot performance stats
png(filename = paste(PLOTS_PATH, 'performance_cross_validation.png', sep=''), height = 400, width = 1000, res = 150)
plot(importer)
dev.off()

# plot states characteristics stats
png(filename = paste(PLOTS_PATH, 'state_stats.png', sep=''), height = 400, width = 1400, res = 150)
plot(importer, type = 'state-stats', inputs = user_names)
dev.off()

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
png(filename = paste(PLOTS_PATH, 'example_state_space.png', sep=''), height = 500, width = 1600, res = 180)
plot(importer, type = 'state-space', inputs = uids_sel)
dev.off()

# temperature & goal plot
df = cbind(weather, 
           Goal = Goal$Goal, 
           Rate = ratevec,
           #Reductions = 10 * c(rep(3, 6)*(0.9 + 0.2*runif(6)), rep(2, 8)*(0.9 + 0.2*runif(8)), rep(5, 4)*(0.9 + 0.2*runif(4)), rep(4,6)*(0.9 + 0.2*runif(6))),
           Time = 1:nrow(weather))
names(df) = c('Temperature T(t) [F]', 
              'Goal g(t) [kWh]', 
              'Penalty q(t) [$/kWh]',
              #paste('Reductions', expression(Delta), '(t) [kWh]'), 
              'Time')
df = melt(df, id.vars = 'Time')

# construct plot
plt = ggplot(df, aes(x = Time, y = value, color = variable)) + 
  geom_point(size = 2.5) + geom_line(size=1.5) + 
  facet_wrap(~variable, nrow =1, scales = 'free')
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
  
png(filename = paste(PLOTS_PATH, 'inputs_profiles.png', sep=''), height = 450, width = 1500, res = 170)
plt
dev.off()
