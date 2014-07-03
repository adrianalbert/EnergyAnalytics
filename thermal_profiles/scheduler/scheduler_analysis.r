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

# source('classes/DataImporter.r', chdir = T)
# source('classes/Scheduler.r', chdir = T)
# 
# PLOTS_PATH = '~/Dropbox/OccupancyStates/plots/scheduling/'
# DATA_PATH  = '~/Dropbox/OccupancyStates/fits/bakersfield/'
# 
# # load in temperature forecast
# weather  = read.csv('~/Dropbox/OccupancyStates/data/weather_stanford_10_29_2013.csv')
# weather  = subset(weather, select = c('Time', 'TemperatureF'))
# weather$Hour = hour(as.POSIXct(weather$Time))
# weather  = aggregate(TemperatureF ~ Hour, data = weather, FUN = mean)
# weather$Hour = NULL
# weather$TemperatureF = weather$TemperatureF + 25
# weather = rbind(weather, weather$TemperatureF[nrow(weather)]*0.95)
# 
# # Goal data
# Goal = data.frame(Goal = 5 * c(rep(3, 6), rep(2, 8), rep(5, 4), rep(4,6)))
# # Goal = data.frame(Goal = rep(1200,23))
# 
# # load in baby names
# baby_names  = read.csv('~/Dropbox/OccupancyStates/data/baby-names.csv')
# 
# # summer rate for PG&E
# ## summer rate from Residential Time-of-Use Schedule E-6 and Rate
# prate  = 0.28719
# pprate = 0.17528
# oprate = 0.10074
# ratevec= c(rep(oprate, 10), rep(prate, 11), rep(pprate, 3))
ratevec0 = rep(1,24)
# 
# # ____________________________________________________
# # Import and format data for scheduler
# 
# # read in data
# importer = new(Class='DataImporter', path = DATA_PATH)
# importer
# 
# # produce inputs for scheduler
# inputs = computeProfiles(importer, x_vec = weather$TemperatureF)
# 
# # give the users some more readable names
# user_names = unique(as.character(baby_names$name))[1:length(inputs)]
# names(user_names) = names(inputs)
# names(inputs) = user_names
# 
# # sample example profiles and assign names to users
# uids_sel        = user_names[c(10, 20, 30)]
# selection       = inputs[uids_sel]
# 
# save.image('~/energy-data/bakersfield_profiles.RData')
load('~/energy-data/bakersfield_profiles.RData')
source('classes/DataImporter.r', chdir = T)
source('classes/Scheduler.r', chdir = T)

# # ____________________________________________________
# # Set up scheduling problem for a given budget
# 
# # initialize scheduler
# setup_info = list(DT            = 5,
#                   goal          = Goal,
#                   tou_rates     = ratevec)
# 
# # initialize scheduler
# scheduler = new(Class = "Scheduler", inputs, setup = setup_info)
# 
# # compute schedules in the most general case (tailored schedule for each user)
# options = list(budget = 10, gamma = 0.01)
# scheduler = solveSchedules(scheduler, options = options)

# # plot example effort profiles
# png(filename = paste(PLOTS_PATH, 'effort-profiles.png', sep=''), height = 1000, width = 2000, res = 180)
# plot(scheduler, type = 'effort-profiles')
# dev.off()
#  
# # plot matching between aggregate profile and goal
# png(filename = paste(PLOTS_PATH, 'goal_match.png', sep=''), height = 600, width = 1000, res = 180)
# plot(scheduler, type = 'goal-match')
# dev.off()

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

res = mclapply(2*(5:15), 
               mc.cores = 5, 
               perform_analysis_cluster)

# ____________________________________________________
# Fit vs no of clusters


# simulation grid



# bg = expand.grid(budget = c(1:8), gamma = c(0))

# dir.create(file.path(paste(PLOTS_PATH, 'simulation/', sep='')))
# 
# source('classes/Scheduler.r')
# res = lapply(1:nrow(bg),
# #               mc.cores = 5,
# #               mc.silent = F,
# #               mc.preschedule = TRUE, 
#               function(i) {
#                 
#                 budget    = bg[i, 'budget']
#                 gamma     = bg[i, 'gamma']                
#                 print(bg[i,])
#                 
#                 # solve scheduling problem
#                 opts      = list(budget = budget, gamma = gamma)
#                 cur_sched = solveSchedulesClusters(scheduler, options = opts)
#                 
#                 # plot current solution
#                 png(filename = paste(PLOTS_PATH, 'simulation/effort-profiles_b_', budget, '_g_', gamma, '.png', sep=''), height = 1000, width = 2000, res = 180)
#                 print(plot(cur_sched, type = 'effort-profiles'))
#                 dev.off()                
#                 png(filename = paste(PLOTS_PATH, 'simulation/goal_match_b-', budget, '_g-', gamma, '.png', sep=''), height = 600, width = 1000, res = 180)
#                 print(plot(cur_sched, type = 'goal-match'))
#                 dev.off()
#                 
#                 # store results
#                 nr = cur_sched@OUTPUT$segments$nr
#                 obj= cur_sched@OUTPUT$segments$EC
#                 
#                 return(list(budget = budget, gamma = gamma, nr = nr, obj = obj))
#               })
# 
# 
# # plots
# obj = sapply(res, function(r) r$obj)
# nrs = sapply(res, function(r) r$nr)
# bdg = sapply(res, function(r) r$budget)
# 
# df  = as.data.frame(t(nrs))
# dh = data.frame(Cost = obj, No.Selected = rowSums(t(nrs)), Budget = bdg)
# names(df) = paste('Segment', 1:ncol(df))
# df$Budget = bdg
# df = melt(df, id.vars = 'Budget')
# 
# # no users selected per class
# plt = ggplot(df, aes(x = Budget, y = value)) + 
#   geom_point(size = 2.5) + geom_line(size=1.5) + facet_wrap(~variable, ncol = 3)
# plt = plt + theme_bw() + 
#   theme(panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.background = element_blank(),
#         strip.text.x     = element_text(size=18),
#         axis.text.y      = element_text(size=18), 
#         axis.text.x      = element_text(size=18),
#         axis.title.y     = element_text(size=18),
#         axis.title.x     = element_text(size=18),
#         plot.title       = element_text(size=20),            
#         legend.text      = element_text(size=18),
#         legend.title     = element_text(size=18),
#         legend.position  = 'none',
#         axis.ticks = element_blank()) + 
#   ggtitle('No. Users Enrolled') + xlab('Effort Budget beta [deg F x 5 / 24hrs]') + ylab('No. Users')
# 
# png(filename = paste(PLOTS_PATH, 'no_users_vs_budget_class.png', sep=''), height = 700, width = 1000, res = 170)
# plt
# dev.off()
# 
# # cost vs budget
# library('pracma')
# png(filename = paste(PLOTS_PATH, 'no_users_vs_budget.png', sep=''), height = 600, width = 1000, res = 150)
# plotyy(dh$Budget, dh$Cost, dh$Budget, dh$No.Selected, gridp = TRUE, box.col = "grey",
#        type = "b", lwd = 3, lty = 1, cex.lab = 1.2, cex.axis = 1.2, cex.main = 1.3, 
#        xlab = "Effort Budget beta [deg F x 5 / 24hrs]", ylab = "Expected Cost", 
#        main = 'DR Program Administration',
#        col.y1 = "navy", col.y2 = "maroon")
# legend(4,990, # places a legend at the appropriate place 
#        c("Expected Deviation Penalty","No. Users Enrolled"), # puts text in the legend               
#        lty=c(1,1), # gives the legend appropriate symbols (lines)       
#        lwd=c(2.5,2.5),col=c("navy","maroon"), cex = 1)
# dev.off()

# # ____________________________________________________
# # Produce example plots: data inputs
# 
# # plot performance stats
# png(filename = paste(PLOTS_PATH, 'performance_cross_validation.png', sep=''), height = 400, width = 1000, res = 150)
# plot(importer)
# dev.off()
# 
# # plot states characteristics stats
# png(filename = paste(PLOTS_PATH, 'state_stats.png', sep=''), height = 400, width = 1400, res = 150)
# plot(importer, type = 'state-stats', inputs = user_names)
# dev.off()
# 
# # plot example profiles
# png(filename = paste(PLOTS_PATH, 'example_profiles.png', sep=''), height = 800, width = 1600, res = 180)
# plot(importer, type = 'profiles', inputs = selection)
# dev.off()
# 
# # plot covariance matrices
# plts = plot(importer, type = 'covmat', inputs = selection)
# png(filename = paste(PLOTS_PATH, 'example_covmat_response.png', sep=''), height = 600, width = 1600, res = 180)
# plts[['response']]
# dev.off()
# png(filename = paste(PLOTS_PATH, 'example_covmat_baseload.png', sep=''), height = 600, width = 1600, res = 180)
# plts[['baseload']]
# dev.off()
# 
# # plot example propagated distributions
# png(filename = paste(PLOTS_PATH, 'example_distribution.png', sep=''), height = 600, width = 2000, res = 180)
# plot(importer, type = 'distributions', inputs = selection)
# dev.off()
# 
# # plot example state space charaterization (a, b, sigma)
# png(filename = paste(PLOTS_PATH, 'example_state_space.png', sep=''), height = 500, width = 1600, res = 180)
# plot(importer, type = 'state-space', inputs = uids_sel)
# dev.off()
# 
# # temperature & goal plot
# df = cbind(weather, 
#            Goal = Goal$Goal, 
#            Rate = ratevec,
#            #Reductions = 10 * c(rep(3, 6)*(0.9 + 0.2*runif(6)), rep(2, 8)*(0.9 + 0.2*runif(8)), rep(5, 4)*(0.9 + 0.2*runif(4)), rep(4,6)*(0.9 + 0.2*runif(6))),
#            Time = 1:nrow(weather))
# names(df) = c('Temperature T(t) [F]', 
#               'Goal g(t) [kWh]', 
#               'Penalty q(t) [$/kWh]',
#               #paste('Reductions', expression(Delta), '(t) [kWh]'), 
#               'Time')
# df = melt(df, id.vars = 'Time')
# 
# # construct plot
# plt = ggplot(df, aes(x = Time, y = value, color = variable)) + 
#   geom_point(size = 2.5) + geom_line(size=1.5) + 
#   facet_wrap(~variable, nrow =1, scales = 'free')
# plt = plt + theme_bw() + 
#   theme(panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.background = element_blank(),
#         strip.text.x     = element_text(size=18),
#         axis.text.y      = element_text(size=18), 
#         axis.text.x      = element_text(size=18),
#         axis.title.y     = element_text(size=18),
#         axis.title.x     = element_text(size=18),
#         plot.title       = element_text(size=20),            
#         legend.text      = element_text(size=18),
#         legend.title     = element_text(size=18),
#         legend.position  = 'none',
#         axis.ticks = element_blank()) + 
#   ggtitle('Input Profiles') + xlab('Hour of Day')
#   
# png(filename = paste(PLOTS_PATH, 'inputs_profiles.png', sep=''), height = 450, width = 1500, res = 170)
# plt
# dev.off()
