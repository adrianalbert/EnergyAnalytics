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
library('parallel')

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
Goal = data.frame(Goal = 5 * c(rep(3, 6), rep(2, 8), rep(5, 4), rep(4,6)))
# Goal = data.frame(Goal = rep(1200,23))

# load in baby names
baby_names  = read.csv('~/Dropbox/OccupancyStates/data/baby-names.csv')

# load in covariance matrix data ("covmat")
# load('~/Dropbox/OccupancyStates/data/bakersfield_covmat.RData')
# covmat.diag = mclapply(1:24, mc.cores = 5,
#                         function(t) {
#                           res = sapply(1:length(covmat), function(i) covmat[[i]][1,t])
#                           res = unlist(res)
#                           res = matrix(res, ncol = sqrt(length(res)), byrow = T)                        
#                           return(res)
#                         })

# summer rate for PG&E
## summer rate from Residential Time-of-Use Schedule E-6 and Rate
prate  = 0.28719
pprate = 0.17528
oprate = 0.10074
ratevec= c(rep(oprate, 10), rep(prate, 11), rep(pprate, 3))

# ____________________________________________________
# Format data for scheduler

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

# # get cross-user covariances for selected users
# idx_sel = which(as.character(idx$i) %in% names(uids_sel) & as.character(idx$j) %in% names(uids_sel) )
# sel_cov = covmat[idx_sel]
# sel_cov = do.call('rbind', sel_cov)
# sel_cov = as.data.frame(sel_cov)
# usr_sel = idx[idx_sel,]
# usr_sel$i = uids_sel[as.character(usr_sel$i)]
# usr_sel$j = uids_sel[as.character(usr_sel$j)]
# sel_cov$UID = paste(usr_sel$i, usr_sel$j, sep='-')
# sel_cov = sel_cov[c(2,3,6),]

# save.image(file = '~/Dropbox/OccupancyStates/data/bakersfield_processed_profiles.RData')
# load('~/Dropbox/OccupancyStates/data/bakersfield_processed_profiles.RData')

# ____________________________________________________
# Set up scheduling problem

# initialize scheduler
setup_info = list(cvx.setup.dir = "/usr/local/MATLAB/cvx/",
                  DT            = 5,
                  noClusters    = c(9, 9),
                  goal          = Goal,
                  max.iter      = 20,
                  tou_rates     = ratevec)

# initialize scheduler
scheduler = new(Class = "Scheduler", inputs, setup = setup_info)

# compute segmentation
scheduler = segmentProfiles(scheduler)

# compute schedules
source('classes/Scheduler.r')
options = list(budget = 10, gamma = 0.01)
scheduler = solveSchedulesClusters(scheduler, options = options)

png(filename = paste(PLOTS_PATH, 'effort-profiles.png', sep=''), height = 1000, width = 2000, res = 180)
plot(scheduler, type = 'effort-profiles')
dev.off()

png(filename = paste(PLOTS_PATH, 'goal_match.png', sep=''), height = 600, width = 1000, res = 180)
plot(scheduler, type = 'goal-match')
dev.off()


# ____________________________________________________
# Simulate solution for different budget/gamma levels

# simulation grid
bg = expand.grid(budget = c(1:8), gamma = c(0))

dir.create(file.path(paste(PLOTS_PATH, 'simulation/', sep='')))

source('classes/Scheduler.r')
res = lapply(1:nrow(bg),
#               mc.cores = 5,
#               mc.silent = F,
#               mc.preschedule = TRUE, 
              function(i) {
                
                budget    = bg[i, 'budget']
                gamma     = bg[i, 'gamma']                
                print(bg[i,])
                
                # solve scheduling problem
                opts      = list(budget = budget, gamma = gamma)
                cur_sched = solveSchedulesClusters(scheduler, options = opts)
                
                # plot current solution
                png(filename = paste(PLOTS_PATH, 'simulation/effort-profiles_b_', budget, '_g_', gamma, '.png', sep=''), height = 1000, width = 2000, res = 180)
                print(plot(cur_sched, type = 'effort-profiles'))
                dev.off()                
                png(filename = paste(PLOTS_PATH, 'simulation/goal_match_b-', budget, '_g-', gamma, '.png', sep=''), height = 600, width = 1000, res = 180)
                print(plot(cur_sched, type = 'goal-match'))
                dev.off()
                
                # store results
                nr = cur_sched@OUTPUT$segments$nr
                obj= cur_sched@OUTPUT$segments$EC
                
                return(list(budget = budget, gamma = gamma, nr = nr, obj = obj))
              })


# plots
obj = sapply(res, function(r) r$obj)
nrs = sapply(res, function(r) r$nr)
bdg = sapply(res, function(r) r$budget)

df  = as.data.frame(t(nrs))
dh = data.frame(Cost = obj, No.Selected = rowSums(t(nrs)), Budget = bdg)
names(df) = paste('Segment', 1:ncol(df))
df$Budget = bdg
df = melt(df, id.vars = 'Budget')

# no users selected per class
plt = ggplot(df, aes(x = Budget, y = value)) + 
  geom_point(size = 2.5) + geom_line(size=1.5) + facet_wrap(~variable, ncol = 3)
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
  ggtitle('No. Users Enrolled') + xlab('Effort Budget beta [deg F x 5 / 24hrs]') + ylab('No. Users')

png(filename = paste(PLOTS_PATH, 'no_users_vs_budget_class.png', sep=''), height = 700, width = 1000, res = 170)
plt
dev.off()

# cost vs budget
library('pracma')
png(filename = paste(PLOTS_PATH, 'no_users_vs_budget.png', sep=''), height = 400, width = 1000, res = 170)
plotyy(dh$Budget, dh$Cost, dh$Budget, dh$No.Selected, gridp = TRUE, box.col = "grey",
       type = "b", lwd = 3, lty = 1, cex.lab = 2, cex.axis = 2, cex.main = 2, 
       xlab = "Effort Budget beta [deg F x 5 / 24hrs]", ylab = "Expected Cost", 
       main = 'DR Program Administration',
       col.y1 = "navy", col.y2 = "maroon")
legend(4.5,950, # places a legend at the appropriate place 
       c("Expected Deviation Penalty","No. Users Enrolled"), # puts text in the legend               
       lty=c(1,1), # gives the legend appropriate symbols (lines)       
       lwd=c(2.5,2.5),col=c("navy","maroon"), cex = 1)
dev.off()



# ____________________________________________________
# Produce example plots 

if (1 == 0) {

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

# # plot cross-user covariance diagonals for selected users
# df = melt(sel_cov, id.vars = 'UID')
# df$variable = as.numeric(df$variable)
# plt = ggplot(df, aes(x = variable, y = value, color = UID)) + 
#   geom_point(size = 2.5) + geom_line(size=1.5) 
#   #facet_wrap (~variable, nrow =1, scales = 'free')
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
#         legend.position  = c(0.3, 0.5), 
#         axis.ticks = element_blank()) + 
#   ggtitle('Input Profiles') + xlab('Hour of Day') + ylab('Covariance Diagonal')
# 
# png(filename = paste(PLOTS_PATH, 'cross_user_covar_example.png', sep=''), height = 500, width = 1000, res = 170)
# plt
# dev.off()

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

# ____________________________________________________
# Produce clustering plots

png(filename = paste(PLOTS_PATH, 'cluster_centers_2.png', sep=''), height = 1000, width = 2000, res = 180)
plot(scheduler, type = 'cluster-centers')
dev.off()

source('classes/Scheduler.r')
png(filename = paste(PLOTS_PATH, 'effort-profiles.png', sep=''), height = 1000, width = 2000, res = 180)
plot(scheduler, type = 'effort-profiles')
dev.off()

}

