# test_scheduler_approx.r
#
# Test for Scheduler class.
# 
# Adrian Albert
# Last modified: July 2014.
# -----------------------------------------------------------------------
  
rm(list = ls())
options(error = recover)

# ____________________________________________________
# Initializations....

library('lubridate')
library('parallel')

setwd('~/EnergyAnalytics/thermal_profiles/scheduler/')

# load input profiles
load('~/energy-data/bakersfield/bakersfield_profiles.RData')

source('classes/Scheduler.r', chdir = T)

PLOTS_PATH = '~/Dropbox/OccupancyStates/plots/scheduling/'

# Goal data
Goal = data.frame(Goal = 5 * c(rep(3, 6), rep(2, 8), rep(5, 4), rep(4,6)))

# summer rate for PG&E
## summer rate from Residential Time-of-Use Schedule E-6 and Rate
prate  = 0.28719
pprate = 0.17528
oprate = 0.10074
ratevec= c(rep(oprate, 10), rep(prate, 11), rep(pprate, 3))
ratevec0 = rep(1,24)

# ____________________________________________________
# Set up scheduling problem 

# initialize scheduler
setup_info = list(DT            = 5,
                  goal          = Goal,
                  tou_rates     = ratevec)

# initialize scheduler
scheduler = new(Class = "Scheduler", inputs, setup = setup_info)

# initialize options
eta   = rep(c(2,5,7,15,17,18,19), 3)
gamma = rep(c(2,3,4), each=7)
options = list(eta = eta, gamma = gamma, beta = 3)

Rprof()
scheduler = solveSchedulesApprox(scheduler, options = options)
Rprof(NULL)
tmp = summaryRprof()
