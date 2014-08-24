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

# options for scheduler
setup_info = list(DT            = 5,
                  goal          = Goal,
                  tou_rates     = ratevec)

# ____________________________________________________
# Vary budget

dir.create(file.path(paste(PLOTS_PATH, 'varying_beta_approx/', sep='')))

# group together all analysis logic
perform_analysis_budget = function(beta) {
  
  cat(paste('beta =', beta, '\n'))
  
  # initialize scheduler
  scheduler.cur = new(Class = "Scheduler", inputs, setup = setup_info)
  
  # initialize options
  eta   = 2*(1:10)
  options = list(eta = eta, beta = beta, verbose = FALSE, NOBJ = 100)
  
  # compute schedules
  scheduler.cur = solveSchedulesApprox(scheduler.cur, options = options)  
  
  # save data for later
  nr        = sum(scheduler.cur@OUTPUT$nr)
  objvec    = scheduler.cur@OUTPUT$OBJVEC
  EC        = scheduler.cur@OUTPUT$EC
  Delta.bar = scheduler.cur@OUTPUT$Delta.bar
  Delta.var = scheduler.cur@OUTPUT$Delta.var
  Utot      = colSums(scheduler.cur@OUTPUT$U)
  
  PLOTS_PATH_CUR = paste(PLOTS_PATH, 'varying_beta_approx/beta_', beta, '/', sep = '')
  dir.create(PLOTS_PATH_CUR, recursive = T)
  
  # plot schedules per classs
  pdf(file = paste(PLOTS_PATH_CUR, 'effort-profiles-', beta, '.pdf', sep=''), height = 3, width = 9)
  print(plot(scheduler.cur, type = 'effort-profiles', selected = c('Robert', 'Joe', 'Lewis')))
  dev.off()
  
  # plot match between goal and aggregate with percentiles usr selected for comparison
  pdf(file = paste(PLOTS_PATH_CUR, 'agg-goal.pdf', sep=''), height = 4, width = 6)
  print(plot(scheduler.cur, type = 'goal-match-approx', compare = c(10, 30, 50, 90)))
  dev.off()
  
  # manually clear memory since R is so bad with this
  rm(list = 'scheduler.cur'); gc()
  
  return(list(nr        = nr,
              beta      = beta, 
              objvec    = objvec,
              Delta.bar = Delta.bar,
              Delta.var = Delta.var,
              EC        = EC,
              Utot      = Utot)) 
}

resbeta = mclapply(1:20, 
                   mc.cores = 5, 
                   perform_analysis_budget)

# get quantities for plots
EC = sapply(resbeta, function(r) r$EC)
nr = sapply(resbeta, function(r) r$nr)
Db = sapply(resbeta, function(r) r$Delta.bar)
Ob = sapply(resbeta, function(r) r$objvec)
beta = sapply(resbeta, function(r) r$beta)
Utot = sapply(resbeta, function(r) r$Utot)

# save data
save.image('~/energy-data/bakersfield/bakersfield_approx_sched.RData')

# ____________________________________________________
# Plots: budget analysis


library(pracma)

# cost vs budget
pdf(file = paste(PLOTS_PATH, 'no_users_vs_budget_approx.pdf', sep=''), height =4, width = 6)
df  = data.frame(Budget = beta, Nr.Selected = nr, Cost = EC)
plotyy(df$Budget, df$Cost, 
       df$Budget, df$Nr.Selected/100, gridp = TRUE, box.col = "grey",
       type = "b", lwd = 3, lty = 1, cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5, 
       xlab = "Effort Budget beta [deg F x 5 / 24hrs]", ylab = "Expected Cost", 
       main = 'DR Program Administration',
       col.y1 = "navy", col.y2 = "maroon")
legend(7,900, # places a legend at the appropriate place 
       c("Expected Penalty Cost","No. Users Enrolled (x100)"), # puts text in the legend               
       lty=c(1,1), # gives the legend appropriate symbols (lines)       
       lwd=c(2.5,2.5),col=c("navy","maroon"), cex = 1)
dev.off()

# cost vs no selected for selected budgets
pdf(file = paste(PLOTS_PATH, 'no_users_vs_cost_approx.pdf', sep=''), height =4, width = 7)
beta_sel = c(1, 3, 5, 10, 15)
cols = c("navy", "red", "black", "green", "magenta")
plot(1:length(Ob[[beta_sel[1]]]), Ob[[beta_sel[1]]],     
     ylim = range(unlist(Ob)),
     type = "b", lwd = 2, lty = 1, cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5, 
     pch = 1, 
     xlab = "No. Consumers Enrolled", ylab = "Expected Cost", 
     main = 'DR Program: Approx. Consumer Selection', col = cols[1])
for (b in 2:length(beta_sel)) {
  points(1:length(Ob[[beta_sel[b]]]), Ob[[beta_sel[b]]],    
         pch = b, 
         type = "b", lwd = 2, lty = 1, cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5, col = cols[b])  
}
legend(1200, 900,# places a legend at the appropriate place 
       sapply(beta_sel, function(b) as.expression(bquote(beta ~ "=" ~ .(b)))), # puts text in the legend     
       pch = 1:5, 
       fill = cols,
       lty=rep(1,length(beta_sel)), # gives the legend appropriate symbols (lines)       
       lwd=rep(2.5,length(beta_sel),col=cols, cex = 1))
dev.off()

# plot total number selected by time of day
beta_sel = c(3, 5, 10)
df = as.data.frame(Utot[,beta_sel]); names(df) = beta_sel; df = cbind(Hour = 1:nrow(df), df)
df = melt(df, id.vars = 'Hour'); names(df)[c(2,3)] = c('Budget', 'No.Selected')

plt = ggplot(df, aes(x = Hour, y = No.Selected, color = Budget, shape = Budget)) + 
  geom_point(size = 2.5) + geom_line(size=1.5)
plt = plt + scale_color_manual(values = cols[2:4])
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
        legend.title      = element_text(size=18),
        legend.position  = c(0.45, 0.75),
        axis.ticks = element_blank()) + 
  ggtitle(paste("Hour-of-day profile of selected customers")) + ylab('# Customers') + xlab('Hour of Day')
pdf(file = paste(PLOTS_PATH, 'no_users_profile.pdf', sep=''), height =4, width = 7)
print(plt)
dev.off()
