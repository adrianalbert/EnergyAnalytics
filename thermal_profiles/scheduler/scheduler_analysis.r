# scheduler_analysis.r
#
# Scheduling analysis.
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

load('~/energy-data/bakersfield/bakersfield_profiles.RData')
load('~/energy-data/bakersfield/scheduler_full.RData')

source('../../clustering/kError.r', chdir = T)
source('classes/DataImporter.r', chdir = T)
source('classes/Scheduler.r', chdir = T)

PLOTS_PATH = '~/Dropbox/OccupancyStates/plots/scheduling/'

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
  print(paste('Clustering for K=',k))
  done = F; att = 0
  while (!done & att < 5) {
    fit = try(kError(Abar, W, k, iter = 100))
    done= class(fit) != 'try-error'
    if (!done) {
      print(paste('Error at K=', k, '... re-trying', att, '/5'))
      att = att + 1
    }
  }
  print(paste('...done! (Clustering for K=',k, ')'))
  
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
  options = list(budget = 3, gamma = 0.01, Nr = Nr, presaved = FALSE)
  scheduler.clust = solveSchedules(scheduler.clust, options = options)
  
  # compute aggregate by applying to each profile its class-based schedule
  r = compute_clust_agg(Abar, W, scheduler.clust@OUTPUT$U, fit$assignment, setup_info$goal, setup_info$tou_rates)
  
  PLOTS_PATH_CUR = paste(PLOTS_PATH, 'varying_K/K_', k, '/', sep = '')
  dir.create(PLOTS_PATH_CUR, recursive = T)
  
  # plot schedules per classs
  pdf(file = paste(PLOTS_PATH_CUR, 'cluster-effort-profiles.pdf', sep=''), height = 6, width = 14)
  print(plot(scheduler.clust, type = 'effort-profiles'))
  dev.off()
  
  # plot match between goal and aggregate
  pdf(file = paste(PLOTS_PATH_CUR, 'cluster-agg-goal.pdf', sep=''), height = 4, width = 6)
  print(plot(scheduler.clust, type = 'goal-match', compare = r))
  dev.off()
  
  # plot comparison of # users selected from each class
  pdf(file = paste(PLOTS_PATH_CUR, 'cluster-selected.pdf', sep=''), height = 4, width = 8)
  print(plot(scheduler.clust, type = 'number-selected'))
  dev.off()
  
  return(list(nr        = scheduler.clust@OUTPUT$nr,
              Delta.bar = r$mu,
              Delta.var = r$var,
              EC        = r$EC)) 
}

res = mclapply(2*(5:15), 
               mc.cores = 5, 
               perform_analysis_cluster)

# ___________________________________________
# fit vs no of clusters

# form data
EC    = sapply(res, function(l) l$EC)
Delta = sapply(res, function(l) l$Delta.bar)
nr    = sapply(res, function(l) sum(l$nr))
ref   = scheduler@OUTPUT$Delta.bar
nr.ref= sum(scheduler@OUTPUT$nr)
df    = as.data.frame(t(Delta)); names(df) = sprintf("%02d", 1:ncol(df))
df    = rbind(df, ref)
df$nr = c(nr, nr.ref); df$R = paste('R=',c(2*(5:15), 0))
df = melt(df, id = c('nr', 'R'))

# plot
plt = ggplot(df, aes(y = value, x = variable, color = R, group = R)) + 
  geom_line(size = 2) + geom_point(size = 4)
plt = plt + theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.text.x     = element_text(size=18),
        axis.text.y      = element_text(size=18), 
        axis.text.x      = element_text(size=18, angle = -20),
        axis.title.y     = element_text(size=18),
        axis.title.x     = element_text(size=18),
        plot.title       = element_text(size=20),            
        legend.text      = element_text(size=18),
        legend.title      = element_text(size=18),
        axis.ticks = element_blank()) + 
  ggtitle(paste("Aggregate and goal: cluster-based schedules")) + ylab('kWh') + xlab('Time of day')

# plot effort distribution
pdf(file = paste(PLOTS_PATH, 'varying_K/aggregate_vs_goal.pdf', sep=''), height = 4, width = 8)
  print(plt)
dev.off()


# ____________________________________________________
# Vary budget

# dir.create(file.path(paste(PLOTS_PATH, 'varying_beta/', sep='')))
# 
# # options 
# setup_info = list(DT            = 5,
#                   goal          = Goal,
#                   tou_rates     = ratevec)
# 
# # prepare data 
# cur_inputs = inputs
# Abar = do.call('rbind', lapply(cur_inputs, function(l) as.numeric(l$a$mu[,1])))
# rownames(Abar)= names(cur_inputs)
# Abar[which(Abar < 0)] = 0
# W = lapply(cur_inputs, function(l) diag(diag(l$a$covmat)))
# 
# # group together all analysis logic
# perform_analysis_budget = function(beta) {
#   
#   cat(paste('beta =', beta, '\n'))
#   
#   # initialize scheduler
#   scheduler.cur = new(Class = "Scheduler", cur_inputs, setup = setup_info)
#   
#   # compute schedules in the most general case (tailored schedule for each user)
#   if (beta == 1) presaved = TRUE else presaved = TRUE
#   options       = list(budget = beta,  presaved = presaved)
#   scheduler.cur = solveSchedules(scheduler.cur, options = options)
#   
#   PLOTS_PATH_CUR = paste(PLOTS_PATH, 'varying_beta/beta_', beta, '/', sep = '')
#   dir.create(PLOTS_PATH_CUR, recursive = T)
#   
#   # plot schedules per classs
#   pdf(file = paste(PLOTS_PATH_CUR, 'effort-profiles-', beta, '.pdf', sep=''), height = 3, width = 9)
#   print(plot(scheduler.cur, type = 'effort-profiles', selected = c('Robert', 'Joe', 'Lewis')))
#   dev.off()
#   
#   # plot match between goal and aggregate
#   pdf(file = paste(PLOTS_PATH_CUR, 'agg-goal.pdf', sep=''), height = 4, width = 6)
#   print(plot(scheduler.cur, type = 'goal-match'))
#   dev.off()
#   
#   # save data for later
#   nr        = sum(scheduler.cur@OUTPUT$nr)
#   EC        = scheduler.cur@OUTPUT$EC
#   Delta.bar = scheduler.cur@OUTPUT$Delta.bar
#   Delta.var = scheduler.cur@OUTPUT$Delta.var
#   
#   # manually clear memory since R is so bad with this
#   rm(list = 'scheduler.cur'); gc()
#   
#   return(list(nr        = nr,
#               Delta.bar = Delta.bar,
#               Delta.var = Delta.var,
#               EC        = EC)) 
# }
# 
# resbeta = lapply(1:20, 
#           #     mc.cores = 3, 
#                perform_analysis_budget)
# 
# # ____________________________________________________
# # Plots: budget analysis
# 
# # plots
# EC = sapply(resbeta, function(r) r$EC)
# nr = sapply(resbeta, function(r) r$nr)
# Db = sapply(resbeta, function(r) r$Delta.bar)
# 
# df  = data.frame(Nr.Selected = nr, Cost = EC, Budget = 1:20)
# 
# # cost vs budget
# library('pracma')
# pdf(file = paste(PLOTS_PATH, 'no_users_vs_budget.pdf', sep=''), height =4, width = 6)
# plotyy(df$Budget, df$Cost, 
#        df$Budget, df$Nr.Selected/100, gridp = TRUE, box.col = "grey",
#        type = "b", lwd = 3, lty = 1, cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5, 
#        xlab = "Effort Budget beta [deg F x 5 / 24hrs]", ylab = "Expected Cost", 
#        main = 'DR Program Administration',
#        col.y1 = "navy", col.y2 = "maroon")
# legend(7,300, # places a legend at the appropriate place 
#        c("Expected Penalty Cost","No. Users Enrolled (x100)"), # puts text in the legend               
#        lty=c(1,1), # gives the legend appropriate symbols (lines)       
#        lwd=c(2.5,2.5),col=c("navy","maroon"), cex = 1)
# dev.off()
# 
