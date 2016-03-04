# describe_pecan.r
#
# Performs a basic statistical description of the Pecan St dataset.
#
# Adrian Albert
# Last modified: June 2014.
# ---------------------------------------------------------

rm(list = ls())
options(error = recover)
library('parallel')
library('lubridate')
library('ggplot2')
library('reshape')

# __________________________________________________
# Initializations...

setwd('~/EnergyAnalytics/batch/pecan/')
source('../../utils/aggregate_data.r')
source('define_categories_pecan.r')

# define some paths
DATA_PATH_SEL = '~/S3L_server/energy-data/pecan_street/usage-select/'
METADATA_PATH = '~/S3L_server/energy-data/pecan_street/metadata/'
PLOTS_PATH    = '~/S3L_server/plots/pecan-street-2'
dir.create(PLOTS_PATH)

# __________________________________________________
# Access usage data

# list all data files for the selected subset of data
# files_01 = list.files(path=paste(DATA_PATH_SEL, '/01min', sep=''), full.names = T, recursive = T)
# files_15 = list.files(path=paste(DATA_PATH_SEL, '/15min', sep=''), full.names = T, recursive = T)
files_60 = list.files(path=paste(DATA_PATH_SEL, '/60min', sep=''), full.names = T, recursive = T)

# extract user IDs 
users_sel = sapply(files_60, function(s) strsplit(tail(strsplit(s, '/')[[1]], 1), '\\.')[[1]][1])

# load IDs and names
usr_name = read.csv(paste(METADATA_PATH, 'user_names_ids.csv', sep = '/'))

# __________________________________________________
# Define plotting functions

# color codes for different usage types
cols = c(total= '#000000', C = '#FE2E2E', H = '#0040FF', user = '#088A08', #nonHVAC = '#5F04B4',
         lights = '#B45F04', always_on = '#01A9DB', scheduled = '#FF00BF', Temperature = '#424242')

# plot ground truth components
plot_user = function(homeData, main = 'minute') {
  
  maxUsage = max(homeData$use, na.rm=T)

  # plot different usages
  plot(homeData$date,homeData$use,type='l',col=cols['total'],
       main=main,ylab='kW',
       ylim=c(0,1.1*maxUsage),xaxt='n',cex=1.5)
  for (var in setdiff(
    names(homeData), c('total', 'date', "nonHVAC"))) {
    points(homeData$date,homeData[,var],type='l',col=cols[var])
  }
  
  # if plotting temperature, adjust values to plot on same axis
  if ('Temperature' %in% names(homeData))
    points(homeData$date, homeData$Temperature/100*maxUsage,type='l',col=cols['Temperature'])
  
  d = homeData$date
  dts = seq(d[1],d[length(d)],by=3600*24) # one per day from one per minute
  axis(1, dts, format(dts, "%a, %m/%d"), cex.axis=1)
  grid(nx=NA,ny=NULL)
  abline(v=dts,col="black",lty=3)  
}

# plot usage by given factor
plot_comps = function(dat, var = 'Hour', title = 'Components', facet = NULL) {
  
  # construct plot
  dat = melt(dat, id.vars = c(var, facet))
  
  p = ggplot(dat, aes_string(x = var, y = 'value'))
  p = p + geom_area(aes(fill = variable, color = variable), position = 'stack')    
  if (!is.null(facet)) 
    p = p + facet_wrap(as.formula(paste("~", facet)), ncol=1, scales = 'fixed')
  p = p + scale_color_manual(values = cols) + scale_fill_manual(values = cols)
  p = p + theme_bw() +
    theme(panel.grid.major  = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          strip.text.x     = element_text(size=18),
          axis.title.x     = element_text(size = 18),
          axis.text.y      = element_text(size=18), 
          axis.text.x      = element_text(size=18),
          axis.title.y     = element_text(size=18),
          axis.title.x     = element_text(size=18),
          plot.title       = element_text(size=20),            
          legend.text      = element_text(size=18),        
          legend.title     = element_blank(),    
          legend.position  = c(0.17, 0.4),
          legend.direction="horizontal",
          axis.ticks = element_blank() ) + 
    theme(plot.title=element_text(family="Times", face="bold", size=20)) + 
    ggtitle(title)  + xlab(var) + ylab('Avg. kWh')
  
  return(p)
}

# # plot on-percentages by different factors
# plot_on = function(dat, var = 'Hour', title = 'ON-State', facet = NULL) {
#   
#   # construct plot
#   dat = melt(dat, id.vars = c(var, facet))
#   
#   p = ggplot(dat, aes_string(x = var, y = 'value', group = 'variable'))
#   p = p + geom_point(aes_string(shape = 'variable'), size=4) + scale_shape(solid=FALSE)
#   p = p + geom_line(size = 2, aes_string(color = 'variable'))
#   if (!is.null(facet)) 
#     p = p + facet_wrap(as.formula(paste("~", facet)), ncol=1, scales = 'free')
#   p = p + theme_bw() +
#     theme(panel.grid.major  = element_blank(),
#           panel.grid.minor = element_blank(),
#           panel.background = element_blank(),
#           strip.text.x     = element_text(size=18),
#           axis.title.x     = element_text(size = 18),
#           axis.text.y      = element_text(size=18), 
#           axis.text.x      = element_text(size=18),
#           axis.title.y     = element_text(size=18),
#           axis.title.x     = element_text(size=18),
#           plot.title       = element_text(size=20),            
#           legend.text      = element_text(size=18),        
#           legend.title     = element_text(size=18),    
#           legend.position  = c(0.1, 0.25),
#           axis.ticks = element_blank() ) + 
#     theme(plot.title=element_text(family="Times", face="bold", size=20)) + 
#     ggtitle(title)  + xlab(var) + ylab('Avg. kWh')
#   
#   return(p)
# }
# 

# __________________________________________________
# Generate plots for each user

i = 138 # Nelson
res = mclapply(1:length(users_sel), 
               mc.cores=4,
               function(i){

# for (i in 1:length(users_sel)) {
  # current user info
  user   = subset(usr_name, ID == users_sel[i])

  # create folder for plots
  dir.create(file.path(PLOTS_PATH, user$ID), recursive = TRUE)
    
  # read in data at different resolutions...
  cat(paste(i, '/', length(users_sel), ':', user$ID, '...\n'))
#   data_01 = read.csv(files_01[i]); data_01$date = as.POSIXct(data_01$date); 
#   data_15 = read.csv(files_15[i]); data_15$date = as.POSIXct(data_15$date); 
  data_60 = read.csv(files_60[i]); data_60$date = as.POSIXct(data_60$date); #data_60 = na.omit(data_60)   
  
  # skip users that don't have enough data, at least a month of hourly data
  if (nrow(data_60) < 30 * 24) {
    cat('Not enough samples (less than one month!) !\n')
    next
  }
    
  # Compute basic stats about end-uses
  # --------------------------------------------
  
  # prepare data...
  df = cbind(data_60, Hour = hour(data_60$date), Month = month(data_60$date, label=T))
  df$Time = 'Night'
  df$Time[which(df$Hour %in% 5:11)]  = 'Morning'
  df$Time[which(df$Hour %in% 12:18)] = 'Afternoon'
  df$Time[which(df$Hour %in% 19:21)] = 'Evening'
  df$Season = 'Winter'
  df$Season[which(df$Month %in% c('Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep'))] = 'Summer'
  # massage components to increase interpretability
  comps = c('user', 'light', 'always_on', 'scheduled')
  if (length(intersect(comps, names(df)))>0) df$nonHVAC = df$nonHVAC - add.columns(df, comps)
    
  # consumption by time of day  
  cur_mu_time = aggregate(. ~ Time + Season, FUN = function(x) mean(x, na.rm=T),
                          data = df[,-which(names(df) %in% c('date', 'Temperature', 'Hour', 'Month', 'total'))])
  cur_mu_hour = aggregate(. ~ Hour + Season, FUN = function(x) mean(x, na.rm=T),
                          data = df[,-which(names(df) %in% c('date', 'Temperature', 'Time', 'Month', 'total'))])
  cur_mu_time = melt(cur_mu_time, id.vars = c('Time', 'Season'))
  cur_mu_time$ID = user$ID; cur_mu_time$name = user$name;

  # are AC and HV on at the same time?
  dfs = df[,c('Season', 'Hour', 'Time')]
  if ('C' %in% names(df) && 'H' %in% names(df)) 
    dfs$C_H_ON = df$C > 0 & df$H > 0
  if ('C' %in% names(df)) {
    dfs$C_ON = df$C > 0    
  }
  if ('H' %in% names(df)) {
    dfs$H_ON = df$H > 0
  }
  cur_on_hour = aggregate(. ~ Hour + Season, FUN = function(x) mean(x, na.rm=T), data = dfs[,-which(names(dfs) == 'Time')])
  cur_on_time = aggregate(. ~ Time + Season, FUN = function(x) mean(x, na.rm=T), data = dfs[,-which(names(dfs) == 'Hour')])
  cur_on_time = melt(cur_on_time, id.vars = c('Time', 'Season'))
  cur_on_time$ID = user$ID; cur_on_time$name = user$name;
  
  # Plot breakdown of avg usage by end use & hr
  # --------------------------------------------
  
  p.hr = plot_comps(cur_mu_hour[intersect(names(cur_mu_hour),
                      c('Hour', 'Season', 'H', 'C', 'user'))], 
                    var = 'Hour', facet = 'Season',
                    title = paste('Avg. Hourly Usage for', user$name))
  
  pdf(file=paste(paste(PLOTS_PATH, user$ID, sep = '/'), paste(user$ID,'_components.pdf',sep=''), sep='/'),width=8,height=6)
  print(p.hr)
  dev.off()
  
#   # Plot breakdown of simultaneous ON end-use
#   # --------------------------------------------
#   
#   p.on = plot_on(cur_on_hour, var = 'Hour', facet = 'Season',
#                     title = paste('Simultaneous ON-time, User', user$name, 'ID', user$ID))
#   
#   pdf(file=paste(paste(PLOTS_PATH, user$ID, sep = '/'), paste(user$ID,'_simultaneous_on.pdf',sep=''), sep='/'),width=8,height=6)
#   print(p.on)
#   dev.off()
  
#   # Plot usage by minute & by hour for each user
#   # --------------------------------------------
#   
#   pdf(file=paste(paste(PLOTS_PATH, user$ID, sep = '/'), paste(user$ID,'_resolutions.pdf',sep=''), sep='/'),width=10,height=6)
#   
#   # define plot layout
#   op <- par(no.readonly = TRUE)
#   m <- matrix(c(1,2,3,4),nrow=4,ncol=1,byrow=T)
#   layout(mat = m,heights = c(0.3,0.3,0.3,0.1))
#   par(oma=c(2,2,2,0),mar=c(2,4,2,1)) # Room for the title
#   
#   # select a random starting point, make sure data is ok for plotting  
#   ok = FALSE
#   while (!ok) {
#     start_date = sample(data_60$date[1:(nrow(data_60)-24*7)], 1)
# #     sel_01 = subset(data_01, date >= start_date & date < start_date + 7*24*3600)
# #     sel_15 = subset(data_15, date >= start_date & date < start_date + 7*24*3600)
#     sel_60 = subset(data_60, date >= start_date & date < start_date + 7*24*3600)
# #   if (nrow(na.omit(sel_01)) > 0 & nrow(na.omit(sel_15)) > 0 & nrow(na.omit(sel_60)) > 0 ) 
#     if (nrow(na.omit(sel_60)) > 0 ) 
#     ok = TRUE
#   }      
#     
# #   # print minute-by-minute data
# #   print(plot_user(sel_01, main = 'minute'))
# #   # print 15 min data
# #   print(plot_user(sel_15, main = '15 minute'))
#   # print hourly data
#   print(plot_user(sel_60, main = 'hourly'))
#   
#   mtext(paste('Pecan Street Experiment User', user$name, ', ID =', user$ID), line=0, font=2, cex=1.2, outer=TRUE)
#   par(mar=c(1,4,1,1))
#   plot.new()
#   vars_sel = setdiff(names(data_60), 'date')
#   legend("center", lty=1,cex=1,lwd=2, legend=vars_sel, col=cols[vars_sel],horiz=T)
# 
#   # save plot to file
#   dev.off()

  # Plot usage by minute & by hour for each user
  # --------------------------------------------

#   p = plot_temperature_profile(data_60, user)
#   
#   nVars = ncol(data_60) - 2
#   pdf(file=paste(file.path(PLOTS_PATH, user$ID), '/', 'temperature_profiles.pdf',sep=''),
#       width=nVars * 18/5, height=4)
#   print(p)
#   dev.off()
  
  
  # types of end uses
  end_uses = setdiff(names(data_60), c('date', 'total', 'Temperature'))
  
  return(list(end_uses = end_uses, cur_mu_time = cur_mu_time, cur_on_time = cur_on_time))
})

# idx_err_res  = which(sapply(res, length) != 3)
# res          = res[-idx_err_res]
# end_uses     = lapply(res, function(l) l[[1]])
# avg_cons_all = lapply(res, function(l) l[[2]])
# on_simul_all = lapply(res, function(l) l[[3]])
# 
# # how many users have different end uses monitored?
# end_uses_tab = table(unlist(end_uses))
# 
# # form analysis dataset for entire population
# avg_cons_all = do.call('rbind', avg_cons_all)
# on_simul_all = do.call('rbind', on_simul_all)
# 
# # save overall data
# cat('Saving stats...')
# save(file = paste(METADATA_PATH, '/stats.RData', sep = ''), 
#      list = c('end_uses', 'avg_cons_all', 'on_simul_all'))
# cat(' done!\n')
# 
# # _______________________________________________________
# # Plot distribution of appliance categories across users
# 
# df = as.data.frame(end_uses_tab / length(end_uses))
# names(df) = c('Type', 'Percentage')
# p = ggplot(df, aes(x = as.factor(Type), y = Percentage, fill = Type))
# p = p + geom_bar()
# p = p + scale_color_manual(values = cols) + scale_fill_manual(values = cols)
# p = p + theme_bw() +
#   theme(panel.grid.major  = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.background = element_blank(),
#         strip.text.x     = element_text(size=18),
#         strip.text.y     = element_text(size=18),
#         axis.title.x     = element_text(size = 18),
#         axis.text.y      = element_text(size=18), 
#         axis.text.x      = element_text(size=18),
#         axis.title.y     = element_text(size=18),
#         axis.title.x     = element_text(size=18),
#         plot.title       = element_text(size=20),            
#         legend.text      = element_text(size=18),        
#         legend.title     = element_text(size=18),    
#         legend.position  = 'none',
#         axis.ticks = element_blank() ) + 
#   theme(plot.title=element_text(family="Times", face="bold", size=20)) + 
#   ggtitle('Consumption by End-Use')  + xlab('End Use Type') + ylab('Percentage of Users Monitored')
# 
# pdf(file=paste(PLOTS_PATH, '/appliances_monitored.pdf',sep=''),width=9.3,height=4)
# print(p)
# dev.off()
# 
# # _______________________________________________________________
# # Plot distribution of consumption by end use and by time of day
# 
# dat = subset(avg_cons_all, variable %in% c('C', 'H', 'nonHVAC', 'lights'))
# dat$Time <- factor(dat$Time, c('Morning', 'Afternoon', 'Evening', 'Night'))
# idx = which(abs(dat$value)>2)
# dat[idx,'value'] = NA
# p = ggplot(dat, aes(x = value, color = Season))
# p = p + geom_density(size = 1.5)
# p = p + geom_density()
# p = p + facet_grid(~variable)
# # p = p + scale_color_manual(values = cols) + scale_fill_manual(values = cols)
# p = p + theme_bw() +
#   theme(panel.grid.major  = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.background = element_blank(),
#         strip.text.x     = element_text(size=18),
#         strip.text.y     = element_text(size=18),
#         axis.title.x     = element_text(size = 18),
#         axis.text.y      = element_text(size=18), 
#         axis.text.x      = element_text(size=18),
#         axis.title.y     = element_text(size=18),
#         axis.title.x     = element_text(size=18),
#         plot.title       = element_text(size=20),            
#         legend.text      = element_text(size=18),        
#         legend.title     = element_text(size=18),    
#         legend.position  = c(0.85, 0.55),
#         axis.ticks = element_blank() ) + 
#   theme(plot.title=element_text(family="Times", face="bold", size=20)) + 
#   ggtitle('Consumption by End-Use')  + xlab('kWh') + ylab('pdf')
# 
# pdf(file=paste(PLOTS_PATH, '/consumption_breakdown_population.pdf',sep=''),width=8,height=3)
# print(p)
# dev.off()
# 
# # _________________________________________________________________
# # Plot distribution of AC+HV on at the same timeand by time of day
# 
# dat = on_simul_all
# dat$Time <- factor(dat$Time, c('Morning', 'Afternoon', 'Evening', 'Night'))
# p = ggplot(dat, aes(x = Time, y = value, color = variable))
# p = p + geom_boxplot()
# p = p + facet_grid(~Season)
# # p = p + scale_color_manual(values = cols) + scale_fill_manual(values = cols)
# p = p + theme_bw() +
#   theme(panel.grid.major  = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.background = element_blank(),
#         strip.text.x     = element_text(size=18),
#         strip.text.y     = element_text(size=18),
#         axis.title.x     = element_text(size = 18),
#         axis.text.y      = element_text(size=18), 
#         axis.text.x      = element_text(size=18),
#         axis.title.y     = element_text(size=18),
#         axis.title.x     = element_text(size=18),
#         plot.title       = element_text(size=20),            
#         legend.text      = element_text(size=18),        
#         legend.title     = element_text(size=18),    
#         axis.ticks = element_blank() ) + 
#   theme(plot.title=element_text(family="Times", face="bold", size=20)) + 
#   ggtitle('Percentage of hours when in use')  + xlab('Time of Day') + ylab('% hours ON')
# 
# pdf(file=paste(PLOTS_PATH, '/simultaneous_on.pdf',sep=''),width=14,height=3)
# print(p)
# dev.off()
# 
