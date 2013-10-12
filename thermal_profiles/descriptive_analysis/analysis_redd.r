#!/usr/bin/Rscript
# analysis_redd.r
#
# HMM analysis of REDD data
# 
# - 
# Adrian Albert
# Last modified: May 2013.

# -----------------------------------------------
# Define constants, load libraries and functions
# -----------------------------------------------

rm(list = ls())
options(error = recover)

# load classes, utils
library('ggplot2')
library('lubridate')

setwd('~/Dropbox/OccupancyStates/')

plots_path_root = 'plots/redd_CA/'
if (!file.exists(plots_path_root)) dir.create(plots_path_root, recursive = TRUE)
data_path  = 'data/redd/redd_pre_processed/'
wthr_path  = 'data/redd/weather_processed/'

# load zip-to-house table
house_zip          = read.csv('data/redd/redd_zips.csv', colClasses = rep('character',4))

# choose some devices of interest
sel.devices = c('mains', 'furnace', 'space_heater', 'floor_heater', 
                'electric_heat', 'furance', 'air_conditioning' )

# choose weather vars of interest
wthr_vars   = c('TemperatureF')

# load data
csv.files  = list.files(path = data_path, pattern = '*.csv', full.names = T)
redd       = list()
for (f in csv.files) {
  # read electricity data
  cur_kwh  = read.csv(f)  
  cur_kwh$date = as.POSIXct(cur_kwh$DateTime)
  cur_kwh$DateTime = NULL
  
  # read in weather data
  idx      = unlist(sapply(as.character(house_zip$house_id), grep, f))
  house_id = names(idx)
  cur_zip  = house_zip$zip[which(house_zip$house_id == house_id)]
  cur_wthr = read.csv(paste(wthr_path, house_id, '_', sprintf('%05s',cur_zip), '_processed.csv', sep=''))
  cur_wthr$date = as.POSIXct(cur_wthr$date)
  
  # align weather and energy data to same time vector
  time_vec = intersect(cur_kwh$date, cur_wthr$date)
  cur_data = subset(cur_kwh, date %in% time_vec)
  sel_wthr = data.frame(subset(cur_wthr, date %in% time_vec, select = wthr_vars))
  names(sel_wthr) = wthr_vars
  cur_data = cbind(cur_data, sel_wthr)
  
  # add types of devices together to obtain functional aggregates
  agg = lapply(sel.devices, function(s) {
    idx = grep(s, names(cur_data))    
    if (length(idx)>0) {
      if (length(idx)>1) df = rowSums(cur_data[,idx]) else df = cur_data[,idx]
      df = data.frame(df)
      names(df) = s
    } else df = NULL
    return(df)
  })
  agg              = agg[!sapply(agg, is.null)]
  agg              = do.call('cbind', agg)
  agg$date         = cur_data$date
  agg              = cbind(agg, sel_wthr)
  
  # basic cleaning
  agg = subset(agg, TemperatureF >= 0)
#   for (i in setdiff(names(cur_data),'date')) {
#     x = cur_data[,i]
#     q = quantile(x, probs = c(0.005, 0.995))
#     x[which(x <= q[1] | x >= q[2])] = NA
#     cur_data[,i] = x
#   }
#   cur_data <- cur_data[,colSums(is.na(cur_data))<nrow(cur_data)]
#   cur_data = na.omit(cur_data)
  
  redd[[house_id]] = agg
}

# house 26 looks like it has the furnace on a separate circuit than the main?
redd[['house_26']]$mains = redd[['house_26']]$mains + redd[['house_26']]$furnace

# --------------------------------------------
# This portion is a hack for MoMo power data
# save everything into a file 
all.data = lapply(redd, function(l) subset(l, select = c('date', 'mains', 'TemperatureF')))
no.hrs   = sapply(redd, nrow)
all.data = do.call('rbind', all.data)
all.data$ZIPCODE = rep(house_zip$zipcode, times = no.hrs)

# read lat/long of zipcodes
geo = read.csv('data/redd/free-zipcode-database-Primary.csv')
geo = subset(geo, select = c('Zipcode', 'Lat', 'Long'))
all.data = merge(all.data, geo, by.x = 'ZIPCODE', by.y = 'Zipcode')

# write to file
write.csv(all.data, '/home/toni/Dropbox/WellDone/WellDone-Data-Simulator-Dashboard/data/power_weather_sample.csv')

# --------------------------------------------

# -----------------------------------------------
# Exploratory analysis of REDD data
# -----------------------------------------------

for (hid in names(redd)) {
  
  # create folder for each test case
  plots_path = paste(plots_path_root, hid, '/', sep='')
  if (!file.exists(plots_path)) dir.create(plots_path, recursive = TRUE)
  
  # format data
  df  = redd[[hid]]
  df1 = melt(df[,-which(names(redd[[hid]]) == 'date')], id.vars = c('TemperatureF'))
  df2 = melt(df, id.vars = c('date'))  
  
  # _____________________________
  # scatterplot with temperature
  
  plt = ggplot(df1, aes(TemperatureF, value, color = variable)) + geom_point(size = 2)
  plt = plt + facet_wrap(~variable, scales = 'free')
  plt = plt + theme_bw() +
    theme(panel.grid.major  = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          strip.text.x     = element_text(size=15),
          axis.title.x     = element_text(size = 18),
          axis.text.y      = element_text(size=15), 
          axis.text.x      = element_text(size=18),
          axis.title.y     = element_text(size=18),
          axis.title.x     = element_text(size=18),
          plot.title       = element_text(size=20),            
          legend.position  = 'none',
          axis.ticks = element_blank() ) + 
    theme(plot.title=element_text(family="Times", face="bold", size=20)) + 
    ggtitle( paste('Temperature Response:', hid )) + ylab('kWh')
  
  png(paste(plots_path,hid,'_temperature_profile.png', sep=''), width = 1600, height = 600, res = 200)
  print(plt)
  dev.off()
  
  # _________________________________
  # Time series view 
  
  plt = ggplot(df2, aes(date, value, color = variable)) + geom_point(size = 1.5) + geom_line()
  plt = plt + facet_wrap(~variable, scales = 'free', ncol = 1)
  plt = plt + theme_bw() +
    theme(panel.grid.major  = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          strip.text.x     = element_text(size=15),
          axis.title.x     = element_text(size = 18),
          axis.text.y      = element_text(size=15), 
          axis.text.x      = element_text(size=18),
          axis.title.y     = element_text(size=18),
          axis.title.x     = element_text(size=18),
          plot.title       = element_text(size=20),            
          legend.position  = 'none',
          axis.ticks = element_blank() ) + 
    theme(plot.title=element_text(family="Times", face="bold", size=20)) + 
    ggtitle( paste('Consumption Readings:', hid) )
  
  png(paste(plots_path,hid,'temperature_time_series.png', sep=''), width = 1800, height = 1400, res = 200)
  print(plt)
  dev.off()
  
}

# -----------------------------------------------
# Perform analysis on individual houses
# -----------------------------------------------

source('code/OccupancyStates.r')
source('code/occupancyAnalysis.r')

simulation = list()
for (i in names(redd)[c(1,2,3,5,6)]) {
  
  # create folder for each test case
  plots_path = paste(plots_path_root, i, '/', sep='')
  if (!file.exists(plots_path)) dir.create(plots_path, recursive = TRUE)

  # add mains to obtain 
  cur_data = redd[[i]][, c('date', 'mains')]
  names(cur_data) = c('date', 'kWh')
  cur_data$date = as.character(cur_data$date)
  cur_wthr = redd[[i]][, c('date', 'TemperatureF')]
  cur_wthr$date = as.character(cur_wthr$date)
  simulation[[i]] = occupancyAnalysis(cur_data, cur_wthr, i, 1, 
                                      Kmin = 2, Kmax = 4, 
                                      verbose = T, 
                                      NOBS_THRESH = 0, 
                                      resp.vars = c('(Intercept)', 'TemperatureF'),
                                      tran.vars = c('(Intercept)', 'TemperatureF'),
                                      addl.vars = c(),
                                      plots_path = plots_path, dump_path = plots_path)  
}

# --------------------------------------------------
# Compute appliance usage per state for each house
# --------------------------------------------------

library('lubridate')
for (i in names(simulation)) {
  
  # create folder for each test case
  plots_path = paste(plots_path_root, i, '/', sep='')
  if (!file.exists(plots_path)) dir.create(plots_path, recursive = TRUE)
  
  # ________________________________
  # Select & format data of interest
  
  states             = simulation[[i]][[2]]$HMM$states$states
  cur_data           = redd[[i]][1:length(states),]
  cur_data           = subset(cur_data, select = c(intersect(names(cur_data), sel.devices), 'date'))
  hvac_vars          = setdiff(intersect(names(cur_data), sel.devices),'mains')
  if (length(hvac_vars)>1)
    kwh_hvac         = rowSums(cur_data[,hvac_vars]) else kwh_hvac = cur_data[,hvac_vars]
  cur_data$non.hvac  = cur_data$mains - kwh_hvac
  cur_data$state     = as.factor(states)
  cur_data$HourOfDay = hour(cur_data$date)
  cur_data$fit.Thermal = simulation[[i]][[2]]$HMM$comps$ts[,'TemperatureF']
  cur_data$fit.Activity= simulation[[i]][[2]]$HMM$comps$ts[,'Activity']
  
  # components 
  # cur_comp           = as.data.frame(abs(simulation[[i]][[2]]$HMM$comps$ts[,'TemperatureF']))
  cur_comp           = data.frame(Hard.Breakpoint = abs(simulation[[i]][[2]]$HMM$comps$hard.resp),
                                  Soft.Breakpoint = abs(simulation[[i]][[2]]$HMM$comps$ts[,'TemperatureF']))
  cur_comp$fit       = NULL
  cur_comp$kWh       = NULL
  cur_comp$HourOfDay = cur_data$HourOfDay
    
  # ________________________________________________________________
  # Plot proportion of HVAC/non-HVAC per hour of day (ground truth)
  
  df       = cur_data[,-which(names(cur_data)=='date')]
  df$mains = NULL
  df1      = df
  df1$state = NULL
  df1$fit.Thermal    = NULL
  df1$fit.Activity   = NULL
  df1 = melt(df1, id.vars = c('HourOfDay'))
  df1 = aggregate(value ~ variable + HourOfDay, data = df1, FUN = sum)
  p1  = ggplot(df1, aes(x = factor(HourOfDay), y = value, fill = variable))
  p1  = p1 + geom_bar(width = 0.9, stat = 'identity') 
  p1  = p1 + theme_bw() +
    theme(panel.grid.major  = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          strip.text.x     = element_text(size=15),
          axis.title.x     = element_text(size = 18),
          axis.text.y      = element_text(size = 18), 
          axis.text.x      = element_text(size=18),
          axis.title.y     = element_text(size=18),
          axis.title.x     = element_text(size=18),
          plot.title       = element_text(size=20),            
          legend.text      = element_text(size=15),        
          legend.title     = element_text(size=15),    
          legend.direction = 'horizontal',
          legend.position  = c(0.8,0.9),
          axis.ticks = element_blank() ) + 
    theme(plot.title=element_text(family="Times", face="bold", size=20)) + 
    ggtitle( paste('Ground Truth:', i) )  + xlab('') + ylab('kWh')
  
  # Plot fit breakdown
  df       = cur_data[,-which(names(cur_data)=='date')]
  df$mains = NULL
  df2      = subset(df, select = c('HourOfDay', 'state'))
  df2      = as.data.frame(table(df2$HourOfDay, df2$state))
  names(df2) = c('HourOfDay', 'state', 'Count')
  p2   = ggplot(df2, aes(x = factor(HourOfDay), y = Count, fill = state))
  p2   = p2 + geom_bar(width = 0.9, stat = 'identity', position = 'stack') 
  p2   = p2 + theme_bw() +
    theme(panel.grid.major  = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          strip.text.x     = element_text(size=15),
          axis.title.x     = element_text(size = 18),
          axis.text.y      = element_text(size = 18), 
          axis.text.x      = element_text(size=18),
          axis.title.y     = element_text(size=18),
          axis.title.x     = element_text(size=18),
          plot.title       = element_text(size=20),            
          legend.text      = element_text(size=15),        
          legend.title      = element_text(size=15),        
          legend.direction = 'horizontal',
          legend.position  = c(0.7,0.9),
          axis.ticks = element_blank() ) + 
    theme(plot.title=element_text(family="Times", face="bold", size=20)) + 
    ggtitle( paste('Fit:', i) )  + xlab('Hour of Day') + ylab('Count')
  
  png(paste(plots_path, i, 'ground_truth_thermal.png', sep=''), width = 1400, height = 1000, res = 150)
  print(multiplot(p1,p2,cols=1))
  dev.off()
  
  # ________________________________________________________________
  # Plot temperature-dependent components
  
  # fit data
  df   = cur_data
  df$Thermal.Hard = cur_comp$Hard.Breakpoint
  df$Thermal.Soft = cur_comp$Soft.Breakpoint
  df$Day   = as.numeric(strftime(as.POSIXct(df$date), format = "%j"))
  df2 = subset(df, select = c('HourOfDay', 'Thermal.Hard', 'Thermal.Soft', hvac_vars, 'Day'))  
  agg.day = aggregate(data = df2, cbind(furnace, Thermal.Hard, Thermal.Soft) ~ Day, FUN = sum)
  agg.day = melt(agg.day, id.vars = 'Day')
  
  # ground truth data
  p3   = ggplot(agg.day, aes(x = Day, y = value, color = variable))
  p3   = p3 + geom_point(size = 2) + geom_line() 
#  p3   = p3 + facet_wrap(~as.factor(HourOfDay))
  p3   = p3 + theme_bw() +
    theme(panel.grid.major  = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          strip.text.x     = element_text(size=15),
          axis.title.x     = element_text(size = 18),
          axis.text.y      = element_text(size = 18), 
          axis.text.x      = element_text(size=18),
          axis.title.y     = element_text(size=18),
          axis.title.x     = element_text(size=18),
          plot.title       = element_text(size=20),            
          legend.text      = element_text(size=15),        
          legend.title      = element_text(size=15),        
          legend.direction = 'horizontal',
          legend.position  = c(0.7,0.9),
          axis.ticks = element_blank() ) + 
    theme(plot.title=element_text(family="Times", face="bold", size=20)) + 
    ggtitle( paste('Ground Truth vs Estimated Thermal Response:', i) )  + xlab('Hour of Day') + ylab('kWh')
  
  png(paste(plots_path, i, 'ground_truth_thermal_comparison.png', sep=''), width = 1400, height = 600, res = 150)
  print(p3)
  dev.off()
  
  # ________________________________________________________________
  # Plot fit time series vs ground
    
  df   = cur_data
  df$Thermal.Hard = cur_comp$Hard.Breakpoint
  df$Thermal.Soft = cur_comp$Soft.Breakpoint
  df$mains = NULL
  df = df[,-which(names(df) %in% c('non.hvac', 'fit.Thermal', 'fit.Activity', 'state', 'HourOfDay'))]
  df = df[200:350,]
  df.mlt = melt(df, id.vars = 'date')
  
  plt = ggplot(df.mlt, aes(x = date, y = value, color = variable)) + geom_line() + geom_point(size = 2)
  plt = plt + theme_bw() +
    theme(panel.grid.major  = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          strip.text.x     = element_text(size=15),
          axis.title.x     = element_text(size = 18),
          axis.text.y      = element_text(size = 18), 
          axis.text.x      = element_text(size=18),
          axis.title.y     = element_text(size=18),
          axis.title.x     = element_text(size=18),
          plot.title       = element_text(size=20),            
          legend.text      = element_text(size=15),        
          legend.title      = element_text(size=15),        
          legend.direction = 'horizontal',
          legend.position  = c(0.7,0.9),
          axis.ticks = element_blank() ) + 
    theme(plot.title=element_text(family="Times", face="bold", size=20)) + 
    ggtitle( paste('Ground Truth vs Estimated Thermal Response:', i) )  + xlab('Time') + ylab('kWh')
  
  png(paste(plots_path, i, 'ground_truth_thermal_time_series.png', sep=''), width = 1400, height = 600, res = 150)
  print(plt)
  dev.off()
    
  # __________________________________________________________________________________________
  # Compute # times soft & hard breakpoints agree w/ ground truth in the inter-quartile range
  
  # perform computations
  qG = quantile(cur_data[,hvac_vars], probs = c(0.25, 0.75))
  S_activity = cur_comp$Soft.Breakpoint >= qG[1]
  H_activity = cur_comp$Hard.Breakpoint >= qG[1] 
  G_activity = cur_data[,hvac_vars] >= qG[1] 
  agreement = data.frame(date = cur_data$date, 
                         HourOfDay = cur_comp$HourOfDay,
                         Ground.Truth    = as.numeric(G_activity), 
                         Hard.Breakpoint = as.numeric(H_activity), 
                         Soft.Breakpoint = as.numeric(S_activity))
  
  # compute agreement by time-of-day
  tab.tod = aggregate(data = agreement, 
                      cbind(Soft.Breakpoint, Hard.Breakpoint, Ground.Truth) ~ HourOfDay, FUN = mean)
  tab.ltx= round(t(tab.tod), digits = 2)
  
  # save to latex
  library(xtable)
  tab.ltx = xtable(tab.ltx, caption = paste('Average thermal kWh assignment (hour of day):', i),
                   label = paste('tab:',i,sep=''), digits = 2)
  print(tab.ltx, file = paste(plots_path, 'validation.tex', sep=''))

  # plot table
  tab.mlt = melt(tab.tod, id.vars = 'HourOfDay')
  p3   = ggplot(tab.mlt, aes(x = HourOfDay, y = value, color = variable))
  p3   = p3 + geom_point(size = 5) + geom_line(size=2)
  #  p3   = p3 + geom_bar(width = 0.9, stat = 'identity', position = 'dodge') 
  #  p3   = p3 + facet_wrap(~as.factor(HourOfDay))
  p3   = p3 + theme_bw() +
    theme(panel.grid.major  = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          strip.text.x     = element_text(size=15),
          axis.title.x     = element_text(size = 18),
          axis.text.y      = element_text(size = 18), 
          axis.text.x      = element_text(size=18),
          axis.title.y     = element_text(size=18),
          axis.title.x     = element_text(size=18),
          plot.title       = element_text(size=20),            
          legend.text      = element_text(size=15),        
          legend.title      = element_text(size=15),        
          legend.direction = 'vertical',
          legend.position  = c(0.2,0.3),
          axis.ticks = element_blank() ) + 
    theme(plot.title=element_text(family="Times", face="bold", size=20)) + 
    ggtitle( paste('Ground Truth vs Estimated Thermal Response:', i) )  + xlab('Hour of Day') + ylab('% Interquartile Match')
  
  png(paste(plots_path, i, 'ground_truth_interquartile_comparison.png', sep=''), width = 1400, height = 600, res = 150)
  print(p3)
  dev.off()
      
}

