# describe_pecan.r
#
# Performs a basic statistical description of the Pecan St dataset.
#
# Adrian Albert
# Last modified: May 2014.
# ---------------------------------------------------------

rm(list = ls())
options(error = recover)
library('lubridate')
library('ggplot2')
library('reshape')

# __________________________________________________
# Initializations...

setwd('~/EnergyAnalytics/batch/pecan/')
source('../../utils/aggregate_data.r')
source('../../batch/pecan/define_categories_pecan.r')

# define some paths
DATA_PATH_SEL = '~/energy-data/pecan_street/usage-select/'
METADATA_PATH = '~/energy-data/pecan_street/metadata/'
PLOTS_PATH    = '~/Dropbox/OccupancyStates/plots/pecan-street'
dir.create(PLOTS_PATH)

# __________________________________________________
# Access usage data

# list all data files for the selected subset of data
files_01 = list.files(path=paste(DATA_PATH_SEL, '/01min', sep=''), full.names = T, recursive = T)
files_15 = list.files(path=paste(DATA_PATH_SEL, '/15min', sep=''), full.names = T, recursive = T)
files_60 = list.files(path=paste(DATA_PATH_SEL, '/60min', sep=''), full.names = T, recursive = T)

# extract user IDs 
users_sel = sapply(files_60, function(s) strsplit(tail(strsplit(s, '/')[[1]], 1), '\\.')[[1]][1])
year_sel  = sapply(files_60, function(s) tail(strsplit(s, '/')[[1]], 2)[1])

# load IDs and names
usr_name = read.csv(paste(METADATA_PATH, 'user_names_ids.csv', sep = '/'))

# __________________________________________________
# Define plotting functions

# color codes for different usage types
cols = c(use= '#000000', AC = '#FE2E2E', HV = '#0040FF', user = '#088A08', nonHVAC = '#5F04B4',
         lights = '#B45F04', always_on = '#01A9DB', scheduled = '#FF00BF', TemperatureF = '#424242')

# plot ground truth components
plot_user = function(homeData, main = 'minute') {
  
  maxUsage = max(homeData$use)

  # plot different usages
  plot(homeData$date,homeData$use,type='l',col=cols['use'],main=main,ylab='kW',ylim=c(0,1.1*maxUsage),xaxt='n')
  for (var in setdiff(names(homeData), c('use', 'date'))) {
    points(homeData$date,homeData[,var],type='l',col=cols[var])
  }
  
  # if plotting temperature, adjust values to plot on same axis
  if ('TemperatureF' %in% names(homeData))
    points(homeData$date, homeData$TemperatureF/100*maxUsage,type='l',col=cols['TemperatureF'])
  
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
    p = p + facet_wrap(as.formula(paste("~", facet)), ncol=1, scales = 'free')
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
          legend.title     = element_text(size=18),    
          legend.position  = c(0.2, 0.5),
          axis.ticks = element_blank() ) + 
    theme(plot.title=element_text(family="Times", face="bold", size=20)) + 
    ggtitle(title)  + xlab(var) + ylab('Avg. kWh')
  
  return(p)
}


# __________________________________________________
# Generate plots for each user
end_uses = list()
avg_cons = list()
for(i in 1:length(users_sel)){

  # current user info
  user   = subset(usr_name, ID == users_sel[i])
  year   = year_sel[i]

  # create folder for plots
  dir.create(file.path(PLOTS_PATH, year, user$ID), recursive = TRUE)
    
  # read in data at different resolutions...
  cat(paste(user$ID, '...\n'))
  data_01 = read.csv(files_01[i]); data_01$date = as.POSIXct(data_01$date); 
  data_15 = read.csv(files_15[i]); data_15$date = as.POSIXct(data_15$date); 
  data_60 = read.csv(files_60[i]); data_60$date = as.POSIXct(data_60$date); data_60 = na.omit(data_60)   
  
  # Compute basic stats about end-uses
  # --------------------------------------------
  
  # types of end uses
  end_uses[[i]] = setdiff(names(data_60), c('date', 'use', 'TemperatureF'))
  
  # consumption by time of day
  df = cbind(data_60, Hour = hour(data_60$date), Month = month(data_60$date, label=T))
  df$Time = 'Night'
  df$Time[which(df$Hour %in% 5:12)]  = 'Morning'
  df$Time[which(df$Hour %in% 13:18)] = 'Afternoon'
  df$Time[which(df$Hour %in% 19:21)] = 'Evening'
  df$Season = 'Fall-Winter'
  df$Season[which(df$Month %in% c('Jun', 'Jul', 'Aug', 'Sep'))] = 'Summer'
  df$Season[which(df$Month %in% c('Feb', 'Mar', 'Apr', 'May'))] = 'Spring'
  cur_mu_hour = aggregate(. ~ Hour + Season, FUN = mean,
                          data = df[,-which(names(df) %in% c('date', 'TemperatureF', 'Time', 'Month', 'use'))])  
  cur_mu_time = aggregate(. ~ Time + Season, FUN = mean,
                          data = df[,-which(names(df) %in% c('date', 'TemperatureF', 'Hour', 'Month', 'use'))])
  cur_mu_time = melt(cur_mu_time, id.vars = c('Time', 'Season'))
  cur_mu_time$ID = user$ID; cur_mu_time$name = user$name; cur_mu_time$year = year
  avg_cons[[i]] = cur_mu_time
  
  # Plot breakdown of avg usage by end use & hr
  # --------------------------------------------
  
  p.hr = plot_comps(cur_mu_hour, var = 'Hour', facet = 'Season',
                    title = paste('Avg. Hourly Usage, User', user$name, 'ID', user$ID))
    
  pdf(file=paste(paste(PLOTS_PATH, year, user$ID, sep = '/'), paste(user$ID,'_components.pdf',sep=''), sep='/'),width=10,height=6)
  print(p.hr)
  dev.off()
  
  # Plot usage by minute & by hour for each user
  # --------------------------------------------
  
  pdf(file=paste(paste(PLOTS_PATH, year, user$ID, sep = '/'), paste(user$ID,'_resolutions.pdf',sep=''), sep='/'),width=10,height=6)
  
  # define plot layout
  op <- par(no.readonly = TRUE)
  m <- matrix(c(1,2,3,4),nrow=4,ncol=1,byrow=T)
  layout(mat = m,heights = c(0.3,0.3,0.3,0.1))
  par(oma=c(2,2,2,0),mar=c(2,4,2,1)) # Room for the title
  
  # select a random starting point, make sure data is ok for plotting  
  ok = FALSE
  while (!ok) {
    start_date = sample(data_60$date[-((nrow(data_60)-24*10):nrow(data_60))], 1)
    sel_01 = subset(data_01, date >= start_date & date < start_date + 10*24*3600)
    sel_15 = subset(data_15, date >= start_date & date < start_date + 10*24*3600)
    sel_60 = subset(data_60, date >= start_date & date < start_date + 10*24*3600)
    if (nrow(na.omit(sel_01)) > 0 & nrow(na.omit(sel_15)) > 0 & nrow(na.omit(sel_60)) > 0 ) 
      ok = TRUE
  }      
    
  # print minute-by-minute data
  print(plot_user(sel_01, main = 'minute'))
  # print 15 min data
  print(plot_user(sel_15, main = '15 minute'))
  # print hourly data
  print(plot_user(sel_60, main = 'hourly'))
  
  mtext(paste('Pecan Street Experiment User', user$name, ', ID =', user$ID), line=0, font=2, cex=1.2, outer=TRUE)
  par(mar=c(1,4,1,1))
  plot.new()
  vars_sel = setdiff(names(data_60), 'date')
  legend("center", lty=1,cex=1,lwd=2, legend=vars_sel, col=cols[vars_sel],horiz=T)

  # save plot to file
  dev.off()
  
}

# how many users have different end uses monitored?
end_uses_tab = table(unlist(end_uses))

# form analysis dataset for entire population
avg_cons_all = do.call('rbind', avg_cons)

# _______________________________________________________
# Plot distribution of appliance categories across users

# _______________________________________________________________
# Plot distribution of consumption by end use and by time of day

