# run_profiler_pecan.r
#
# Applies HMM decoding on Pecan Street data. 
#
# Adrian Albert
# Last modified: May 2014.
# ---------------------------------------------------------

rm(list = ls())
options(error = recover)
library('segmented')

# __________________________________________________
# Initializations...

setwd('~/EnergyAnalytics/thermal_profiles/profiler/')
source('stateProcessorWrapper.r')
source('stateVisualizerWrapper.r')

DATA_PATH  = '~/energy-data/pecan_street/usage-processed/2013/'
PLOTS_PATH = '~/Dropbox/OccupancyStates/plots/pecan-street'
dir.create(PLOTS_PATH)

# load weather data
weather.hourly = read.csv('~/energy-data/pecan_street/weather/weather_hourly.csv')
weather.15mins = read.csv('~/energy-data/pecan_street/weather/weather_15mins.csv')
weather.hourly$date = as.POSIXct(as.character(weather.hourly$date))
weather.15mins$date = as.POSIXct(as.character(weather.15mins$date))

# load baby names
baby_names  = read.csv('~/Dropbox/OccupancyStates/data/baby-names.csv')

# __________________________________________________
# Select appliances of interest

HVAC   = c('AIR1','AIR2','FURNACE1','FURNACE2')
HV     = c('FURNACE1','FURNACE2')
AC     = c('AIR1','AIR2')
lights = c('LIGHTS_PLUGS1','LIGHTS_PLUGS2','LIGHTS_PLUGS3','LIGHTS_PLUGS4')
total   = c('USE')

# __________________________________________________
# Load up user data

# list all data files for 2013
files    = list.files(path=DATA_PATH, full.names = T)
files_01 = files[grep('minute',files)]
files_15 = files[grep('15mins',files)]
files_60 = files[grep('hourly', files)]

# __________________________________________________
# Plot usage by minute & by hour for each user

# aggregate data from columns defined by labels
usage = function(data,labels) {
  sub = match(labels,names(data))
  sub = sub[!is.na(sub)]
  return(apply(as.matrix(data[,sub]),1,sum))
}

# plot ground truth components
plot_user = function(homeData, main = 'minute') {
  
  names(homeData) = toupper(names(homeData))
  
  # aggregate components
  AC_kwh        = usage(homeData,AC)
  HV_kwh        = usage(homeData,HV)
  total_kwh     = usage(homeData,total)
  occupancy_kwh = total_kwh - AC_kwh - HV_kwh
  
  maxUsage = max(total_kwh)
  plot(homeData$DATE,total_kwh,type='l',col=cols[1],main=main,ylab='kW',ylim=c(0,1.1*maxUsage),xaxt='n')
  points(homeData$DATE,AC_kwh,type='l',col=cols[3])
  points(homeData$DATE,HV_kwh,type='l',col=cols[2])
  points(homeData$DATE,occupancy_kwh,type='l',col=cols[4])
  if ('TEMPERATUREF' %in% names(homeData))
    points(homeData$DATE,homeData$TEMPERATUREF/100*maxUsage,type='l',col=cols[5])
  d = homeData$DATE
  dts = seq(d[1],d[length(d)],by=3600*24) # one per day from one per minute
  axis(1, dts, format(dts, "%a, %m/%d"), cex.axis=1)
  grid(nx=NA,ny=NULL)
  abline(v=dts,col="black",lty=3)
  
}

all_data = list()
for(i in 1:length(files_01)){

  # some gymnastics to get user ID
  userID = strsplit(rev(strsplit(files_01[i], '/')[[1]])[1], '_')[[1]][1]
  
  # read in data at different resolutions...
  cat(paste(userID, '...\n'))
  data_01 = read.csv(files_01[i]); data_01$localminute = as.POSIXct(data_01$localminute); names(data_01)[1] = 'date';
  data_15 = read.csv(files_15[i]); data_15$date = as.POSIXct(data_15$date);
  data_60 = read.csv(files_60[i]); data_60$date = as.POSIXct(data_60$date);
  
  # add in weather (temperature) data
  data_15 = merge(data_15, subset(weather.15mins, select = c('date', 'TemperatureF')), by = 'date')
  data_60 = merge(data_60, subset(weather.hourly, select = c('date', 'TemperatureF')), by = 'date')  
  
  # skip users that don't have enough data
  if (nrow(data_60) < 2 * 30 * 24) {
    cat('Not enough data!\n')
    next
  }
  
  dir.create(file.path(PLOTS_PATH, userID))
  pdf(file=paste(paste(PLOTS_PATH, userID, sep = '/'), paste(userID,'.pdf',sep=''), sep='/'),width=10,height=6)
  
  # define plot parameters
  op <- par(no.readonly = TRUE)
  m <- matrix(c(1,2,3,4),nrow=4,ncol=1,byrow=T)
  layout(mat = m,heights = c(0.3,0.3,0.3,0.1))
  par(oma=c(2,2,2,0),mar=c(2,4,2,1)) # Room for the title
  
  cols = c('black','#FE2E2E','#0040FF', '#088A08', '#424242')
  # print minute-by-minute data
  print(plot_user(data_01[1:(7*24*60),], main = 'minute'))
  # print 15 min data
  print(plot_user(data_15[1:(7*24*4),], main = '15 minute'))
  # print hourly data
  print(plot_user(data_60[1:(7*24),], main = 'hourly'))
  
  mtext(paste('Pecan Street Experiment User', userID), line=0, font=2, cex=1.2, outer=TRUE)
  par(mar=c(1,4,1,1))
  plot.new()
  legend("center", lty=1,cex=1,lwd=2,
         legend=c('total','HV', 'AC', 'occupant', 'temperature (F)'), 
         col=cols,horiz=T)

  # save plot to file
  dev.off()
  
  # store data for later processsing
  all_data[[userID]] = list(min_15 = data_15, min_60 = data_60)
}

# save data to RData file
save(list = c('all_data'), file = '~/energy-data/pecan_street/')