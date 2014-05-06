# prepare_pecan.r
#
# Constructs an analysis dataset for pecan st validation. 
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
source('../../utils/aggregate_data.r')

# use data from 2013
DATA_PATH     = '~/energy-data/pecan_street/usage-processed/2012/'
DATA_PATH_RAW = '~/energy-data/pecan_street/usage-orig/2012/'
PLOTS_PATH    = '~/Dropbox/OccupancyStates/plots/pecan-street'
dir.create(PLOTS_PATH)

# load weather data
weather.hourly = read.csv('~/energy-data/pecan_street/weather/weather_hourly.csv')
weather.15mins = read.csv('~/energy-data/pecan_street/weather/weather_15mins.csv')
weather.hourly$date = as.POSIXct(as.character(weather.hourly$date))
weather.15mins$date = as.POSIXct(as.character(weather.15mins$date))

# load baby names
# we're going to name each user for later easiness of use
baby_names  = read.csv('~/Dropbox/OccupancyStates/data/baby-names.csv')

# __________________________________________________
# Define appliance categories

# some interesting components
appliances   = as.character(read.csv('~/energy-data/pecan_street/metadata/appliances.csv')$Appliance)
select.keep  = c('dataid', 'localminute', 'use')
select.AC    = c("air1", "air2", "air3", "airwindowunit1", "housefan1")
select.HV    = c("furnace1", "furnace2", "heater1", "housefan1")
select.light = c("lights_plugs1", "lights_plugs2", "lights_plugs3", "lights_plugs4", "lights_plugs5", "lights_plugs6",
                 "outsidelights_plugs1", "outsidelights_plugs2")
select.alwOn = c('refridgerator1', 'refridgerator2', 'winecooler1', 'aquarium1',
                 "freezer1")
select.sched = c("pool1", "pool2", 'sprinkler1', "poolpump1", "pump1")
select.total = c('use')
select.dhw   = c('waterheater1', 'waterheater2')
select.user  = c("bathroom1", "bathroom2", "bedroom1", "bedroom2", "bedroom3", "bedroom4", "bedroom5",
                 "clotheswasher1", "clotheswasher_dryg1", "diningroom1", "diningroom2", "dishwasher1",
                 "disposal1", "drye1", "dryg1", "garage1", "garage2", "icemaker1", "jacuzzi1", 
                 "kitchenapp1", "kitchenapp2", "lights_plugs1", "lights_plugs2", "lights_plugs3",
                 "lights_plugs4", "lights_plugs5", "lights_plugs6", "livingroom1", "livingroom2", 
                 "microwave1", "office1", "outsidelights_plugs1", "outsidelights_plugs2", "oven1", 
                 "poollight1",  "range1", "security1", "shed1", "utilityroom1", "venthood1")
select.solar = c('gen')
select.ev    = c('car1')

# __________________________________________________
# Load up user data

# list all data files for the selected subset of data
files    = list.files(path=DATA_PATH, full.names = T)
files_01 = list.files(path=DATA_PATH_RAW, full.names = T)
files_15 = files[grep('15mins',files)]
files_60 = files[grep('hourly', files)]

# select those users for which enough data is available


# __________________________________________________
# Plot usage by minute & by hour for each user

# plot ground truth components
plot_user = function(homeData, main = 'minute') {
  
  names(homeData) = toupper(names(homeData))
  
  # aggregate components
  AC_kwh        = add.columns(homeData,AC)
  HV_kwh        = add.columns(homeData,HV)
  total_kwh     = add.columns(homeData,total)
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