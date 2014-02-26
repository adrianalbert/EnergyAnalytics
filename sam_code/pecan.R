# pecan.R
# R functions for working with Pecan Street sample data by loading data 
# straight from xlsx files into data frames.
#
# Written by Sam Borgeson (sborgeson@berkeley.edu)
# 
# You are free to re-use and modify for non-commercial purposes.
# Contact Sam with questions.
# Last modified 10/7/2013

# Usage:
# 1) make a directory called 'pecan' off your working directory 
# 2) copy Pecan Street '1-minute interval data' directory into it
# 3) double check you are in the working directory that contains the 'pecan'
#    directory.
# setwd('/path/to/working/directory')
#
# See the example at the bottom of this file for more on usage.
require('reshape2')
require('stringr')
require('xlsx') # R packages that uses Java's POI xls interface to read files

source(file.path(getwd(),'wunderground.R'))

PECAN_DATA_DIR = 'pecan/1-minute interval data'
PECAN_DATA_FILE = file.path(PECAN_DATA_DIR,'pecan.RData')
PECAN_WEATHER_FILE = file.path(PECAN_DATA_DIR,'weather.RData')
PECAN_COMBINED_FILE = file.path(PECAN_DATA_DIR,'pecan_plus_weather.RData')

# default station selected by hand for decent weather history Bouldin-South Austin, Austin, TX
zip5 = 78723 # Autsin Tx for weather data...
PECAN_STATION = 'KTXAUSTI90'
tz = 'CST6CDT'

pecan = function(forceReload=F) {
  if(file.exists(PECAN_DATA_FILE) && ! forceReload) {
    load(PECAN_DATA_FILE)
  } else {
    pecanData = list()
    dayData = c()
    for(homeDir in list.dirs(PECAN_DATA_DIR,full.name=F,recursive=F)) {
      parts = strsplit(homeDir,'[/\\]')[[1]]
      home = str_trim(parts[length(parts)]) # just the name, with whitespace removed.
      print(home)
      homeData = data.frame(c()) 
      for(df in list.files(homeDir,pattern='^[^~].*xlsx')) {
        print(df)
        dayData <- read.xlsx2(file.path(PECAN_DATA_DIR,home,df),1,colClasses=rep('numeric',100))
        # DAMMIT Excel!!! Get your act together on dates!
        # The two day difference is explained here:
        # http://r.789695.n4.nabble.com/Convert-number-to-Date-td1691251.html
        # Not sure where the 5 hour difference comes from...
        dates = round(as.POSIXlt(as.numeric(dayData$Date...Time)*24*3600+5*3600,tz=tz,origin='1899-12-30'),units='mins')
        dayData$date...Time <- as.POSIXct(dates) # the Time part will be stripped later
        print(dim(dayData))
        homeData <- rbind(homeData,dayData)
      }
      home <-  str_replace(home,' ','')         # eliminate spaces
      home <-  strsplit(home,'-')[[1]][1]       # Truncate at first '-' (if any)
      home <-  str_trim(home)                   # whitespace removed.
      parts <- strsplit(names(homeData),'\\.+') # now fix the column names
      names(homeData) = sapply(parts,function(x) {
         out = c()
         if(x[length(x)] == 'kVA'){
           out = paste(paste(x[-length(x)],collapse='_'),x[length(x)],sep='_' )
         } else {
           out = paste(x[-length(x)],collapse='_')
         }
         return(out) })
      class(homeData) = c('data.frame','pecan')
      pecanData[[home]] <- homeData
      rm(homeData)
      gc()
    }
  }
  save(list=c('pecanData'),file=PECAN_DATA_FILE)
  print(names(pecanData))
  return(pecanData)
}

usage = function(data,labels) {
  sub = match(labels,names(data))
  sub = sub[!is.na(sub)]
  return(apply(as.matrix(data[,sub]),1,sum))
}

toHourly = function(data,dateCol='date'){
  data$dayhr = strftime(data[,dateCol],'%Y-%m-%d.%H',tz='CST6CDT')
  hourly = aggregate(. ~ dayhr,data=data, mean, na.rm=T )
  hourly$date = strptime(hourly$dayhr,'%Y-%m-%d.%H')
  return(hourly)
}

# get all weather in the range of the pecan street data
# yes, there is a lot of hard coded stuff in here that could be more parameterized
pecanWeatherData = function(forceReload=F) {
  pecanWeather = c()
  if(file.exists(PECAN_WEATHER_FILE) && ! forceReload) {
    load(PECAN_WEATHER_FILE)
  } else {
    pecanWeather = wuWeather(PECAN_STATION,'2012-08-01','2012-10-01',tz='CST6CDT')
    save(list='pecanWeather',file=PECAN_WEATHER_FILE)
  }
  return(pecanWeather)
}


pecanPlusWeather = function(forceReload=F) {
  pecanPlusWeather = list()
  if(file.exists(PECAN_COMBINED_FILE) && ! forceReload) {
    load(PECAN_COMBINED_FILE)
  } else {
    pecanData = pecan(forceReload)
    # append matching weather data to each pecan data.frame
    source('wunderground.R')
    pecanWeat = pecanWeatherData(forceReload)
    for(homeName in names(pecanData)) {
      print(homeName)
      homeData       = pecanData[[homeName]]
      matchedWeather = interpolateTime(pecanWeat,homeData$date)
      pecanPlusWeather[[homeName]] <- merge(homeData,matchedWeather,by.x='date',by.y='Time')
    }
    save(list=c('pecanPlusWeather'),file=PECAN_COMBINED_FILE) # update the data on disk
  }
  return(pecanPlusWeather)
}


HVAC   = c('AIR1','AIR2','FURNACE1','FURNACE2')
fridge = c('REFRIGERATOR1')
dhw    = c('WATERHEATER1')
aux    = c('SPRINKLER1')
user   = c('DININGROOM1','DISHWASHER1','DRYG1','KITCHEN1','KITCHEN2','KSAC1','KSAC2','SECURITY1',
         'LIGHTING_PLUGS1','LIGHTING_PLUGS2','LIGHTING_PLUGS3','LIGHTING_PLUGS4',
         'LIVINGROOM1','MASTERBATH1','MICROWAVE1','OVEN1','RANGE1','WASHER1','WASHINGMACHINE1',
         'BATHROOM1','BEDROOM1','BEDROOM2','DISPOSAL1','DRYE1','FAMROOM1','BATH1','DRYER2','GARAGE','GARAGE1',
         'GENLIGHT1','SMALLAPPLIANCE1','SMALLAPPLIANCE2','SMALLAPPLIANCE3','COOKTOP1',
         'BACKYARD1','OFFICE1','TVROOM1','THEATER1','MASTERBED1')
total   = c('use')
solar   = c('gen')
ev      = c('CAR1')
net     = c('Grid')

runExample = F
if(runExample) {
  ppw = pecanPlusWeather()
  for(homeName in names(ppw)){
    print(homeName)
    homeData = ppw[[homeName]]
    #print(names(homeData))
    print(dim(homeData))
    
    op <- par(no.readonly = TRUE)
    cols = c('black','#ff9999','blue','#bbbbbb')
    m <- matrix(c(1,2,3),nrow=3,ncol=1,byrow=T)
    layout(mat = m,heights = c(0.4,0.4,0.2))
    par(oma=c(2,2,2,0),mar=c(2,4,2,1)) # Room for the title
    maxUsage = max(usage(homeData,total))
    plot(homeData$date,usage(homeData,total),type='l',col=cols[1],main='minute',ylab='kW',ylim=c(0,1.1*maxUsage),xaxt='n')
    points(homeData$date,usage(homeData,HVAC),type='l',col=cols[2])
    points(homeData$date,usage(homeData,user),type='l',col=cols[3])
    points(homeData$date,homeData$TemperatureF/100*maxUsage,type='l',col=cols[4])
    d = homeData$date
    dts = seq(d[1],d[length(d)],by=3600*24) # one per day from one per minute
    axis(1, dts, format(dts, "%a, %m/%d"), cex.axis=1)
    grid(nx=NA,ny=NULL)
    abline(v=dts,col="black",lty=3)
    
    hrData = toHourly(homeData)
    maxUsage = max(usage(hrData,total))
    plot(hrData$date,usage(hrData,total),type='l',col=cols[1],main='hourly',ylim=c(0,1.1*maxUsage),ylab='kW',xaxt='n')
    points(hrData$date,usage(hrData,HVAC),type='l',col=cols[2])
    points(hrData$date,usage(hrData,user),type='l',col=cols[3])
    points(hrData$date,hrData$TemperatureF/100*maxUsage,type='l',col=cols[4])
    d = as.POSIXct(hrData$date)
    dts = seq(d[1],d[length(d)],by=3600*24) # one per day from one per hour
    axis(1, dts, format(dts, "%a, %m/%d"), cex.axis=1)
    grid(nx=NA,ny=NULL)
    abline(v=dts,col="black",lty=3)
    mtext(paste('Example','kW demand'), line=0, font=2, cex=1.2, outer=TRUE)
    par(mar=c(1,4,1,1))
    plot.new()
    legend("center", lty=1,cex=1,lwd=2,
           legend=c('total','HVAC','user','temperature'), 
           col=cols,horiz=T)
    par(op)
    #break
    dev.copy2pdf(file=file.path(getwd(),PECAN_DATA_DIR,paste(homeName,'.pdf',sep='')),width=10,height=6)
  }
}

#plot(homeData$date,usage(homeData,'THEATER1'),type='l')
#points(homeData$date,usage(homeData,fridge),type='l',col='green')

