# make a directory called pecan off your working directory 
# and copy the 1-minute... directory of pecan data into it
# double check you are in the working directory that contains the pecan data

PECAN_DATA_DIR = 'pecan/1-minute interval data'
PECAN_DATA_FILE = file.path(PECAN_DATA_DIR,'pecan.RData')
PECAN_WEATHER_FILE = file.path(PECAN_DATA_DIR,'weather.RData')
PECAN_COMBINED_FILE = file.path(PECAN_DATA_DIR,'pecan_plus_weather.RData')

require('reshape2')
require('stringr')

zip5 = 78723 # TODO: get weather data...
station = 'KTXAUSTI204'
tz = 'CST6CDT'

pecan = function(forceReload=F) {
  if(file.exists(PECAN_DATA_FILE) && ! forceReload) {
    load(PECAN_DATA_FILE)
  } else {
    require('xlsx') # R packages that uses Java's POI xls interface to read files
    pecanData = list()
    dayData = c()
    for(homeDir in list.dirs(PECAN_DATA_DIR,full.name=F,recursive=F)) {
      parts = strsplit(homeDir,'[/\\]')[[1]]
      home = parts[length(parts)]
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
      home <- str_replace(home,' ','')   # eliminate spaces
      home <- strsplit(home,'-')[[1]][1] # Truncate at first '-' (if any)
      parts = strsplit(names(homeData),'\\.+')
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

pecanCombinedFile = c()
if(file.exists(PECAN_COMBINED_FILE) && ! forceReload) {
  load(PECAN_COMBINED_FILE)
} else {
  pecanData = pecan()
  # append matching weather data to each pecan data.frame
  source('wunderground.R')
  pecanWeat = pecanWeatherData()
  pecanCombinedData = list()
  for(homeName in names(pecanData)) {
    print(homeName)
    homeData       = pecanData[[homeName]]
    matchedWeather = interpolateTime(pecanWeat,homeData$date)
    pecanCombinedData[[homeName]] <- merge(homeData,matchedWeather,by.x='date',by.y='Time')
  }
  save(list=c('pecanCombinedData'),file=PECAN_COMBINED_FILE) # update the data on disk
}


HVAC = c('AIR1','AIR2','FURNACE1','FURNACE2')
fridge = c('REFRIGERATOR1')
dhw   = c('WATERHEATER1')
aux = c('SPRINKLER1')
user = c('DININGROOM1','DISHWASHER1','DRYG1','KITCHEN1','KITCHEN2','KSAC1','KSAC2','SECURITY1',
         'LIGHTING_PLUGS1','LIGHTING_PLUGS2','LIGHTING_PLUGS3','LIGHTING_PLUGS4',
         'LIVINGROOM1','MASTERBATH1','MICROWAVE1','OVEN1','RANGE1','WASHER1','WASHINGMACHINE1',
         'BATHROOM1','BEDROOM1','BEDROOM2','DISPOSAL1','DRYE1','FAMROOM1','BATH1','DRYER2','GARAGE','GARAGE1',
         'GENLIGHT1','SMALLAPPLIANCE1','SMALLAPPLIANCE2','SMALLAPPLIANCE3','COOKTOP1',
         'BACKYARD1','OFFICE1','TVROOM1','THEATER1','MASTERBED1')
total = c('use')
solar = c('gen')
ev    = c('CAR1')
net   = c('Grid')

for(homeName in names(pecanCombinedData)){
  print(homeName)
  homeData = pecanCombinedData[[homeName]]
  #print(names(homeData))
  print(dim(homeData))
  
  op <- par(no.readonly = TRUE)
  par( mfrow=c(2,1), oma=c(2,2,3,0),mar=c(2,2,2,2)) # Room for the title
  maxUsage = max(usage(homeData,total))
  plot(homeData$date,usage(homeData,total),type='l',col='gray',main='minute',ylim=c(0,maxUsage))
  points(homeData$date,usage(homeData,HVAC),type='l',col='#ff9999')
  points(homeData$date,usage(homeData,user),type='l',col='blue')
  points(homeData$date,homeData$TemperatureF/*120/maxUsage,type='l',col='#99ff99')
  
  hrData = toHourly(homeData)
  maxUsage = max(usage(hrData,total))
  plot(hrData$date,usage(hrData,total),type='l',col='gray',main='hourly',ylim=c(0,maxUsage))
  points(hrData$date,usage(hrData,HVAC),type='l',col='#ff9999')
  points(hrData$date,usage(hrData,user),type='l',col='blue')
  points(hrData$date,hrData$TemperatureF*120/maxUsage,type='l',col='#99ff99')
  par(op)
}


#points(homeData$date,usage(homeData,fridge),type='l',col='green')
#plot(homeData$date,usage(homeData,'THEATER1'),type='l')

