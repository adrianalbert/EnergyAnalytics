# wunderground.R
# R functions for working with Weather Underground station data in R
# Given a stantion and date range, this code will download and parse
# hourly weather data into a data.frame for use in R.
#
# Adapted from Dylan Beaudette's blog post on using R to retrieve WU data: 
# http://casoilresource.lawr.ucdavis.edu/drupal/node/991
#
# Written by Sam Borgeson (sborgeson@berkeley.edu)
# 
# You are free to re-use and modify for non-commercial purposes.
# Contact Sam with questions.
# Last modified 10/7/2013
#
# See test code example at the bottom of this file for more on usage.
require(plyr)

# i.e. wuWeatherDay('KTXAUSTI90','2012-10-01',tz='CST6CDT')
wuWeatherDay <- function(station, dateStr,tz=NULL) {
  base_url <- 'http://www.wunderground.com/weatherstation/WXDailyHistory.asp?'
  
  day = as.POSIXlt(dateStr)
  # parse date
  m <- day$mon + 1     # mon is zero based
  d <- day$mday
  y <- day$year + 1900 # year is since 1900
  
  # compose final url
  final_url <- paste(base_url,
                     'ID=',     station,
                     '&month=', m,
                     '&day=',   d,
                     '&year=',  y,
                     '&format=1', sep='')
  
  #print(final_url)
  u <- url(final_url)
  urlData <- readLines(u)   # reading in as raw lines from the web server
  close(u)                  # contains <br> tags on every other line
  # only keep records with more than 5 rows of data
  if(length(urlData) <= 5 ) {
    print('[wuWeatherDay]:No usable data found.')
    print(urlData)
  } else {
    # remove the first and last lines
    urlData <- urlData[-c(1, length(urlData))]
    
    # remove odd numbers starting from 3 --> end
    urlData <- urlData[-seq(3, length(urlData), by=2)]
    
    # extract header and cleanup
    head <- urlData[1]
    head <- make.names(strsplit(head, ',')[[1]])
    
    # convert to CSV, without header
    tC <- textConnection(paste(urlData, collapse='\n'))
    weatherData <- read.csv(tC, as.is=TRUE, row.names=NULL, header=FALSE, skip=1)
    close(tC)
    
    # remove the last column, created by trailing comma
    weatherData <- weatherData[, -ncol(weatherData)]
    
    # assign column names
    names(weatherData) <- head
    
    # convert Time column into properly encoded date time
    if(!is.null(tz)) {
      weatherData$Time <- as.POSIXct(strptime(weatherData$Time, format='%Y-%m-%d %H:%M:%S'),tz=tz)
    } else { 
      weatherData$Time <- as.POSIXct(strptime(weatherData$Time, format='%Y-%m-%d %H:%M:%S'))
    }
    
    # remove UTC and software type columns
    weatherData$DateUTC.br. <- NULL
    weatherData$SoftwareType <- NULL
    
    # sort and fix rownames
    weatherData <- weatherData[order(weatherData$Time), ]
    row.names(weatherData) <- 1:nrow(weatherData)
    
    # done
    return(weatherData)
  }
}

# be sure to load the function from above
# get a single day's worth of (hourly) data
#w <- wuWeather(station, as.Date('2012-09-05'))

wuWeather = function(station,startDate,endDate=NULL,tz=NULL) {
  if(is.null(endDate)){endDate = startDate}
  # get data for a range of dates
  library(plyr)
  date.range <- seq.Date(from=as.Date(startDate), to=as.Date(endDate), by='1 day')

  l <- vector(mode='list', length=length(date.range)) # pre-allocate list
  
  # loop over dates, and fetch data
  for(i in seq_along(date.range))
  {
    print(date.range[i])
    l[[i]] <- wuWeatherDay(station, date.range[i],tz=tz)
  }
  
  # stack elements of list into DF, filling missing columns with NA
  df <- ldply(l) # from plyr
  return(df)
}

interpolateTime = function(data,newTimes,dateCol='Time') {
  weatherCols = c(  'TemperatureF',
		    'DewpointF',
		    'PressureIn',
		    'Humidity',
		    'SolarRadiationWatts.m.2',
		    'WindSpeedMPH' )
  t = data$Time
  # run approx of the data to interpolate values for the times passed in
  newData = data.frame(sapply(weatherCols,function(x) approx(t,data[,x],xout=newTimes,rule=2)$y))
  newData[[dateCol]] = newTimes
  return(newData)
}

test=F
if(test){
  station='KTXAUSTI90' # this is ausitn Tx
  
  # get single day
  #wuWeather(station,'2009-01-25')
  
  weather = wuWeather(station,'2012-08-01','2012-9-01',tz='CST6CDT')
  # useful to save the weather data at this point:
  # save(weather,file=paste(station,'data.RData',sep='')
  
  # interpolate data to arbitrary times
  t = weather$Time
  minuteTimes = seq.POSIXt(from=trunc(t[1],'day'),to=round(t[length(t)],'day'),by='min')  
  minuteData = interpolateTime(weather,minuteTimes)
  
  # time average to hourly data
  minuteData$dayhr = strftime(minuteData[,'Time'],'%Y-%m-%d.%H',tz='CST6CDT')
  hourData = aggregate(. ~ dayhr,data=minuteData, mean, na.rm=T )
  hourData$date = strptime(hourData$dayhr,'%Y-%m-%d.%H')
  print(names(hourData))
  dt = hourData$date
  plot(dt,hourData$TemperatureF,type='l',
       ylim=c(50,185),
       main=paste('Weather data for',station,dt[1],'to',dt[length(time)]),
       xlab='Date',ylab='deg F + arbitrary other units')
  abline(v=as.POSIXct(dt)[c(1,which(dt$hour==0))],lty=3,col='#dddddd') # label day starts
  points(dt,
         105 +(hourData$WindSpeedMPH / max(hourData$WindSpeedMPH))*20,
         type='l',lty=1,col='#cccccc')
  points(dt,
         125 +(hourData$Humidity / max(hourData$Humidity))*20,
         type='l',lty=1,col='#ffcccc')
  points(dt,
         145 +(hourData$SolarRadiationWatts.m.2 / max(hourData$SolarRadiationWatts.m.2))*20,
         type='l',lty=1,col='#ccffcc')
  points(dt,
         165 +(hourData$PressureIn / max(hourData$PressureIn))*20,
         type='l',lty=1,col='#ccccff')
  points(dt,
         hourData$DewpointF,
         type='l',lty=1,col='#ffccff')
  
}

