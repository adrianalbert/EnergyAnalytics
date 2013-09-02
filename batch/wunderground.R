

# adapted from blog post: http://casoilresource.lawr.ucdavis.edu/drupal/node/991
wuWeatherDay <- function(station, dateStr) {
  base_url <- 'http://www.wunderground.com/weatherstation/WXDailyHistory.asp?'
  
  day = as.POSIXlt(dateStr)
  # parse date
  m <- day$mon + 1     # mon is zero based
  d <- day$mday
  y <- day$year + 1900 # year is since 1900
  
  # compose final url
  final_url <- paste(base_url,
                     'ID=', station,
                     '&month=', m,
                     '&day=', d,
                     '&year=', y,
                     '&format=1', sep='')
  
  
  u <- url(final_url)
  urlData <- readLines(u)   # reading in as raw lines from the web server
  close(u)                  # contains <br> tags on every other line
  
  # only keep records with more than 5 rows of data
  if(length(urlData) <= 5 ) {
    print('No usable data found.')
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
    weatherData$Time <- as.POSIXct(strptime(weatherData$Time, format='%Y-%m-%d %H:%M:%S'))
    
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
#w <- wuWeatherDay(station, as.Date('2012-09-05'))

wuWeather = function(station,startDate,endDate=NULL) {
  if(is.null(endDate)){endDate = startDate}
  # get data for a range of dates
  library(plyr)
  date.range <- seq.Date(from=as.Date(startDate), to=as.Date(endDate), by='1 day')

  l <- vector(mode='list', length=length(date.range)) # pre-allocate list
  
  # loop over dates, and fetch data
  for(i in seq_along(date.range))
  {
    print(date.range[i])
    l[[i]] <- wuWeatherDay(station, date.range[i])
  }
  
  # stack elements of list into DF, filling missing columns with NA
  df <- ldply(l)
  return(df)
}

interpolateTime = function(data,newTimes,dateCol='Time') {
  weatherCols = c('TemperatureF','DewpointF','PressureIn','Humidity','SolarRadiationWatts.m.2','WindSpeedMPH')
  t = data$Time
  newData = data.frame(sapply(weatherCols,function(x) approx(t,data[,x],xout=newTimes,rule=2)$y))
  newData[[dateCol]] = newTimes
  return(newData)
}

# get single day
#wuWeather(station,'2009-01-25')
# get all weather in the range of the pecan street data
# yes, there ais a lot of hard coded stuff in here that could be more parameterized
pecanWeatherData = function(forceReload = F) {
  # default station selected by hand for decent weather history Bouldin-South Austin, Austin, TX
  station='KTXAUSTI90'
  pecanWeather = c()
  if(file.exists(PECAN_WEATHER_FILE) && ! forceReload) {
    load(PECAN_WEATHER_FILE)
  } else {
    pecanWeather = wuWeather(station,'2012-08-01','2012-10-01')
    save(list='pecanWeather',file=PECAN_WEATHER_FILE)
  }
  return(pecanWeather)
}

test=F
if(test){
  pecanWeather = pecanWeatherData()
  
  # interpolate data to arbitrary times
  t = pecanWeather$Time
  minuteTimes = seq.POSIXt(from=trunc(t[1],'day'),to=round(t[length(t)],'day'),by='min')  
  minuteData = interpolateTime(pecanWeather,minuteTimes)
  
  # average data to hourly intervals
  hourly = toHourly(pecanWeather[,c('Time',weatherCols)],'Time')

}

