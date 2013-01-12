require(RMySQL)

run.query = function(query,db='pge_res') {
  #print(query)
  data <- c()
  tryCatch({
    con  <- conf.dbCon(db)
    res  <- dbGetQuery(con, query)
    if(length(res)>0) data  <- res
    
  },
           error = function(e) {print(e)},
           finally = {
             # close the results set if necessary
             resultSet <- dbListResults(con)
             if(length(resultSet)>0) { 
               dbClearResult(resultSet[[1]])
               rm(resultSet)
             }
             dbDisconnect(con)
             rm(con)
           } )
  return(data)
}

db.getZips = function() {
  query    <- paste("select distinct zip5 from",conf.weatherTable(),'order by zip5')
  return(run.query(query,conf.weatherDB())[[1]])
}

db.getSPs = function(zip=NA) {
  if(is.na(zip)) { query <- paste("select distinct sp_id from",conf.meterTable()) }
  else {           query  <- paste("select distinct sp_id from",conf.meterTable(zip),"where zip5=",zip) }
  return(run.query(query,conf.meterDB())[[1]])
}

db.getAllData = function(zip=NULL) {
  query = paste(
         'SELECT 
         sp_id, zip5, DATE,
         hkw1, hkw2, hkw3, hkw4, hkw5, hkw6, hkw7, hkw8, hkw9, hkw10,hkw11,hkw12,
         hkw13,hkw14,hkw15,hkw16,hkw17,hkw18,hkw19,hkw20,hkw21,hkw22,hkw23,hkw24 
         FROM',conf.meterTable(zip),'ORDER BY sp_id, DATE')
  return(run.query(query,conf.meterDB()))
}

# class structure based on example from
# http://bryer.org/2012/object-oriented-programming-in-r
RowWeatherClass = function(zipcode){
  query = paste(
    'SELECT `date`,
    htempF1, htempF2, htempF3, htempF4, htempF5, htempF6, htempF7, htempF8, htempF9, htempF10,htempF11,htempF12,
    htempF13,htempF14,htempF15,htempF16,htempF17,htempF18,htempF19,htempF20,htempF21,htempF22,htempF23,htempF24
    FROM',conf.temperatureTable(),'where zip5 =',zipcode,'ORDER BY DATE')
  raw = run.query(query,conf.weatherDB())
  
  days = as.POSIXct(raw[,1],tz="PST8PDT", '%Y-%m-%d')
  # create a row of hourly values for each day
  daySteps = 24
  dtDay = daySteps/24 * 60 * 60 # in seconds
  # sapply returns an array of numeric epoch seconds (for origin '1970-01-01')
  dateMat = sapply(days,FUN=function(x) x + (0:(daySteps-1) * dtDay))
  
  # flatten into a vector and re-convert into date objects
  dates = as.POSIXlt(as.vector(dateMat),origin='1970-01-01')
  
  toutMat = raw[,2:25]
  # reshape the tout readings into a vector matching the dates
  tout = as.vector(t(toutMat))
  
  # TODO: do we need to do anything about the NA values?
  
  obj = list (
    days    = days,
    dateMat = dateMat,
    dates   = dates,
    tout    = tout,
    toutMat = toutMat,
    get     = function(x) obj[[x]],
    # Not sure why <<- is used here
    # <<- searches parent environments before assignment
    # http://stat.ethz.ch/R-manual/R-patched/library/base/html/assignOps.html
    set     = function(x, value) obj[[x]] <<- value,
    props   = list()
  )
  
  obj$add = function(name, value) {
    obj[[name]] = value
  }
  
  # note how list manipulation requires the use of assign
  # not sure why values can't be set in place, but it 
  # appears to have to do with variable scoping
  obj$addProp = function(name, value) {
    p <- obj$props
    p[[name]] <- value
    assign('props', p, envir=obj)
  }
  
  obj$resample = function(newDates) {
    # approx returns both the newDates and the interpreted values
    # but we only need the values
    a = approx(obj$dates, obj$tout, newDates, method="linear")[[2]]
    #b = a[2]
    if (all(is.na(a))){
      stop("No weather data available")
    }
    return(a)
  }
  
  #obj <- list2env(obj)
  class(obj) = "RowWeatherClass"
  return(obj)
}

# class structure based on example from
# http://bryer.org/2012/object-oriented-programming-in-r
WeatherClass = function(zipcode){
  query = paste(
    'SELECT `date`, TemperatureF, Pressure, HourlyPrecip
    FROM',conf.weatherTable(),'where zip5 =',zipcode,'ORDER BY DATE')
  raw = run.query(query,conf.weatherDB())
  if(length(raw)==0) stop(paste('No data found for zipcode',zipcode))
  
   
  rawData = data.frame(
    dates = as.POSIXlt(raw[,1],tz="PST8PDT",'%Y-%m-%d %H:%M:%S'),
    tout = raw[,'TemperatureF'],
    pout = raw[,'Pressure'],
    rain = raw[,'HourlyPrecip']
  )
  
  days = unique(as.Date(rawData$dates))
  
  # FYI, spring forward causes NA dates to find these:
  # which(is.na(dates))
  
  # TODO: do we need to do anything about the NA values?
  
  obj = list (
    days    = days,
    dates   = rawData$dates,
    tout    = rawData$tout,
    rawData = rawData,
    get     = function(x) obj[[x]],
    # Not sure why <<- is used here
    # <<- searches parent environments before assignment
    # http://stat.ethz.ch/R-manual/R-patched/library/base/html/assignOps.html
    set     = function(x, value) obj[[x]] <<- value,
    props   = list()
  )
  
  obj$add = function(name, value) {
    obj[[name]] = value
  }
  
  # note how list manipulation requires the use of assign
  # not sure why values can't be set in place, but it 
  # appears to have to do with variable scoping
  obj$addProp = function(name, value) {
    p <- obj$props
    p[[name]] <- value
    assign('props', p, envir=obj)
  }
  
  obj$resample = function(newDates,name='tout') {
    # approx returns both the newDates and the interpreted values
    # but we only need the values
    a = approx(obj$dates, obj$rawData[,name], newDates, method="linear")[[2]]
    #b = a[2]
    if (all(is.na(a))){ stop("No weather data available") }
    return(a)
  }
  #obj <- list2env(obj)
  class(obj) = "WeatherClass"
  return(obj)
}

ResDataClass = function(sp_id,zip=NULL,weather=NULL,data=NULL,db='pge_res'){
  if(is.null(data) || length(data) == 0) {
    query = paste(
      'SELECT 
        sp_id,zip5,DATE,
      hkw1, hkw2, hkw3, hkw4, hkw5, hkw6, hkw7, hkw8, hkw9, hkw10,hkw11,hkw12,
      hkw13,hkw14,hkw15,hkw16,hkw17,hkw18,hkw19,hkw20,hkw21,hkw22,hkw23,hkw24 
      FROM',conf.meterTable(zip),'WHERE sp_id =',sp_id,'ORDER BY DATE')
    data = run.query(query,conf.meterDB())
  }
  if(length(data)==0) stop(paste('No data found for sp_id',sp_id))
  
  zipcode = data[1,'zip5']
  kwMat = data[,4:27]
  # reshape the kW readings into a vector matching the dates
  kw    = as.vector(t(kwMat))
  
  days = as.POSIXct(data[,'DATE'],tz="PST8PDT", '%Y-%m-%d')
  # create a row of hourly values for each day
  daySteps = 24
  dtDay = daySteps/24 * 60 * 60 # in seconds
  # sapply returns an array of numeric epoch seconds (for origin '1970-01-01')
  dateMat = sapply(days,FUN=function(x) x + (0:(daySteps-1) * dtDay))
  # flatten into a vector and re-convert into date objects
  dates = as.POSIXlt(as.vector(dateMat),origin='1970-01-01')
  
  if (is.null(weather)) weather = WeatherClass(zipcode) # todo: pass in the dates to interpolate ,dates)
  
  # TODO: clear out obviously bad readings
  #keepers   = which(kw > 0)
  
  obj = list (
    dates = dates,
    kw  = kw,
    kwMat = kwMat,
    days = days,
    zipcode = zipcode,
    weather = weather,
    tout = weather$resample(dates,'tout'),
    get = function(x) obj[[x]],
    # Not sure why <<- is used here
    # <<- searches parent environments before assignment
    # http://stat.ethz.ch/R-manual/R-patched/library/base/html/assignOps.html
    set = function(x, value) obj[[x]] <<- value,
    props = list()
  )
  
  obj$w = function(name='tout') {
    return( obj$weather$resample(dates,name) )
  }
  
  obj$norm = function(data) {
    # divide data by the 97th %ile
    return (data/quantile(data,0.97,na.rm=TRUE))
  }
  
  obj$add = function(name, value) {
    obj[[name]] = value
  }
  
  # note how list manipulation requires the use of assign
  # not sure why values can't be set in place, but it 
  # appears to have to do with variable scoping
  obj$addProp = function(name, value) {
    p <- obj$props
    p[[name]] <- value
    assign('props', p, envir=obj)
  }
  
  obj$matchDates = function(newDates) {
    a = approx(obj$dates, obj$tout, newDates, method="linear" )[[2]]
    return(a)
  }
  
  #obj <- list2env(obj)
  class(obj) = "ResDataClass"
  return(obj)
}

