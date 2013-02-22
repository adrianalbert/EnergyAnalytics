# this file contains all the major knowledge about the structure of the PGE data
# in the database and the methods and classes required to access and manipulate
# that data.
# It requires that localConf.R is properly setup and has been run prior to invocation
# of the functions it contains.
# The two main classes are WeatherClass and ResDataClass, whcih provide a simple
# object interface to the data in the database.
require(RMySQL)

# utility function that 
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

# return the available summary information for every zipcode
# including sp_id count, climate zone, weather station, and income stats
db.getZipData = function(zip=NULL) {
  where = ''
  if(!is.null(zip)) { where = paste('where zip5=',zip) }
  query = paste(
    'SELECT zip5, COUNT(DISTINCT sp_id), cecclmzn, climate, GCOUNTY, WTHRSTN,
    median_income, median_income_quantiles
    FROM', conf.accountTable(), where, 'GROUP BY zip5')
  print(query)
  return(run.query(query,conf.meterDB()))
}

# class structure based on example from
# http://bryer.org/2012/object-oriented-programming-in-r
WeatherClass = function(zipcode){
  query = paste(
    'SELECT `date`, TemperatureF, Pressure, DewpointF, HourlyPrecip
    FROM',conf.weatherTable(),'where zip5 =',zipcode,'ORDER BY DATE')
  raw = run.query(query,conf.weatherDB())
  if(length(raw)==0) stop(paste('No data found for zipcode',zipcode))
  
   
  rawData = data.frame(
    dates = as.POSIXlt(raw[,1],tz="PST8PDT",'%Y-%m-%d %H:%M:%S'),
    tout = raw[,'TemperatureF'],
    pout = raw[,'Pressure'],
    rain = raw[,'HourlyPrecip'],
    dp   = raw[,'DewpointF']
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
  
  # returns relative humidity as decimal from 0 to 1 given temperature and dewpoint
  # using August-Roche-Magnus approximation: http://andrew.rsmas.miami.edu/bmcnoldy/humidity_conversions.pdf
  obj$rh = function(tout,dp) {
    a = 17.271
    b = 237.7
    tout = (tout - 32) * 5/9
    dp   = (dp   - 32) * 5/9
    rh = exp(a*dp/(b + dp)) / exp(a*tout/(b + tout))
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
  
  obj$resample = function(newDates,name='tout') {
    # approx returns both the newDates and the interpreted values
    # but we only need the values
    a = approx(obj$dates, obj$rawData[,name], newDates, method="linear")[[2]]
    #b = a[2]
    if (all(is.na(a))){ 
      print(paste(obj$dates[1],obj$dates[-1]))
      print(paste(newDates[1],newDates[-1]))
      stop("No weather data available") 
    }
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
  tout = weather$resample(dates,'tout')
  
  # TODO: clear out obviously bad readings
  #keepers   = which(kw > 0)
  
  obj = list (
    id = sp_id,
    dates = dates,
    kw  = kw,
    kwMat = kwMat,
    days = days,
    zipcode = zipcode,
    weather = weather,
    tout = tout,
    toutMat = matrix(tout,ncol=24,byrow=TRUE),
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
  
  obj$df = function() {
    return( data.frame(kw=obj$kw,
                       tout=obj$tout,
                       dates=obj$dates ) )
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

validateRes = function(r) {
  issues = data.frame(id=r$id)
  timeDiffs = diff(r$dates)
  units(timeDiffs) <- "hours"
  maxtd = max(timeDiffs) / 24
  span = difftime(tail(r$dates, n=1),r$dates[1],units='days')
  zerospct = sum((r$kw == 0)*1,na.rm=TRUE) / length(r$kw)
  kwmean = mean(r$kw,na.rm=TRUE)
  daylen = length(r$days)
  if( daylen < 180)     issues$days180    = daylen # less than 180 days (could be non-consecutive)
  #if( span < 270 )      issues$span270    = span # spanning less than 270 days total
  #if( maxtd > 60 )      issues$bigdiff    = maxtd # more than 2 months of missing data
  if( kwmean < 0.110 )  issues$lowmean    = kwmean # mean less than 150W is almost always empty or bad readings
  if( zerospct > 0.15 ) issues$zerospct15 = zerospct # over 15% of readings are zero
  return(issues)
}

mapColors = function(data,colorMap=NA,log=FALSE) {
  if(length(colorMap) < 2) { colorMap = heat.colors(100) }
  mn = min(data,na.rm=TRUE)
  data_mn = data - mn + 0.001
  if(log) {
    data_mn = log(data_mn + 1) # nothing below zero after we take the log
  }
  idx = ceiling(data_mn / max(data_mn,na.rm=TRUE) * length(colorMap))
  return(colorMap[idx])
}

hmap = function(data,colorMap=NA,yvals=NA,xvals=NA,log=FALSE,...) {
  n = dim(data)[1]
  m = dim(data)[2]
  # defailt values
  if(length(colorMap) < 2) { colorMap = heat.colors(100) }
  if(length(xvals)    < 2) { xvals=1:m }
  if(length(yvals)    < 2) { yvals=1:n }
  
  cols = rep(xvals,n) # duplicate column position values across all rows
  rows = rep(yvals,each=m) # duplicate y values across a whole row of data
  vals = as.vector(t(as.matrix(data))) # linearize the matrix of data
  plot(cols,rows,col=mapColors(vals,colorMap,log=log),
       ylim=c(max(yvals),min(yvals)),
       xlim=c(min(xvals),max(xvals)),
       axes=F,
       cex=5,pch=15,
       xlab='Hour of day',ylab='Date',...)
  axis(1,at=(0:23)+0.5,labels=(1:24),mgp=c(1,0,0),tcl=0.5,tick=F) # 1 = xaxis, labels in the center of each range
  axis(1,at=(0:24),labels=F,tick=T,mgp=c(1,0,0),tcl=0.5)          # ticks at boundaries
  #axis(2,pretty(c(min(yvals),max(yvals))),mgp=c(1,0,0),tcl=0.5) # yaxis
}

plot.ResDataClass = function(r,colorMap=NA,main=NA,type='summary') {
  # needs a list, called r with:
  # r$id unique identifier (just for the title)
  # r$zip zipcode for the title
  # r$days (1 date per day of data)
  # r$kw (vector of kw readings)
  # r$kwMat (matrix of kw readings with 24 columns)
  # r$toutMat (matrix of Tout readings with 24 columns)
  if(type=='summary') {
    if(is.na(main)) { main <- paste(r$id,' (',r$zip,') summary info',sep='') }
    op <- par( mfrow=c(2,2), oma=c(2,0,3,0),mar=c(2,2,2,2))# Room for the title
    #plot(r$kw,xlab='Date',ylab='kWh/h',main='Raw usage')
    hmap(r$kwMat,yvals=r$days,colorMap=colorMap,log=TRUE,main='Heatmap',mgp=c(1,0,0),tcl=0.5) # axis label on row 1, axis and ticks on 0, with ticks facing in
    end = min(length(r$kw),240)
    plot(r$dates[1:end],r$kw[1:end],type='l',xlab='Date',ylab='kWh/h',main='Raw usage zoom',mgp=c(1,0,0),tcl=0.5)
    plot(r$days,rowMeans(r$toutMat),col='grey',axes=F,ylab='',xlab='',,mgp=c(1,0,0),tcl=0.5)
    axis(4, pretty(c(0, 1.1*rowMeans(r$toutMat)),n=5), col='grey',col.axis='grey',mgp=c(1,0,0),tcl=0.5)
    mtext("T out (F)", side=4, line=1, cex=0.9, col='grey')
    par(new=T) # plot the next plot call on the same figure as previous
    plot(r$days,rowSums(r$kwMat),ylab='kWh/day',xlab='Day',main='kWh/day',mgp=c(1,0,0),tcl=0.5)
    
    plot(rowMeans(r$toutMat),rowSums(r$kwMat),main='kWh/day vs mean outside temp (F)',
         xlab='mean T (degs F)',ylab='kWh/day',mgp=c(1,0,0),tcl=0.5)
    #par(op) # Leave the last plot
    mtext(main, line=0, font=2, cex=1.2,outer=TRUE)
    par(op)
    #heatmap(as.matrix(r$kwMat),Rowv=NA,Colv=NA,labRow=NA,labCol=NA)
  } else if(type=='temp') {
    if(is.na(main)) { main <- paste(r$id,' (',r$zip,') temperature info',sep='') }
    colors = rainbow(24)
    plot(r$tout,r$kw,xlab='Tout',ylab='kW',main=main,type='p',col=mapColors(r$dates$hour,colors)) # rainbow(n, start=2/6, end=1)
    legend('right', paste("hr ", 1:24), fill=colors, ncol = 1, cex = 0.8)
  } else if(type=='hourly') {
    if(is.na(main)) { main <- paste(r$id,' (',r$zip,') hourly info',sep='') }
    grid = cbind(matrix(c(1:24),nrow=4,ncol=6,byrow=TRUE),c(25,25,25,25))
    layout(grid,widths=c(rep(2,6),1))
    op <- par(oma=c(2,3,3,0),mar=c(1,0,1,0))# Room for the title
    pallete = rainbow(12) #,start=2/6, end=1)
    colvals = r$dates$mon
    #colvals = r$dates$wday
    #colvals = r$dates$wday == 0 | r$dates$wday == 6 # 0 = Sun, 6 = Sat
    colors = mapColors(colvals,pallete)
    xlm = c(min(r$tout,na.rm=TRUE), max(r$tout,na.rm=TRUE))
    ylm = c(min(r$kw,  na.rm=TRUE), max(r$kw,  na.rm=TRUE))
    for(i in 0:23) {
      yax = 'n' # supress axis plotting
      xax = 'n'
      if(i %%  6 == 0) yax = 's' # standard y-axis for first of row
      if(i %/% 6 == 3) xax = 's' # standard x-axis for bottom row
      sub = r$dates$hour==i
      plot(subset(r$tout,sub),subset(r$kw,sub),
           col=subset(colors,sub),
           xlim=xlm,ylim=ylm,yaxt=yax,xaxt=xax )
      grid()
      text(mean(xlm),ylm[2] * 0.95,paste('hr',i),font=2, cex=1.2)
    }
    mtext(main, line=0, font=2, cex=1.2,outer=TRUE)
    #par(xpd=TRUE)
    #plot.new()
    plot(1, type = "n", axes=FALSE, xlab="", ylab="")
    legend('center', month.abb, fill=pallete, ncol = 1, cex = 1)
    #par(xpd=FALSE)
    
    par(op)
  }
  
}
