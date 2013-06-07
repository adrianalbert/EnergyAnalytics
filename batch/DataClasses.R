# this file contains all the major knowledge about the structure of the PGE data
# in the database and the methods and classes required to access and manipulate
# that data.
# It requires that localConf.R is properly setup and has been run prior to invocation
# of the functions it contains.
# The two main classes are WeatherClass and ResDataClass, whcih provide a simple
# object interface to the data in the database.
require(RMySQL)
library(RColorBrewer)

source(file.path(getwd(),'solaRUtil.R'))

# utility function that returns a list of all the zipcodes in the data set
db.getZips = function(useCache=F,forceRefresh=F) {
  query    <- paste("select distinct zip5 from",conf.weatherTable(),'order by zip5')
  cacheFile = NULL
  if(useCache) { cacheFile='zipList.RData' }
  return(run.query(query,conf.weatherDB(),cacheFile=cacheFile,forceRefresh=forceRefresh)[[1]])
}

db.getSPs = function(zip=NA,useCache=F,forceRefresh=F) {
  cacheFile = NULL
  if(is.na(zip)) { query <- paste("select distinct sp_id from",conf.meterTable()) }
  else {
    if(useCache) {
      cacheFile = paste('spids_',zip,'.RData',sep='') # only cache sp list for individual zips
    }
    query  <- paste("select distinct sp_id from",conf.meterTable(zip),"where zip5=",zip) 
  }
  return(run.query(query,conf.meterDB(),cacheFile=cacheFile,forceRefresh=forceRefresh)[[1]])
}

db.getZipCounts = function() { 
  query = 'SELECT zip5, COUNT(DISTINCT sp_id) as count FROM pge_res_final GROUP BY zip5'
  return(run.query(query,conf.meterDB()))
}

db.getAllData = function(zip=NULL,useCache=F,forceRefresh=F) {
  cacheFile = NULL
  if(useCache) { cacheFile=paste('meterData_',zip,'.RData',sep='') }
  query = paste(
    'SELECT 
         sp_id, zip5, DATE,
         hkw1, hkw2, hkw3, hkw4, hkw5, hkw6, hkw7, hkw8, hkw9, hkw10,hkw11,hkw12,
         hkw13,hkw14,hkw15,hkw16,hkw17,hkw18,hkw19,hkw20,hkw21,hkw22,hkw23,hkw24 
         FROM',conf.meterTable(zip),'ORDER BY sp_id, DATE')
  zipData = run.query(query,conf.meterDB(),cacheFile=cacheFile,forceRefresh=forceRefresh)
  return(zipData)
}

# return the available summary information for every zipcode
# including sp_id count, climate zone, weather station, and income stats
ZIP_DATA = NULL
db.getZipData = function(zip=NULL,useCache=F) {
  if(is.null(ZIP_DATA)) {
    query = paste(
      'SELECT zip5, COUNT(DISTINCT sp_id), cecclmzn, climate, GCOUNTY, WTHRSTN,
      median_income, median_income_quantiles
      FROM', conf.accountTable(), 'GROUP BY zip5')
    # has to go into the global env to persist as a 'cache'
    cacheFile=NULL
    if(useCache) { cacheFile='zipData.RData' }
    assign('ZIP_DATA', run.query(query,conf.resDB(),cacheFile=cacheFile), envir = .GlobalEnv) 
  }
  #else { print('Using zip data cache') }
  if(is.null(zip)) { out = ZIP_DATA                       }
  else             { out = ZIP_DATA[ZIP_DATA$zip5 == zip,]}
  return(out)
}

sumDay = function(df) {
  sums = 24 * colMeans(df,na.rm=T) # calculate the sums for all columns
}

meanDay = function(df) {
  means = colMeans(df,na.rm=T) # calculate the means for all columns
  means['rain'] = 24*means['rain'] # rain should be the total, not mean
  return(means)
}

minDay = function(df) {
  # see: http://stackoverflow.com/questions/13676878/fastest-way-to-get-min-from-every-column-in-a-matrix
  #mins = do.call(pmin,c(lapply(1:nrow(df), function(i)df[i,]),na.rm=T)) # calculate the minsacross all rows
  mins = suppressWarnings(apply(df, 2, min, na.rm=T))
  mins[mins == Inf] = NA
  return(mins)
}

maxDay = function(df) {
  #maxs = do.call(pmax,c(lapply(1:nrow(df), function(i)df[i,]),na.rm=T)) # calculate the maxs across all rows
  maxs = suppressWarnings(apply(df, 2, max, na.rm=T))
  maxs[maxs == -Inf] = NA
  return(maxs)
}

dailySummary = function(rawData,fn) {
  days = factor(as.Date(rawData$dates,tz="PST8PDT"))
  dayMeans = do.call(rbind,as.list(by(rawData[-1],days,fn))) # -1 to get rid of the dates
  dayMeans[is.nan(dayMeans)] = NA
  dayMeans = data.frame(dayMeans)
  dayMeans$day = as.Date(rownames(dayMeans))
  rownames(dayMeans) <- c()
  return(dayMeans)
}

dailyMeans = function(rawData) {
  return(dailySummary(rawData,meanDay))
}

dailyMins = function(rawData) {
  return(dailySummary(rawData,minDay))
}

dailyMaxs = function(rawData) {
  return(dailySummary(rawData,maxDay))
}

dailySums = function(rawData) {
  return(dailySummary(rawData,sumDay))
}

# class structure based on example from
# http://bryer.org/2012/object-oriented-programming-in-r
WeatherClass = function(zipcode,doMeans=T,useCache=F,doSG=F){
  query = paste(
    'SELECT `date`, TemperatureF, Pressure, DewpointF, HourlyPrecip
    FROM',conf.weatherTable(),'where zip5 =',zipcode,'ORDER BY DATE')
  cacheFile=NULL
  if(useCache) {
    cacheFile=paste('weather_',zipcode,'.RData',sep='')
  }
  raw = run.query(query,conf.weatherDB(),cacheFile=cacheFile)
  if(length(raw)==0) stop(paste('No data found for zipcode',zipcode))
  
   
  rawData = data.frame(
    dates = as.POSIXlt(raw[,1],tz="PST8PDT",'%Y-%m-%d %H:%M:%S'),
    tout = raw[,'TemperatureF'],
    pout = raw[,'Pressure'],
    rain = raw[,'HourlyPrecip'],
    dp   = raw[,'DewpointF']
  )
  
  days = unique(as.Date(rawData$dates))
  
  sg = list()
  if(doSG) {
    sg = solarGeom(rawData$dates,zip=zipcode)
  }
  
  
  # FYI, spring forward causes NA dates to find these:
  # which(is.na(dates))
  
  # TODO: do we need to do anything about the NA values?
  
  dayMeans    = c()
  dayMins     = c()
  dayMaxs     = c()
  dayLengths  = c()
  if(doMeans) {
    dayMeans  = dailyMeans(rawData)
    dayMins   = dailyMins(rawData)
    dayMaxs   = dailyMaxs(rawData)
    if(doSG) {
      dayLengths = dailySums(sg[,c('dates','daylight')])
    }
  }
  obj = list (
    zip      = zipcode,
    days     = days,
    dates    = rawData$dates,
    tout     = rawData$tout,
    sg       = sg,
    daylight = sg$daylight,
    rawData  = rawData,
    dayMeans = dayMeans,
    dayMins  = dayMins,
    dayMaxs  = dayMaxs,
    dayLengths = dayLengths,
    get      = function(x) obj[[x]],
    # Not sure why <<- is used here
    # <<- searches parent environments before assignment
    # http://stat.ethz.ch/R-manual/R-patched/library/base/html/assignOps.html
    set      = function(x, value) obj[[x]] <<- value,
    props    = list()
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
    if(all(is.na(obj$rawData[,name]))) {
      return( rep(NA,length(newDates)))
    }
    a = approx(obj$dates, obj$rawData[,name], newDates, method="linear")[[2]]
    #b = a[2]
    if (all(is.na(a)) & name == 'tout'){ 
      print(paste(obj$dates[1],obj$dates[length(obj$dates)]))
      print(paste(newDates[1], newDates[length(newDates)]))
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
  
  days = as.POSIXct(data[,'DATE'],tz="PST8PDT", '%Y-%m-%d')
  # some days are duplicated in the DB. Remove all but the first one.
  dup = which(diff(days) == 0)
  if(any(dup > 0)) {
    data = data[-dup,] # delete the row
    days = as.POSIXct(data[,'DATE'],tz="PST8PDT", '%Y-%m-%d')
  }
  zipcode = data[1,'zip5']
  kwMat = data[,4:27]
  # reshape the kW readings into a vector matching the dates
  kw    = as.vector(t(kwMat))
  
  
  # create a row of hourly values for each day
  daySteps = 24
  dtDay = daySteps/24 * 60 * 60 # in seconds
  # sapply returns an array of numeric epoch seconds (for origin '1970-01-01')
  dateMat = sapply(days,FUN=function(x) x + (0:(daySteps-1) * dtDay))
  # flatten into a vector and re-convert into date objects
  dates = as.POSIXlt(as.vector(dateMat),origin='1970-01-01')
  
  if (is.null(weather)) weather = WeatherClass(zipcode,doSG=T)
  tout = weather$resample(dates,'tout')
  pout = weather$resample(dates,'pout')
  rain = weather$resample(dates,'rain')
  dp   = weather$resample(dates,'dp')
  rh   = weather$rh(tout,dp)
  
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
    pout = pout,
    rain = rain,
    dp   = dp,
    rh   = rh,
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

  obj$daily = function(var='kw',fun=mean) {
    numDays = dim(obj$kwMat)[1]
    daily = sapply(1:numDays, function(x) fun(obj[[var]][(24*(x-1)+1):(24*x)],na.rm=T))
    return(daily)
  }
  
  #obj <- list2env(obj)
  class(obj) = "ResDataClass"
  return(obj)
}

validateRes = function(r) {
  issues = data.frame(id=r$id)
  #timeDiffs = diff(r$dates)
  #units(timeDiffs) <- "hours"
  #maxtd = max(timeDiffs) / 24
  #span = difftime(tail(r$dates, n=1),r$dates[1],units='days')
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
  if(is.numeric(data)) {
    mn = min(data,na.rm=TRUE)
    data_mn = data - mn + 0.001
    if(log) {
      data_mn = log(data_mn + 1) # nothing below zero after we take the log
    }
    idx = ceiling(data_mn / max(data_mn,na.rm=TRUE) * length(colorMap))
  }
  if(is.character(data)) { data = factor(data) }
  if(is.factor(data)) { idx = data }
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

# quickly find the estimates for a simple change point model with one cp
quickEst = function(cp,const,b,tlim=c(0,100)) {
  tlim=c(floor(tlim[1]),ceiling(tlim[2]))
  tRng = tlim[1]:tlim[2]
  print(cp)
  print(b)
  pieces = regressor.piecewise(tRng,cp)
  X = cbind(1,pieces)
  y = X %*% c(const,b)
  #print(y)
  #y = c(c(const + tlim[1]:cp * lower),c(const + cp * lower + ((cp+1):tlim[2] -cp) * upper))
  return(cbind(tRng,y))
}

save.png.plot = function(r,path,issues=NULL) {
  tryCatch( {
    png(path)
    colorMap = rev(colorRampPalette(brewer.pal(11,"RdBu"))(100))
    issueTxt = ''
    if(length(issues)>1) { issueTxt = paste(colnames(issues)[-1],collapse=', ') }
    plot( r, colorMap=colorMap, main=paste(r$zip, r$id),issueTxt=issueTxt,estimates=toutChangePointFast(df=rDFA(r)) )
  }, 
  error = function(e) { print(e) },
  finally = { dev.off() } )
}

plot.ResDataClass = function(r,colorMap=NA,main=NULL,issueTxt='',type='summary',colorBy='hours',estimates=NULL,extraPoints=NULL) {
  # needs a list, called r with:
  # r$id unique identifier (just for the title)
  # r$zip zipcode for the title
  # r$days (1 date per day of data)
  # r$kw (vector of kw readings)
  # r$kwMat (matrix of kw readings with 24 columns)
  # r$toutMat (matrix of Tout readings with 24 columns)
  if(type=='summary') {
    if(is.null(main)) { main <- paste(r$id,' (',r$zip,') summary info',sep='') }
    if(length(colorMap) < 2) { colorMap = rev(colorRampPalette(brewer.pal(11,"RdBu"))(100)) } #colorMap = heat.colors(100)
    op <- par(no.readonly = TRUE)
    par( mfrow=c(2,2), oma=c(2,2,3,0),mar=c(2,2,2,2))# Room for the title
    #plot(r$kw,xlab='Date',ylab='kWh/h',main='Raw usage')
    
    # image is messed up. we need to reverse the rows, convert to a matrix and transpose the data
    # to get the right orientation!
    image(t(as.matrix(r$kwMat[rev(1:dim(r$kwMat)[1]),])),col=colorMap,axes=F,main='kW')
    axis(1, at = seq(0, 1, by = 1/6),labels=0:6 * 4,mgp=c(1,0,0),tcl=0.5)
    if(length(r$days) > 16) {
      axis(2, at = seq(1,0, by = -1/15),labels=format(r$days[seq(1/16, 1, by = 1/16) * length(r$days)],'%m/%d/%y'),las=1,mgp=c(1,0,0),tcl=0.5)
    } else {
      axis(2, at = seq(1,0, by = -1/(length(r$days)-1)),labels=format(r$days,'%m/%d/%y'),las=1,mgp=c(1,0,0),tcl=0.5)
    }
    
    #hmap(,yvals=r$days,colorMap=colorMap,log=TRUE,main='Heatmap',mgp=c(1,0,0),tcl=0.5) # axis label on row 1, axis and ticks on 0, with ticks facing in
    end = min(length(r$kw),240)
    plot(r$dates[1:end],r$kw[1:end],type='l',xlab='Date',ylab='kWh/h',main='Raw usage zoom',mgp=c(1,0,0),tcl=0.5)
    toutMeans = rowMeans(r$toutMat)
    kWh       = rowSums(r$kwMat)
    plot(r$days,toutMeans,col='grey',axes=F,ylab='',xlab='',,mgp=c(1,0,0),tcl=0.5)
    axis(4, pretty(c(0, 1.1*toutMeans),n=5), col='grey',col.axis='grey',mgp=c(1,0,0),tcl=0.5)
    mtext("T out (F)", side=4, line=1, cex=0.9, col='grey')
    par(new=T) # plot the next plot call on the same figure as previous
    plot(r$days,kWh,ylab='kWh/day',xlab='Day',main='kWh/day',mgp=c(1,0,0),tcl=0.5)
    
    xlm = c(min(toutMeans,na.rm=TRUE), max(toutMeans,na.rm=TRUE))
    ylm = c(min(kWh,  na.rm=TRUE), max(kWh,  na.rm=TRUE))
    plot(toutMeans,rowSums(r$kwMat),main='kWh/day vs mean outside temp (F)',
         xlab='mean T (degs F)',ylab='kWh/day',
         xlim=xlm,ylim=ylm, # use a well defined set of ranges so they can be matched by any estimates below
         mgp=c(1,0,0),tcl=0.5)
    txtCol = 'black'
    if(nchar(issueTxt) > 0) { 
      main = paste(main,issueTxt)             
      txtCol = 'red'
    }
    mtext(paste(main,issueTxt), line=0, font=2, col=txtCol, cex=1.2, outer=TRUE)
    if(length(estimates) > 0) {
      #print(estimates)
      par(new=T)
      fit = estimates
      color = 'blue'
      if (fit['AIC_0'] < fit['AIC_cp']) color = 'gray'
      if (fit['nullModelTest'] > 0.1)   color = 'red'
      cp = fit[grep('^cp',names(fit))]
      middle = fit[grep('^middle',names(fit))]
      qe = quickEst(cp,fit['(Intercept)'],c(fit['lower'],middle,fit['upper']),
                  c(min(toutMeans,na.rm=T),max(toutMeans,na.rm=T)))
      #print(qe[,1])
      plot(qe[,1],qe[,2],
           type='l',
           xlim=xlm,ylim=ylm,
           axes=F,lwd=2,col=color )
    }
    par(new=F)
    par(op)
    #heatmap(as.matrix(r$kwMat),Rowv=NA,Colv=NA,labRow=NA,labCol=NA)
  } else if(type=='temp') {
    if(is.null(main)) { main <- paste(r$id,' (',r$zip,') temperature info',sep='') }
    #colors = rainbow(24)
    colorBy = r$dates$hour
    legendLabs = paste("hr ", 1:24)
    colors = colorRampPalette(brewer.pal(11,"PRGn"))(length(unique(colorBy)))
    colors = colors[(1:length(colors) + 8) %% length(colors)] # shift values backwards by 8 hours (so color discontinuity is in the afternoon)
    plot(r$tout,r$kw,xlab='Tout',ylab='kW',main=main,type='p',col=mapColors(colorBy,colors),cex=0.8,pch=19) # rainbow(n, start=2/6, end=1)
    legend('right', legendLabs, fill=colors, ncol = 1, cex = 0.5)
    if(!is.null(extraPoints)) {
      color = 'blue'
      points(extraPoints$x,extraPoints$y,col=color )
           #lwd=2 )
    }
  } else if(type=='hourly') {
    if(is.null(main)) { main <- paste(r$id,' (',r$zip,') hourly info',sep='') }
    op <- par(no.readonly = TRUE)
    grid = cbind(matrix(c(1:24),nrow=4,ncol=6,byrow=TRUE),c(25,25,25,25))
    layout(grid,widths=c(rep(2,6),3))
    par(oma=c(2,3,3,0),mar=c(1,0,1,0))# Room for the title
    
    pallete = colorRampPalette(brewer.pal(11,"Spectral"))(12) #rainbow(12) #,start=2/6, end=1)
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
           #col='black',
           pch=subset(colvals,sub),
           xlim=xlm,ylim=ylm,yaxt=yax,xaxt=xax,cex=1.0
           )
      if(length(estimates) > 1) {
        if(length(estimates[,i+1]) > 1) {
          fit = estimates[,i+1]
          color = 'black'
          if (fit['AIC_0'] < fit['AIC_cp']) color = 'gray'
          if (fit['nullModelTest'] > 0.1)   color = 'red'
          par(new=T)
          qe = quickEst(fit['cp'],fit['(Intercept)'],fit['lower'],fit['upper'],
                   c(min(subset(r$tout,sub),na.rm=T),max(subset(r$tout,sub),na.rm=T)))
          #print(qe[,1])
          plot(qe[,1],qe[,2],
                type='l',
                xlim=xlm,ylim=ylm,axes=F,lwd=2,col=color )
          par(new=F)
        }
      }
      grid()
      text(mean(xlm),ylm[2] * 0.95,paste('hr',i),font=2, cex=1.2)
    }
    mtext(main, line=0, font=2, cex=1.2,outer=TRUE)
    #par(xpd=TRUE)
    #plot.new()
    plot(1, type = "n", axes=FALSE, xlab="", ylab="")
    legend('center', month.abb, col=pallete, pch=1:24, ncol = 1, cex = 1.3)
    #par(xpd=FALSE)
    par(op)
  }
  
}