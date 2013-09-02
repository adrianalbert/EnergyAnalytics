getCAISO = function() {
  caisoFile = 'caiso/caiso.RData'
  if(file.exists(caisoFile)) {
    load(caisoFile)
  } else {
    caisoRaw = read.csv('caiso/caiso.csv')
    dates = strptime(caisoRaw$date,'%m/%d/%Y')
    caisoRaw$year = dates$year + 1900
    caisoRaw$mon = dates$mon
    caisoRaw$mday = dates$mday
    cam = as.matrix(caisoRaw[,9:32])
    
    dateMat = t(apply(as.matrix(as.POSIXct(dates)),1,FUN=function(d) { d + 3600 * 0:23 }))
    dateTime = as.POSIXlt(as.vector(t(dateMat)),origin='1970-01-01')
    dt = as.POSIXct(dateTime)
    kw = as.vector(t(cam))
    caiso = data.frame(kw=kw,date=as.POSIXct(dateTime),year=dateTime$year+1900,month=dateTime$mon+1,day=dateTime$mday,hour=dateTime$hour)
    save(list=c('caiso'),file=caisoFile)
  }
  return(caiso)
}

annualTopDemand = function(n=5,years=2009:2012) {
  caiso = getCAISO()
  out = data.frame(c())
  for(yr in years){
    yrData = caiso[caiso$year == yr,]
    tops = match(sort(yrData$kw,decreasing=T)[1:n],yrData$kw)
    out = rbind(out,yrData[tops,])
  }
  return(out)
}

annualPctDemand = function(pct=0.99,years=2009:2012) {
  caiso = getCAISO()
  out = data.frame(c())
  for(yr in years){
    yrData = caiso[caiso$year == yr,]
    threshold = quantile(yrData$kw,pct,na.rm=T)
    tops = which(yrData$kw >= threshold)
    out = rbind(out,yrData[tops,])
  }
  return(out)
}

peakCoincidentReadings = function(peakHrs=1,peakPct=NULL) {
  allZips  <- DATA_SOURCE$getZips(useCache=T)
  if(is.null(peakPct)) {
    topLoad = annualTopDemand(peakHrs)
  } else {
    topLoad = annualPctDemand(peakPct)
  }
  topData = data.frame(c())
  topDays = sort(unique(strptime(paste(topLoad$year,topLoad$month,topLoad$day,sep='-'),'%Y-%m-%d')))
  n = length(allZips)
  z = 0
  tic('peakCoincidentReadings')
  for (zip in allZips) {
    z = z + 1
    print(paste(z,'of',n))
    zipData <- DATA_SOURCE$getAllData(zip,useCache=T)
    days = strptime(zipData$DATE,'%Y-%m-%d') # convert to foolproof day comparison format
    # todo: this would be faster if it grew a pre-allocated list and 
    # then ran as.data.frame(do.call(rbind,list)) after...
    topData = rbind(topData,zipData[days %in% topDays,]) 
  }
  toc('peakCoincidentReadings')
  kwMat = NULL
  if(dim(topData)[2] > 40) { # average 15 minute data into hourly.
    kwMat96 = as.matrix(topData[,4:99])/1000
    kwMat = sapply(t(1:24),function(x) rowMeans(kwMat96[,(4*(x-1)+1):(4*x)]))
    topData = cbind(topData[,1:3],kwMat)
    names(topData)[-1:-2] <- c('day',paste('h',1:24,sep=''))
  }
  coincidentReadings = data.frame(c())
  for(td in topLoad$date){
    td = as.POSIXlt(td,origin='1970-01-01')
    print(td)
    day = strptime(td,'%Y-%m-%d')
    h = td$hour
    newReadings = topData[topData$day == day,c(1,2,3,4+h)] # h is zero based, with the first hour in the 4th column
    names(newReadings)[4] <- 'kw'
    print(names(newReadings))
    coincidentReadings = rbind(coincidentReadings,newReadings) 
  }
  return(coincidentReadings)
}

test = F
if(test) {
  caiso = getCAISO()
  onePct = annualPctDemand(0.99)
  top10 = annualTopDemand(10)
  plot(caiso$date,caiso$kw,type='l',main='CAISO hourly system demand',ylab='kW',xlab='date')
  points(onePct$date,onePct$kw,col='red')
  points(top10$date,top10$kw,col='blue')
  
  coincidentReadings = peakCoincidentReadings(peakHrs=1)
  plot((1:length(coincidentReadings$kw))/length(coincidentReadings$kw)*100,
       cumsum(sort(coincidentReadings$kw,decreasing=T))/sum(coincidentReadings$kw,na.rm=T),
       xlab='percentile',ylab='fraction of load',
       main=paste('Cumulative load at peak hour of system demand (N=',length(unique(coincidentReadings$sp_id)),')'))
  zipContribution = dcast(coincidentReadings,zip5 ~ .,value.var='kw',fun.aggregate=mean,na.rm=T)
  names(zipContribution)[2] <- 'mean.kw'
  calMap(df=zipContribution,plotCol='mean.kw',intervalStyle='fixed',intervalBreaks=c(0,1,2,3,4,5,10),main='Mean kW demand at peak load by zip code')

}