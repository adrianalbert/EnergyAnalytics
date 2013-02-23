weatherFeatures = function(w){ # r is an instance of WeatherClass
  wMeans = colMeans(w$rawData[,-1],na.rm=T)
  features = c(zip5=w$zip,wMeans)
  return(features)
}

getWeatherSummary = function(w) {
  summaryFile = paste(getwd(),'/weatherSummary.RData',sep='')
  if(file.exists(summaryFile)) {
    load(summaryFile)
  }
  else {
    zips = db.getZips()
    i = 0
    n = length(zips)
    wList = as.list(rep(NA,length(zips)))
    for(zip in zips) {
      i = i+1
      print(paste(zip,'(',i,'/',n,')'))
      wList[[i]] = weatherFeatures(WeatherClass(zip,doMeans=F))
    }
    weatherSummary = do.call(rbind,wList)
    save(weatherSummary,file=summaryFile)
  }
  return(weatherSummary)
}