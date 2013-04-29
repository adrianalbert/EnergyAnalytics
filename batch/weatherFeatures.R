weatherFeatures = function(w){ # r is an instance of WeatherClass
  summerMon = 4:8 # May through Sept - zero based
  summerSubset = as.POSIXlt(w$rawData$dates)$mon %in% summerMon
  wMeans = colMeans(w$rawData[,-1],na.rm=T)
  sMeans = colMeans(subset(w$rawData[,-1],subset=summerSubset),na.rm=T)
  names(sMeans) = paste('s.',names(sMeans),sep='')
  features = c(zip5=w$zip,wMeans,sMeans)
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
    weatherSummary = data.frame(do.call(rbind,wList))
    colnames(weatherSummary)[1] <- c('zip5')
    for(col in colnames(weatherSummary)) {
      if(is.factor(weatherSummary[[col]])) {
        # as.numeric(levels(f))[f] is the preferred way to convert factors to numerics
        weatherSummary[[col]] = as.numeric(levels(weatherSummary[[col]] ))[weatherSummary[[col]]]
      }
    }
    save(weatherSummary,file=summaryFile)
  }
  return(weatherSummary)
}

#print(weatherFeatures(WeatherClass(94610)))