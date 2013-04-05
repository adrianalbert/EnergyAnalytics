weatherFeatures = function(w){ # r is an instance of WeatherClass
  wMeans = colMeans(w$rawData[,-1],na.rm=T)
  features = c(w$zip,wMeans)
  return(features)
}

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