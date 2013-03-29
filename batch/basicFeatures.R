basicFeatures = function(r){ # r is an instance of ResDataClass
  data <- as.matrix(r$kwMat)
  id   <- r$id

  # todo: apply/quantile is slow. Is it worth softening min,max this way?
  qtle     <- apply(data,1,FUN=quantile,c(0.01,0.99),na.rm=T)
  dMin     <- qtle[1,]
  dMax     <- qtle[2,]
  
  dMean    <- rowMeans(data,dims=1,na.rm=T) # rowMeans is faster than apply with FUN=mean ?rowMeans for details
  
  dRange   <- dMax - dMin
  dHalfway <- dMin + (dMax - dMin) / 2
  dHighD   <- rowSums(data > dHalfway,dims=1)
  dMn2mx   <- dMin / dMax
  dN2d     <- data[,4] / data[,16] # 4am comp to 4pm
  # we need dates for this one...
  #wkend = something to do with the DOW in the dates
  #wkdays = ! wkend
  #wknd2wk = mean(wkend) / mean(wkdays)

  soft <- quantile(r$kw,c(0.03,0.97),na.rm=T)
  max  <- soft[2]
  min  <- soft[1]
  nObs <- prod(dim(data))

  kw.mean     = mean(r$kw,na.rm=T)
  
  kw.tout.cor        = cor(r$kw,r$tout,               use='complete.obs') # correlation between tout and kw
  kw.pout.cor        = cor(r$kw,r$w('pout'),          use='complete.obs') # correlation between pressure and kw
  kw.var             = var(r$kw/kw.mean,              use='complete.obs') # normed by mean of kW
  daily.kw.var       = var(rowMeans(r$kwMat)/kw.mean, use='complete.obs') # normed by mean of kW
  daily.kw.min.var   = var(dMin/kw.mean,              use='complete.obs') # normed by mean of kW
  daily.kw.max.var   = var(dMax/kw.mean,              use='complete.obs') # normed by mean of kW
  
  lags = 0:24
  lag.cor = apply(as.matrix(lags),  1,function(x) cor(r$kw,lag(r$tout,x),use='complete.obs')) 
  lag.ma  = apply(as.matrix(lags[-1]),1,function(x) cor(r$kw,ma(r$tout,x), use='complete.obs')) # no width 0
  names(lag.cor) <- c(paste('lag',lags,  sep=''))
  names(lag.ma)  <- c(paste('ma', lags[-1],sep=''))

  basics = c(id=id,
             nObs=nObs,
             kw.mean=kw.mean,
             max=max,  # will be named 'max.97%' due to quantile origin
             min=min,  # will be called 'min.3%' due to quantile origin
             kw.var=kw.var,
             daily.kw.var=daily.kw.var,
             daily.kw.min.var=daily.kw.min.var,
             daily.kw.max.var=daily.kw.max.var,
             kw.pout.cor=kw.pout.cor,
             kw.tout.cor=kw.tout.cor,
             lag.cor,
             lag.ma)
  
  dailyFeatures <- cbind(dMax,dMin,dMean,dRange,dHighD,dMn2mx,dN2d)
  dailyFeatures[dailyFeatures == Inf] = NA # these will be caught by the na.rm in the mean fn
  colnames(dailyFeatures) <- c('max','min','mean','range','dur','mn2mx','n2d')
  
  # return the means across all days
  ret <- c(basics,apply(dailyFeatures,2,FUN=mean,na.rm=T))
  #print(ret)
  return(ret)
}