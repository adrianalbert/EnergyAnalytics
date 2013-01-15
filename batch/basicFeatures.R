basicFeatures = function(data,id=NULL){
  dMax     <- apply(data,1,FUN=quantile,.99,na.rm=TRUE)
  dMin     <- apply(data,1,FUN=min,.1,na.rm=TRUE)
  dMean    <- rowMeans(data,dims=1,na.rm=TRUE) # rowMeans is faster than apply with FUN=mean ?rowMeans for details
  dRange   <- dMax - dMin
  dHalfway <- dMin + (dMax - dMin) / 2
  dHighD   <- rowSums(data > dHalfway,dims=1)
  dMn2mx   <- dMin / dMax
  dN2d     <- data[,3] / data[,16] # 3am comp to 4pm
  # we need dates for this one...
  #wkend = something to do with the DOW in the dates
  #wkdays = ! wkend
  #wknd2wk = mean(wkend) / mean(wkdays)
  
  softMax  <- quantile(data,0.97,na.rm=TRUE)
  softMin  <- quantile(data,0.03,na.rm=TRUE)
  nObs     <- prod(dim(data))
  variance <- var(as.vector(t(as.matrix(data))),na.rm=TRUE)
  # todo: as dataframe?
  dailyFeatures <- cbind(dMax,dMin,dMean,dRange,dHighD,dMn2mx,dN2d)
  dailyFeatures[dailyFeatures == Inf] = NA # these will be caught by the na.rm in the mean fn
  colnames(dailyFeatures) <- c('max','min','mean','range','dur','mn2mx','n2d')
  # return the means across all days
  ret <- c(apply(dailyFeatures,2,FUN=mean,na.rm=TRUE),softMax,softMin,nObs=nObs,id=id)
  return(ret)
}