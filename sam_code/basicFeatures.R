Mode <- function(x) {
  ux <- unique(x)
  tab = tabulate(match(x, ux))
  maxCount = max(tab)
  
  modes = ux[tab == maxCount] # there could be more than one
  sample(rep(modes,2),1) # random choice breaks ties; single mode is deterministic
}

basicFeatures = function(r){ # r is an instance of ResDataClass
  data <- as.matrix(r$kwMat)
  id   <- r$id
  
  monthTotals = sapply(split(r$kw,factor(format(r$dates,'%b'),levels=month.abb)),mean,na.rm=T) * 30 * 24
  names(monthTotals) <- paste('kw.total.',names(monthTotals),sep='')
  
  summerMon = 4:8 # May through Sept - zero based
  summerSubset = as.POSIXlt(r$dates)$mon %in% summerMon
  
  # todo: apply/quantile is slow. Is it worth softening min,max this way?
  qtle     <- apply(data,1,FUN=quantile,c(0.01,0.99),na.rm=T)
  dMin     <- qtle[1,]
  dMax     <- qtle[2,]
  
  hMean    <- colMeans(data,na.rm=T)
  dMean    <- rowMeans(data,dims=1,na.rm=T) # rowMeans is faster than apply with FUN=mean ?rowMeans for details
  tMean    <- rowMeans(as.matrix(r$toutMat),dims=1,na.rm=T)
  # index of the hottest 10% of days
  t90Idx   <- which(tMean > quantile(tMean,0.9,na.rm=T))
  # index (hour) of the max demand on the hottest days
  t90kwhr = apply(data[t90Idx,],1,FUN=which.max)
  # index (hour) of the max temperature on the hottest days
  t90thr = apply(r$toutMat[t90Idx,],1,FUN=which.max)
  
  # index of the coldest 10% of days
  t10Idx   <- which(tMean > quantile(tMean,0.1,na.rm=T))
  # index (hour) of the max demand on the coldest days
  t10kwhr = apply(data[t10Idx,],1,FUN=which.max)
  # index (hour) of the max temperature on the coldest days
  t10thr = apply(r$toutMat[t10Idx,],1,FUN=which.max)
  
  # index of the top hours of consumption
  kw90Idx = which(r$kw > quantile(r$kw,0.9,na.rm=T))
  kw90hr = r$dates[kw90Idx]$hour
  
  maxHOD <- which.max(hMean)
  kw90   <- Mode(kw90hr)
  t90kw  <- Mode(t90kwhr)
  t90t   <- Mode(t90thr)
  t10kw  <- Mode(t10kwhr)
  t10t   <- Mode(t10thr)
  
  dRange   <- dMax - dMin
  dHalfway <- dMin + (dMax - dMin) / 2
  dHighD   <- rowSums(data > dHalfway,dims=1)
  dMn2mx   <- dMin / dMax
  dN2d     <- rowMeans(data[,2:5],na.rm=T) / rowMeans(data[,16:19],na.rm=T) # 2-5am comp to 4-7pm
  nv2dv    <- var(as.vector(data[,2:5]),na.rm=T) / var(as.vector(data[,16:19]),na.rm=T) # var for night 2-5am comp to evening 4-7pm
  # we need dates for this one...
  #wkend = something to do with the DOW in the dates
  #wkdays = ! wkend
  #wknd2wk = mean(wkend) / mean(wkdays)

  soft <- quantile(r$kw,c(0.03,0.97),na.rm=T)
  max  <- soft[2]
  min  <- soft[1]
  nObs <- prod(dim(data))

  kw.mean        = mean(r$kw,na.rm=T)
  kw.mean.summer = mean(r$kw[ summerSubset],na.rm=T)
  kw.mean.winter = mean(r$kw[!summerSubset],na.rm=T)

  therm.mean = NA
  therm.mean.summer = NA
  therm.mean.winter = NA
  therm.min = NA
  
  if(length(r$therms > 0)) {
    therm.mean     = mean(r$therms,na.rm=T)
    gasSummer      = as.POSIXlt(r$gasDays)$mon %in% summerMon
    therm.mean.summer = mean(r$therms[ gasSummer],na.rm=T)
    therm.mean.winter = mean(r$therms[!gasSummer],na.rm=T)
    thermdf = data.frame(therm=r$therms,day=r$gasDays,month=as.POSIXlt(r$gasDays)$mon)
    thermMonth = dcast(thermdf,month ~ .,mean,value.var='therm')
    names(thermMonth) <- c('month','therm')
    therm.min = mean(thermMonth$therm[thermMonth$month %in% 5:7])
  }
  maxIdx = which.max(r$kw)
  minIdx = which.min(r$kw)
  
  max.hr.kw      = r$kw[maxIdx]
  min.hr.kw      = r$kw[minIdx]
  
  max.hr.tout = r$tout[maxIdx]
  min.hr.tout = r$tout[minIdx]
  
  #print(r$dates[maxIdx])
  max.hr.date = as.POSIXct(r$dates[maxIdx])
  min.hr.date = as.POSIXct(r$dates[minIdx])
  
  maxIdx = which.max(dMean)
  minIdx = which.min(dMean)
  
  max.day.kw      = dMean[maxIdx] * 24
  min.day.kw      = dMean[minIdx] * 24
  
  max.day.tout = tMean[maxIdx]
  min.day.tout = tMean[minIdx]
  
  max.day.date = r$days[maxIdx]
  min.day.date = r$days[minIdx]
  
  # calculate percentiles of touts during max/min days
  #print(max.day.tout)
  #print(min.day.tout)
  max.day.pct = mean(tMean <= max.day.tout)
  min.day.pct = mean(tMean <= min.day.tout)
  
  kw.tout.cor        = cor(r$kw,r$tout,               use='complete.obs') # correlation between tout and kw
  kw.pout.cor = NA
  #if(any(! is.na(r$w('pout')))) {
  #  kw.pout.cor = cor(r$kw,r$w('pout'), use='complete.obs') # correlation between pressure and kw
  #}
  kw.var             = var(r$kw/kw.mean,                use='complete.obs') # normed by mean of kW
  kw.var.summer      = var(r$kw[summerSubset]/kw.mean,na.rm=T)   # normed by mean of kW
  kw.var.winter      = var(r$kw[!summerSubset]/kw.mean,na.rm=T)  # normed by mean of kW
  
  daily.kw.var       = var(rowMeans(r$kwMat)/kw.mean,   use='complete.obs') # normed by mean of kW
  daily.kw.min.var   = var(dMin/kw.mean,                use='complete.obs') # normed by mean of kW
  daily.kw.max.var   = var(dMax/kw.mean,                use='complete.obs') # normed by mean of kW
  
  lags = 0:24
  #lag.cor = apply(as.matrix(lags),  1,function(x) cor(r$kw,lag(r$tout,x),use='complete.obs')) 
  lag.ma  = apply(as.matrix(lags[-1]),1,function(x) cor(r$kw,ma(r$tout,x), use='complete.obs')) # no width 0
  #names(lag.cor) <- c(paste('lag',lags,  sep=''))
  #names(lag.ma)  <- c(paste('ma', lags[-1],sep=''))
  max.MA = which.max(lag.ma)
  basics = c(id=id,
             nObs=nObs,
             kw.mean=kw.mean,kw.mean.summer=kw.mean.summer,kw.mean.winter=kw.mean.winter,
             kw.total=kw.mean * 365 * 24,
             therm.mean=therm.mean,therm.mean.summer=therm.mean.summer,therm.mean.winter=therm.mean.winter,
             therm.total=therm.mean * 365 * 24,
             therm.min=therm.min,
             max=max,  # will be named 'max.97.' due to quantile origin
             min=min,  # will be called 'min.3.' due to quantile origin
             kw.var=kw.var,kw.var.summer=kw.var.summer,kw.var.winter=kw.var.winter,
             max.hr.kw=max.hr.kw,min.hr.kw=min.hr.kw,max.hr.tout=max.hr.tout,min.hr.tout=min.hr.tout,max.hr.date=max.hr.date,max.hr.date=max.hr.date,
             max.day.kw=max.day.kw,min.day.kw=min.day.kw,max.day.tout=max.day.tout,min.day.tout=min.day.tout,max.day.date=max.day.date,max.day.date=max.day.date,
             max.day.pct=max.day.pct,min.day.pct=min.day.pct,
             t90kw=t90kw,t90t=t90t,t10kw=t10kw,t10t=t10t,maxHOD=maxHOD,kw90=kw90,
             daily.kw.var=daily.kw.var,
             daily.kw.min.var=daily.kw.min.var,
             daily.kw.max.var=daily.kw.max.var,
             kw.pout.cor=kw.pout.cor,
             kw.tout.cor=kw.tout.cor,
             max.MA=max.MA,
             monthTotals,
             nv2dv=nv2dv
             #lag.cor,
             #lag.ma
             )
  
  dailyFeatures <- cbind(dMax,dMin,dMean,dRange,dHighD,dMn2mx,dN2d)
  dailyFeatures[dailyFeatures == Inf] = NA # these will be caught by the na.rm in the mean fn
  colnames(dailyFeatures) <- c('max','min','mean','range','dur','mn2mx','n2d')
  monthFeatures = t(sapply(split(data.frame(dailyFeatures),factor(format(r$days,'%b'),levels=month.abb)),colMeans,na.rm=T))
  # return the means across all days
  augF = monthFeatures[c('Aug'),]
  names(augF) <- paste('Aug_',names(augF),sep='')
  janF = monthFeatures[c('Jan'),]
  names(janF) <- paste('Jan_',names(janF),sep='')
  ret <- c(basics,apply(dailyFeatures,2,FUN=mean,na.rm=T),augF,janF)
  #print(ret)
  return(ret)
}