quantileOutlierDates = function(e,dates,q=0.95) {
  if(length(e) != length(dates)) {
    stop(paste('[occupancy.quantileOutlierDates] Errors and dates have different lengths:',length(e),length(dates)))
  }
  threshold = quantile(e,q,na.rm=T)
  outliers = which(e > threshold) # which strips NAs
  return(dates[outliers])
}

monthlyQuantileOutlierDates = function(e,dates,q=0.95) {
  if(length(e) != length(dates)) {
    stop(paste('[occupancy.monthlyQuantileOutlierDates] Errors and dates have different lengths:',length(e),length(dates)))
  }
  outDates = data.frame(c())
  #monthProfiles = c()
  if('POSIXct' %in% class(dates)) { dates = as.POSIXlt(dates) }
  for(m in sort(unique(dates$mon))) {
    mnth = which(dates$mon == m)
    threshold = quantile(e[mnth],q,na.rm=T)
    outliers = which(e[mnth] > threshold) # which strips NAs
    odts = as.POSIXct(dates[mnth][outliers])
    #monthProfiles = rbind(monthProfiles,eventDensities(odts)$hr)
    outDates = rbind(outDates,data.frame(d=odts))
  }
  #matplot(t(as.matrix(monthProfiles[,1:24])))
  return(outDates$d)
}

# this supports combining the dates from several runs by adding a count for each matching date
overlappingOutlierDates = function(summariesDF) { # find the date of outliers from multiple model runs and count the overlapping dates
  if(! 'data.frame' %in% class(summariesDF)) { summariesDF = data.frame(summariesDF) }
  outDates = data.frame(c())
  residualsL = summariesDF$residuals
  datesL     = summariesDF$dates
  for(i in 1:length(residualsL)) {
    residuals = residualsL[[i]]
    dates     = datesL[[i]]
    outDates = rbind(outDates,data.frame(dates=quantileOutlierDates(residuals,dates,q=0.95)))
  }
  #print(sort(outDates$dates))
  dateCounts = data.frame(table(outDates)) # names will be outDates, Freq
  names(dateCounts)[1] <- 'dates'
  #print(names(dateCounts))
  #print(levels(dateCounts$dates))
  dateCounts$dates = as.POSIXct(levels(dateCounts$dates))[dateCounts$dates]
  
  return(dateCounts)
}

eventDensities = function(outDates) {
  if('POSIXct' %in% class(outDates)) {
    outDates = as.POSIXlt(outDates)
  }
  wday = as.array(table(factor(outDates$wday,levels=0:6)))
  names(wday) = c('Su','Mo','Tu','We','Th','Fr','Sa')[as.numeric(names(wday))+1]
  wday = c(wday/sc,scale=sum(wday))
  
  hr= as.array(table(factor(outDates$hour,levels=0:23)))
  hr = c(hr/sc,scale=sum(hr))
  
  wknd= as.array(table(factor(outDates[outDates$wday %in% c(0,6)]$hour,levels=0:23)))
  wknd = c(wknd/sc*7/2,scale=sum(wknd))
  
  wkdy= as.array(table(factor(outDates[outDates$wday %in% c(1:5)]$hour,levels=0:23)))
  wkdy = c(wkdy/sc*7/5,scale=sum(wkdy))
  
  mon = as.array(table(factor(outDates$mon,levels=0:11)))
  names(mon) = month.abb
  mon = c(mon/sc,scale=sum(mon))
  
  how = as.array(table(factor(outDates$wday*24+outDates$hour,levels=0:(24*7-1))))
  how = c(how/sc,scale=sum(how))
  
  return(list(wday=wday,hr=hr,wknd=wknd,wkdy=wkdy,mon=mon,how=how))
}

quantileDensities = function(e,dates,q=0.95,monthly=F) {
  qd=NA
  if(monthly) {
    qd = monthlyQuantileOutlierDates(e,dates,q)
  } else {
    qd = quantileOutlierDates(e,dates,q)
  }
  return(eventDensities(qd))
}

distance = function(centers,membership,actual) {
  for(ci in sort(unique(membership))) {
    members = actual[membership == ci,]
    center = centers[ci,]
    print(mean(rowSums(abs(members - center/4))))
  }
}

clusterSumOfSquares = function(data,countSeq=1:15,iter.max=10) {
  # Determine number of clusters
  max_err <- c()
  for (i in countSeq){
    km <- kmeans(data,centers=i,iter.max=iter.max)
    counts = table(sort(km$cluster))
    worst = c()
    for(j in 1:dim(km$centers)[1]){
      centerL2 = sqrt(sum(km$centers[j,]^2)) # magnitude of the center
      diffL2 = sqrt(rowSums((data[km$cluster == j,] - km$centers[rep(j,counts[j]),])^2)) # difference of each member and the center squared
      worst[j] = mean(diffL2 / centerL2)
    }
    max_err[i] <- mean(worst)
  }
  return(max_err)
}


showClusters = function(km,obs,ncol=4,xlabs=NULL,main=NULL,ylab=NULL,xlab=NULL,...) {
  op <- par(no.readonly = TRUE)
  par( mfrow=c(ceiling(nrow(km$centers)/ncol),ncol), oma=c(2,2,3,0),mar=c(2,1,2,1)) # Room for the title
  for(c in sort(unique(km$cluster))){
    hmap(data=as.matrix(obs[km$cluster == c,]),main=paste(c,' (n=',sum(km$cluster==c),')',sep=''),xlabs=xlabs,...)
  }
  if(!is.null(main)) { mtext(main, line=0, font=2, cex=1.2, outer=TRUE) }
  if(!is.null(ylab)) { mtext(ylab, line=0, font=2, cex=0.8, side=2, outer=TRUE) }
  if(!is.null(xlab)) { mtext(xlab, line=0, font=2, cex=0.8, side=1, outer=TRUE) }
  par(op)
}

testData = function() {
  setwd('f:/dev/pge_collab/EnergyAnalytics/batch')
  sam = read.csv('data/sam_home_data.csv')
  #plot(sam$reading)
  sam$dates = as.POSIXlt(sam$date,tz="PST8PDT",'%Y-%m-%d %H:%M')
  sam$kw = sam$reading / 1000
  
  # hack in fake weather data so we can use canned regressorDF function
  sam$tout = sam$kw
  sam$pout = sam$kw
  sam$rain = sam$kw
  sam$dp = sam$kw
  sam$rh = sam$kw
  return(sam)
}

test = F
if(test) {
  sam= testData()
  df = regressorDF(sam)
  m = lm(kw ~ HOW -1,data=df)
  e = residuals(m)
  #plot(e)
  
  dens = quantileDensities(e,sam$dates)
}

doFigures = F
if(doFigures) {
  occupancyDir = 'C:/Users/Sam/Dropbox/writing/occupancy_paper/'
  load(file.path(getwd(),'w_results_occ','occupancy.RData'))
  #dg = DescriptorGenerator(name='toutPieces',genImpl=toutPieces24Generator,subset=list(all="TRUE"),terms='+ HOW - 1',breaks=c(65),diverge=T)
  #dg = DescriptorGenerator(name='toutPieces',genImpl=toutPieces24Generator,subset=list(all="TRUE"),terms='+ HOD + DOW - 1',breaks=c(65),diverge=T)
  #dg = DescriptorGenerator(name='toutPieces',genImpl=toutPieces24Generator,subset=list(all="TRUE"),terms='+ WKND:HOD - 1',breaks=c(65),diverge=T)
  dgWindow = DescriptorGenerator(name='toutPieces',genImpl=toutPieces24Generator,subset='idxSteps',terms='+ HOW - 1',breaks=c(60,70),diverge=T)
  dg       = DescriptorGenerator(name='toutPieces',genImpl=toutPieces24Generator,subset=list(all="TRUE"),terms='+ HOW - 1',breaks=c(60,70),diverge=T)
  
  # warning...this uses a Stanford home, but the rest is baed on Wharton data
  r = ResDataClass(553991005,93304,useCache=T,doSG=F)
  #r = ResDataClass(554622151,93304,useCache=T,doSG=F)
  #r = ResDataClass(1064423310,93304,useCache=T,doSG=F)
  reg = dg$run(r)$summaries # standards regression model run with piecewise temeprature interacted with HOD and HOW fixed effects
  edReg = data.frame(data.frame(e=reg[[1,'residuals']], d=reg[[1,'dates']], p=reg[[1,'prediction']])) # same regression over sliding window
  
  window = dgWindow$run(r)$summaries
  ed = data.frame(c())
  for(i in 1:dim(window)[1]){
    ed = rbind(ed,data.frame(e=window[[i,'residuals']], d=window[[i,'dates']], p=window[[i,'prediction']]))
  }
  outliers = overlappingOutlierDates(window) # dates, Freq
  
  
  #plot(ed$d,ed$e)
  #outIdx = which(ed$d %in% outliers$dates)
  #plot(ed$d,ed$e)
  #points(ed$d[outIdx],ed$e[outIdx],col='red',cex=0.4,pch=15)
  
  cols = c('#000000','#dddddd','#FED976','#FD8D3C','#B10026')
  labels = c('Observed consumption','Outside temp. (scaled)','1 of 3 matches','2 of 3 matches','3 of 3 matches')
  plot(r$dates,r$kw,ylab='kW',xlab='date',ylim=c(0,max(r$kw,na.rm=T)* 1.1),main='Observations and modeled outliers',col=cols[1])
  
  points(r$dates,r$tout / 21,col=cols[2],type='l')
  points(r$dates,r$kw,col=cols[1])
  for(hitCount in 1:3){
    outIdx = which(as.POSIXct(r$dates) %in% outliers$dates[outliers$Freq==hitCount])
    points(r$dates[outIdx],r$kw[outIdx],col=cols[hitCount+2],cex=0.7,pch=16)
  }
  legend('top',inset=c(0.0,0.05),lty=c(NA,1,NA,NA,NA),
         legend=c('Observed consumption','Outside temp. (scaled)','1 of 3 matches','2 of 3 matches','3 of 3 matches'),
         col=cols,cex=0.7,pch=c(1,NA,16,16,16))
  dev.copy2pdf(file = paste(occupancyDir,'window_matches.pdf',sep=''),width=12,height=8)
  
  plot(r$tout,r$kw,col=cols[1],ylab='kW',xlab='Outside temp.',main='Temperature vs. kW with outliers')
  #points(r$dates,r$tout / 25,col='grey',type='l')
  #points(r$dates,r$kw)
  for(hitCount in 1:3){
    outIdx = which(as.POSIXct(r$dates) %in% outliers$dates[outliers$Freq==hitCount])
    points(r$tout[outIdx],r$kw[outIdx],col=cols[hitCount+2],cex=0.7,pch=16)
  }
  legend('topleft',inset=c(0.02,0.05),
         legend=labels[-2],
         col=cols[-2],cex=0.7,pch=c(1,16,16,16))
  dev.copy2pdf(file = paste(occupancyDir,'window_matches_scatter.pdf',sep=''),width=12,height=8)
  
  dens3 = eventDensities(outliers$dates[outliers$Freq<4])
  dens2 = eventDensities(outliers$dates[outliers$Freq<3])
  dens1 = eventDensities(outliers$dates[outliers$Freq<2])
  
  plot(dens3$mon[1:12],type='l')
  points(dens2$mon[1:12],type='l',col='blue')
  points(dens1$mon[1:12],type='l',col='green')
  plot(dens3$wkdy[1:24]/meankw)
  
  qod  = edReg$d %in% quantileOutlierDates(edReg$e,edReg$d)
  mqod = edReg$d %in% monthlyQuantileOutlierDates(edReg$e,edReg$d)
  op <- par(no.readonly = TRUE)
  par( mfrow=c(2,1), oma=c(2,2,1,0),mar=c(2,4,2,2))
  plot(r$dates,r$kw,col='#444444',main='Observations (grey) & model (green)',ylab='kW')
  points(edReg$d,edReg$p, col='#22ff22',cex=0.4,pch=15)
  grid()
  legend('top',inset=c(0,0.03),legend=c('observed','predicted'),col=c('#444444','#22ff22'),pch=c(1,16),horiz=T,cex=0.6)
  plot(edReg$d,edReg$e,xlab='observation #',col='#444444',main='Residuals (grey) & monthly top 5% (red)',ylab='kW')
  #points(qod,e[qod],col='blue',ylab='kW',cex=0.4)
  points(edReg$d[mqod],edReg$e[mqod],col='red',cex=0.4,pch=16)
  #mtext('Model predictions, residuals, and outliers', line=0, font=2, cex=1.4, outer=TRUE)
  grid()
  legend('top',inset=c(0,0.03),legend=c('model error','monthly top 5%'),col=c('#444444','red'),pch=c(1,16),horiz=T,cex=0.6)
  par(op)
  
  dev.copy2pdf(file = paste(occupancyDir,'single_model.pdf',sep=''),width=12,height=8)
  
  meankw = c()
  for(h in sort(unique(r$dates$hour))) {
    hr = which(r$dates$hour == h)
    meankw = c(meankw,mean(r$kw[hr],na.rm=T))
  }
  meankw = meankw / sum(meankw)
  
  meane = c()
  for(h in sort(unique(r$dates$hour))) {
    hr = which(r$dates$hour == h)
    meane = c(meane,mean(e[hr],na.rm=T))
  }
  meane = meane / sum(meane)
  
  dens = quantileDensities(edReg$e,edReg$d)
  densM = quantileDensities(edReg$e,edReg$d,monthly=T)
  plot(densM$hr[1:24]*densM$hr[25],type='o',main='Count of occupancy events by hour',ylab='count',xlab='Hour of day',pch=16,col='red',xaxt='n',ylim=c(0,max(densM$hr[1:24]*densM$hr[25])))
  #points(densM$hr[1:24]*densM$hr[25],type='o',col='blue')
  points(meankw*dens$hr[25],col='#999999',type='l',lty=2)
  points(dens3$hr[1:24]*densM$hr[25],col='blue',type='o')
  points(dens2$hr[1:24]*densM$hr[25],col='#6666ff',type='o')
  points(dens1$hr[1:24]*densM$hr[25],col='#aaaaff',type='o')
  grid()
  legend('topleft',inset=c(0.01,0.03),c('Count of top 5%','Rescaled mean load'),lty=c(1,2),pch=c(16,-1),cex=0.8,col=c('red','#999999'))
  xlabs=0:23
  axis(1, at=1:length(xlabs), labels=xlabs)
  #points(meane,col='green')
  dev.copy2pdf(file = paste(occupancyDir,'p_occ_event.pdf',sep=''),width=12,height=6)
    
  fits = clusterSumOfSquares(as.matrix(occ.hr.m[,2:25]),seq(2,60,4),20)
  vals = which(! is.na(fits))
  plot(vals,fits[vals],type='o',ylab='error',xlab='# of clusters',main='Fit error vs. number of clusters')
  dev.copy2pdf(file = paste(occupancyDir,'clutser_err.pdf',sep=''))
  
  fits = clusterSumOfSquares(as.matrix(occ.wday.m[,2:8]),seq(2,20,1),20)
  vals = which(! is.na(fits))
  plot(vals[-1],fits[vals][-1],type='o',ylab='error',xlab='# of clusters',main='Fit error vs. number of clusters')
  dev.copy2pdf(file = paste(occupancyDir,'wday_clutser_err.pdf',sep=''))
  
  kwd = kmeans(as.matrix(occ.wday.m[,2:8],ncol=7),12,iter.max=20)
  showClusters(kwd,occ.wday.m[,2:8],ncol=4,xlabs=c('Su','Mo','Tu','We','Th','Fr','Sa'),main='Members of HOD occupant event clusters',ylab='members',xlab='Day of week')
  dev.copy2pdf(file = paste(occupancyDir,'wday_clutser_members.pdf',sep=''))
  
  matplot(t(kwd$centers),type='l',main='Day of week K-means cluster centers',ylab='probability density',xlab='Day of week',xaxt='n')
  xlabs=c('Su','Mo','Tu','We','Th','Fr','Sa')
  axis(1, at=1:length(xlabs), labels=xlabs)
  dev.copy2pdf(file = paste(occupancyDir,'wday_clutser_centers.pdf',sep=''))
  
  km = kmeans(as.matrix(occ.hr.m[,2:25],ncol=24),12,iter.max=20)
  showClusters(km,occ.hr.m[,2:25],ncol=4,xlabs=c(0:23),main='Members of hour of day occupant event clusters',ylab='members',xlab='Hour of day')
  dev.copy2pdf(file = paste(occupancyDir,'clutser_members.pdf',sep=''))
  
  #op <- par(no.readonly = TRUE)
  #par( mfrow=c(2,1), oma=c(2,2,3,0),mar=c(2,2,2,2)) # Room for the title
  #matplot(t(km$centers),type='l',main='K=Means cluster centers (all days)',ylab='density',xlab='hour of day')  
  #kmwd = kmeans(as.matrix(occ.wkdy[,2:25],ncol=24),10,iter.max=20)
  #matplot(t(kmwd$centers),type='l',main='K=Means cluster centers (weekdays only)',ylab='density',xlab='hour of day')
  #par(op)
  #dev.copy2pdf(file = paste(occupancyDir,'clutser_centers_orig.pdf',sep=''))

  op <- par(no.readonly = TRUE)
  par( mfrow=c(2,1), oma=c(2,3,1,0),mar=c(2,1,2,1)) # Room for the title
  kmwd = kmeans(as.matrix(occ.wkdy.m[,2:25],ncol=24),10,iter.max=20)
  matplot(t(kmwd$centers),type='l',main='Weekday K-means cluster centers (by monthly percentile)',ylab='density',xlab='hour of day',xaxt='n')  
  axis(1, at=(1:length(xlabs)), labels=xlabs,cex=0.6)
  grid()
  kmnd = kmeans(as.matrix(occ.wknd.m[,2:25],ncol=24),10,iter.max=20)
  matplot(t(kmnd$centers),type='l',main='Weekend K-means cluster centers (by monthly percentile)',ylab='density',xlab='hour of day',xaxt='n')
  xlabs = 0:23
  axis(1, at=(1:length(xlabs)), labels=xlabs,cex=0.6)
  grid()
  mtext('Hour of day', line=0, font=2, cex=1, side=1, outer=TRUE)
  mtext('probability density', line=2, font=2, cex=1, side=2, outer=TRUE)
  par(op)
  
  dev.copy2pdf(file = paste(occupancyDir,'clutser_centers.pdf',sep=''),width=12,height=8)
}

