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
  outDates = c()
  #monthProfiles = c()
  if('POSIXct' %in% class(dates)) { dates = as.POSIXlt(dates) }
  for(m in sort(unique(dates$mon))) {
    mnth = which(dates$mon == m)
    threshold = quantile(e[mnth],q,na.rm=T)
    outliers = which(e[mnth] > threshold) # which strips NAs
    odts = as.POSIXct(dates[mnth][outliers])
    #monthProfiles = rbind(monthProfiles,densities(odts)$hr)
    outDates = c(outDates,odts)
  }
  #matplot(t(as.matrix(monthProfiles[,1:24])))
  return(as.POSIXlt(outDates,origin='1970-01-01'))
}

densities = function(outDates) {
  if('POSIXct' %in% class(outDates)) {
    outDates = as.POSIXlt(outDates)
  }
  wday = as.array(table(factor(outDates$wday,levels=0:6)))
  names(wday) = c('Su','Mo','Tu','We','Th','Fr','Sa')[as.numeric(names(wday))+1]
  sc = sum(wday)
  wday = c(wday/sc,scale=sc)
  
  hr= as.array(table(factor(outDates$hour,levels=0:23)))
  sc = sum(hr)
  hr = c(hr/sc,scale=sc)
  
  wknd= as.array(table(factor(outDates[outDates$wday %in% c(0,6)]$hour,levels=0:23)))
  wknd = c(wknd/sc*7/2,scale=sc)
  
  wkdy= as.array(table(factor(outDates[outDates$wday %in% c(1:5)]$hour,levels=0:23)))
  wkdy = c(wkdy/sc*7/5,scale=sc)
  
  
  mon = as.array(table(factor(outDates$mon,levels=0:11)))
  names(mon) = month.abb
  sc = sum(mon)
  mon = c(mon/sc,scale=sc)
  
  how = as.array(table(factor(outDates$wday*24+outDates$hour,levels=0:(24*7-1))))
  sc = sum(how)
  how = c(how/sc,scale=sc)
  
  return(list(wday=wday,hr=hr,wknd=wknd,wkdy=wkdy,mon=mon,how=how))
}

quantileDensities = function(e,dates,q=0.95,monthly=F) {
  qd=NA
  if(monthly) {
    qd = monthlyQuantileOutlierDates(e,dates,q)
  } else {
    qd = quantileOutlierDates(e,dates,q)
  }
  return(densities(qd))
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
  wss <- (nrow(data)-1)*mean(apply(data,2,var))
  #wss <- sum(apply(data,2,var))
  for (i in countSeq){
    km <- kmeans(data,centers=i,iter.max=iter.max)
    
    counts = apply(as.matrix(1:i),1,function(x) sum(km$cluster==x))
    wss[i] <- mean(km$withinss) #(km$withinss/counts) / rowSums(km$centers)
  }
  print(wss/nrow(data))
  #print(wss/nrow(data))
  return(wss/nrow(data))
}


showClusters = function(km,obs,ncol=4,xlabs=NULL,main=NULL,ylab=NULL,xlab=NULL,...) {
  op <- par(no.readonly = TRUE)
  par( mfrow=c(ceiling(nrow(kwd$centers)/ncol),ncol), oma=c(2,2,3,0),mar=c(2,1,2,1)) # Room for the title
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
  dg = DescriptorGenerator(name='toutPieces',genImpl=toutPieces24Generator,subset=list(all="TRUE"))
  
  r = ResDataClass(553991005,93304,useCache=T,doSG=F)
  #r = ResDataClass(554622151,93304,useCache=T,doSG=F)
  #r = ResDataClass(1064423310,93304,useCache=T,doSG=F)
  pieces = dg$run(r)
  e = pieces$summaries[[9]]
  
  qod  = match(quantileOutlierDates(e,r$dates),r$dates)
  mqod = match(monthlyQuantileOutlierDates(e,r$dates),r$dates)

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
  
  op <- par(no.readonly = TRUE)
  par( mfrow=c(2,1), oma=c(2,2,1,0),mar=c(2,4,2,2))
  plot(r$dates,r$kw,col='#444444',main='Observations (grey) & model (green)',ylab='kW')
  points(r$dates,pieces$summaries[,'prediction'][[1]], col='#22ff22',cex=0.4,pch=15)
  grid()
  plot(r$dates,pieces$summaries[,'residuals'][[1]],xlab='observation #',col='#444444',main='Residuals (grey) & monthly top 5% (red)',ylab='kW')
  #points(qod,e[qod],col='blue',ylab='kW',cex=0.4)
  points(r$dates[mqod],e[mqod],col='red',cex=0.4,pch=15)
  #mtext('Model predictions, residuals, and outliers', line=0, font=2, cex=1.4, outer=TRUE)
  grid()
  par(op)
  
  dev.copy2pdf(file = paste(occupancyDir,'single_model.pdf',sep=''),width=12,height=8)
  
  dens = quantileDensities(e,r$dates)
  densM = quantileDensities(e,r$dates,monthly=T)
  plot(densM$hr[1:24]*densM$hr[25],type='o',main='Count of occupancy events by hour',ylab='count',xlab='Hour of day',pch=16,col='red',xaxt='n',ylim=c(0,max(densM$hr[1:24]*densM$hr[25])))
  #points(densM$hr[1:24]*densM$hr[25],type='o',col='blue')
  points(meankw*dens$hr[25],col='#999999',type='l',lty=2)
  grid()
  legend(37,c('Count of top 5%','Rescaled mean load'),lty=c(1,2),pch=c(16,-1),col=c('red','#999999'))
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

