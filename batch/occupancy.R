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
    wss[i] <- (km$withinss/counts) / rowSums(km$centers)
  }
  print(wss/nrow(data))
  #print(wss/nrow(data))
  return(wss/nrow(data))
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
  outDir = 'occupancy/'
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
  par( mfrow=c(2,1), oma=c(2,2,3,0),mar=c(2,2,2,2)) # Room for the title
  plot(r$kw,ylab='kW',col='#444444')
  points(pieces$summaries[[26]],type='l', col='blue')
  plot(pieces$summaries[[9]],ylab='kW',xlab='observation #',col='#444444')
  points(qod,e[qod],col='blue',ylab='kW')
  points(mqod,e[mqod],col='red',ylab='kW')
  
  mtext('Model predictions, residuals, and outliers', line=0, font=2, cex=1.2, outer=TRUE)
  par(op)
  dev.copy2pdf(file = paste(outDir,'single_model.pdf',sep=''))
  
  dens = quantileDensities(e,r$dates)
  densM = quantileDensities(e,r$dates,monthly=T)
  plot(dens$hr[1:24],type='l',main='p(Occupancy Event)',ylab='density',xlab='hour of day')
  points(densM$hr[1:24],type='o',col='blue')
  points(meankw,col='red')
  points(meane,col='green')
  dev.copy2pdf(file = paste(outDir,'p_occ_event.pdf',sep=''))
    
  fits = clusterSumOfSquares(as.matrix(occ.hr[,2:25]),seq(2,60,4),20)
  vals = which(! is.na(fits))
  plot(vals,fits[vals],type='o',ylab='error',xlab='# of clusters',main='Fit error vs. number of clusters')
  dev.copy2pdf(file = paste(outDir,'clutser_err.pdf',sep=''))
  
  
  km = kmeans(as.matrix(occ.hr.m[,2:25],ncol=24),10,iter.max=20)
  op <- par(no.readonly = TRUE)
  par( mfrow=c(3,4), oma=c(2,2,3,0),mar=c(2,2,2,2)) # Room for the title
  for(c in sort(unique(km$cluster))){
    hmap(data=as.matrix(occ.hr[km$cluster == c,2:25]),main=paste(c,' (n=',sum(km$cluster==c),')',sep=''))
  }
  par(op)
  
  dev.copy2pdf(file = paste(outDir,'clutser_members.pdf',sep=''))
  
  op <- par(no.readonly = TRUE)
  par( mfrow=c(2,1), oma=c(2,2,3,0),mar=c(2,2,2,2)) # Room for the title
  matplot(t(km$centers),type='l',main='K=Means cluster centers (all days)',ylab='density',xlab='hour of day')  
  kmwd = kmeans(as.matrix(occ.wkdy[,2:25],ncol=24),10,iter.max=20)
  matplot(t(kmwd$centers),type='l',main='K=Means cluster centers (weekdays only)',ylab='density',xlab='hour of day')
  par(op)
  
  dev.copy2pdf(file = paste(outDir,'clutser_centers_orig.pdf',sep=''))

  op <- par(no.readonly = TRUE)
  par( mfrow=c(2,1), oma=c(2,2,3,0),mar=c(2,2,2,2)) # Room for the title
  kmwd = kmeans(as.matrix(occ.wkdy.m[,2:25],ncol=24),10,iter.max=20)
  matplot(t(kmwd$centers),type='l',main='K=Means cluster centers (annual percentile)',ylab='density',xlab='hour of day')  
  kmnd = kmeans(as.matrix(occ.wknd.m[,2:25],ncol=24),10,iter.max=20)
  matplot(t(kmnd$centers),type='l',main='K=Means cluster centers (monthly percentile)',ylab='density',xlab='hour of day')
  par(op)
  
  dev.copy2pdf(file = paste(outDir,'clutser_centers.pdf',sep=''))
}
