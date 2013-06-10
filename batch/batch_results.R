#.libPaths('~/R/library') # use my local R library even from the command line

require(ggplot2)
require(gtools)
require(gridExtra)
require(reshape2)
require(hexbin)
require(plyr)
require(scales) # for muted

# todo: is there a better way to detect the current directory?
conf.basePath = file.path('~/EnergyAnalytics/batch')
if(Sys.info()['sysname'] == 'Windows' & Sys.info()['user'] == 'Sam') {
  conf.basePath = file.path('f:/dev/pge_collab/EnergyAnalytics/batch')
}
setwd(conf.basePath)
source(file.path(getwd(),'resultAnalysis.R'))
source(file.path(getwd(),'dbUtil.R'))
source(file.path(getwd(),'DataClasses.R'))
source(file.path(getwd(),'weatherFeatures.R'))
source(file.path(getwd(),'zipMap.R'))

resultsDir = 'results_daily'       # daily models for all homes
#resultsDir = 'results_daily_flex'  # 2 change point models
resultsDir = 'results_daily_standard'
resultsDir = 'results_daily_standard2'
resultsDir = 'results_daily_nestedCP'
resultsDir = 'results_daily_full'
dirZips = do.call(rbind,strsplit(list.files(file.path(getwd(),resultsDir),pattern='modelResults.RData'),'_'))[,1]
allZips = dirZips

resultsDir = 'results_daily_DL'
dirZips = do.call(rbind,strsplit(list.files(file.path(getwd(),resultsDir),pattern='modelResults.RData'),'_'))[,1]
allZips = dirZips

resultsDir = 'results_basics'
dirZips = do.call(rbind,strsplit(list.files(file.path(getwd(),resultsDir),pattern='modelResults.RData'),'_'))[,1]
allZips = dirZips

#allZips = c(94923,94503,94574,94559,94028,94539,94564,94702,94704,94085,
#            95035,94041,95112,95113,95765,95648,95901,94531,94585,95205,
#            95202,93619,93614,93304,93701,95631,95726,95223,95666)
#summary = combineSummaries(allZips,resultType='d_summaries') # fails from too much data?!

#allZips = c(94923,94503,94574,94559)
#allZips = c(93304)

# get a list of all the data filees in the dir and extract their zips

resultBasicsFile  = file.path(getwd(),resultsDir,'resultBasics.RData')
resultScalarsFile = file.path(getwd(),resultsDir,'resultScalars.RData')
cpDataFile        = file.path(getwd(),resultsDir,'cpData.RData')
coeffDataFile     = file.path(getwd(),resultsDir,'coeffData.RData')
bestModelDataFile = file.path(getwd(),resultsDir,'bestModelData.RData')
invalidIdsFile    = file.path(getwd(),resultsDir,'invalidIds.RData')

#loadData = function(zip) {
#  dataFile = file.path(getwd(),resultsDir,paste(zip,'_modelResults.RData',sep=''))
#  load(dataFile)
#}

if(file.exists(invalidIdsFile)) {
  load(invalidIdsFile)
} else {
  invalids = combine(allZips, resultType='invalid.ids',
                    fun=function(x) { return(x$id) }, appendZipData=F)
  colnames(invalids) = c('id','zip')
  save(list=c('invalids'),file=invalidIdsFile)
}
gc()

if(file.exists(resultScalarsFile)) {
  load(resultScalarsFile)
} else {
  sclrs      = combine(allZips,resultType='d_summaries',fun=scalars,appendZipData=T)
  save(list=c('sclrs'),file=resultScalarsFile)
}
gc()

if(file.exists(resultBasicsFile)) {
  load(resultBasicsFile)
} else {
  basics     = combine(allZips,resultType='features.basic',fun=as.matrix,appendZipData=T)
  basicMeans = combine(allZips,resultType='features.basic',fun=function(x) { t(colMeans(x)) },appendZipData=T )
  save(list=c('basics','basicMeans'),file=resultBasicsFile)
}
gc()

if(file.exists(coeffDataFile)) {
  load(coeffDataFile)
} else {
  results.cfs    = list()
  results.stde   = list()
  results.stdeNW = list()
  results.pvs    = list()
  results.pvsNW  = list()
  results.tvs    = list()
  models = unique(sclrs$model.name)
  for(model in models) {
    
    print(paste(model,'from',paste(models,collapse=',')))
    results.cfs[[model]]  = combine(allZips,resultType='d_summaries',fun=cf,         model.name=model,appendZipData=F)
    results.stde[[model]] = combine(allZips,resultType='d_summaries',fun=stderrs,    model.name=model,appendZipData=F)
    results.pvs[[model]]  = combine(allZips,resultType='d_summaries',fun=pvals,      model.name=model,appendZipData=F)
    results.tvs[[model]]  = combine(allZips,resultType='d_summaries',fun=tvals,      model.name=model,appendZipData=F)
    results.stdeNW[[model]] = combine(allZips,resultType='d_summaries',fun=stderrsNW,model.name=model,appendZipData=F)
    results.pvsNW[[model]]  = combine(allZips,resultType='d_summaries',fun=pvalsNW,  model.name=model,appendZipData=F)
    
  }
  save(list=c('results.cfs','results.stde','results.pvs','results.tvs','results.pvsNW','results.stdeNW'),file=coeffDataFile)
  #save(list=c('results.cfs','results.stde','results.pvs','results.tvs'),file=coeffDataFile)
}
gc()


ws = getWeatherSummary()
if(file.exists(cpDataFile)) {
  load(cpDataFile)
} else {
  cp1Data = combine(allZips,
                   resultType='d_others',
                   subResultType='DOW_toutCP_DL_l1',
                   fun=function(x) {
                     cpMatrix = t(apply(x, 1, 
                     function(y) { 
                       out = c(id=y$id,t(y$data[rownames(y$data),]))
                       names(out) <- c('id',rownames(y$data));
                       return(out)
                       }));
                     return(data.frame(cpMatrix))
                   },
                   appendZipData=T)
  cp1Data = merge(ws,cp1Data,by.x='zip5',by.y='zip5')
  cp1Data = merge(basics[,c('id','nObs','kw.mean')],cp1Data,by.x='id',by.y='id')
  cp2Data = combine(allZips,
                    resultType='d_others',
                    subResultType='DOW_tout2CP_DL_l1',
                    fun=function(x) {
                      cpMatrix = t(apply(x, 1, 
                       function(y) { 
                         out = c(id=y$id,t(y$data[rownames(y$data),]))
                         names(out) <- c('id',rownames(y$data));
                         return(out)
                       }));
                      return(data.frame(cpMatrix))
                    },
                    appendZipData=T)
  cp2Data = merge(ws,cp2Data,by.x='zip5',by.y='zip5')
  cp2Data = merge(basics[,c('id','nObs','kw.mean')],cp2Data,by.x='id',by.y='id')
  save(list=c('cp1Data','cp2Data'),file=cpDataFile)
}
gc()

# find the best r.squared result for each sp_id
if(file.exists(bestModelDataFile)) {
  load(bestModelDataFile)
} else {
  bestModels = do.call(rbind,by(sclrs,sclrs$id,function(df) df[which.max(df$r.squared),]))
  save(list=c('bestModels'),file=bestModelDataFile)
}

sclrsPlus = merge(sclrs,basics[c('id','kw.var','kw.mean')],by.x='id', by.y='id')

ggplot(subset(sclrsPlus,
              model.name %in% 
                c('tout','DOW_tout','DOW_tout_DL','DOW_tout_DL_l1','DOW_tout.min_DL', 'DOW_tout.max_DL','DOW_DD_DL','DOW_toutCP_DL_l1DailyCP')),
              aes(x=adj.r.squared,color=model.name)) + geom_density() + 
        xlim(0,60) + labs(title='Adj R2 sequence for various thermal properties')

# note that sigma^2 = 1/(n-p) Sum(w[i] R[i]^2)
ggplot(subset(sclrsPlus,
              model.name %in% 
                c('tout','DOW_tout','DOW_tout_DL','DOW_tout_DL_l1','DOW_tout.min_DL', 'DOW_tout.max_DL','DOW_DD_DL','DOW_toutCP_DL_l1DailyCP')),
              aes(x=sigma/(kw.mean*24)*100,color=model.name)) + geom_density() + 
        xlim(0,60)  + labs(title='Coeff of Variation for various thermal properties')

# heating cumsum
ord = order(cp1Data$lower)
plot(-1 * cp1Data$lower[ord][cp1Data$lower[ord] < 0],pch=20,main='Magnitude of heating demand (kWh/HDD)',ylab='kWh/HDD',xlab='Count of residences',ylim=c(0,-1*min(cp1Data$lower)))
grid()

plot(cumsum(-1 * cp1Data$lower[ord][cp1Data$lower[ord] < 0]),pch=20,main='Magnitude of heating demand (kWh/HDD)',ylab='kWh/HDD',xlab='Count of residences')
grid()

# cooling cumsum
ord = rev(order(cp1Data$upper))
plot(cp1Data$upper[ord][cp1Data$upper[ord] > 0],pch=20,main='Magnitude of cooling demand (kWh/CDD)',ylab='kWh/CDD',xlab='Count of residences',ylim=c(0,max(cp1Data$upper)))
grid()

cooling = cp1Data$upper*pmax(0,90 - cp1Data$cp)
ord = rev(order(cooling))
plot(cumsum( cooling[ord][cooling[ord] > 0] )/1000,type='l',,lty=1,main='Cumlative cooling demand (kWh/day)',ylab='MWh/day',xlab='Count of residences')

ts = c(100,90,80,70,60,50)
for(i in 1:length(ts)) {
  t = ts[i]
  cooling = cp1Data$upper*pmax(0,t - cp1Data$cp)
  ord = rev(order(cooling))
  cum = cumsum( cooling[ord][cooling[ord] > 0] )/1000
  if(i==1) {
    plot(cum,type='l',lty=i,main='Cumlative cooling demand (kWh/day)',ylab='MWh/day',xlab='Count of residences')
  }
  else {
    points(cum,pch=20,type='l',lty=i)
  }
}
grid()
# add the legend
# xapx and yaxp return the lowest and highest tick mark for the plot
# usr returns the outer dimensions of the current plot as c(xmin, xmax, ymin, ymax)
legend(par('xaxp')[1],par('usr')[4]*0.95, paste(ts,'F'),lty=1:length(ts), title="Temperature",cex=0.8)

ts = c(30,40,50,60,70)
for(i in 1:length(ts)) {
  t = ts[i]
  heating = cp1Data$lower*pmin(0,t - cp1Data$cp) # negative values only multiply negative slopes to give kWh
  ord = rev(order(heating))
  if(i==1) {
    plot(cumsum( heating[ord][heating[ord] > 0] )/1000,type='l',lty=i,main='Cumlative heating demand (kWh/day)',ylab='MWh/day',xlab='Count of residences')
  }
  else {
    points(cumsum( heating[ord][heating[ord] > 0] )/1000,pch=20,type='l',lty=i)
  }
}
grid()
# add the legend
# xapx and yaxp return the lowest and highest tick mark for the plot
# usr returns the outer dimensions of the current plot as c(xmin, xmax, ymin, ymax)
legend(par('xaxp')[1],par('usr')[4]*0.95, paste(ts,'F'),lty=1:length(ts), title="Temperature",cex=0.8)

# kurtosis of change point model runs
ggplot(sclrs,aes(x=kurtosis,color=model.name)) + geom_density() + xlim(0,10) + labs(title='Kurtosis of residuals across model runs by model type',x='kurtosis')

# R2 of change point model runs
ggplot(sclrs,aes(x=r.squared,color=model.name)) + geom_density() + xlim(0,1) + labs(title='R2 for model runs by model type',x='R2')

# rmse of change point model runs
ggplot(sclrs,aes(x=cv.rmse,color=model.name)) + geom_density() + xlim(0,20) + labs(title='Cross validated RMSE for model runs by model type',x='RMSE (cross validated)')

# rmse normed by mean kw for model runs
ggplot(sclrsPlus,aes(x=cv.rmse/(kw.mean*24),color=model.name)) + 
  geom_density() + xlim(0,1) + 
  labs(title='Cross validated CV for model runs by model type',x='RMSE/mean (cross validated)')


ggplot(cp2Data,aes(x=cp2-cp1)) + geom_histogram(binwidth=1) + xlim(0,30) + labs(title='Distance between CPs', x='Degs F', y='count')

sclrs2 = subset(sclrs,model.name=='tout2CPDailyFlexCP') # 2 change points
sclrs1 = subset(sclrs,model.name=='tout1CPDailyCP')     # 1 change point
sclrs0 = subset(sclrs,model.name=='dailyTout')          # no change point 

trueCP1 = cp1Data$upper > 0 & cp1Data$lower <= 0 & cp1Data$pval.upper < 0.01
trueCP2 = cp2Data$upper > 0 & cp2Data$lower <= 0 & cp2Data$pval.upper < 0.01 & cp2Data$middle_1 > -0.6 & cp2Data$middle_1 < 0.6

trueCP1Data = cp1Data[trueCP1,]
trueCP2Data = cp2Data[trueCP2,]
# histograms of cp1 and cp2 overlaid together
ggplot(trueCP2Data,aes(x=cp1)) + 
  geom_histogram(binwidth=1,fill='#0000ff',alpha=0.5) + 
  geom_histogram(aes(x=cp2),binwidth=1,fill='#ff0000',alpha=0.5) +
  labs(title='Histogram of upper and lower change point',
       x='change points (lower and upper)')

ggplot(cp2Data,aes(x=cp1)) + 
  geom_histogram(binwidth=1,fill='#0000ff',alpha=0.5) + 
  geom_histogram(aes(x=cp2),binwidth=1,fill='#ff0000',alpha=0.5) +
  labs(title='Histogram of upper and lower change point',
       x='change points (lower and upper)')

# filter by "realistic" gap between CPs and "realistic" lower bound ot upper CP
cpCool = subset(cp2Data,subset=cp2-cp1 > 5 & cp2 > 60 & cp1 < 70)
ggplot(cpCool,aes(x=cp1)) + 
  geom_histogram(binwidth=1,fill='#0000ff',alpha=0.5) + 
  geom_histogram(aes(x=cp2),binwidth=1,fill='#ff0000',alpha=0.5) +
  labs(title='Histogram of "realistic" upper and lower change point',
       x='change points (lower and upper)')

# for daylight duration model, plot the impacts of a marginal hour of daylight
ggplot(results.cfs$tout_mean_WKND,aes(x=day.length*1000)) + 
    geom_histogram(binwidth=50) + 
    labs(title='Impact of hours of daylight on power consumption',
         x='delta W during daylight', y='density') + 
    scale_x_continuous(limits=c(-4000,3000),breaks=seq(-4000,3000,by = 500))
mean(results.cfs$tout_mean_WKND[,'day.length'])


calMap(db.getZipCounts(),'count','zip5',main='Meter count by zip code',colorMap=brewer.pal(9,"Blues") )
calMap(db.getZipData(),'cecclmzn','zip5',main='CEC climate zones',colorMap=brewer.pal(12,"Paired")[c(-1,-9,-11)] )
calMap(db.getZipData(),'climate' ,'zip5',main='PGE climate zones',colorMap=brewer.pal(12,"Paired")[c(-1,-9,-11)] )
ws = getWeatherSummary()
ws$rain = ws$rain * 365 * 24
ws$rain[ws$rain > 120] = NA # there is junk rain data (suprise!!)
calMap(ws,'tout','zip5',main='Mean annual temperature',colorMap=rev(brewer.pal(9,"RdBu")) )
calMap(ws,'rain','zip5',main='Total annual rain (inches)',colorMap=brewer.pal(9,"Blues") )
calMap(ws,'dp','zip5',main='Mean annual dew point (F)',colorMap=brewer.pal(9,"Purples") )

basics$kw.total = basics$kw.mean * 365 * 24
basicMeans$kw.total = basicMeans$kw.mean * 365 * 24

# maps of various metrics 
# see zipMap.R for the calMap implementation and other plots
calMap(basicMeans,'kw.total',main='Mean annual electricity consumption (kWh)',colorMap=brewer.pal(9,"Reds") )
calMap(basicMeans,'lag0',main='Simple correlation between Tout and kW denamd: corr(Tout,kW)',colorMap=rev(brewer.pal(9,"RdBu")) )
basicMeansWS = merge(basicMeans,ws,by.x='zip5',by.y='zip5')
basicMeansWS$rankDiff = rank(basicMeansWS$kw.total) - rank(basicMeansWS$tout)
calMap(basicMeansWS,'rankDiff',main='Rank difference between kW demand and temperature',colorMap=rev(brewer.pal(9,"RdBu")) )
ggplot(basicMeansWS,aes(x=tout,y=kw.total)) + geom_point() + labs(title='Annual energy vs. Tout',x='Annual average Tout',y='Annual energy (kWh)') + ylim(0,20000)

# histograms of annual energy demand
ggplot(basics,aes(x=kw.total,color=cecclmzn)) + geom_density() + xlim(0,40000) + labs(title='Prob. density of annual kWh by climate zone',x='Annual kWh',y='Count')

basics$minLevels = cut(basics$min * 1000,seq(0,1000,100))
ggplot(basics,aes(x=kw.total,fill=minLevels)) + geom_histogram(binwidth=100) + xlim(0,40000) + labs(title='Count of homes by annual kWh',x='Annual kWh',y='Count')
ggplot(basics,aes(x=kw.total,color=minLevels)) + geom_density() + xlim(0,40000) + labs(title='Count of homes by annual kWh',x='Annual kWh',y='Count')
ggplot(basics,aes(x=min,y=kw.total)) + geom_point() + ylim(0,40000) + labs(title='Constant loads vs. total annual energy',x='Average daily min demand (kW)',y='annual energy (kWh)')
baseFraction = basics$min*365*24 / basics$kw.total
plot(sort(baseFraction),type='l',ylab='Fraction of energy',main='Fraction of annual kWh from 24x7 loads')


basics$maxLevels = cut(basics$max * 1000,seq(0,4000,200))
ggplot(basics,aes(x=kw.total,fill=maxLevels)) + geom_histogram(binwidth=100) + xlim(0,40000) + labs(title='Count of homes by annual kWh',x='Annual kWh',y='Count')
ggplot(basics,aes(x=kw.total,color=maxLevels)) + geom_density() + xlim(0,40000) + labs(title='Count of homes by annual kWh',x='Annual kWh',y='Count')


# income vs energy use
ggplot(basicMeans,aes(x=median_income,y=kw.total)) + geom_point(colour="black") + ylim(0,20000) + labs(x='',y='')
ggplot(basicMeans,aes(x=median_income,y=kw.total,color=cecclmzn)) + geom_point() + ylim(0,20000) + labs(title='Zip code median income vs. mean annual energy (kWh)',y='Annual kWh',x='Median income')


MAm = as.matrix(basics[,grep("^ma[0-9]",colnames(basics))])
basics$maxMA = as.numeric(apply(MAm,1,which.max))

# example daily scatter plot types
toutScatter = function(r,main='',type='daily',useMean=F,xlim=c(40,95),ylim=NULL) {
  xLabel = 'mean T (degs F)'
  yLabel = 'kW (hourly average)'
  x = r$tout
  y = r$kw
  if(type=='daily') {
    dfa = rDFA(r)
    xLabel = 'mean T (degs F)'
    
    x = dfa$tout.mean
    if(useMean) {
      y = dfa$kw.mean
      yLabel = 'kW (daily average)'
    } else {
      y = dfa$kwh
      yLabel = 'kWh/day'
    }
  }
  if(type=='monthly') {
    dfa = regressorDFAggregated(r)$monthly
    xLabel = 'mean T (degs F)'
    yLabel = 'kW (monthly average)'
    x = dfa$tout.mean
    if(useMean) {
      y = dfa$kw.mean
    } else {
      y = dfa$kwh
    }
  }
  
  plot(x,y,
       xlab=xLabel,ylab=yLabel,main=main,xlim=xlim,ylim=ylim,
       mgp=c(1,0,0),tcl=0.5)
}

op <- par(no.readonly = TRUE)
par( mfrow=c(3,3), oma=c(2,0,3,0),mar=c(2,2,2,2))# Room for the title
toutScatter(ResDataClass(1064429605,93304),'1') # _/  good fit
toutScatter(ResDataClass(1418575405,93304),'2') # \/  good fit
toutScatter(ResDataClass(6481381805,93304),'3') # \_/ good fit
toutScatter(ResDataClass(6502182810,93304),'4') # \_/ good fit
toutScatter(ResDataClass(6533703010,93304),'5') # WTF bad fit
toutScatter(ResDataClass(3044473610,94704),'6') # \_  ok fit
toutScatter(ResDataClass(3074319310,94704),'7') # \_  ok fit
toutScatter(ResDataClass(3064027805,94704),'8') # __  bad fit
toutScatter(ResDataClass(3074465310,94704),'9') # __  ok fit
mtext("Typical daily kWh vs. mean Tout scatters", line=0, font=2, cex=1.2,outer=TRUE)
par(op)


op <- par(no.readonly = TRUE)
par( mfrow=c(1,2), oma=c(2,0,3,0),mar=c(2,2,2,2))# Room for the title
r = ResDataClass(1064429605,93304)
toutScatter(r,'Hourly averages',type='hourly',  useMean=T,xlim=c(30,105),ylim=c(0,8)) # _/ hourly
toutScatter(r,'Daily averages',type='daily',    useMean=T,xlim=c(30,105),ylim=c(0,8)) # _/ daily
#toutScatter(r,'Monthly averages',type='monthly',useMean=T,xlim=c(30,105),ylim=c(0,8)) # _/ monthly
mtext("Impact of time averaging on electricity consumption data", line=0, font=2, cex=1.2,outer=TRUE)
par(op)

r = ResDataClass(1064429605,93304);  plot(r,estimates=toutDoubleChangePoint(rDFA(r))) # _/  good fit
r = ResDataClass(1418575405,93304);  plot(r,estimates=toutDoubleChangePoint(rDFA(r))) # \/  good fit
r = ResDataClass(6481381805,93304);  plot(r,estimates=toutDoubleChangePoint(rDFA(r))) # \_/ good fit
r = ResDataClass(6502182810,93304);  plot(r,estimates=toutDoubleChangePoint(rDFA(r))) # \_/ good fit
r = ResDataClass(6533703010,93304);  plot(r,estimates=toutDoubleChangePoint(rDFA(r))) # WTF bad fit
r = ResDataClass(3044473610,94704);  plot(r,estimates=toutDoubleChangePoint(rDFA(r))) # \_  ok fit
r = ResDataClass(3074319310,94704);  plot(r,estimates=toutDoubleChangePoint(rDFA(r))) # \_  ok fit
r = ResDataClass(3064027805,94704);  plot(r,estimates=toutDoubleChangePoint(rDFA(r))) # __  bad fit
r = ResDataClass(3074465310,94704);  plot(r,estimates=toutDoubleChangePoint(rDFA(r))) # __  ok fit


ggplot(aes(x=cp,color=climate),data=cp1Data) + geom_density()
ggplot(aes(x=cp,color=cecclmzn),data=cp1Data) + geom_density() +
  labs(title='Change point distribution by Climate',x='change point')

g = ggplot(trueCP1Data,aes(x=cp,y=tout))
g = ggplot(cp1Data,aes(x=cp,y=tout))
g + geom_point(aes(color=cecclmzn))
g + geom_hex()
g + stat_density2d(geom="tile", aes(fill = ..density..), contour = FALSE) + 
  ylim(50,70) + 
  scale_fill_gradient(limits=c(1e-5,0.02)) + 
  labs(title='Density of change point vs annual average Tout',y='tout',x='change point')

g + stat_density2d(aes(fill = ..level..), geom="polygon",contour=T) +
  labs(title='Density of change point vs annual average Tout',y='tout',x='change point')


g = ggplot(cp1Data,aes(x=cp,y=cut(tout, breaks=c(seq(40,80,1),Inf))))
g + stat_bin(aes(fill=..density..), geom="tile", binwidth=3, position="identity")

g + stat_bin(aes(fill=..count..), geom="tile", binwidth=3, position="identity") + 
  labs(title='Annual mean temperature vs. change point',
           x='Change point (F)',
           y='tout mean (F)')

# box plot for change points as a function of summer mean temp
g = ggplot(trueCP1Data,aes(y=cp,x=cut(summer.tout, breaks=c(seq(40,80,3),Inf))))
g + geom_boxplot() + ylim(10,90) +
  labs(title='Change points as a function of mean summer temperatures',
       x='Summer temperature bin (deg F)',y='Change Point') +
  coord_flip() 

g = ggplot(trueCP2Data,aes(y=cp2,x=cut(summer.tout, breaks=c(seq(40,80,3),Inf))))
g + geom_boxplot() + ylim(10,90) +
  labs(title='Upper change points as a function of mean summer temperatures',
       x='Summer temperature bin (deg F)',y='Change Point') +
  coord_flip() 

g = ggplot(cp2Data,aes(y=cp2,x=cut(winter.tout, breaks=c(seq(10,80,3),Inf))))
g + geom_boxplot() + ylim(10,90) +
  labs(title='Lower change points as a function of mean winter temperatures',
       x='Winter temperature bin (deg F)',y='Change Point') +
  coord_flip() 

g = ggplot(cp2Data,aes(y=cp2,x=cecclmzn))
g + geom_boxplot() + ylim(10,90) +
  labs(title='Upper change points as a function of climate zone',
       x='Climate zone',y='Upper change point') +
  coord_flip() 

# histograms for each temperature range
g = ggplot(cp1Data,aes(x=cp,color=cut(summer.tout, breaks=c(seq(55,90,2),Inf)))) + geom_density() + xlim(55,90)
g

ggplot(aes(x=model.name),data=bestModels) + geom_bar() # counts of best models.

ab = merge(basics[c('id','kw.var')],bestModels,by.x='id',by.y='id',all=F) # add the model variance so sigma can be normalized
ggplot(aes(y=sigma,x=r.squared),data=bestModels) + geom_point() + ylim(0,50) # relationship between our two metrics

ggplot(aes(y=sigma/kw.var^2,x=r.squared),data=ab) + geom_point() + ylim(0,1000) # relationship between our two metrics

bestCP = merge(bestModels,cp1Data,by.x='id',by.y='id')
ggplot(bestCP,aes(y=cut(cp,seq(30,90,5)),x=pval.upper)) + geom_point()
sbcp = subset(bestCP,subset=bestCP$pval.upper < 0.1)
ggplot(sbcp,aes(x=cp,color=model.name)) + geom_density()

hists(sclrs,metric='r.squared',zip='all data (200k homes)')
zipHists(sclrs,metric='r.squared',model.name='toutDailyCP')
hists(sclrs,metric='sigma',xlim=c(0,20),zip='all data (200k homes)')
zipHists(sclrs,metric='sigma',model.name='toutDailyCP',xlim=c(0,20))

zip='94923'
#d_summary = modelResults$d_summary #combineSummaries(zips,'d_summaries')
load(file.path(getwd(),resultsDir,paste(zip,'_modelResults.RData',sep='')))
summary   = clean(modelResults$d_summaries)
cp = modelResults$changePoints
cpgrid = t(apply(cp,1,function(X) X[['changePoints']]['cp',]))
changes = ksc(cpgrid,8)
plotClusterCenters(changes)
showClusters(changes)
clusterHist(changes)
p = clusterHeatMaps(cpgrid,changes)
do.call(grid.arrange,p)

hists(summary,metric='sigma',zip=zips[2])
hists(summary,metric='r.squared',zip=zips[2])
hists(summary,metric='adj.r.squared',zip=zips[2])

hists(d_summary,metric='sigma',zip=zips[1])
hists(d_summary,metric='r.squared',zip=zips[1])
hists(d_summary,metric='adj.r.squared',zip=zips[1])

estimates = cf(summary,'toutPIECES24','all')
pvals     = cf(summary,'toutPIECES24','all',col='Pr(>|t|)')

slopeCols = grep('tout75_Inf',colnames(estimates),value=T)[-c(1:8)]
slope75   = estimates[,c('id',slopeCols)]
slope75p  = pvals[,c('id',slopeCols)]
sig       = slope75p[,-1] >= 0.0
slope75[cbind(rep(FALSE,dim(sig)[1]),!sig)] <- NA
slope75m = melt(slope75,id.vars='id')
ggplot(data=slope75m,aes(x=variable, y=value)) + scale_x_discrete(labels=c('12am',1:11,'12pm',1:11)) + geom_line(aes(group=id),alpha = 0.05)

a = as.matrix(slope75[,-1])



hodCols   = grep('^HODH[0-9]+$',colnames(estimates),value=T)
hodConst  = estimates[,c('id',hodCols)]
hodConstp = pvals[,c('id',hodCols)]
sig       = hodConstp[,-1] > 0.05
#hodConst[cbind(rep(FALSE,dim(sig)[1]),sig)] <- NA
hodConstm = melt(hodConst,id.vars='id')
ggplot(data=hodConstm,aes(x=variable, y=value)) + scale_x_discrete(labels=c('12am',1:11,'12pm',1:11)) + geom_line(aes(group=id))

# create estimates for temperatures 0:100 by doing some linear algbra
pmCols    = grep('HODH15',colnames(estimates),value=T)
pmSlopes  = estimates[,c('id',pmCols)]
#pmConst   = estimates[,c('id',pmCols[1])]
pmSlopesp = pvals[,c('id',pmCols)]
X = regressor.piecewise(0:100,c(55,65,75)) # piecewise tout values from 0:100
X = cbind(1,X)                             # constant term
# one col of values per id
est = apply(as.matrix(pmSlopes[,-1]),MARGIN=1,FUN=function(row) t(X %*% as.matrix(row)))
colnames(est) <- paste('id_',pmSlopes[,1],sep='') # restore id names
estm = melt(est)
ggplot(data=estm,aes(x=X1, y=value,color=X2)) + geom_line(aes(group=X2))


zip = 94610
load(file.path(getwd(),resultsDir,paste(zip,'_modelResults.RData',sep='')))
summary   = clean(modelResults$summaries)
estimates.s = cf(summary,'toutTOD_WKND','summer')
estimates.w = cf(summary,'toutTOD_WKND','winter')
pvals.w     = cf(summary,'toutTOD_WKND','winter',col='Pr(>|t|)')
pvals.s     = cf(summary,'toutTOD_WKND','summer',col='Pr(>|t|)')

quad_chart(estimates.s,title=paste(zip,model,' (coast/moderate income', 'summer', 'only)'))


plot(colMeans(pvals.s[,grep("^HODWKWK", colnames(pvals.s), value=TRUE)]),main=paste(zip,' p-values'))

g1 = hists(scalars(summary,subset.name='summer'),'sigma',zip)
  
zip = 95223
load(file.path(getwd(),resultsDir,paste(zip,'_modelResults.RData',sep='')))
summary   = clean(modelResults$summaries)
estimates.s = cf(summary,model,'summer')
estimates.w = cf(summary,model,'winter')
pvals.w     = cf(summary,model,'winter',col='Pr(>|t|)')
pvals.s     = cf(summary,model,'summer',col='Pr(>|t|)')

# eliminate values not different from zero
pmax = 0.4 # maximum acceptable p-value in (0,1), with 0.05 being 95% confidence
estimates.s[cbind(FALSE,pvals.s[,-1] > pmax)] = NA
estimates.w[cbind(FALSE,pvals.w[,-1] > pmax)] = NA

# demean each set of coefficients
estimates.s = demean(estimates.s,"^tout65.HODWKWK")
estimates.w = demean(estimates.w,"^tout65.HODWKWK")
estimates.s = demean(estimates.s,"^tout65.HODWKND")
estimates.w = demean(estimates.w,"^tout65.HODWKND")
estimates.s = demean(estimates.s,"^HODWKWK")
estimates.w = demean(estimates.w,"^HODWKWK")
estimates.s = demean(estimates.s,"^HODWKND")
estimates.w = demean(estimates.w,"^HODWKND")

c_summary(estimates.s,pvals.s,0.05,pattern=c('^HODWKWK','^tout65.HODWKWK'))

p_summary(pvals.s,0.05,pattern=c('^HODWKWK','^tout65.HODWKWK'))


d_summary = clean(modelResults$d_summaries)

g1 = hists(scalars(d_summary,subset.name='all'),'r.squared',zip)
g2 = hists(scalars(d_summary,subset.name='all'),'sigma',zip)
grid.arrange(g1,g2)

# eliminate values not different from zero
#pmax = 0.05 # maximum acceptable p-value in (0,1), with 0.05 being 95% confidence
#estimates.s[cbind(FALSE,pvals.s[,-1] > pmax)] = NA
#estimates.w[cbind(FALSE,pvals.w[,-1] > pmax)] = NA

quad_chart(estimates.w,title=paste(zip,model,'(inland/lower income', 'winter', 'only)'))


par(mfrow=c(1,2))
plot(colMeans(pvals.w[,grep("^HODWKWK", colnames(pvals), value=TRUE)]),
     col='blue',
     ylim=c(0,1),
     main=paste(zip,'HOD p-values'),
     xlab='Hour of day',
     ylab='p-value (significant at < 0.05)')
points(colMeans(pvals.s[,grep("^HODWKWK", colnames(pvals), value=TRUE)]) / 2,pch=2,col='red')

plot(colMeans(pvals.w[,grep("^tout65.HODWKWK", colnames(pvals), value=TRUE)]),
     col='blue',
     ylim=c(0,1),
     main=paste(zip,'tout65xHOD p-values'),
     xlab='Hour of day',
     ylab='p-value (significant at < 0.05)')
points(colMeans(pvals.s[,grep("^tout65.HODWKWK", colnames(pvals), value=TRUE)]) / 2,pch=2,col='red')


g2 = hists(scalars(summary,subset.name='summer'),'sigma',zip)
grid.arrange(g2)
  
g2 = hists(scalars(summary,subset.name='all'),'r.squared',zip) + xlab('kWh/day')
grid.arrange(g2)


b = subset(modelResults$summaries,model.name == 'toutTOD' & subset.name == 'all',select=c(id,coefficients))
c = adply(b,.margins=1,.fun = function(row) {row$coefficients[[1]][1:10,'Estimate']})
c$coefficients = c()
c = c[, c('id',grep("^tout", colnames(c), value=TRUE))]
q = quantile(c$tout70_80,c(0.9))
c = subset(c,tout70_80 > q)
d = melt(c,id.vars=c('id'))
d$id = as.numeric(d$id) # somehow d$ids is a list of single entry lists?!
ggplot(d, aes(variable, value, group = id)) + geom_path(alpha = 0.1) + ylim(-0.25,0.5)

b = subset(modelResults$summaries,model.name == 'toutTOD' & subset.name == 'all',select=c(id,coefficients))
cnames = grep('^HOD',rownames(b[[1,'coefficients']]),value=TRUE)
c = adply(b,.margins=1,.fun = function(row) {row$coefficients[[1]][cnames,'Estimate']})
c$coefficients = c()
d = melt(c,id.vars=c('id'))
d$id = as.numeric(d$id) # somehow d$ids is a list of single entry lists?!
ggplot(d, aes(variable, value, group = id)) + geom_path(alpha = 0.05) #+ ylim(-0.1,0.15)

b = subset(modelResults$summaries,subset.name == 'all',select=c(id,model.name,rmse))
g = ggplot(b, aes(sigma, color = model.name)) + geom_density(size=1, alpha=0.2) + xlim(0,2.1)

a = subset(modelResults$features.basic,select=c(id,mean,min))
a$kWh = a$mean * 24 * 365
a$base = a$min * 24 * 365
b = subset(modelResults$summaries,model.name == 'toutPIECES' & subset.name == 'all',select=c(id,coefficients))
c = adply(b,.margins=1,.fun = function(row) {row$coefficients[[1]][,'Estimate']})
c$coefficients = c()
c = c[, c('id',grep("^tout", colnames(c), value=TRUE))]
q = quantile(c$tout70_80,c(0.75))
c$lower <- c$tout70_80 < q
c$quartile = paste('q',ceiling(rank(c$tout70_80) / (max(rank(c$tout70_80))/4)),seq='')
d = merge(a,c,by.y=c('id')) # get the kWh for the high Tout slope residence
d$rank = rank(d$kWh)
g = ggplot(a, aes(kWh)) + geom_density(size=1, alpha=0.2)
g = ggplot(d, aes(kWh,color=lower,fill=lower)) + geom_density(size=1, alpha=0.2) + geom_density(aes(kWh), size=1, alpha=0.2, color='black')
g = ggplot(d, aes(kWh,color=lower,fill=lower)) + geom_histogram(binwidth=500, alpha=.5, position="stack")
g = ggplot(d, aes(kWh,color=quartile,fill=quartile)) + geom_histogram(binwidth=500, alpha=.5, position="stack")
g = ggplot(d, aes(y=kWh,x=rank)) + geom_line()

resc[['z_94610']]$hourly$rmse[c('id','toutPIECES')]
gt = resc[['z_94610']]$hourly$rmse[c('id','toutPIECES')] < resc[['z_94610']]$hourly$rmse[c('id','standard')]
improved = gt[,2]     

                                                                                                    

