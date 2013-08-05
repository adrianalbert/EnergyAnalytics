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
source(file.path(getwd(),'census.R'))

lm_eqn = function(m) {
  
  l <- list(a = format(coef(m)[1], digits = 2),
            b = format(abs(coef(m)[2]), digits = 2),
            r2 = format(summary(m)$r.squared, digits=2));
  
  if (coef(m)[2] >= 0)  {
    eq <- substitute(italic(kWh) == a + b %.% italic(Tout)*","~~italic(r)^2~"="~r2,l)
  } else {
    eq <- substitute(italic(kWh) == a - b %.% italic(Tout)*","~~italic(r)^2~"="~r2,l)    
  }
  
  as.character(as.expression(eq));                 
}

resultsDir = 'results_daily'       # daily models for all homes
#resultsDir = 'results_daily_flex'  # 2 change point models
resultsDir = 'results_daily_standard'
resultsDir = 'results_daily_standard2'
resultsDir = 'results_daily_nestedCP'
resultsDir = 'results_daily_full'
#resultsDir = 'results_daily_residuals'
# load the results data for hte residuals
#allRes = t(sapply(s$residuals,
#                FUN=function(res) return(c(res,rep(NA,1200-length(res)))),simplify=T))
#allIds = matrix(sapply(s$id, FUN=function(id) return(c(id)),simplify=T),ncol=1)
#allRes = cbind(allIds,allRes)
#write.csv(allRes,file='93304_residuals.csv')

dirZips = do.call(rbind,strsplit(list.files(file.path(getwd(),resultsDir),pattern='modelResults.RData'),'_'))[,1]
allZips = dirZips

resultsDir = 'results_daily_DL'
dirZips = do.call(rbind,strsplit(list.files(file.path(getwd(),resultsDir),pattern='modelResults.RData'),'_'))[,1]
allZips = dirZips


resultsDir = 'w_results_basics'
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


if(file.exists(resultBasicsFile)) {
  load(resultBasicsFile)
} else {
  basics     = combine(allZips,resultType='features.basic',fun=as.matrix,appendZipData=T)
  basics$idZip = paste(basics$id,basics$zip5,sep='.')
  
  basicMeans = combine(allZips,resultType='features.basic',fun=function(x) { t(colMeans(x)) },appendZipData=T )
  
  basics$kw.total = basics$kw.mean * 365 * 24
  basicMeans$kw.total = basicMeans$kw.mean * 365 * 24
  
  rm(MAm)
  save(list=c('basics','basicMeans'),file=resultBasicsFile)
}
gc()

ws = getWeatherSummary()
#wbasics = basics
#wbasicMeans = basicMeans
#wbasicMeansWS = merge(wbasicMeans,ws,by.x='zip5',by.y='zip5')
#wbasicsWS = merge(wbasics,ws,by.x='zip5',by.y='zip5')
#wbasicMeansWS$rankDiff = rank(wbasicMeansWS$kw.total) - rank(wbasicMeansWS$tout)


basicMeansWS = merge(basicMeans,ws,by.x='zip5',by.y='zip5')
basicsWS = merge(basics,ws,by.x='zip5',by.y='zip5')
basicMeansWS$rankDiff = rank(basicMeansWS$kw.total) - rank(basicMeansWS$tout)


if(file.exists(resultScalarsFile)) {
  load(resultScalarsFile)
} else {
  sclrs       = combine(allZips,resultType='d_summaries',fun=scalars,appendZipData=T)
  sclrs$idZip = paste(sclrs$id,sclrs$zip5,sep='.')
  save(list=c('sclrs'),file=resultScalarsFile)
}
gc()

sclrsBasicsWS   = merge(sclrs,basicsWS[c('idZip','kw.var','kw.mean','maxMA')],by.x='idZip', by.y='idZip',all.x=T,all.y=F)
gc()


if(file.exists(coeffDataFile)) {
  load(coeffDataFile)
} else {
  results.cfs    = list()
  results.stde   = list()
  results.stdeNW = list()
  results.pvs    = list()
  results.pvsNW  = list()
  #results.tvs    = list()
  models = unique(sclrs$model.name)
  for(model in models) {
    
    print(paste(model,'from',paste(models,collapse=',')))
    results.cfs[[model]]    = combine(allZips,resultType='d_summaries',fun=cf,         model.name=model,appendZipData=F)
    #results.stde[[model]]   = combine(allZips,resultType='d_summaries',fun=stderrs,    model.name=model,appendZipData=F)
    results.pvs[[model]]    = combine(allZips,resultType='d_summaries',fun=pvals,      model.name=model,appendZipData=F)
    #results.tvs[[model]]    = combine(allZips,resultType='d_summaries',fun=tvals,      model.name=model,appendZipData=F)
    #results.stdeNW[[model]] = combine(allZips,resultType='d_summaries',fun=stderrsNW,model.name=model,appendZipData=F)
    results.pvsNW[[model]]  = combine(allZips,resultType='d_summaries',fun=pvalsNW,  model.name=model,appendZipData=F)
    
  }
  save(list=c('results.cfs','results.stde','results.pvs','results.pvsNW','results.stdeNW'),file=coeffDataFile)
  #save(list=c('results.cfs','results.stde','results.pvs','results.tvs'),file=coeffDataFile)
}
gc()


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


metricArray = function(metric='r.squared') {
  mNames = unique(sclrs$model.name)
  metrics = data.frame(idZip=unique(sclrs$idZip))
  for(mName in mNames) {
    print(mName)  
    res = subset(sclrs[,c('idZip',metric)],sclrs$model.name==mName)
    colnames(res) <- c('idZip',paste(mName,'.',metric,sep=''))
    #print(dim(res))
    #print(length(unique(res$idZip)))
    #print(dim(metrics))
    #print(table(res$idZip)[table(res$idZip) > 1])
    metrics = merge(metrics,res,by.x='idZip',by.y='idZip',all.x=T,all.y=F)
    gc()
  }
  return(metrics)
}

# find the best r.squared result for each sp_id
if(file.exists(bestModelDataFile)) {
  load(bestModelDataFile)
} else {
  #bestModels = do.call(rbind,by(sclrs,sclrs$id,function(df) df[which.max(df$r.squared),]))
  metrics.cv.rmse = metricArray('cv.rmse')
  metrics.r.squared = metricArray('r.squared')
  metrics.adj.r.squared = metricArray('adj.r.squared')
  metrics.AIC = metricArray('AIC')
  metrics.pacf.1 = metricArray('pacf.1')
  metrics.DW = metricArray('DW')
  
  mNames = factor(unique(sclrs$model.name))
  best = data.frame(idZip=metrics.cv.rmse$idZip)
  best$cv.rmse =       mNames[apply(metrics.cv.rmse[,-1],1,which.min)]
  best$r.squared =     mNames[apply(metrics.r.squared[,-1],1,which.max)]
  best$adj.r.squared = mNames[apply(metrics.adj.r.squared[,-1],1,which.max)]
  best$AIC =           mNames[apply(metrics.AIC[,-1],1,which.min)]
  best$pacf.1 =        mNames[apply(metrics.pacf.1[,-1],1,which.min)]
  best$DW =            mNames[apply(2 - abs(metrics.DW[,-1]),1,which.min)]
  
  save(list=c('best','metrics.cv.rmse','metrics.r.squared',
              'metrics.adj.r.squared','metrics.AIC',
              'metrics.pacf.1','metrics.DW'),file=bestModelDataFile)
}




bestm = melt(best,id.vars='idZip',variable.name='metric',value.name='model_name')
ggplot(bestm,aes(x=metric,color=model_name,fill=model_name)) + geom_bar(position='dodge') + 
  labs(title='Count of best fit model by metric of assessment',x='metric',y='count')
rm(bestm)

png(file.path(getwd(),'figures','density_ar2.png'),width=800,height=600)
ggplot(subset(sclrsBasicsWS,
              model.name %in% 
                c('tout','DOW_tout','DOW_tout_DL','DOW_tout_DL_l1',
                  'DOW_tout.min_DL','DOW_tout.max_DL','DOW_DD_DL',
                  'DOW_toutCP_DL_l1DailyCP')),
              aes(x=adj.r.squared,color=model.name)) + geom_density() + 
              xlim(0,1.0) + 
              labs(title='Adj R2 sequence for various thermal properties')
dev.off()



png(file.path(getwd(),'figures','density_ar2.png'),width=800,height=600)
ggplot(sclrs,
  aes(x=adj.r.squared,color=model.name)) + geom_density() + 
  ylim(0,5) + xlim(0,1.0) + 
  labs(title='Adj R2 sequence for various thermal properties')
dev.off()


# note that sigma^2 = 1/(n-p) Sum(w[i] R[i]^2)
ggplot(subset(sclrsBasicsWS,
              model.name %in% 
                c('tout','DOW_tout','DOW_tout_DL','DOW_tout_DL_l1',
                  'DOW_tout.min_DL', 'DOW_tout.max_DL','DOW_DD_DL',
                  'DOW_toutCP_DL_l1DailyCP','DOW_tout2CP_DL_l1DailyFlexCP')),
              aes(x=cv.rmse/(kw.mean*24)*100,color=model.name)) + geom_density() + 
              xlim(0,60)  + 
              labs(title='Coeff of Variation for various thermal properties')

ggplot(subset(sclrsBasicsWS,
              model.name %in% 
                c('DOW_tout_DL_65','DOW_tout_DL_CP65','DOW_toutNP_DL_l1DailyCP','DOW_toutCP_DLDailyCP',
                  'DOW_toutCP_DL_l1DailyCP','DOW_tout2CP_DL_l1DailyFlexCP')),
              aes(x=adj.r.squared,color=model.name)) + geom_density() + 
              xlim(0,1)  + 
              labs(title='Adjusted R2 for change point models')

ggplot(subset(sclrsBasicsWS,
              model.name %in% c('DOW_toutCP_DLDailyCP')),
       aes(x=adj.r.squared,color=cecclmzn)) + geom_density() + 
       xlim(0,1)  + 
       labs(title='Adjusted R2 for single change point model across CZs')

ggplot(subset(sclrsBasicsWS,
              zip5 %in% c('93304')),
       aes(x=adj.r.squared,color=model.name)) + geom_density() + 
  xlim(0,1)  + ylim(0,4) +
  labs(title='Adjusted R2 for Bakersfield zip')

# kurtosis of change point model runs
ggplot(subset(sclrs,
              model.name %in% 
                c('tout','DOW_tout','DOW_tout_DL','DOW_tout_DL_l1',
                  'DOW_tout.min_DL', 'DOW_tout.max_DL','DOW_DD_DL',
                  'DOW_toutCP_DL_l1DailyCP','DOW_tout2CP_DL_l1DailyFlexCP')),
       aes(x=kurtosis,color=model.name)) + geom_density() + 
  xlim(0,10) + 
  labs(title='Kurtosis of residuals across model runs by model type',
       x='kurtosis')

ggplot(subset(sclrs,
              model.name %in% 
                c('DOW_tout_DL','DOW_tout_DL_l1',
                  'DOW_toutCP_DL_l1DailyCP','DOW_tout2CP_DL_l1DailyFlexCP')),
       aes(x=pacf.1,color=model.name)) + geom_density() + 
  xlim(-0.4,1) + 
  labs(title='Partial auto correlation across model runs by model type',
       x='pacf(1)')

ggplot(subset(sclrs,
              model.name %in% 
                c('DOW_toutCP_DL_l1DailyCP')),aes(y=r.squared,x=tout)) + 
  geom_point(alpha = 0.03) + labs(title='R2 vs. Tout',x='Tout',y='r.squared')


ggplot(subset(sclrsBasicsWS,
              model.name %in% 
                c('DOW_toutCP_DL_l1DailyCP')),aes(y=r.squared,x=kw.mean)) + 
  geom_point(alpha = 0.1) + labs(title='R2 vs. kwh',x='kwh',y='r.squared')

ggplot(subset(sclrs,
              model.name %in% 
                c('DOW_toutCP_DL_l1DailyCP')),aes(y=pacf.1,x=tout)) + 
  geom_point(alpha = 0.03) + labs(title='pacf.1 vs. Tout',x='Tout',y='pacf.1')

ggplot(subset(sclrs,
              model.name %in% 
                c('DOW_toutCP_DL_l1DailyCP')),aes(y=pacf.1,x=r.squared)) + 
  geom_point(alpha = 0.05) + labs(title='pacf.1 vs. r2',x='R2',y='pacf.1')

ggplot(subset(sclrs,
              model.name %in% 
                c('tout','DOW_tout','DOW_tout_DL','DOW_tout_DL_l1',
                  'DOW_tout.min_DL', 'DOW_tout.max_DL','DOW_DD_DL',
                  'DOW_toutCP_DL_l1DailyCP','DOW_tout2CP_DL_l1DailyFlexCP')),
       aes(x=DW,color=model.name)) + geom_density() + 
  labs(title='Durbin-Watson test statistic across model runs by model type',
       x='DW')


# heating cumsum
ord = order(cp1Data$lower)
plot(-1 * cp1Data$lower[ord][cp1Data$lower[ord] < 0],pch=20,ylim=c(0,10),main='Magnitude of heating demand (kWh/HDD)',ylab='kWh/HDD',xlab='Count of residences')
grid()

plot(cumsum(-1 * cp1Data$lower[ord][cp1Data$lower[ord] < 0]),pch=20,main='Magnitude of heating demand (kWh/HDD)',ylab='kWh/HDD',xlab='Count of residences')
grid()

# cooling cumsum
ord = rev(order(cp1Data$upper))
plot(cp1Data$upper[ord][cp1Data$upper[ord] > 0],pch=20,main='Magnitude of cooling demand (kWh/CDD)',ylab='kWh/CDD',xlab='Count of residences',ylim=c(0,10))
grid()

cooling = cp1Data$upper*pmax(0,90 - cp1Data$cp)
ord = rev(order(cooling))
plot(cumsum( cooling[ord][cooling[ord] > 0] )/1000,type='l',lty=1,main='Cumlative cooling demand (kWh/day)',ylab='MWh/day',xlab='Count of residences')

ts = c(100,90,80,70,60,50)
for(i in 1:length(ts)) {
  t = ts[i]
  cooling = cp1Data$upper*pmax(0,t - cp1Data$cp)
  ord = rev(order(cooling))
  cum = cumsum( cooling[ord][cooling[ord] > 0] )/1000
  if(i==1) {
    plot(cum,type='l',lty=i,main='Cumlative cooling demand (kWh/day)',ylim=c(0,5000),ylab='MWh/day',xlab='Count of residences')
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
    plot(cumsum( heating[ord][heating[ord] > 0] )/1000,type='l',ylim=c(0,5000),lty=i,main='Cumlative heating demand (kWh/day)',ylab='MWh/day',xlab='Count of residences')
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



#extract legend
#https://github.com/hadley/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs
getLegend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

loadShapeBreaks = c('aug','jan','year')
loadShapeColors = c('red','blue','grey30')
loadShapePlots = list()
loadShapePlots[['max']] = ggplot(wbasicsWS,aes(x=max,color='year')) + geom_density() + 
  geom_density(data=wbasicsWS,aes(x=Jan_max,color='jan'))  +
  geom_density(data=wbasicsWS,aes(x=Aug_max,color='aug'))   + 
  theme_bw() + theme(text=element_text(size=16)) + 
  labs(title='max',x='max (kWh/hr)') + 
  #theme(legend.position="none") + 
  scale_colour_manual("timeframe", breaks=loadShapeBreaks, values=loadShapeColors) +
  xlim(0,7.5)

loadShapePlots[['min']] = ggplot(wbasicsWS,aes(x=min)) + geom_density(color='grey10') + 
  geom_density(data=wbasicsWS,aes(x=Jan_min),color='blue')  +
  geom_density(data=wbasicsWS,aes(x=Aug_min),color='red')   +  
  labs(title='min',x='min (kWh/hr)') + 
  theme_bw() + theme(text=element_text(size=16)) + 
  xlim(0,1.5)


loadShapePlots[['mean']] = ggplot(wbasicsWS,aes(x=mean)) + geom_density(color='grey10') + 
  geom_density(data=wbasicsWS,aes(x=Jan_mean),color='blue')  +
  geom_density(data=wbasicsWS,aes(x=Aug_mean),color='red')   +  
  labs(title='mean',x='mean (kWh/hr)') + 
  theme_bw() + theme(text=element_text(size=16)) + 
  xlim(0,3)

loadShapePlots[['range']] = ggplot(wbasicsWS,aes(x=range)) + geom_density(color='grey10') + 
  geom_density(data=wbasicsWS,aes(x=Jan_range),color='blue')  +
  geom_density(data=wbasicsWS,aes(x=Aug_range),color='red')   +  
  labs(title='range',x='range (kWh/hr)') + 
  theme_bw() + theme(text=element_text(size=16)) + 
  xlim(0,7)

loadShapePlots[['mn2mx']] = ggplot(wbasicsWS,aes(x=mn2mx)) + geom_density(color='grey10') + 
  geom_density(data=wbasicsWS,aes(x=Jan_mn2mx),color='blue')  +
  geom_density(data=wbasicsWS,aes(x=Aug_mn2mx),color='red')   +  
  labs(title='min/max',x='min/max') + 
  theme_bw() + theme(text=element_text(size=16)) + 
  xlim(0,0.75)

loadShapePlots[['n2d']] = ggplot(wbasicsWS,aes(x=n2d)) + geom_density(color='grey10') + 
  geom_density(data=wbasicsWS,aes(x=Jan_n2d),color='blue')  +
  geom_density(data=wbasicsWS,aes(x=Aug_n2d),color='red')   +  
  labs(title='night/day',x='night/day') + 
  theme_bw() + theme(text=element_text(size=16)) + 
  xlim(0,3)

mylegend <- getLegend(loadShapePlots[['max']])
arg.list = c(loadShapePlots,list(
  #main='Daily load shape metrics averaged over Jan., Aug. and all year',
  ncol=3))
do.call(grid.arrange,arg.list)   



sclrs2 = subset(sclrs,model.name=='DOW_tout2CP_DL_l1DailyFlexCP') # 2 change points
sclrs1 = subset(sclrs,model.name=='DOW_toutCP_DL_l1DailyCP')     # 1 change point
sclrs0 = subset(sclrs,model.name=='DOW_tout_DL_l1')          # no change point 

trueCP1 = cp1Data$upper > 0 & cp1Data$lower <= 0 & cp1Data$pval.upper < 0.01
trueCP2 = cp2Data$upper > 0 & cp2Data$lower <= 0 & cp2Data$pval.upper < 0.01 & cp2Data$middle_1 > -0.6 & cp2Data$middle_1 < 0.6

trueCP1Data = cp1Data[trueCP1,]
trueCP2Data = cp2Data[trueCP2,]
# histograms of cp1 and cp2 overlaid together
ggplot(trueCP2Data,aes(x=cp1)) + 
  geom_histogram(binwidth=1,fill='#0000ff',alpha=0.5) + 
  geom_histogram(aes(x=cp2),binwidth=1,fill='#ff0000',alpha=0.5) +
  labs(title=paste('Histogram of upper and lower change point (n=',dim(trueCP2Data)[1],')',sep=''),
       x='change points (lower and upper)')

ggplot(trueCP2Data,aes(x=cp2-cp1)) + geom_histogram(binwidth=1) + xlim(0,30) + 
  labs(title='Distance between CPs', x='Degs F', y='count')

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
ggplot(data.frame(results.cfs$DOW_toutCP_DL_l1DailyCP),aes(x=day.length*1000)) + 
    #geom_histogram(binwidth=50) + 
    geom_density() + 
    labs(title='Impact of hours of daylight on power consumption',
         x='delta W during daylight', y='density') + 
    scale_x_continuous(limits=c(-4000,3000),breaks=seq(-4000,3000,by = 500))
mean(results.cfs$DOW_tout_DL_l1[,'day.length'])

zipData = DATA_SOURCE$getZipData()
calMap(zipData,'counts','zip5',main='Meter count by zip code',colorMap=brewer.pal(9,"Blues"),legend.title='# per zip code' )
calMap(StanfordData()$getZipData(),'cecclmzn','zip5',main='CEC climate zones',colorMap=brewer.pal(12,"Paired")[c(-1,-9,-11)] )
calMap(zipData,'climate' ,'zip5',main='PGE climate zones',colorMap=brewer.pal(12,"Paired")[c(-1,-9,-11)] )
calMap(zipData,'geography' ,'zip5',main='PGE sampling zones (10,000 accounts each)',colorMap=brewer.pal(3,"Dark2") )


calMap(ws,'tout','zip5',main='Mean annual temperature (F)',colorMap=rev(brewer.pal(9,"RdBu")) )
calMap(ws,'rain','zip5',main='Total annual rain (inches)',colorMap=brewer.pal(9,"Blues") )
calMap(ws,'dp','zip5',main='Mean annual dew point (F)',colorMap=brewer.pal(9,"Purples") )


# maps of various metrics 
# see zipMap.R for the calMap implementation and other plots
calMap(wbasicMeans,'kw.total',main='Mean annual electricity consumption (kWh)',colorMap=brewer.pal(9,"Reds") )
calMap(basicMeans,'lag0',main='Simple correlation between Tout and kW denamd: corr(Tout,kW)',colorMap=rev(brewer.pal(9,"RdBu")) )


calMap(wbasicMeansWS,'rankDiff',main='Rank difference between annual kWh demand and temperature',colorMap=rev(brewer.pal(9,"RdBu")) )
ggplot(wbasicMeansWS,aes(x=toutC,y=kw.total)) + 
  geom_point() + 
  labs(title='Annual home energy vs. Tout (means by zip code)',
       x='Annual average Tout (C)',
       y='Annual home energy (kWh)') + 
  ylim(0,20000) +
  theme_bw() + theme(text=element_text(size=16)) + 
  scale_x_continuous(breaks=seq(5,30,by=5)) + 
  geom_smooth(method = "lm") + 
  annotate("text",x=11,y=18000,
           label=lm_eqn(lm(kw.total ~ toutC,wbasicMeansWS)),
           colour='black',size=5,parse=T)


ggplot(wbasicMeansWS,aes(x=summer.toutC,y=kw.total)) + 
  geom_point() + 
  labs(title='Annual energy vs. Tout (means by zip code)',x='Summer average Tout (C)',y='Annual energy (kWh)') + 
  ylim(0,20000) +
  theme_bw() + theme(text=element_text(size=16)) + 
  scale_x_continuous(breaks=seq(10,30,by=5)) + 
  geom_smooth(method = "lm") + 
  annotate("text",x=16.5,y=18000,
           label=lm_eqn(lm(kw.total ~ summer.toutC,wbasicMeansWS)),
           colour='black',size=5,parse=T)


wbasicMeansWS$var.winter = wbasicMeansWS$kw.var.winter * wbasicMeansWS$kw.mean
wbasicMeansWS$var.summer = wbasicMeansWS$kw.var.summer * wbasicMeansWS$kw.mean

ggplot(melt(wbasicMeansWS[,c('kw.mean','var.summer','var.winter')],id=c('kw.mean')),aes(x=value*kw.mean,color=variable,fill=variable)) + 
  geom_density(size=1,alpha=0.2) + 
  xlim(0,3) + 
  theme_bw() + theme(text=element_text(size=16)) + 
  labs(title='Seasonal hourly kW variance',x='Hourly variance (kW)') + 
  #scale_colour_brewer(palette="Set1") + 
  scale_fill_brewer(name="season",
                      breaks=c('var.summer', 'var.winter'),
                      labels=c('summer','winter'),palette='Set1' ) +
  scale_color_brewer(name="season",
                      breaks=c('var.summer', 'var.winter'),
                      labels=c('summer','winter'),palette='Set1' )

wbm = melt(wbasicMeansWS[,c('kw.total','median_hh_income_val','owner_hh_size_val','tout','owner_occupied_pct','pct_below_poverty_pct','median_home_value_val','median_pop_age_val','pop_above_65_val','median_rooms_val')],id.vars=c('kw.total'))
ggplot(wbm,aes(x=value,y=kw.total)) + geom_point() + 
  facet_wrap(~ variable,scales='free_x',nrow=2) + 
  ylim(0,30000) +
  theme_bw() + theme(text=element_text(size=16))

# histograms of annual energy demand
ggplot(basics,aes(x=kw.total,color=cecclmzn)) + geom_density() + xlim(0,40000) + 
  labs(title='Prob. density of annual kWh by climate zone',x='Annual kWh',y='density')

ggplot(wbasics,aes(x=kw.total)) + geom_histogram(binwidth=100,colour="gray50",fill='gray50') + xlim(0,30000) + 
  theme_bw() + theme(text=element_text(size=16)) +
  labs(title=paste('Histogram of annual kWh (N=',round_any(dim(wbasics)[1],100,floor),')',sep=''),x='Annual kWh',y='count')  +
  geom_density(aes(y=100*..count..), colour="firebrick",size=1,alpha=0.5)

s   = sum(wbasics$kw.total)
q   = quantile(wbasics$kw.total,c(0.01,0.05,0.1,0.25,0.5,0.75,0.9,0.95,0.99))
s90 = sum(wbasics$kw.total[wbasics$kw.total > q['90%']])

plot(1:length(wbasics$kw.total)/length(wbasics$kw.total)*100,cumsum(sort(wbasics$kw.total,decreasing=T))/s,type='l',
ylab='fraction of total annual energy',xlab=paste('percentile of home (N=',round_any(length(wbasics$kw.total),100,floor),')',sep=''),
main='Cumulative distribution of annual energy', axes=F)
axis(side = 2, at = seq(0,1,0.1))
axis(side = 1, at = seq(0,100,10))
abline(h=seq(0,1,0.1),v=seq(0,100,10),col='gray',lty=3)

ggplot(wbasics,aes(x=log(kw.total))) + geom_histogram(binwidth=0.02,colour="gray50",fill='gray50') + 
  theme_bw() + theme(text=element_text(size=16)) +
  labs(title=paste('Histogram of annual kWh (N=',round_any(dim(wbasics)[1],100,floor),')',sep=''),x='Annual kWh',y='count')

# we have to re-order the sampling regions to get the logical layout for this plot
ggplot(transform(wbasics, geography=factor(geography,levels=c("Coast","Inland Hills","Central Valley"))), aes(x = kw.total)) +
  stat_density(aes(ymax = ..density..,  ymin = -..density..),
               fill = "grey50", colour = "grey50",
               geom = "ribbon", position = "identity") +
  facet_grid(. ~ geography) + 
  theme_bw() + theme(text=element_text(size=16)) +
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank()) +
  coord_flip() + xlim(0,20000) + labs(title='Annual energy by climate',x='Annual kWh')

p1 = ggplot(wbasics,aes(x=min)) + geom_density(size=1) + theme_bw() + 
  labs(title='min (kW)',x='min (kW)') + xlim(0,1)
p2 = ggplot(wbasics,aes(x=mean)) + geom_density(size=1) + theme_bw() + 
  labs(title='mean (kW)',x='mean (kW)') + xlim(0,3)
p3 = ggplot(wbasics,aes(x=max)) + geom_density(size=1) + theme_bw() + 
  labs(title='max (kW)',x='max (kW)') + xlim(0,9)
p4 = ggplot(wbasics,aes(x=range)) + geom_density(size=1) + theme_bw() + 
  labs(title='range (kW)',x='range (kW)') + xlim(0,7)
p5 = ggplot(wbasics,aes(x=mn2mx)) + geom_density(size=1) + theme_bw() + 
  labs(title='min/max',x='min/max') + xlim(0,1)
p6 = ggplot(wbasics,aes(x=n2d)) + geom_density(size=1) + theme_bw() + 
  labs(title='night/day',x='night/day') + xlim(0,3)
grid.arrange(p1,p2,p3,p4,p5,p6,ncol=3,main=paste('Distributions of load metrics (N=',round_any(dim(wbasics)[1],100,floor),')',sep=''))

#Central Valley
#1,550,138
#Inland Hills
#1,781,481
#Coast
#1,139,145


ggplot(wbasics,aes(x=(max.hr.tout-32)*5/9)) + geom_density(size=1) + 
  theme_bw() + theme(text=element_text(size=16)) + 
  labs(title='Temperature at peak hour of demand',x='temp. (C) at peak hour of demand') + 
  xlim(-5,45)

ggplot(basics,aes(x=(max.day.tout-32)*5/9)) + geom_density(size=1) + 
  theme_bw() + theme(text=element_text(size=16)) + 
  labs(title='Temperature for peak day of demand',x='temp. (C) for peak day of demand') + 
  xlim(-5,45)

ggplot(wbasics,aes(x=(min.day.tout-32)*5/9)) + geom_density(size=1) + 
  theme_bw() + theme(text=element_text(size=16)) + 
  labs(title='Temperature for minimum day of demand',x='temp. (C) for minimum day of demand') + 
  xlim(-5,35)

ggplot(wbasics,aes(x=t90kw)) + geom_density(size=1) +
  theme_bw() + theme(text=element_text(size=16)) + 
  labs(title='Average timing of top 10% of demand hours',x='Average hour for top 10% of demand') + 
  scale_x_continuous(breaks=1:24)


ggplot(basics,aes(x=kw90)) + geom_histogram(binwidth=1) +
  theme_bw() + theme(text=element_text(size=16)) + 
  labs(title='Modal timing of top 10% of demand hours',x='Modal hour for top 10% of demand') + 
  scale_x_continuous(breaks=1:24)

ggplot(basics,aes(x=maxHOD)) + geom_histogram(binwidth=1) +
  theme_bw() + theme(text=element_text(size=16)) + 
  labs(title='Max hourly mean',x='Hour with highest hourly mean') + 
  scale_x_continuous(breaks=1:24)



ggplot(wbasics,aes(x=month(as.POSIXlt(max.day.date,origin='1970-01-01')))) + 
  geom_histogram(binwidth=1) + scale_x_continuous(breaks=1:12) + 
  theme_bw() + theme(text=element_text(size=16)) +
  labs(title='Month of peak daily usage',x='Month number')

ggplot(wbasics,aes(x=as.POSIXlt(max.day.date,origin='1970-01-01')$wday)) + 
  geom_histogram(binwidth=1,color='white') + scale_x_continuous(breaks=0:6,labels=c('Su','Mo','Tu','We','Th','Fr','Sa')) + 
  theme_bw() + theme(text=element_text(size=16)) +
  labs(title='Day of peak daily usage',x='Day of week')

ggplot(wbasics,aes(x=as.POSIXlt(max.hr.date,origin='1970-01-01')$wday)) + 
  geom_histogram(binwidth=1,color='white') + scale_x_continuous(breaks=0:6,labels=c('Su','Mo','Tu','We','Th','Fr','Sa')) + 
  theme_bw() + theme(text=element_text(size=16)) +
  labs(title='Day of peak hourly usage',x='Day of week')

ggplot(wbasics,aes(x=hour(as.POSIXlt(max.hr.date,origin='1970-01-01')))) + 
  geom_histogram(binwidth=1) + scale_x_continuous(breaks=0:24) + 
  theme_bw() + theme(text=element_text(size=16)) +
  labs(title='Time of peak hour of usage',x='Hour of day')

# add test for peak due to cooling
wbasics$cooling = wbasics$max.day.pct > 0.65

ggplot(wbasics,aes(x=hour(as.POSIXlt(max.hr.date,origin='1970-01-01')))) + 
  geom_histogram(binwidth=1) + scale_x_continuous(breaks=0:24) + 
  theme_bw() + theme(text=element_text(size=16)) +
  labs(title='Time of peak hour of usage',x='Hour of day') + 
  facet_grid(. ~ cooling,labeller=label_both)

ggplot(wbasics,aes(x=kw.total,color=cooling)) + 
  geom_density(size=1) + xlim(0,30000) +
  theme_bw() + theme(text=element_text(size=16)) +
  labs(title='Annual kWh conditional on high temp. peak demand',x='Annual energy (kWh)')

ggplot(wbasics,aes(x=kw.total,fill=cooling)) + 
  geom_histogram(binwidth=100) + xlim(0,30000) +
  theme_bw() + theme(text=element_text(size=16)) +
  labs(title='Annual kWh conditional on high temp. peak demand',x='Annual energy (kWh)')


basics$minLevels = cut(basics$min * 1000,seq(0,1000,100))
ggplot(basics,aes(x=kw.total,fill=minLevels)) + geom_histogram(binwidth=100) + xlim(0,40000) + labs(title='Count of homes by annual kWh',x='Annual kWh',y='Count')
ggplot(basics,aes(x=kw.total,color=minLevels)) + geom_density() + xlim(0,40000) + labs(title='Density of annual energy use by daily minimum',x='Annual kWh',y='Count')
ggplot(basics,aes(x=min,y=kw.total)) + geom_point() + ylim(0,40000) + labs(title='Constant loads vs. total annual energy',x='Average daily min demand (kW)',y='annual energy (kWh)')

# fraction of annual energy from daily min loads * 24
baseFraction = wbasics$min*365*24 / wbasics$kw.total
plot(1:length(baseFraction)/length(baseFraction)*100,sort(baseFraction),type='l',
     ylab='fraction of annual energy',xlab=paste('percentile of home (N=',round_any(length(baseFraction),100,floor),')',sep=''),
     main='Fraction of annual kWh from minimum loads', axes = FALSE)
axis(side = 2, at = seq(0,1,0.1))
axis(side = 1, at = seq(0,100,10))
abline(h=seq(0,1,0.1),v=seq(0,100,10),col='gray',lty=3)

basics$maxLevels = cut(basics$max * 1000,seq(0,4000,600))
ggplot(basics,aes(x=kw.total,fill=maxLevels)) + geom_histogram(binwidth=100) + xlim(0,40000) + labs(title='Count of homes by annual kWh',x='Annual kWh',y='Count')
ggplot(basics,aes(x=kw.total,color=maxLevels)) + geom_density() + xlim(0,40000) + labs(title='Count of homes by annual kWh',x='Annual kWh',y='Count')


# income vs energy use
ggplot(basicMeans,aes(x=median_income,y=kw.total)) + geom_point(colour="black") + ylim(0,20000) + labs(x='',y='')
ggplot(basicMeans,aes(x=median_income,y=kw.total,color=cecclmzn)) + geom_point() + ylim(0,20000) + labs(title='Zip code median income vs. mean annual energy (kWh)',y='Annual kWh',x='Median income')

# best fit moving average width
ggplot(sclrsBasicsWS,aes(x=maxMA,color=cecclmzn,fill=cecclmzn)) + geom_histogram(binwidth=1) +
    labs(title='Best fit moving average width by climate zone',x='moving average width',y='count')
# todo: use MA results with good fit from thermal models.
# i.e. define a class of high thermal and low thermal response and locate them within 
# the MA histogram 

dSE = results.stde$DOW_tout2CP_DL_l1DailyFlexCP - results.stdeNW$DOW_tout2CP_DL_l1DailyFlexCP
par(oma=c(0,6,0,0))
barplot(colMeans(dSE[,-c(1,2,dim(dSE)[2])]),
        las=1,horiz=T,
        main='Difference in SE after Newey-West')


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

g = ggplot(trueCP1Data,aes(y=cp,x=cut(tout, breaks=c(seq(40,80,3),Inf))))
g + geom_boxplot() + ylim(10,90) +
  labs(title='Change points as a function of mean annual temperatures',
       x='annual temperature bin (deg F)',y='change point') +
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

g = ggplot(trueCP1Data,aes(y=cp,x=cecclmzn))
g + geom_boxplot() + ylim(10,90) +
  labs(title='Change points as a function of climate zone',
       x='Climate zone',y='change point') +
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

hists(   sclrs,metric='r.squared',zip='all data (200k homes)')
zipHists(sclrs,metric='r.squared',model.name='toutDailyCP')
hists(   sclrs,metric='sigma',xlim=c(0,20),zip='all data (200k homes)')
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

                                                                                                    

