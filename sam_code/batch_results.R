#.libPaths('~/R/library') # use my local R library even from the command line

require(ggplot2)
require(gtools)
require(gridExtra)
require(reshape2)
require(hexbin)
require(plyr)
require(xtable)
require(scales) # for muted
library(fitdistrplus) # for fitdist for log normal
library(extrafont) # allows the use of system fonts
font_import(pattern="[T/t]imes",prompt=F) # Times new roman
font_import(pattern="[A/a]rial",prompt=F) # Arial
#fonts()
# fonttable()

setMethod("display", signature = className("xtable", "xtable"),definition = xtable::display)


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
source(file.path(getwd(),'occupancy.R'))
source(file.path(getwd(),'xtable_printbold.R'))

theme_thesis = function(size=16) {
  return(theme_bw() + theme(text=element_text(size=size,family="Times New Roman")))
}

theme_pub = function(size=16) {
  #return(theme_bw() + theme(text=element_text(size=size,family="Arial")))
  return(theme_bw() + theme(text=element_text(size=size,family="Times")))
}



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
empDir = 'C:/Users/Sam/Dropbox/writing/empirical_paper'
dailyDir = 'C:/Users/Sam/Dropbox/writing/daily_regression_paper'

resultsDir = 'results_pixel_test' 

resultsDir = 'results_cp1_energy'
resultsDir = 'results_daily_full'       # daily models for all homes
resultsDir = 'results_daily_full_2'
#resultsDir = 'results_daily_flex'  # 2 change point models
resultsDir = 'results_daily_standard'
resultsDir = 'results_daily_nestedCP'
resultsDir = 'w_results_occ'
resultsDir = 'w_results_occ_segments'
resultsDir = 'w_results_occ_segments_CP24'
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
resultEnergyFile  = file.path(getwd(),resultsDir,'resultEnergy.RData')
cpDataFile        = file.path(getwd(),resultsDir,'cpData.RData')
cp1File           = file.path(getwd(),resultsDir,'cp1.RData')
occupancyFile     = file.path(getwd(),resultsDir,'occupancy.RData')
coeffDataFile     = file.path(getwd(),resultsDir,'coeffData.RData')
bestModelDataFile = file.path(getwd(),resultsDir,'bestModelData.RData')
invalidIdsFile    = file.path(getwd(),resultsDir,'invalidIds.RData')

#loadData = function(zip) {
#  dataFile = file.path(getwd(),resultsDir,paste(zip,'_modelResults.RData',sep=''))
#  load(dataFile)
#}

fixNames = function(n) {
  n <- gsub('DailyFlexCP','',n)
  n <- gsub('DailyCP','',n)
  n <- gsub('tout_DL_65$','toutCP65_DL',n)
  n <- gsub('l1','L1',n)
  n <- gsub('DD','CDH',n)
  n <- gsub('toutNP','toutFCP',n)
  n <- gsub('_','+',n)
  return(n)
  
}

if(file.exists(invalidIdsFile)) {
  load(invalidIdsFile)
} else {
  invalids = combine(allZips, resultType='invalid.ids',
                     fun=function(x) { return(x$id) }, appendZipData=F)
  colnames(invalids) = c('id','zip')
  save(list=c('invalids'),file=invalidIdsFile)
}
gc() 

segmentedOut = function(x,lowestCount=0) {
  out = data.frame(c())
  x$sid = unlist(x$id)
  for(sp_id in unique(x$sid)) {
    print(sp_id)
    sub = x[x$sid == sp_id,]
    dce = overlappingOutlierDates(sub)
    dens = eventDensities(list(d=dce$dateCounts$dates[dce$dateCounts$Freq >= lowestCount],e=dce$meane))
    densRow = data.frame(t(unlist(dens))) # convert to df with sensibly named columns
    densRow$id = sp_id
    #print(names(densRow))
    out = rbind(out,densRow)
  }
  #print(eventDensities(dateCounts$dates)) # all
  #print(eventDensities(dateCounts$dates[dateCounts$Freq > 1])) # counts of 2 or 3
  #print(eventDensities(dateCounts$dates[dateCounts$Freq > 2])) # counts of 2 or 3
  #print(names(out))
  return(out)
}

renorm = function(x) { # this fixes a normalization error made during the runs and is harmless if the data is normed properly
  #hw = dim(x)
  #width = hw[2]
  #print(hw)
  #counts = x[,1:(width-1)] * x[,width]
  #sc = rowSums(counts)
  #x[,1:(width-1)] = counts / sc
  #x[,width] = sc
  return(x)
}
# lump all the values from the sliding window regressions into a single data frame
segmentedOutliers = combine(allZips,resultType='summaries',fun=segmentedOut,appendZipData=F,lowestCount=3)
save(list='segmentedOutliers',file=file.path(getwd(),resultsDir,'segmentedOutliers.RData'))

load(file.path(getwd(),resultsDir,'segmentedOutliers.RData'))

# split the data frame into the various parts that are of interest
occ.wkdy = renorm(segmentedOutliers[,grep('wkdy',names(segmentedOutliers))])
occ.wknd = renorm(segmentedOutliers[,grep('wknd',names(segmentedOutliers))])
occ.wday = renorm(segmentedOutliers[,grep('wday',names(segmentedOutliers))])
occ.mon  = renorm(segmentedOutliers[,grep('mon' ,names(segmentedOutliers))])
occ.how  = renorm(segmentedOutliers[,grep('how' ,names(segmentedOutliers))])
occ.hr   = renorm(segmentedOutliers[,grep('hr'  ,names(segmentedOutliers))])

occ.hr$id = segmentedOutliers$id

toutHourlyCoeff = function(x) {
  out = data.frame(c())
  x$sid = unlist(x$id)
  #print(x$coefficients[,'Estimate'])
  b = lapply(x$coefficients,function(y){data.frame(t(y[,'Estimate']))} )
  #print(names(b))
  c = do.call(rbind.fill,b)
  out = data.frame(id=x$sid,c)
  #print(class(out))
  #out = cf(x,unique(x$model.name),fill=T)
  return(out)
}
hourlyCoeff = combine(allZips,resultType='summaries',fun=toutHourlyCoeff,appendZipData=F,fill=T)
save(list=c('occ.wkdy','occ.wknd','occ.wday','occ.mon','occ.how','occ.hr','hourlyCoeff'),file=occupancyFile)

load(file=occupancyFile)
#tSlopes = hourlyCoeff[,c('id',sort(names(hourlyCoeff)[grep('tout60_70',names(hourlyCoeff))]))]
tSlopes = hourlyCoeff[,c('id',names(hourlyCoeff)[grep('tout_upper',names(hourlyCoeff))])]


coverage = apply(tSlopes,1,function(x) sum(! is.na(x[-1]))) # count the number of hours with slope estimates
fullCoverage = which(coverage > 23)
fullSlopes= tSlopes[fullCoverage,]
coolinghrs = data.frame(fullSlopes[-1] > 0.010)
occ.hr.sub = occ.hr[occ.hr$id %in% fullSlopes$id,]
coolingOA = occ.hr.sub[,1:24] * coolinghrs
hmap(coolingOA > 0 & coolingOA < 0.02,colorMap=c('#ffffff','#000000'),main='Hours with p(OA) < 2% and tslope > 10w/deg F')
dev.copy2pdf(file = paste(occupancyDir,'pOA_lt_2_tslope_gt_10w.pdf',sep=''),width=10,height=6)

hmap(coolingOA > 0 & coolingOA < 0.04,colorMap=c('#ffffff','#000000'),main='Hours with p(OA) < 4% and tslope > 10w/deg F')
dev.copy2pdf(file = paste(occupancyDir,'pOA_lt_4_tslope_gt_10w.pdf',sep=''),width=10,height=6)

hrs = apply(coolingOA > 0 & coolingOA < 0.02,1,sum)
hist(hrs,breaks=12,main='Hist of count of hrs p(OA) < 2% and tslope > 10W/deg F')
dev.copy2pdf(file = paste(occupancyDir,'unocc_2pct_cooling_hist.pdf',sep=''),width=10,height=6)

hrs = apply(coolingOA > 0 & coolingOA < 0.04,1,sum)
hist(hrs,breaks=12,main='Hist of count of hrs p(OA) < 4% and tslope > 10W/deg F')
dev.copy2pdf(file = paste(occupancyDir,'unocc_4pct_cooling_hist.pdf',sep=''),width=10,height=6)


matplot(t(fullSlopes[,-1]),type='l')


if(file.exists(resultBasicsFile)) {
  load(resultBasicsFile)
} else {
  basics     = combine(allZips,resultType='features.basic',fun=as.matrix,appendZipData=T)
  basics$idZip = paste(basics$id,basics$zip5,sep='.')
  basics$therm.heating = (basics$therm.mean - basics$therm.min) * 365
  
  basicMeans = combine(allZips,resultType='features.basic',fun=function(x) { t(colMeans(x,na.rm=T)) },appendZipData=T )
  
  basics$kw.total = basics$kw.mean * 365 * 24
  basicMeans$kw.total = basicMeans$kw.mean * 365 * 24
  
  save(list=c('basics','basicMeans'),file=resultBasicsFile)
}
gc()

ws = getWeatherSummary()

wbasics = basics
wbasicMeans = basicMeans
wbasicMeansWS = merge(wbasicMeans,ws,by.x='zip5',by.y='zip5')
wbasicsWS = merge(wbasics,ws,by.x='zip5',by.y='zip5')
wbasicMeansWS$rankDiff = rank(wbasicMeansWS$kw.total) - rank(wbasicMeansWS$tout)


basicMeansWS = merge(basicMeans,ws,by.x='zip5',by.y='zip5')
basicsWS = merge(basics,ws,by.x='zip5',by.y='zip5')
basicMeansWS$rankDiff = rank(basicMeansWS$kw.total) - rank(basicMeansWS$tout)

#resType = 'summaries'
resType = 'd_summaries'
if(file.exists(resultScalarsFile)) {
  load(resultScalarsFile)
} else {
  sclrs       = combine(allZips,resultType='d_summaries',fun=scalars,appendZipData=T)
  sclrs$idZip = paste(sclrs$id,sclrs$zip5,sep='.')
  # fix up the naming for later display / accuracy
  sclrs$model.name <- factor(sclrs$model.name)
  levels(sclrs$model.name) <- fixNames(levels(sclrs$model.name))
  #levels(sclrs$model.name) <- gsub('DailyFlexCP','',levels(sclrs$model.name))
  #levels(sclrs$model.name) <- gsub('DailyCP','',levels(sclrs$model.name))
  #levels(sclrs$model.name) <- gsub('tout_DL_65$','toutCP65_DL',levels(sclrs$model.name))
  #levels(sclrs$model.name) <- gsub('l1','L1',levels(sclrs$model.name))
  #levels(sclrs$model.name) <- gsub('DD','CDH',levels(sclrs$model.name))
  #levels(sclrs$model.name) <- gsub('toutNP','toutFCP',levels(sclrs$model.name))
  #levels(sclrs$model.name) <- gsub('_','+',levels(sclrs$model.name))
  
  save(list=c('sclrs'),file=resultScalarsFile)
}
gc()

sclrsBasicsWS   = merge(sclrs,basicsWS[c('idZip','kw.var','kw.mean','maxMA')],by.x='idZip', by.y='idZip',all.x=T,all.y=F)
gc()

#resType = 'summaries'
resType = 'd_summaries'
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
    results.cfs[[model]]    = combine(allZips,resultType=resType,fun=cf,     model.name=model,appendZipData=F,fill=T)
    #results.stde[[model]]   = combine(allZips,resultType=resType,fun=stderrs,model.name=model,appendZipData=F,fill=T)
    #results.pvs[[model]]    = combine(allZips,resultType=resType,fun=pvals,  model.name=model,appendZipData=F,fill=T)
    #results.tvs[[model]]    = combine(allZips,resultType=resType,fun=tvals,      model.name=model,appendZipData=F)
    #results.stdeNW[[model]] = combine(allZips,resultType=resType,fun=stderrsNW,model.name=model,appendZipData=F)
    #results.pvsNW[[model]]  = combine(allZips,resultType=resType,fun=pvalsNW,  model.name=model,appendZipData=F)
    
  }
  save(list=c('results.cfs','results.stde','results.pvs','results.pvsNW','results.stdeNW'),file=coeffDataFile)
  #save(list=c('results.cfs','results.stde','results.pvs','results.tvs'),file=coeffDataFile)
}
gc()
#names(results.cfs) <- fixNames(names(results.cfs))
#names(results.stde) <- fixNames(names(results.stde))
#names(results.stdeNW) <- fixNames(names(results.stdeNW))
#names(results.pvs) <- fixNames(names(results.pvs))

toutvskwhData = combine(allZips,
                  resultType='d_others',
                  fun=function(x) {
                    dataMatrix = t(sapply(x,
                    function(y) {
                      #print(names(y[,'data']$data))
                      out = c(y[,'id'],t(y[,'data']$data))
                      #names(out) <- c('id',rownames(y$data));
                      return(out)
                    }));
                    return(data.frame(dataMatrix))
                  },
                  appendZipData=F)
badRows = apply(toutvskwhData,1,function(x) any(is.na(x)))
toutvskwhData = toutvskwhData[!badRows,]
plotCenters(kmeans(toutvskwhData[,c(-1,-227)],12),ncol=4,colorMap=colorRampPalette(c('#ffffff','#000000'))(100),main='Cluster centers for scatter density')
dev.copy2pdf(file = file.path(empDir,'scatter_centers_12.pdf'),width=10,height=8)

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
                    appendZipData=T,fill=T)
  cp2Data = merge(ws,cp2Data,by.x='zip5',by.y='zip5')
  cp2Data = merge(basics[,c('id','nObs','kw.mean')],cp2Data,by.x='id',by.y='id')
  save(list=c('cp1Data','cp2Data'),file=cpDataFile)
}
gc()


if(file.exists(cp1File)) {
  load(cp1File)
} else {
  cp1_fit = list()
  cp1_fit$cfs = data.frame(results.cfs$'DOW+toutCP+DL+L1')
  cp1_fit$pvs = data.frame(results.pvs$'DOW+toutCP+DL+L1')
  cp1_fit$pvsNW = data.frame(results.pvsNW$'DOW+toutCP+DL+L1')
  
  names(cp1_fit$cfs) = c('id',paste('cf.',names(cp1_fit$cfs)[-1],sep=''))
  names(cp1_fit$pvs) = c('id',paste('pv.',names(cp1_fit$pvs)[-1],sep=''))
  names(cp1_fit$pvsNW) = c('id',paste('pvNW.',names(cp1_fit$pvsNW)[-1],sep=''))
  
  cp1 = subset(sclrs,subset=sclrs$model.name=='DOW+toutCP+DL+L1')
  cp1 = merge(cp1,cp1_fit$cfs,all.y=F)
  cp1 = merge(cp1,cp1_fit$pvs,all.y=F)
  cp1 = merge(cp1,cp1_fit$pvsNW,all.y=F)
  
  save(list=c('cp1'),file=cp1File)
  rm(cp1_fit)
}
gc()

energyContrib = function(x) {
  out = data.frame(c())
  x$sid = unlist(x$id)
  #print(names(x))
  #print(x$contribution)
  b = lapply(x$contribution,function(y){
      df = data.frame(t(y))
      names(df) <- gsub('tou[.]m','tout.m',names(df))
      return(df)
    } )
  c = do.call(rbind.fill,b)
  out = data.frame(id=x$sid,zip5=unlist(x$zip),model.name=unlist(x$model.name),c)
  #print(names(out))
  return(out)
}

if(file.exists(resultEnergyFile)) {
  load(resultEnergyFile)
} else {
  energyContribution= combine( allZips,
                               resultType='d_summaries',
                               fun=energyContrib,
                               appendZipData=F,fill=T)
  levels(energyContribution$model.name) <- fixNames(levels(energyContribution$model.name))
  save(list=c('energyContribution'),file=resultEnergyFile)
}
gc()

ec = energyContribution[energyContribution$model.name=='DOW+toutCP+DL+L1',]
ggplot(ec,aes(x=tout.mean_upper)) + geom_density() + xlim(-2000,7500)
ggplot(ec,aes(x=tout.mean_lower)) + geom_density() + xlim(-1000,5000)
ecna = ec
ecna$tout.mean_upper[ecna$tout.mean_upper<0 | ecna$tout.mean_upper>15000] = 0
mean(ecna$tout.mean_upper,na.rm=T)

zipCooling = aggregate(tout.mean_upper ~ zip5,data=ecna,FUN=mean)
zipCooling$tout.mean_upper[zipCooling$tout.mean_upper < 0] <- 0
calMap(zipCooling,'tout.mean_upper',intervalStyle='equal',
       main='Mean modeled annual cooling energy by zip code',
       family='Times New Roman',colorMap=brewer.pal(7,"Blues"),legend.title='kWh/yr')
dev.copy2pdf(file = file.path(dailyDir,'cooling_energy_map.pdf'),family='serif',width=8,height=8)

echna = ec
echna$tout.mean_lower[echna$tout.mean_lower<0 | echna$tout.mean_lower>15000] = 0
mean(echna$tout.mean_lower,na.rm=T)

zipHeating = aggregate(tout.mean_lower ~ zip5,data=echna,FUN=mean)
#zipHeating$tout.mean_lower[zipHeating$tout.mean_lower < 0] <- 0
calMap(zipHeating,'tout.mean_lower',intervalStyle='equal',
       main='Mean modeled annual heating energy by zip code',
       family='Times New Roman',colorMap=brewer.pal(7,"Reds"),legend.title='kWh/yr')
dev.copy2pdf(file = file.path(dailyDir,'heating_energy_map.pdf'),width=8,height=8)


ecBasics = merge(energyContribution,basics)
ecRatio = ecBasics$tout.mean_upper / (ecBasics$nObs * ecBasics$kw.mean)
plot(sort(ecRatio[ecBasics$tout.mean_lower > 0 & ecBasics$tout.mean_upper > 0 & ecBasics$X.Intercept. > 0]),type='l',main='Fraction of annual kWh from cooling loads')

cp1 = subset(sclrs,subset=sclrs$model.name=='DOW+toutCP+DL+L1')
cp1 = merge(energyContribution,cp1,all.y=F)
occ.wday.m= combine( allZips,
                     resultType='summaries',
                     #subResultType='DOW_toutCP_DL_l1DailyCP',
                     fun=function(x) {
                       hrMatrix = t(apply(x, 1, 
                                          function(y) {
                                            #print(y)
                                            out = c(id=y$id,y$occ.wday.m)
                                            #names(out) <- c('id',rownames(y$data));
                                            return(out)
                                          }));
                       return(data.frame(hrMatrix))
                     },
                     appendZipData=F)


occ.wknd.m= combine( allZips,
                     resultType='summaries',
                     #subResultType='DOW_toutCP_DL_l1DailyCP',
                     fun=function(x) {
                       hrMatrix = t(apply(x, 1, 
                                          function(y) {
                                            #print(y)
                                            out = c(id=y$id,y$occ.wknd.m)
                                            #names(out) <- c('id',rownames(y$data));
                                            return(out)
                                          }));
                       return(data.frame(hrMatrix))
                     },
                     appendZipData=F)

occ.wkdy.m= combine( allZips,
                     resultType='summaries',
                     #subResultType='DOW_toutCP_DL_l1DailyCP',
                     fun=function(x) {
                       hrMatrix = t(apply(x, 1, 
                                          function(y) {
                                            #print(y)
                                            out = c(id=y$id,y$occ.wkdy.m)
                                            #names(out) <- c('id',rownames(y$data));
                                            return(out)
                                          }));
                       return(data.frame(hrMatrix))
                     },
                     appendZipData=F)


occ.hr.m= combine( allZips,
                   resultType='summaries',
                   #subResultType='DOW_toutCP_DL_l1DailyCP',
                   fun=function(x) {
                     hrMatrix = t(apply(x, 1, 
                                        function(y) {
                                          #print(y)
                                          out = c(id=y$id,y$occ.hr.m)
                                          #names(out) <- c('id',rownames(y$data));
                                          return(out)
                                        }));
                     return(data.frame(hrMatrix))
                   },
                   appendZipData=F)

occ.hr= combine( allZips,
                 resultType='summaries',
                 #subResultType='DOW_toutCP_DL_l1DailyCP',
                 fun=function(x) {
                   hrMatrix = t(apply(x, 1, 
                                      function(y) {
                                        #print(y)
                                        out = c(id=y$id,y$occ.hr)
                                        #names(out) <- c('id',rownames(y$data));
                                        return(out)
                                      }));
                   return(data.frame(hrMatrix))
                 },
                 appendZipData=F)

occ.wkdy= combine( allZips,
                   resultType='summaries',
                   #subResultType='DOW_toutCP_DL_l1DailyCP',
                   fun=function(x) {
                     hrMatrix = t(apply(x, 1, 
                                        function(y) {
                                          #print(y)
                                          out = c(id=y$id,y$occ.wkdy)
                                          #names(out) <- c('id',rownames(y$data));
                                          return(out)
                                        }));
                     return(data.frame(hrMatrix))
                   },
                   appendZipData=F)

occ.mon= combine( allZips,
                  resultType='summaries',
                  #subResultType='DOW_toutCP_DL_l1DailyCP',
                  fun=function(x) {
                    hrMatrix = t(apply(x, 1, 
                                       function(y) {
                                         #print(y)
                                         out = c(id=y$id,y$occ.mon)
                                         #names(out) <- c('id',rownames(y$data));
                                         return(out)
                                       }));
                    return(data.frame(hrMatrix))
                  },
                  appendZipData=F)


occ.how= combine( allZips,
                     resultType='summaries',
                     #subResultType='DOW_toutCP_DL_l1DailyCP',
                     fun=function(x) {
                       hrMatrix = t(apply(x, 1, 
                                          function(y) {
                                            #print(y)
                                            out = c(id=y$id,y$occ.how)
                                            #names(out) <- c('id',rownames(y$data));
                                            return(out)
                                          }));
                       return(data.frame(hrMatrix))
                     },
                     appendZipData=F)



# todo: what is wrong with occ.mon?!
save(list=c('occ.hr.m','occ.wkdy.m','occ.wday.m','occ.wknd.m'),file=occupancyFile)

load(occupancyFile)
k = ksc(as.matrix(occ.hr[,2:25],ncol=24),8,max.iter=30)


zips = unique(cp1$zip5)
c = 0
n = length(zips)
for(z in zips) {
  c = c + 1
  print(paste(z,' (',c,'/',n,')',sep=''))
  w = WeatherClass(z,doMeans=T,useCache=T,doSG=T)
  zipData <- DATA_SOURCE$getAllData(z,useCache=T)
  print('zip data loaded')
  cpSub = subset(cp1,subset=cp1$zip5==z)
  d = dim(cpSub)
  print('subset of cp1 complete')
  cd = 0
  for(sp_id in cpSub$id) {
    cd=cd+1
    tic()
    print(paste(sp_id,' (',cd,'/',d[1],')',sep=''))
    resData = zipData[zipData[,'sp_id']== sp_id,] 
    r = ResDataClass(sp_id,z,weather=w,data=resData,useCache=T)
    dfd = rDFA(r)
    print(names(dfd))
    toc()
  }
  break
  gc()
}


metricArray = function(metric='r.squared',include=NULL,exclude=c(),sclr=NULL) {
  if(is.null(sclr)) { sclr = sclrs } # allow overrige of sclrs to support getting metrics form past runs to fill gaps
  if(! is.null(include)) {
    mNames = include
  } else{
    mNames = setdiff(unique(sclrs$model.name),exclude)
  }
  mNames = mNames[!mNames %in% 'DOW+tout+DL+CP65']
  
  metrics = data.frame(idZip=unique(sclr$idZip))
  for(mName in mNames) {
    print(mName)  
    res = subset(sclr[,c('idZip',metric)],sclr$model.name==mName)
    res = res[!duplicated(res['idZip']),] # there is a bug with duplicates for "DOW+toutCP65+DL"
    #print(dim(res))
    colnames(res) <- c('idZip',paste(mName,'.',metric,sep=''))
    metrics = merge(metrics,res,by='idZip')
    #print(dim(metrics))
    gc()
  }
  return(metrics)
}

# todo:
# r2 vs cv r2?
# f-test for nested models?
# clustered errors?

meanRank = function(metricArray,mult=1) {
  rnks = t(apply(metricArray[,-1],1,FUN=function(x) rank(x*mult)))
  colnames(rnks) <- unlist(strsplit(colnames(rnks),'[.][^.]+$'))
  return(colMeans(rnks))
}

metricDiff = function(metricArray,fun=min,noDiff=F,...) {
  if(noDiff){
    md = metricArray
  } else {
    bestMetric = apply(metricArray,1,fun,...)
    md = metricArray - bestMetric
  }
  names(md) = gsub('tout[.]m','tout_m',names(md))
  names(md) = gsub('[.].+$','',names(md))
  out = data.frame(t(abs(apply(md,2,mean,na.rm=T))))
  names(out) = gsub('[.]','+',names(out))
  names(out) = gsub('_','.',names(out))
  return(out)
}

# find the best r.squared result for each sp_id
if(file.exists(bestModelDataFile)) {
  load(bestModelDataFile)
} else {
  #bestModels = do.call(rbind,by(sclrs,sclrs$id,function(df) df[which.max(df$r.squared),]))
  crapModels = c('DOW+toutCP65+DL')#'DOW+toutFCP+DL+L1','DOW+tout2CP+DL+L1')
  
  # get old metrics from previous run because we need DW and pacf.1
  load(file.path(getwd(),'results_daily_full','resultScalars.RData'))
  oldSclrs = sclrs
  load(file.path(getwd(),'results_daily_full_2','resultScalars.RData'))
  # overwrite the good/new metrics
  metrics.rmse          = metricArray('sigma',        exclude=crapModels)
  metrics.cv.rmse       = metricArray('cv.rmse',      exclude=crapModels)
  metrics.cv.mape       = metricArray('cv.mape',      exclude=crapModels)
  metrics.r.squared     = metricArray('r.squared',    exclude=crapModels)
  metrics.adj.r.squared = metricArray('adj.r.squared',exclude=crapModels)
  metrics.AIC           = metricArray('AIC',          exclude=crapModels)
  metrics.kurtosis      = metricArray('kurtosis',     exclude=crapModels)
  
  # fix the length of the old metrics to match the new
  metrics.lag1cor       = metricArray('pacf.1', exclude=crapModels,sclr=oldSclrs)
  metrics.DW            = metricArray('DW',     exclude=crapModels,sclr=oldSclrs)
  metrics.lag1cor       = merge(metrics.rmse[,c(1,2)],metrics.lag1cor,all.x=T)[,-2]
  metrics.DW            = merge(metrics.rmse[,c(1,2)],metrics.DW,all.x=T)[,-2]
  #metrics.lag1cor       = metricArray('pacf.1',       exclude=crapModels) # this is the correlation between lag0 and lag1 errors
  #metrics.DW            = metricArray('DW',           exclude=crapModels)
  
  naKeeper = function(x){ # long story, but unlist drops NAs
    if(length(x)==0) {
      return(NA)
    } else { return(x) }
  }
  
  
  mNames = factor(setdiff(unique(sclrs$model.name),crapModels))
  best = data.frame(idZip=metrics.cv.rmse$idZip)
  best$rmse =          mNames[laply(apply(metrics.rmse[,-1],1,which.min),naKeeper)]
  best$cv.rmse =       mNames[laply(apply(metrics.cv.rmse[,-1],1,which.min),naKeeper)]
  best$cv.mape =       mNames[laply(apply(metrics.cv.mape[,-1],1,which.min),naKeeper)]
  best$r.squared =     mNames[laply(apply(metrics.r.squared[,-1],1,which.max),naKeeper)]
  best$adj.r.squared = mNames[laply(apply(metrics.adj.r.squared[,-1],1,which.max),naKeeper)]
  best$AIC =           mNames[laply(apply(metrics.AIC[,-1],1,which.min),naKeeper)]
  best$lag1cor =       mNames[laply(apply(metrics.lag1cor[,-1],1,which.min),naKeeper)]
  best$DW =            mNames[laply(apply(abs(2-metrics.DW[,-1]),1,which.min),naKeeper)]
  best$kurtosis =      mNames[laply(apply(abs(metrics.kurtosis[,-1]),1,which.min),naKeeper)]
  
  rnk = data.frame(rmse=meanRank(metrics.rmse[,-1]))
  rnk$cv.rmse = meanRank(metrics.cv.rmse[,-1])
  rnk$cv.mape = meanRank(metrics.cv.mape[,-1])
  rnk$r.squared = meanRank(metrics.r.squared[,-1],mult=-1)
  rnk$adj.r.squared = meanRank(metrics.adj.r.squared[,-1],mult=-1)
  rnk$AIC = meanRank(metrics.AIC[,-1])
  rnk$lag1cor = meanRank(metrics.lag1cor[,-1])
  rnk$DW = meanRank(2 - metrics.DW[,-1])
  rnk$kurtosis = meanRank(abs(metrics.kurtosis[,-1]))
  
  metricDiffs = metricDiff(metrics.r.squared[,-1],max,na.rm=T)
  metricDiffs = rbind(metricDiffs,metricDiff(metrics.adj.r.squared[,-1],max,na.rm=T))
  metricDiffs = rbind(metricDiffs,metricDiff(metrics.rmse[,-1],min,na.rm=T))
  metricDiffs = rbind(metricDiffs,metricDiff(metrics.cv.rmse[,-1],min,na.rm=T))
  metricDiffs = rbind(metricDiffs,metricDiff(metrics.cv.mape[,-1],min,na.rm=T))
  metricDiffs = rbind(metricDiffs,metricDiff(metrics.lag1cor[,-1],min,na.rm=T))
  metricDiffs = rbind(metricDiffs,metricDiff(abs(2-metrics.DW[,-1]),min,na.rm=T))
  metricDiffs = rbind(metricDiffs,metricDiff(abs(metrics.kurtosis[,-1]),min,na.rm=T))
  
  metricMeans = metricDiff(metrics.r.squared[,-1],na.rm=T,noDiff=T)
  metricMeans = rbind(metricMeans,metricDiff(metrics.adj.r.squared[,-1],na.rm=T,noDiff=T))
  metricMeans = rbind(metricMeans,metricDiff(metrics.rmse[,-1],na.rm=T,noDiff=T))
  metricMeans = rbind(metricMeans,metricDiff(metrics.cv.rmse[,-1],na.rm=T,noDiff=T))
  metricMeans = rbind(metricMeans,metricDiff(metrics.cv.mape[,-1],na.rm=T,noDiff=T))
  metricMeans = rbind(metricMeans,metricDiff(metrics.lag1cor[,-1],na.rm=T,noDiff=T))
  metricMeans = rbind(metricMeans,metricDiff(abs(2-metrics.DW[,-1]),na.rm=T,noDiff=T))
  metricMeans = rbind(metricMeans,metricDiff(abs(metrics.kurtosis[,-1]),na.rm=T,noDiff=T))
  
  rowMod <- list()
  rowMod$pos <- list()
  rowMod$pos[[1]] <- c(0,1)
  rowMod$command <- c('\\vspace{2mm} \n')
    
  save(list=c('best','rnk','metricMeans','metricDiff','metrics.rmse','metrics.cv.rmse','metrics.r.squared','metrics.cv.mape',
              'metrics.adj.r.squared','metrics.AIC','metrics.kurtosis',
              'metrics.lag1cor','metrics.DW'),file=bestModelDataFile)
}

tbl = metricDiffs 
rownames(tbl) <- c('R$^2$','adjR$^2$','RMSE','cvRMSE','cvMAPE','lag1cor','DW','kurtosis')
texTable = xtable(t(tbl))
digits(texTable)[-1] = 3
#print(texTable,floating=F,rotate.colnames=F,latex.environments="center",sanitize.text.function=function(x){x})
printbold(texTable,floating=F,rotate.colnames=F,each="column",max=c(F,F,F,F,F,F,F,F),
          latex.environments="center",sanitize.colnames.function=function(x){x},
          file=file.path(dailyDir,'metric_diff_best.tex'))

tbl = metricMeans 
rownames(tbl) <- c('R$^2$','adjR$^2$','RMSE','cvRMSE','cvMAPE','lag1cor','DW','kurtosis')
texTable = xtable(t(tbl))
digits(texTable)[-1] = 3
#print(texTable,floating=F,rotate.colnames=F,latex.environments="center",sanitize.text.function=function(x){x})
printbold(texTable,floating=F,rotate.colnames=F,each="column",max=c(T,T,F,F,F,F,F,F),
          latex.environments="center",sanitize.colnames.function=function(x){x},
          file=file.path(dailyDir,'metric_means.tex'))


# this code can be altered to look at the AIC of nested models.
nested=c('WKND','DOW', 'tout','DOW+tout', 'DOW+tout+DL', 'DOW+toutCP+DL+L1','DOW+tout2CP+DL+L1')#, 'DOW+toutCP+DL+L1')
nested=c('WKND','DOW','DOW+tout')#,'DOW+tout+DL')#, 'DOW+toutCP+DL+L1')
nested.AIC = metricArray('AIC',include=nested)
bestAIC = nested[apply(nested.AIC[,-1],1,which.min)]
table(bestAIC)
aicdf = data.frame(t(as.matrix(table(bestAIC))))
aicdfm = melt(aicdf)
aicdfm$variable <- factor(aicdfm$variable)
levels(aicdfm$variable) <- gsub('[.]','+',levels(aicdfm$variable))
aicdfm$variable <- factor(aicdfm$variable,levels=nested)
ggplot(aicdfm,aes(x=variable,y=value)) + geom_bar(stat='identity')
nested.AIC$best = bestAIC
nested.AIC$zip5 = sapply(strsplit(as.character(nested.AIC$idZip),'[.]'),function(x){x[2]} )

zipData = DATA_SOURCE$getZipData()
cpfans = nested.AIC[nested.AIC$best != 'DOW+tout',]
cpZipCount = aggregate(best ~ zip5,cpfans,FUN=length)
cpZipCount = merge(cpZipCount,zipData[,c('zip5','COUNT(DISTINCT sp_id)')],by='zip5')
names(cpZipCount)[3] <- 'count'
cpZipCount$pct = cpZipCount$best / cpZipCount$count
calMap(cpZipCount,'pct')

#sapply(sclrs,function(x) return(x$id))
cvrmse = metrics.cv.rmse
names(cvrmse) <- gsub('.cv.rmse','',names(cvrmse))
# names(cvrmse)
#[1] "idZip"             "tout"              "WKND"              "DOW"               "DOW+tout"          "DOW+tout+DL"       "DOW+tout+DL+65"   
#[8] "DOW+tout+DL+CP65"  "DOW+tout+DL+L1"    "DOW+tout.min+DL"   "DOW+tout.max+DL"   "DOW+CDH+DL"         "DOW+tout+DL+vac"   "DOW+toutCP+DL"    
#[15] "DOW+toutCP+DL+L1"  "DOW+toutNP+DL+L1"  "DOW+tout2CP+DL+L1"
plotmatrix(cvrmse[sample(1:nrow(cvrmse),10000),c("WKND","DOW","tout","DOW+tout")],
           mapping=aes(alpha=1/10)) + 
  xlim(0,50) + ylim(0,50) + geom_abline(intercept=0,slope=1,color='red') +
  theme_pub() + 
  theme(legend.position = "none") +
  labs(x='row model cross validated RMSE (kWh/day)',y='column model cross validated RMSE (kWh/day)') +
  theme(legend.position = "none")
#dev.copy2pdf(file = file.path(dailyDir,'rmse_comp_primary_drivers.pdf'),width=10,height=10,compress=T)
dev.copy(png,file = file.path(dailyDir,'rmse_comp_primary_drivers.png'),width=10,height=10,units='in',res=300); dev.off()

daylen = plotmatrix(cvrmse[sample(1:nrow(cvrmse),10000),c("DOW+tout","DOW+tout+DL")],
                    mapping=aes(alpha=1/10)) + 
  xlim(0,20) + ylim(0,20) + geom_abline(intercept=0,slope=1,color='red') +
  theme_pub() +
  theme(legend.position = "none") +
  labs(title='Day length term') +
  theme(legend.position = "none")
#dev.copy2pdf(file = file.path(dailyDir,'rmse_comp_day_length.pdf'),width=10,height=10)
linear  = plotmatrix(cvrmse[sample(1:nrow(cvrmse),10000),c("DOW+tout+DL","DOW+tout+DL+L1")],
                     mapping=aes(alpha=1/10)) + 
  xlim(0,20) + ylim(0,20) + geom_abline(intercept=0,slope=1,color='red') +
  theme_pub() +
  theme(legend.position = "none") +
  labs(title='Temperature lag term') +
  theme(legend.position = "none")
breakpt = plotmatrix(cvrmse[sample(1:nrow(cvrmse),10000),c("DOW+tout+DL","DOW+toutCP+DL")],
                     mapping=aes(alpha=1/10)) + 
  xlim(0,20) + ylim(0,20) + geom_abline(intercept=0,slope=1,color='red') +
  theme_pub() +
  theme(legend.position = "none") +
  labs(y='column model cross validated RMSE (kWh/day)',x='row model cross validated RMSE (kWh/day)',title='Linear vs. breakpoint temperature') +
  theme(legend.position = "none")
breakptlag = plotmatrix(cvrmse[sample(1:nrow(cvrmse),10000),c("DOW+toutCP+DL","DOW+toutCP+DL+L1")],
                        mapping=aes(alpha=1/10)) + 
  xlim(0,20) + ylim(0,20) + geom_abline(intercept=0,slope=1,color='red') +
  theme_pub() +
  theme(legend.position = "none") +
  labs(title='breakpoint temperature with/without lag') +
  theme(legend.position = "none")
grid.arrange(daylen,linear,breakpt,breakptlag,ncol=2)#,main=paste('Comparisons',sep=''))
#dev.copy2pdf(file = file.path(dailyDir,'rmse_comp_various.pdf'),width=10,height=10)
dev.copy(png,file = file.path(dailyDir,'rmse_comp_various.png'),width=10,height=10,units='in',res=300); dev.off()

plotmatrix(cvrmse[sample(1:nrow(cvrmse),10000),
                  c("DOW+tout+DL","DOW+tout.min+DL","DOW+tout.max+DL","DOW+CDH+DL")],
           mapping=aes(alpha=1/10)) + 
  xlim(0,40) + ylim(0,40) + geom_abline(intercept=0,slope=1,color='red') +
  theme_pub() +
  theme(legend.position = "none") +
  labs(x='row model cross validated RMSE (kWh/day)',y='column model cross validated RMSE (kWh/day)')
#dev.copy2pdf(file = file.path(dailyDir,'rmse_comp_tout_aggregation_method.pdf'),width=10,height=10)
dev.copy(png,file = file.path(dailyDir,'rmse_comp_tout_aggregation_method.png'),width=10,height=10,units='in',res=300); dev.off()


plotmatrix(cvrmse[sample(1:nrow(cvrmse),40000),
                  c("DOW+CDH+DL","DOW+toutCP65+DL","DOW+toutCP+DL","DOW+tout2CP+DL+L1")],
           mapping=aes(alpha=1/10)) + 
  xlim(0,30) + ylim(0,30) + 
  geom_abline(intercept=0,slope=1,color='red') +
  theme_pub() + 
  theme(legend.position = "none") +
  labs(x='row model cross validated RMSE (kWh/day)',y='column model cross validated RMSE (kWh/day)')
#dev.copy2pdf(file = file.path(dailyDir,'rmse_comp_tout_break_method.pdf'),width=10,height=10)
dev.copy(png,file = file.path(dailyDir,'rmse_comp_tout_break_method.png'),width=10,height=10,units='in',res=300); dev.off()


#plotmatrix(cvrmse[sample(1:nrow(cvrmse),40000),c("DOW+toutCP+DL","DOW+toutCP+DL+L1","DOW+tout2CP+DL+L1")]) + 
#  xlim(0,30) + ylim(0,30) + geom_abline(intercept=0,slope=1,color='red') +
#  theme_bw() + theme(text=element_text(size=16)) + 
#  labs(x='row model cross validated RMSE (kWh/day)',y='column model cross validated RMSE (kWh/day)')

comp = c()
choice = data.frame(c())
metrics = setdiff(names(best),c('idZip','AIC'))
for(nm in names(best[,metrics])){
  sameBest = apply(best[,metrics],2,function(x) best[,nm] == x)
  comp = rbind(comp,colSums(sameBest,na.rm=T)*100 / nrow(sameBest))
  choice = rbind.fill(choice,data.frame(t(as.matrix(table(best[,nm])))))
}

#choice = choice / nrow(best)
metnm = c('RMSE','cvRMSE','cvMAPE','R2','adjR2','lag1cor','DW','kurtosis')
metord = c('RMSE','R2','adjR2','cvRMSE','cvMAPE','lag1cor','DW','kurtosis')

rownames(choice) = metnm #names(best[,metrics])

colnames(comp) <- metnm
rownames(comp) <- colnames(comp)
#hmap(comp,xlabs=colnames(comp),ylabs=rownames(comp),
#     main='Correlation of best fit model choice across metrics',
#     colorMap=colorRampPalette(c('#ffffff','#000000'))(100))
compm = melt(comp)
compm$Var1 = factor(compm$Var1,levels=metord)
compm$Var2 = factor(compm$Var2,levels=rev(metord))
ggplot(compm,aes(x=Var1,y=Var2,fill=value)) + 
  theme_thesis() +
  scale_fill_gradient(name='percent',low='#ffffff',high='#000000',limits=c(0,100)) + 
  labs(x='',y='',title='Correlation of best fit model choice across metrics') + 
  geom_tile()
dev.copy2pdf(file = file.path(dailyDir,'best_fit_metric_comparison.pdf'),width=10,height=6)

choicem = melt(as.matrix(choice))
choicem$Var1 = factor(choicem$Var1,levels=metord)
choicem$Var2 = factor(choicem$Var2) #,levels=rev(metord))
levels(choicem$Var2) <- gsub('tout[.]m','tout!m',levels(choicem$Var2))
levels(choicem$Var2) <- gsub('[.]','+',levels(choicem$Var2))
levels(choicem$Var2) <- gsub('!','.',levels(choicem$Var2))
choicem$Var2 <- factor(choicem$Var2,levels=modelord)
ggplot(choicem,aes(x=Var1,y=Var2,fill=value)) +
theme_thesis() +
  scale_fill_gradient(name='count',
                      low='#000000',high='#ffffff',
                      #limits=c(0,max(choice)),
                      trans='log',
                      breaks=c(10,100,1000,10000,100000),
                      labels = trans_format(round,function(x) { x })) + 
  labs(x='',y='',title='Count of best fit models by metric') + 
  geom_tile()
dev.copy2pdf(file = file.path(dailyDir,'best_fit_counts.pdf'),width=10,height=6)


rnka = rnk[,setdiff(names(rnk),c('AIC'))]
colnames(rnka) <- metnm
rnkm = melt(as.matrix(rnka))
rnkm$Var2 <- factor(rnkm$Var2,levels=metord)
rnkm$Var1 <- factor(rnkm$Var1,levels=modelord)
ggplot(rnkm,aes(x=Var2,y=Var1,fill=value)) +
  theme_thesis() +
  scale_fill_gradient(name='mean rank',
                      low='#ffffff',high='#000000',
                      breaks=seq(1,max(rnk),2),
                      limits=c(1,max(rnk))) +
                      #limits=c(0,max(choice)),
  labs(x='',y='',title='Average rank of model by metric') + 
  geom_tile()# +
  #theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.copy2pdf(file = file.path(dailyDir,'fit_ranks.pdf'),width=10,height=6)


sameBest = apply(best[,-1],2,function(x) best[,'cv.rmse'] == x)
bestOverlap = colSums(sameBest[,-2) / nrow(sameBest)
bom = melt(bestOverlap*100)
bom$variable = names(bestOverlap)
ggplot(bom,aes(y=value,x=variable)) + geom_bar(stat='identity') + 
  theme_thesis() +
  ylim(0,100) +
  labs(x='metric of evaluation',y='percent',title='Percentage of homes with same top choice model for CV RMSE and other metrics')
dev.copy2pdf(file = file.path(dailyDir,'metric_comparison.pdf'),width=10,height=6)


jet <-
  colorRampPalette(rev(c("#00007F", "blue", "#007FFF", "cyan",
                     "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000")))
bestm = melt(best[,c('idZip','r.squared','rmse','cv.rmse','cv.mape','lag1cor','kurtosis')],id.vars='idZip',variable.name='metric',value.name='model_name')
modelord = rev(c("DOW+tout2CP+DL+L1",
                 "DOW+toutCP+DL+L1",
                 "DOW+toutFCP+DL+L1",
                 "DOW+toutCP+DL",
                 "DOW+toutCP65+DL",
                 "DOW+tout+DL+65",
                 "DOW+tout+DL+L1",
                 "DOW+tout+DL+vac",
                 "DOW+CDH+DL",
                 "DOW+tout.max+DL",
                 "DOW+tout.min+DL",
                 "DOW+tout+DL",
                 "DOW+tout",
                 "tout",
                 "DOW",
                 "WKND"))
bestm$model_name = factor(bestm$model_name,levels =modelord)
options(scipen=6)
ggplot(bestm,aes(x=metric,color=NA,fill=model_name)) + geom_bar(position='stack') + 
  labs(title='Count of best fit model by metric of assessment',x='metric',y='count') +
  theme_pub() + 
  scale_fill_manual(values=jet(length(levels(bestm$model_name)))) +
  #scale_y_continuous(labels = comma) + 
  scale_y_reverse(labels=comma) # this ensures that the legend colors have the same stacking order as the figure
dev.copy2pdf(file = file.path(dailyDir,'best_fit_counts_bars.pdf'),width=10,height=6)

rm(bestm)

png(file.path(getwd(),'figures','density_ar2.png'),width=800,height=600)
ggplot(subset(sclrsBasicsWS,
              model.name %in% 
                c('tout','DOW+tout','DOW+tout+DL','DOW+tout+DL+L1',
                  'DOW+tout.min+DL','DOW+tout.max+DL','DOW+CDH+DL',
                  'DOW+toutCP+DL+L1')),
              aes(x=adj.r.squared,color=model.name,linetype=model.name)) + geom_density(size=0.7) + 
              xlim(-0.05,1.0) + 
              theme_thesis(size=14) + 
              labs(title='Adjusted R2 sequence for various thermal properties')
dev.copy2pdf(file = file.path(dailyDir,'adj_R2.pdf'),width=10,height=5)
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
                c(
                  #'tout',
                  'DOW',
                  'DOW+tout','DOW+tout+DL','DOW+tout+DL+L1',
                  'DOW+toutCP+DL+L1','DOW+tout2CP+DL+L1')),
              aes(x=cv.rmse/(kw.mean*24)*100,color=model.name,linetype=model.name)) + geom_density(size=0.7) + 
              xlim(0,120)  + 
              theme_thesis(size=14) + 
              labs(title='Coefficient of Variation for various thermal properties',x='Variation (%)')
dev.copy2pdf(file = file.path(dailyDir,'coeff_var.pdf'),width=10,height=5)

ggplot(subset(sclrsBasicsWS,
              model.name %in% 
                c('DOW+tout+DL','DOW+tout.min+DL','DOW+tout.max+DL','DOW+CDH+DL','DOW+tout+DL+CP65','DOW+toutCP+DL')),
       aes(x=cv.rmse/(kw.mean*24)*100,color=model.name,linetype=model.name)) + geom_density(size=0.7) + 
  xlim(0,120)  + 
  theme_thesis(size=14) + 
  labs(title='Coefficient of Variation for various thermal properties',x='Variation (%)')
dev.copy2pdf(file = file.path(dailyDir,'coeff_var_tout_aggregation.pdf'),width=10,height=5)

ggplot(subset(sclrs,
              model.name %in% 
                c(
                  #'tout',
                  'DOW',
                  'DOW+tout','DOW+tout+DL','DOW+tout+DL+L1',
                  'DOW+toutCP+DL+L1','DOW+tout2CP+DL+L1')),
       aes(x=cv.rmse,color=model.name,linetype=model.name)) + geom_density(size=0.7) + 
  xlim(0,25)  + 
  theme_thesis(size=14) + 
  labs(title='Cross validated RMSE for various thermal models',x='cv.RMSE (kWh/day)')
dev.copy2pdf(file = file.path(dailyDir,'CV_RMSE_all_models.pdf'),width=10,height=5)

ggplot(subset(sclrs,
              model.name %in% 
                c(
                  #'tout',
                  'DOW',
                  'DOW+tout','DOW+tout+DL','DOW+tout+DL+L1',
                  'DOW+toutCP+DL+L1','DOW+tout2CP+DL+L1')),
       aes(x=adj.r.squared,color=model.name,linetype=model.name)) + geom_density(size=0.7) + 
  xlim(-0.05,1)  + 
  theme_thesis(size=14) + 
  labs(title='Adjusted R2 for different Tout aggregation methods')
dev.copy2pdf(file = file.path(dailyDir,'adj_R2_tout_method.pdf'),width=10,height=5)


ggplot(subset(sclrs,
              model.name %in% 
                c('DOW+tout+DL+65','DOW+tout+DL+CP65','DOW+toutNP+DL+L1','DOW+toutCP+DL',
                  'DOW+toutCP+DL+L1','DOW+tout2CP+DL+L1')),
              aes(x=adj.r.squared,color=model.name,linetype=model.name)) + geom_density(size=0.7) + 
              xlim(-0.05,1)  + 
  theme_thesis(size=14) + 
              labs(title='Adjusted R2 for change point models')
dev.copy2pdf(file = file.path(dailyDir,'adj_R2_CP_models.pdf'),width=10,height=5)

ggplot(subset(sclrs,
              model.name %in% c('DOW+toutCP+DL+L1')),
       aes(x=adj.r.squared,color=cecclmzn,linetype=cecclmzn)) + geom_density(size=0.7) + 
       xlim(-0.2,1)  + 
       theme_thesis(size=14) + 
       labs(title='Adjusted R2 for single change point model across CZs')
dev.copy2pdf(file = file.path(dailyDir,'adj_R2_by_CZ.pdf'),width=10,height=5)

ggplot(subset(sclrs,
              zip5 %in% c('93304')),
       aes(x=adj.r.squared,color=model.name)) + geom_density() + 
  xlim(0,1)  + ylim(0,4) +
  labs(title='Adjusted R2 for Bakersfield zip')

# kurtosis of change point model runs
ggplot(subset(sclrs,
              model.name %in% 
                c(
                  #'tout',
                  'DOW',
                  'DOW+tout','DOW+tout+DL','DOW+tout+DL+L1',
                  'DOW+toutCP+DL+L1','DOW+tout2CP+DL+L1')),
       aes(x=kurtosis,color=model.name,linetype=model.name)) + geom_freqpoly(size=0.7,binwidth=0.1,fill=NA) + 
  #xlim(-2,7) + 
  coord_cartesian(xlim=c(-2, 7)) +
  theme_thesis(size=14) + 
  labs(title='Kurtosis of residuals across model runs by model type',
       x='kurtosis')
dev.copy2pdf(file = file.path(dailyDir,'kurtosis.pdf'),width=10,height=5)

ggplot(subset(sclrs,
              model.name %in% 
                c(
                  #'tout',
                  'DOW',
                  'DOW+tout','DOW+tout+DL','DOW+tout+DL+L1',
                  'DOW+toutCP+DL+L1','DOW+tout2CP+DL+L1')),
       aes(x=pacf.1,color=model.name,linetype=model.name)) + geom_density(size=0.7) + 
  xlim(-0.4,1) + 
  theme_thesis(size=14) + 
  labs(title='Lag 1 correlation across model runs by model type',
       x='Corr(t,t-1)')
dev.copy2pdf(file = file.path(dailyDir,'correlation_t_t-1.pdf'),width=10,height=5)

ggplot(subset(sclrs,
              model.name %in% 
                c('DOW+toutCP+DL+L1')),aes(y=r.squared,x=tout)) + 
  geom_point(alpha = 0.03) + labs(title='R2 vs. Tout',x='Tout',y='r.squared')


ggplot(subset(sclrsBasicsWS,
              model.name %in% 
                c('DOW+toutCP+DL+L1')),aes(y=r.squared,x=kw.mean)) + 
  geom_point(alpha = 0.1) + labs(title='R2 vs. kwh',x='kwh',y='r.squared')

ggplot(subset(sclrs,
              model.name %in% 
                c('DOW+toutCP+DL+L1')),aes(y=pacf.1,x=tout)) + 
  geom_point(alpha = 0.03) + labs(title='pacf.1 vs. Tout',x='Tout',y='pacf.1')

ggplot(subset(sclrs,
              model.name %in% 
                c('DOW+toutCP+DL+L1')),aes(y=pacf.1,x=r.squared)) + 
  geom_point(alpha = 0.05) + labs(title='pacf.1 vs. r2',x='R2',y='pacf.1')

ggplot(subset(sclrs,
              model.name %in% 
                c('tout','DOW+tout','DOW+tout+DL','DOW+tout+DL+L1',
                  'DOW+tout.min+DL', 'DOW+tout.max+DL','DOW+CDH+DL',
                  'DOW+toutCP+DL+L1','DOW+tout2CP+DL+L1')),
       aes(x=DW,color=model.name,linetype=model.name)) + geom_density(size=0.7) + 
  theme_pub(size=14) +  
  labs(title='Durbin-Watson test statistic across model runs by model type',
       x='d')
dev.copy2pdf(file = file.path(dailyDir,'durbin-watson.pdf'),width=10,height=5)

# NG data required for these ones
ggplot(basics,aes(x=therm.mean)) + geom_histogram(binwidth=0.1) + xlim(0,6) + 
  labs(x='mean therms per day',y='household count',title='Mean therms per day')
dev.copy2pdf(file = file.path(empDir,'daily_therms_hist.pdf'),width=6,height=6)

ggplot(basics,aes(x=therm.min)) + geom_histogram(binwidth=0.1) + xlim(0,6) + 
  labs(x='base therms per day',y='household count',title='Mean base therms per day')
dev.copy2pdf(file = file.path(empDir,'base_therms.pdf'),width=6,height=6)

#basics$therm.heating = (basics$therm.mean - basics$therm.min) * 365

ggplot(basics,aes(x=therm.heating)) + geom_histogram(binwidth=10) + xlim(-200,1000) + 
  labs(x='Heating therms per year',y='household count',title='Estimated heating therms per year')
dev.copy2pdf(file = file.path(empDir,'heating_therms.pdf'),width=6,height=6)


ggplot(basics,aes(x=therm.heating,color=climate)) + geom_density() + xlim(-200,1000) + 
  labs(x='Heating therms per year',y='density',title='Estimated heating therms per year')
dev.copy2pdf(file = file.path(empDir,'heating_therms_by_CZ.pdf'),width=6,height=6)


ngCounts = aggregate(therm.mean ~ zip5,basics,FUN=function(x) {sum(!is.na(x))}) # count non-null NG usage by zip
heating  = aggregate(therm.heating ~ zip5,basics,FUN=function(x) {mean(x,na.rm=T)})
calMap(ngCounts,'therm.mean')

calMap(basicMeans,'therm.mean') # all NG
calMap(heating,'therm.heating')

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

# add the population weight term to the cp1Data
load('zip_weights.RData')
cp1Data = merge(cp1Data,zip_weights[,c('zip5','weight')],by='zip5',all.x=T)

heatl = list()
countl = list()
ts = c(30,40,50,60,70)
for(i in 1:length(ts)) {
  t = ts[i]
  
  heating = cp1Data$lower*pmin(0,t - cp1Data$cp) # negative values only multiply negative slopes to give kWh
  #print(heating)
  print(sum(heating))
  ord = rev(order(heating))
  #heating = heating * cp1Data$weight
  print(sum(heating))
  heating[which(heating < 0)] <- 0
  
  heatl[[paste(t,'F')]] = cumsum( heating[ord] )/1000
  countl[[paste(t,'F')]] = cumsum( cp1Data$weight[ord] )
}
heatdf = as.data.frame(heatl)
hcountdf = as.data.frame(countl)
names(heatdf) <- names(heatl)
names(hcountdf) <- names(countl)
heatdf$idx = 1:nrow(heatdf)
hcountdf$idx = 1:nrow(hcountdf)

cooll  = list()
countl = list()
ts = c(100,90,80,70,60,50)
for(i in 1:length(ts)) {
  t = ts[i]
  cooling = cp1Data$upper*pmax(0,t - cp1Data$cp)
  ord = rev(order(cooling)) # order before the population weighting
  cooling = cooling * cp1Data$weight
  cooling[which(cooling < 0)] <- 0
  cooll[[paste(t,'F')]] = cumsum( cooling[ord] )/1000
  countl[[paste(t,'F')]] = cumsum( cp1Data$weight[ord] )
}
cooldf = as.data.frame(cooll)
ccountdf = as.data.frame(countl)
names(cooldf) <- names(cooll)
names(ccountdf) <- paste('count',names(countl))
cooldf$idx = 1:nrow(cooldf)
ccountdf$idx = 1:nrow(ccountdf)

heatdfm   = melt(heatdf[seq(1,nrow(heatdf),100),],id.vars='idx')
hcountdfm = melt(hcountdf[seq(1,nrow(hcountdf),100),],id.vars='idx')
heatdfm$type = 'heating'
cooldfm = melt(cooldf[seq(1,nrow(cooldf),100),],id.vars='idx')
ccountdfm = melt(ccountdf[seq(1,nrow(ccountdf),100),],id.vars='idx')
cooldfm$type = 'cooling'
hcdfm = rbind(cooldfm,heatdfm)
counts = rbind(cooldfm,heatdfm)
hcdfm$count = counts$value

#rm(list=c('heatl','cooll','heatdfm','cooldfm'))

hcdfm$type = factor(hcdfm$type,levels=c('cooling','heating'))
hcdfm$variable = factor(hcdfm$variable,levels=rev(paste((3:10)*10,'F')))
levels(hcdfm$variable) <- paste(levels(hcdfm$variable),paste(round((10:3 * 10 - 32)*5/9),'C'),sep='/')
ggplot(hcdfm,aes(x=idx,y=value,linetype=variable,color=variable)) + 
  facet_grid(. ~ type) +
  geom_line(size=1) + 
  theme_thesis() + 
  scale_linetype_discrete("Outside temeprature") + 
  scale_colour_discrete("Outside temeprature") +
  scale_colour_manual("Outside temeprature",values=colorRampPalette(colors=c('#B2182B','#4D4D4D','#2166AC'))(8)) +
  theme(legend.key.width=unit(2,"cm")) +
  labs(title='Cumulative conditioning demand',y='MWh/day',x='Count of residences')
dev.copy2pdf(file = file.path(dailyDir,'cumulative_conditioning.pdf'),width=10,height=6)

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
  theme_pub(size=16) +  
  labs(title='max',x='max (kWh/hr)') + 
  #theme(legend.position="none") + 
  scale_colour_manual("timeframe", breaks=loadShapeBreaks, values=loadShapeColors) +
  xlim(0,7.5)

loadShapePlots[['min']] = ggplot(wbasicsWS,aes(x=min)) + geom_density(color='grey10') + 
  geom_density(data=wbasicsWS,aes(x=Jan_min),color='blue')  +
  geom_density(data=wbasicsWS,aes(x=Aug_min),color='red')   +  
  labs(title='min',x='min (kWh/hr)') + 
  theme_pub(size=16) +  
  xlim(0,1.5)


loadShapePlots[['mean']] = ggplot(wbasicsWS,aes(x=mean)) + geom_density(color='grey10') + 
  geom_density(data=wbasicsWS,aes(x=Jan_mean),color='blue')  +
  geom_density(data=wbasicsWS,aes(x=Aug_mean),color='red')   +  
  labs(title='mean',x='mean (kWh/hr)') + 
  theme_pub(size=16) + 
  xlim(0,3)

loadShapePlots[['range']] = ggplot(wbasicsWS,aes(x=range)) + geom_density(color='grey10') + 
  geom_density(data=wbasicsWS,aes(x=Jan_range),color='blue')  +
  geom_density(data=wbasicsWS,aes(x=Aug_range),color='red')   +  
  labs(title='range',x='range (kWh/hr)') + 
  theme_pub(size=16) +  
  xlim(0,7)

loadShapePlots[['mn2mx']] = ggplot(wbasicsWS,aes(x=mn2mx)) + geom_density(color='grey10') + 
  geom_density(data=wbasicsWS,aes(x=Jan_mn2mx),color='blue')  +
  geom_density(data=wbasicsWS,aes(x=Aug_mn2mx),color='red')   +  
  labs(title='min/max',x='min/max') + 
  theme_pub(size=16) +  
  xlim(0,0.75)

loadShapePlots[['n2d']] = ggplot(wbasicsWS,aes(x=n2d)) + geom_density(color='grey10') + 
  geom_density(data=wbasicsWS,aes(x=Jan_n2d),color='blue')  +
  geom_density(data=wbasicsWS,aes(x=Aug_n2d),color='red')   +  
  labs(title='night/day',x='night/day') + 
  theme_pub(size=16) + 
  xlim(0,3)

mylegend <- getLegend(loadShapePlots[['max']])
arg.list = c(loadShapePlots,list(
  #main='Daily load shape metrics averaged over Jan., Aug. and all year',
  ncol=3))
do.call(grid.arrange,arg.list)   



sclrs2 = subset(sclrs,model.name=='DOW+tout2CP+DL+L1') # 2 change points
sclrs1 = subset(sclrs,model.name=='DOW+toutCP+DL+L1')     # 1 change point
sclrs0 = subset(sclrs,model.name=='DOW+tout+DL+L1')          # no change point 

trueCP1 = cp1Data$upper > 0 & cp1Data$lower <= 0 & cp1Data$pval.upper < 0.01
trueCP2 = cp2Data$upper > 0 & cp2Data$lower <= 0 & cp2Data$pval.upper < 0.01 & cp2Data$middle_1 > -0.6 & cp2Data$middle_1 < 0.6

trueCP1Data = cp1Data[trueCP1,]
trueCP2Data = cp2Data[trueCP2,]
# histograms of cp1 and cp2 overlaid together
ggplot(trueCP2Data,aes(x=cp1)) + 
  geom_histogram(binwidth=1,fill='#0000ff',alpha=0.5) + 
  geom_histogram(aes(x=cp2),binwidth=1,fill='#ff0000',alpha=0.5) +
  theme_thesis() + 
  labs(title=paste('Histogram of upper and lower change point (n=',dim(trueCP2Data)[1],')',sep=''),
       x='change points (lower and upper)')
dev.copy2pdf(file = file.path(dailyDir,'2CP_distros.pdf'),width=8,height=4)

ggplot(trueCP2Data,aes(x=cp2-cp1)) + geom_histogram(binwidth=1) + xlim(0,30) + 
  labs(title='Distance between CPs', x='Degs F', y='count')
dev.copy2pdf(file = file.path(dailyDir,'2CP_distance.pdf'),width=8,height=4)

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
mnlight = mean(results.cfs$'DOW+toutCP+DL+L1'[,'day.length'])
ggplot(data.frame(results.cfs$'DOW+toutCP+DL+L1'),aes(x=day.length)) + 
    #geom_histogram(binwidth=50) + 
    geom_density() + 
    labs(title='Distribution of day length coefficient values',
         x='delta kW during daylight', y='density') + 
    theme_thesis() + 
    scale_x_continuous(limits=c(-4,3),breaks=seq(-4,3,by = 0.5)) +
    geom_vline(xintercept=mnlight,linetype="dashed",color='#555555') + 
    geom_vline(xintercept=0,color='#555555') + 
    annotate('text',
           x=mnlight-0.180, y=0.2,
           angle=90,
           color='#555555',
           label=paste('mean',sprintf('%3.0f',mnlight*1000),'W'))
dev.copy2pdf(file = file.path(dailyDir,'day_length.pdf'),width=8,height=4)


zipData = DATA_SOURCE$getZipData()
calMap(zipData,'counts','zip5',main='Meter count by zip code',colorMap=brewer.pal(9,"Blues"),legend.title='# per zip code',nIntervals=7,intervalStyle='fixed',intervalBreaks=c(0,5,40,80,120,160,200,240,302 ))
dev.copy2pdf(file = file.path(empDir,'meter_count_map.pdf'),width=6,height=6)

calMap(StanfordData()$getZipData(),'cecclmzn','zip5',main='CEC climate zones',colorMap=brewer.pal(12,"Paired")[c(-1,-9,-11)] )
calMap(zipData,'climate' ,'zip5',main='PGE climate zones',colorMap=brewer.pal(12,"Paired")[c(-1,-9,-11)] )

calMap(zipData,'geography' ,'zip5',main='PGE sampling zones (10,000 accounts each)',colorMap=c("#5AAE61","#D53E4F","#4EB3D3")) #brewer.pal(3,"Spectral") )
dev.copy2pdf(file = file.path(empDir,'PGE_sample_zone_map.pdf'),width=6,height=6)


calMap(ws,'tout','zip5',main='Mean annual temperature (F)',colorMap=rev(brewer.pal(9,"RdBu")),nIntervals=9 )
dev.copy2pdf(file = file.path(empDir,'mean_temp_map.pdf'),width=6,height=6)

calMap(ws,'rain','zip5',main='Total annual rain (inches)',colorMap=brewer.pal(9,"Blues") )
calMap(ws,'dp','zip5',main='Mean annual dew point (F)',colorMap=brewer.pal(9,"Purples") )


# maps of various metrics 
# see zipMap.R for the calMap implementation and other plots
calMap(wbasicMeans,'kw.total',main='Mean annual electricity consumption (kWh)',colorMap=brewer.pal(9,"Reds"),nIntervals=7,intervalStyle='fixed',intervalBreaks=c(0,2000,6000,10000,14000,18000,22000,max(wbasicMeans$kw.total,na.rm=T)+1) )
dev.copy2pdf(file = file.path(empDir,'mean_annual_kwh_map.pdf'),width=6,height=6)

calMap(basicMeans,'lag0',main='Simple correlation between Tout and kW denamd: corr(Tout,kW)',colorMap=rev(brewer.pal(9,"RdBu")) )


calMap(wbasicMeansWS,'rankDiff',main='Rank difference between annual kWh demand and temperature',colorMap=rev(brewer.pal(9,"RdBu")) )

ggplot(wbasicMeansWS,aes(x=toutC,y=kw.total)) + 
  geom_point() + 
  labs(title='Annual kWh vs. Tout (by zip code)',
       x='Annual average Tout (C)',
       y='Annual home energy (kWh)') + 
  ylim(0,20000) +
  theme_pub(size=16) +  
  scale_x_continuous(breaks=seq(5,30,by=5)) + 
  geom_smooth(method = "lm") + 
  annotate("text",x=11,y=18000,
           label=lm_eqn(lm(kw.total ~ toutC,wbasicMeansWS)),
           colour='black',size=5,parse=T)
dev.copy2pdf(file = file.path(empDir,'annual_energy_vs_tout.pdf'),width=6,height=6)

ggplot(wbasicMeansWS,aes(x=summer.toutC,y=kw.total)) + 
  geom_point() + 
  labs(title='Annual kWh vs. summer Tout (by zip code)',x='Summer average Tout (C)',y='Annual energy (kWh)') + 
  ylim(0,20000) +
  theme_pub(size=16) + 
  scale_x_continuous(breaks=seq(10,30,by=5)) + 
  geom_smooth(method = "lm") + 
  annotate("text",x=16.5,y=18000,
           label=lm_eqn(lm(kw.total ~ summer.toutC,wbasicMeansWS)),
           colour='black',size=5,parse=T)
dev.copy2pdf(file = file.path(empDir,'annual_energy_vs_tout_summer.pdf'),width=6,height=6)

wbasicMeansWS$var.winter = wbasicMeansWS$kw.var.winter * wbasicMeansWS$kw.mean
wbasicMeansWS$var.summer = wbasicMeansWS$kw.var.summer * wbasicMeansWS$kw.mean

ggplot(melt(wbasicMeansWS[,c('kw.mean','var.summer','var.winter')],id=c('kw.mean')),aes(x=value*kw.mean,color=variable,fill=variable)) + 
  geom_density(size=1,alpha=0.2) + 
  xlim(0,3) + 
  theme_pub(size=16) + 
  labs(title='Seasonal hourly kW variance',x='Hourly variance (kW)') + 
  #scale_colour_brewer(palette="Set1") + 
  scale_fill_brewer(name="season",
                      breaks=c('var.summer', 'var.winter'),
                      labels=c('summer','winter'),palette='Set1' ) +
  scale_color_brewer(name="season",
                      breaks=c('var.summer', 'var.winter'),
                      labels=c('summer','winter'),palette='Set1' )
dev.copy2pdf(file = file.path(empDir,'summer_vs_winter_var.pdf'),width=8,height=5)

wbasicMeansWS$median_hh_income_k = wbasicMeansWS$median_hh_income_val / 1000
wbasicMeansWS$median_home_value_k = wbasicMeansWS$median_home_value_val / 1000
wbm = melt(wbasicMeansWS[,c('kw.total','kw.var','n2d','min','max','range','summer.tout','median_hh_income_k','owner_hh_size_val','tout','owner_occupied_pct','pct_below_poverty_pct','median_home_value_k','median_pop_age_val','median_rooms_val')],id.vars=c('kw.total','min','max','range','kw.var','n2d','summer.tout'))
ggplot(wbm,aes(x=value,y=kw.total)) + geom_point() + 
  facet_wrap(~ variable,scales='free_x',nrow=2) + 
  ylim(0,30000) +
  theme_pub(size=16)
dev.copy2pdf(file = file.path(empDir,'kwh_vs_demographics_faceted.pdf'),width=10,height=5)

ggplot(wbm,aes(x=value,y=min)) + geom_point() + 
  facet_wrap(~ variable,scales='free_x',nrow=2) + 
  ylim(0,1) +
  theme_pub(size=16)

ggplot(wbm,aes(x=value,y=max)) + geom_point() + 
  facet_wrap(~ variable,scales='free_x',nrow=2) + 
  ylim(0,7) +
  theme_pub(size=16)

ggplot(wbm,aes(x=value,y=range)) + geom_point() + 
  facet_wrap(~ variable,scales='free_x',nrow=2) + 
  ylim(0,5) +
  theme_pub(size=16)

ggplot(wbm,aes(x=value,y=kw.var)) + geom_point() + 
  facet_wrap(~ variable,scales='free_x',nrow=2) + 
  ylim(0,5) +
  theme_pub(size=16)

ggplot(wbm,aes(x=value,y=n2d)) + geom_point() + 
  facet_wrap(~ variable,scales='free_x',nrow=2) + 
  ylim(0,4) +
  theme_pub(size=16)

ggplot(wbm,aes(x=value,y=summer.tout)) + geom_point() + 
  facet_wrap(~ variable,scales='free_x',nrow=2) + 
  ylim(40,90) +
  theme_pub(size=16)

# histograms of annual energy demand
ggplot(basics,aes(x=kw.total,color=cecclmzn)) + geom_density() + xlim(0,40000) + 
  labs(title='Prob. density of annual kWh by climate zone',x='Annual kWh',y='density')

slnorm = function(scale,...) { return(scale*dlnorm(...))} # scale log normal function by arbitrary amount
lnm = fitdist(wbasics$kw.total,'lnorm')
ggplot(wbasics,aes(x=kw.total)) + geom_histogram(binwidth=100,colour="gray50",fill='gray50') + xlim(0,30000) + 
  theme_pub(size=16) + 
  labs(title=paste('Histogram of annual kWh (N=',round_any(dim(wbasics)[1],100,floor),')',sep=''),x='Annual kWh',y='count')  +
  #geom_density(aes(y=100*..count..), colour="firebrick",size=1,alpha=0.5) + 
  stat_function(fun=slnorm,args=c(scale=length(wbasics$kw.total)*100,lnm$estimate),color='firebrick',size=1,alpha=0.9)
  #annotate("text",label=lnm$estimate, x =15000,y=250, size = 6, colour = "black")
dev.copy2pdf(file = file.path(empDir,'annual_energy_hist.pdf'),width=10,height=5)

ggplot(wbasics,aes(x=kw.total)) + geom_histogram(binwidth=200) 


s   = sum(wbasics$kw.total)
q   = quantile(wbasics$kw.total,c(0.01,0.05,0.1,0.25,0.5,0.75,0.9,0.95,0.99))
s90 = sum(wbasics$kw.total[wbasics$kw.total > q['90%']])

plot(1:length(wbasics$kw.total)/length(wbasics$kw.total)*100,cumsum(sort(wbasics$kw.total,decreasing=T))/s,type='l',
ylab='fraction of total annual energy',xlab=paste('percentile of home (N=',round_any(length(wbasics$kw.total),100,floor),')',sep=''),
main='Cumulative distribution of annual energy', axes=F)
axis(side = 2, at = seq(0,1,0.1))
axis(side = 1, at = seq(0,100,10))
abline(h=seq(0,1,0.1),v=seq(0,100,10),col='gray',lty=3)
dev.copy2pdf(file = file.path(empDir,'cum_annual_energy.pdf'),width=8,height=5)


ggplot(wbasics,aes(x=log(kw.total))) + geom_histogram(binwidth=0.02,colour="gray50",fill='gray50') + 
  theme_pub(size=16) + 
  labs(title=paste('Histogram of annual kWh (N=',round_any(dim(wbasics)[1],100,floor),')',sep=''),x='Annual kWh',y='count')

# we have to re-order the sampling regions to get the logical layout for this plot
ggplot(transform(wbasics, geography=factor(geography,levels=c("Coast","Inland Hills","Central Valley"))), aes(x = kw.total)) +
  stat_density(aes(ymax = ..density..,  ymin = -..density..),
               fill = "grey50", colour = "grey50",
               geom = "ribbon", position = "identity") +
  facet_grid(. ~ geography) + 
  theme_pub(size=16) + 
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank()) +
  coord_flip() + xlim(0,20000) + labs(title='Annual energy by climate',x='Annual kWh')
dev.copy2pdf(file = file.path(empDir,'annual_energy_violin.pdf'),width=8,height=5)

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

dev.copy2pdf(file = file.path(empDir,'load_metric_distros.pdf'),width=10,height=5)


#Central Valley
#1,550,138
#Inland Hills
#1,781,481
#Coast
#1,139,145


ggplot(wbasics,aes(x=(max.hr.tout-32)*5/9)) + geom_density(size=1) + 
  theme_pub(size=16) + 
  labs(title='Temperature at peak hour of demand',x='temp. (C) at peak hour of demand') + 
  xlim(-5,45)
dev.copy2pdf(file = file.path(empDir,'tout_at_peak_kw.pdf'),width=6,height=6)

ggplot(basics,aes(x=(max.day.tout-32)*5/9)) + geom_density(size=1) + 
  theme_pub(size=16) +  
  labs(title='Temperature for peak day of demand',x='temp. (C) for peak day of demand') + 
  xlim(-5,45)

ggplot(basics,aes(x=max.day.pct*100)) + geom_histogram(binwidth=1) + 
  theme_pub(size=16) +  
  labs(title='Temperature percentile for peak day of demand',x='temp. percentile') + 
  xlim(0,100)
dev.copy2pdf(file = file.path(empDir,'tout_pct_peak_day.pdf'),width=6,height=6)

basics$max.day.pct
ggplot(wbasics,aes(x=(min.day.tout-32)*5/9)) + geom_density(size=1) + 
  theme_pub(size=16) + 
  labs(title='Temperature for minimum day of demand',x='temp. (C) for minimum day of demand') + 
  xlim(-5,35)

ggplot(wbasics,aes(x=t90kw)) + geom_density(size=1) +
  theme_pub(size=16) + 
  labs(title='Average timing of top 10% of demand hours',x='Average hour for top 10% of demand') + 
  scale_x_continuous(breaks=1:24)


ggplot(basics,aes(x=kw90)) + geom_histogram(binwidth=1) +
  theme_pub(size=16) + 
  labs(title='Modal timing of top 10% of demand hours',x='Modal hour for top 10% of demand') + 
  scale_x_continuous(breaks=1:24)
dev.copy2pdf(file = file.path(empDir,'modal_top_10_hrs.pdf'),width=6,height=6)

ggplot(basics,aes(x=maxHOD)) + geom_histogram(binwidth=1) +
  theme_pub(size=16) +  
  labs(title='Max hourly mean',x='Hour with highest hourly mean') + 
  scale_x_continuous(breaks=1:24)



ggplot(wbasics,aes(x=month(as.POSIXlt(max.day.date,origin='1970-01-01')))) + 
  geom_histogram(binwidth=1) + scale_x_continuous(breaks=1:12) + 
  theme_pub(size=16) + 
  labs(title='Month of peak daily usage',x='Month number')
dev.copy2pdf(file = file.path(empDir,'month_of_peak_day.pdf'),width=6,height=6)

ggplot(wbasics,aes(x=as.POSIXlt(max.day.date,origin='1970-01-01')$wday)) + 
  geom_histogram(binwidth=1,color='white') + scale_x_continuous(breaks=0:6,labels=c('Su','Mo','Tu','We','Th','Fr','Sa')) + 
  theme_pub(size=16) + 
  labs(title='Day of peak daily usage',x='Day of week')
dev.copy2pdf(file = file.path(empDir,'day_of_peak_day.pdf'),width=6,height=6)

ggplot(wbasics,aes(x=as.POSIXlt(max.hr.date,origin='1970-01-01')$wday)) + 
  geom_histogram(binwidth=1,color='white') + scale_x_continuous(breaks=0:6,labels=c('Su','Mo','Tu','We','Th','Fr','Sa')) + 
  theme_pub(size=16) + 
  labs(title='Day of peak hourly usage',x='Day of week')
dev.copy2pdf(file = file.path(empDir,'day_of_peak_hour.pdf'),width=6,height=6)

ggplot(wbasics,aes(x=hour(as.POSIXlt(max.hr.date,origin='1970-01-01')))) + 
  geom_histogram(binwidth=1) + scale_x_continuous(breaks=0:24) + 
  theme_pub(size=16) + 
  labs(title='Time of peak hour of usage',x='Hour of day')
dev.copy2pdf(file = file.path(empDir,'time_of_peak.pdf'),width=6,height=6)


# add test for peak due to cooling
wbasics$cooling = wbasics$max.day.pct > 0.65
wbasics$cooling[is.na(wbasics$cooling)] <- F

ggplot(wbasics,aes(x=hour(as.POSIXlt(max.hr.date,origin='1970-01-01')))) + 
  geom_histogram(binwidth=1) + scale_x_continuous(breaks=0:24) + 
  theme_pub(size=16) + 
  labs(title='Time of peak hour of usage',x='Hour of day') + 
  facet_grid(. ~ cooling,labeller=label_both)
dev.copy2pdf(file = file.path(empDir,'time_of_peak_facet.pdf'),width=10,height=5)

ggplot(wbasics,aes(x=kw.total,color=cooling)) + 
  geom_density(size=1) + xlim(0,30000) +
  theme_pub(size=16) + 
  labs(title='Annual kWh conditional on high temp. peak demand',x='Annual energy (kWh)')

ggplot(wbasics,aes(x=kw.total,fill=cooling)) + 
  geom_histogram(binwidth=100) + xlim(0,30000) +
  theme_pub(size=16) + 
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

dev.copy2pdf(file = file.path(empDir,'fraction_min.pdf'),width=10,height=5)

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

dSE = results.stde$'DOW+toutCP+DL+L1' - results.stdeNW$'DOW+toutCP+DL+L1'
nwDiffs = colMeans(dSE[,-c(1,2,dim(dSE)[2])])
nwm = melt(nwDiffs)
nwm$names = rownames(nwm)
barNames =  c('tout < CP','tout > CP','Mo','Tu','We','Th','Fr','Sa','day length','lag(tout65)')
nwm$names = factor(barNames,levels=barNames[c(8,7,6,5,4,3,9,10,1,2)]) # control the bar plot order
ggplot(nwm,aes(x=names,y=value)) + geom_bar(stat='identity')+
  coord_flip() + theme_thesis() + labs(title='Change in SE after Newey-West')
dev.copy2pdf(file = file.path(dailyDir,'NW_SE.pdf'),width=8,height=4)

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
toutScatter(ResDataClass(1064429605,93304,useCache=T,doSG=F),'1') # _/  good fit
toutScatter(ResDataClass(1418575405,93304,useCache=T,doSG=F),'2') # \/  good fit
toutScatter(ResDataClass(6481381805,93304,useCache=T,doSG=F),'3') # \_/ good fit
toutScatter(ResDataClass(6502182810,93304,useCache=T,doSG=F),'4') # \_/ good fit
toutScatter(ResDataClass(6533703010,93304,useCache=T,doSG=F),'5') # WTF bad fit
toutScatter(ResDataClass(3044473610,94704,useCache=T,doSG=F),'6') # \_  ok fit
toutScatter(ResDataClass(3074319310,94704,useCache=T,doSG=F),'7') # \_  ok fit
toutScatter(ResDataClass(3064027805,94704,useCache=T,doSG=F),'8') # __  bad fit
toutScatter(ResDataClass(3074465310,94704,useCache=T,doSG=F),'9') # __  ok fit
mtext("Typical daily kWh vs. mean Tout scatters", line=0, font=2, cex=1.2,outer=TRUE)
par(op)
dev.copy2pdf(file = file.path(dailyDir,'typical_daily_scatters.pdf'),width=8,height=7)


op <- par(no.readonly = TRUE)
par( mfrow=c(1,2), oma=c(2,0,3,0),mar=c(2,2,2,2))# Room for the title
r = ResDataClass(1064429605,93304,useCache=T,doSG=F)
toutScatter(r,'Hourly averages',type='hourly',  useMean=T,xlim=c(30,105),ylim=c(0,8)) # _/ hourly
toutScatter(r,'Daily averages',type='daily',    useMean=T,xlim=c(30,105),ylim=c(0,8)) # _/ daily
#toutScatter(r,'Monthly averages',type='monthly',useMean=T,xlim=c(30,105),ylim=c(0,8)) # _/ monthly
mtext("Impact of time averaging on electricity consumption data", line=0, font=2, cex=1.2,outer=TRUE)
par(op)
dev.copy2pdf(file = file.path(dailyDir,'hr_vs_daily_scatters.pdf'),width=8,height=5)

r = ResDataClass(1064429605,93304,useCache=T,doSG=F);  plot(r,estimates=toutDoubleChangePoint(rDFA(r))) # _/  good fit
r = ResDataClass(1418575405,93304,useCache=T,doSG=F);  plot(r,estimates=toutDoubleChangePoint(rDFA(r))) # \/  good fit
r = ResDataClass(6481381805,93304,useCache=T,doSG=F);  plot(r,estimates=toutDoubleChangePoint(rDFA(r))) # \_/ good fit
r = ResDataClass(6502182810,93304,useCache=T,doSG=F);  plot(r,estimates=toutDoubleChangePoint(rDFA(r))) # \_/ good fit
r = ResDataClass(6533703010,93304,useCache=T,doSG=F);  plot(r,estimates=toutDoubleChangePoint(rDFA(r))) # WTF bad fit
r = ResDataClass(3044473610,94704,useCache=T,doSG=F);  plot(r,estimates=toutDoubleChangePoint(rDFA(r))) # \_  ok fit
r = ResDataClass(3074319310,94704,useCache=T,doSG=F);  plot(r,estimates=toutDoubleChangePoint(rDFA(r))) # \_  ok fit
r = ResDataClass(3064027805,94704,useCache=T,doSG=F);  plot(r,estimates=toutDoubleChangePoint(rDFA(r))) # __  bad fit
r = ResDataClass(3074465310,94704,useCache=T,doSG=F);  plot(r,estimates=toutDoubleChangePoint(rDFA(r))) # __  ok fit


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
  theme_thesis() +
  scale_y_continuous(breaks=seq(20,90,5)) +
  coord_flip()
dev.copy2pdf(file = file.path(dailyDir,'CP_summer_distros.pdf'),width=8,height=6)

g = ggplot(trueCP1Data,aes(y=cp,x=cut(tout, breaks=c(seq(40,80,3),Inf))))
g + geom_boxplot() + ylim(10,90) +
  labs(title='Change points as a function of mean annual temperatures',
       x='annual temperature bin (deg F)',y='change point') +
  theme_thesis() +
  coord_flip() 
dev.copy2pdf(file = file.path(dailyDir,'CP_annual_distros.pdf'),width=8,height=6)

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
zipHists(sclrs,metric='r.squared',model.name='tout')
hists(   sclrs,metric='sigma',xlim=c(0,20),zip='all data (200k homes)')
zipHists(sclrs,metric='sigma',model.name='tout',xlim=c(0,20))

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

                                                                                                    

