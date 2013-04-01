# define model classes:
# 1. Hourly CP (x WKND + lag + solar)
# 2. Pieces CP (x WKND + lag + solar)
# 3. TOD (x DOW)

# test residuals for:
# serial correlation
# RMSE (+ penalty)
# kurtosis

conf.basePath = file.path('~/EnergyAnalytics/batch')
if(Sys.info()['sysname'] == 'Windows') {
  conf.basePath = file.path('f:/dev/pge_collab/EnergyAnalytics/batch')
} else {
  .libPaths('~/R/library') # use my local R library even from the comand line
}
setwd(conf.basePath)
# run 'source' on all includes to load them 
source(file.path(getwd(),'localConf.R'))         # Local computer specific configuration, especially db account info 
source(file.path(getwd(),'dbUtil.R'))            # generic database support functions for things like connection management
source(file.path(getwd(),'DataClasses.R'))       # Object code for getting meter and weather data 
source(file.path(getwd(),'basicFeatures.R'))     # typical max, min, mean, range
source(file.path(getwd(),'regressionSupport.R')) # mostly regressor manipulation
source(file.path(getwd(),'solaRUtil.R'))         # solar geometry
source(file.path(getwd(),'timer.R'))             # adds tic() and toc() functions

library(reshape)
library(timeDate)
library(RColorBrewer)


testModelRun = function(cfg) {
  print(cfg)
  #r = ResDataClass(553991005,93304);
  r = ResDataClass(2547072505,94610);
  df = regressorDF(r)
  summaries   = c()
  d_summaries = c()
  others      = list()
  d_others    = list()
  
  # convert string formulas to ModelDescriptors
  for(mdName in names(cfg$models.hourly)) {
    md = cfg$models.hourly[[mdName]]
    if(class(md) == 'character') {
      cfg$models.hourly[[mdName]] = ModelDescriptor(formula=md,name=mdName)
      print('converted string to ModelDescriptor')
      print(cfg$models.hourly[[mdName]])
    }
  }
  for(mdName in names(cfg$models.daily)) {
    md = cfg$models.daily[[mdName]]
    if(class(md) == 'character') {
      cfg$models.daily[[mdName]] = ModelDescriptor(formula=md,name=mdName)
      print('converted string to ModelDescriptor')
      print(cfg$models.daily[[mdName]])
    }
  }
  
  for(mdName in names(cfg$models.hourly)) {
    md = cfg$models.hourly[[mdName]]
    print(md$name)
    if( length(others[[md$name]]) == 0 ) { others[md$name] = c() }
    runOut = md$run(r,df)
    summaries = rbind(summaries,runOut$summaries)
    if(! empty(runOut$other)) {
        others[[md$name]] = rbind(others[[md$name]],list(id=r$id,data=runOut$other))
    }
  }
  #dfl = regressorDFAggregated(r,norm=FALSE)
  dfd = rDFA(r)
  # daily regressions
  models = cfg$models.daily
  for(mdName in names(cfg$models.daily)) {
    md = cfg$models.daily[[mdName]]
    print(md$name)
    runOut = md$run(r,dfd)
    d_summaries = rbind(d_summaries,runOut$summaries)
    if(! empty(runOut$other)) {
      d_others[[md$name]] = rbind(d_others[[md$name]],list(id=r$id,data=runOut$other))
    }
  }
  return(list(summaries=summaries,others=others,d_summaries=d_summaries,d_others=d_others))
}
  
cfg = list()
cfg$subset = list(
  summer= paste("MOY %in% c(",paste("'M",5:10,"'",sep='',collapse=','),")"),
  winter= paste("MOY %in% c(",paste("'M",c(11:12,1:4),"'",sep='',collapse=','),")"),
  day   = paste("HOD %in% c(",paste("'H",8:19,"'",sep='',collapse=','),")")
  
)
cfg$subset$summer_day=paste(cfg$subset$summer,'&',cfg$subset$day)

cfg$models.hourly = list(
  #foo          = DescriptorGenerator(name='foo',genImpl=cp24Generator,subset=list(all="TRUE",summer=cfg$subset$summer)), 
  #lag          = DescriptorGenerator(name='geomDecay',genImpl=geometricLagGenerator,subset=list(all="TRUE")),
  #cp24         = DescriptorGenerator(name='cp',genImpl=cp24Generator,subset=list(all="TRUE"))
  #lagPieces    = DescriptorGenerator(name='toutPiecesL',genImpl=toutPieces24LagGenerator,subset=list(all="TRUE")),
  #maPieces     = DescriptorGenerator(name='toutPiecesMA',genImpl=toutPieces24MAGenerator,subset=list(all="TRUE")),
  #pieces       = DescriptorGenerator(name='toutPieces',genImpl=toutPieces24Generator,subset=list(all="TRUE")),
  
  #wea          = ModelDescriptor(name='wea',formula="kw ~ tout   + pout + rh + HOW + MOY",subset=list(all="TRUE",summer=cfg$subset$summer)), 
  #wea65        = ModelDescriptor(name='wea65',formula="kw ~ tout65 + pout + rh + HOW + MOY",subset=list(all="TRUE")),
  #HOW          = "kw ~ tout + HOW"
  #toutTOD_WKND = ModelDescriptor(name='toutTOD_WKND',formula="kw ~ 0 + tout65:HODWK + HODWK",subset=list(all="TRUE"))
)
cfg$models.daily = list(
  #tout           = "kwh ~ tout.mean",
  DOW            = "kwh ~ DOW",
  tout_mean      = "kwh ~ tout.mean + DOW",
  tout_mean_WKND = "kwh ~ tout.mean + WKND",
  #tout_mean_vac  = "kwh ~ tout.mean + WKND + vac",
  tout_max       = "kwh ~ tout.max  + DOW",
  #tout_CDD       = "kwh ~ CDD + HDD + DOW",
  tout_CDD_WKND  = "kwh ~ CDD + HDD + WKND",
  wea_mean       = "kwh ~ tout.mean + pout.mean + rh.mean + WKND + vac",
  dailyCPFixed   = DescriptorGenerator(name='toutFixed',genImpl=toutDailyFixedCPGenerator,subset=list(all="TRUE")),
  dailyCP        = DescriptorGenerator(name='tout',genImpl=toutDailyCPGenerator,subset=list(all="TRUE"))
)

runOut = testModelRun(cfg)
summaries = runOut$summaries
others = runOut$others
d_summaries = runOut$d_summaries
d_others = runOut$d_others
#rm(runOut)
print(others$toutPiecesL[,'data']) # all 'other' data is accessed this way
r = ResDataClass(553991005,93304);plot(r,estimates= toutChangePoint(df=rDFA(r),trange=c(50:85),reweight=F))
