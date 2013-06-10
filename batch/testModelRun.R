#!/usr/bin/Rscript

# todo: is there a better way to detect the current directory?
conf.basePath = file.path('~/EnergyAnalytics/batch')
if(Sys.info()['sysname'] == 'Windows') {
  conf.basePath = file.path('f:/dev/pge_collab/EnergyAnalytics/batch')
} else {
  .libPaths('~/R/library') # use my local R library even from the comand line
}
setwd(conf.basePath)

library(reshape)
library(timeDate)
library(RColorBrewer)

# run 'source' on all includes to load them 
source(file.path(getwd(),'localConf.R'))         # Local computer specific configuration, especially db account info 
source(file.path(getwd(),'dbUtil.R'))            # generic database support functions for things like connection management
source(file.path(getwd(),'DataClasses.R'))       # Object code for getting meter and weather data 
source(file.path(getwd(),'basicFeatures.R'))     # typical max, min, mean, range
source(file.path(getwd(),'regressionSupport.R')) # mostly regressor manipulation
source(file.path(getwd(),'solaRUtil.R'))         # solar geometry
source(file.path(getwd(),'timer.R'))             # adds tic() and toc() functions

testModelRun = function(cfg,r=NULL) {
  print(cfg)
  #r = ResDataClass(553991005,93304);
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
  #parts       = DescriptorGenerator(name='parts',genImpl=partsGenerator,subset=list(all="TRUE"))
)
cfg$models.daily = list(
  tout              = ModelDescriptor(    name='tout',             formula="kwh ~ tout.mean",cvReps=4),
  WKND              = ModelDescriptor(    name='WKND',             formula="kwh ~ WKND",cvReps=4),
  DOW               = ModelDescriptor(    name='DOW',              formula="kwh ~ DOW",cvReps=4),
  DOW_tout          = ModelDescriptor(    name='DOW_tout',         formula="kwh ~ DOW + tout.mean",cvReps=4),
  DOW_tout_DL       = ModelDescriptor(    name='DOW_tout_DL',      formula="kwh ~ DOW + tout.mean + day.length",cvReps=4),
  DOW_tout_DL_65    = ModelDescriptor(    name='DOW_tout_DL_65',   formula="kwh ~ DOW + tout.mean + day.length + tout.mean.65",cvReps=4),
  DOW_tout_DL_CP65  = ModelDescriptor(    name='DOW_tout_DL_CP65', formula="kwh ~ DOW + tout.mean.65lower + tout.mean.65upper + day.length",cvReps=4),
  DOW_tout_DL_l1    = ModelDescriptor(    name='DOW_tout_DL_l1',   formula="kwh ~ DOW + tout.mean + day.length + tout.mean.65.l1",cvReps=4),
  DOW_tout.min_DL   = ModelDescriptor(    name='DOW_tout.min_DL',  formula="kwh ~ DOW + tout.min  + day.length",cvReps=4),
  DOW_tout.max_DL   = ModelDescriptor(    name='DOW_tout.max_DL',  formula="kwh ~ DOW + tout.max  + day.length",cvReps=4),
  DOW_DD_DL         = ModelDescriptor(    name='DOW_DD_DL',        formula="kwh ~ DOW + CDH + day.length",cvReps=4),
  DOW_tout_DL_vac   = ModelDescriptor(    name='DOW_tout_DL_vac',  formula="kwh ~ DOW + tout.mean + day.length + vac",cvReps=4),
  DOW_toutCP_DL     = DescriptorGenerator(name='DOW_toutCP_DL',    genImpl=toutDailyCPGenerator, terms='+ DOW + day.length',subset=list(all="TRUE"),cvReps=1), # 1 CP
  DOW_toutCP_DL_l1  = DescriptorGenerator(name='DOW_toutCP_DL_l1', genImpl=toutDailyCPGenerator, terms='+ DOW + day.length + tout.mean.65.l1',subset=list(all="TRUE"),cvReps=4), # 1 CP
  DOW_toutNP_DL_l1  = DescriptorGenerator(name='DOW_toutNP_DL_l1', genImpl=toutDailyNPCPGenerator, terms='+ DOW + day.length + tout.mean.65.l1',subset=list(all="TRUE"),cvReps=4), # non parametric change point model fixed at 55,65,75
  DOW_tout2CP_DL_l1 = DescriptorGenerator(name='DOW_tout2CP_DL_l1',genImpl=toutDailyFlexCPGenerator, terms='+ DOW + day.length + tout.mean.65.l1',subset=list(all="TRUE"),cvReps=4)  # 2 CPs
  #tout_mean_WKND = "kwh ~ tout.mean + WKND + day.length",
  #tout_DL        = "kwh ~ day.length"
  ##tout_mean_vac  = "kwh ~ tout.mean + WKND + vac",
  #tout_max       = "kwh ~ tout.max  + DOW",
  ##tout_CDD       = "kwh ~ CDD + HDD + DOW",
  #tout_CDD_WKND  = "kwh ~ CDD + HDD + WKND",
  #wea_mean       = "kwh ~ tout.mean + pout.mean + rh.mean + WKND + vac",
  #dailyCPFixed   = DescriptorGenerator(name='toutFixed',genImpl=toutDailyFixedCPGenerator,subset=list(all="TRUE")),
  
  #dailyWKND      = ModelDescriptor(    name='dailyWKDN',formula="kwh ~ WKND - 1",subset=list(all="TRUE"),cvReps=8), # no CP
  #dailyDOW       = ModelDescriptor(    name='dailyDOW', formula="kwh ~ DOW - 1",subset=list(all="TRUE"),cvReps=8), # no CP
  #dailyTout      = ModelDescriptor(    name='dailyTout',formula="kwh ~ tout.mean + DOW - 1",subset=list(all="TRUE"),cvReps=8), # no CP
  #dailyDL        = ModelDescriptor(    name='dailyDL',  formula="kwh ~ tout.mean + DOW + day.length",subset=list(all="TRUE"),cvReps=4, step=T) # no CP
  
  #  dailyTout      = ModelDescriptor(    name='dailyTout',formula="kwh ~ tout.mean + DOW - 1",subset=list(all="TRUE"),cvReps=50), # no CP
  #  dailyCP        = DescriptorGenerator(name='tout',     genImpl=toutDailyCPGenerator,       subset=list(all="TRUE"),cvReps=4) # 1 CP
  #  dailyFlexCP    = DescriptorGenerator(name='tout',     genImpl=toutDailyFlexCPGenerator,   subset=list(all="TRUE"),cvReps=50)  # 2 CPs
)
init = F
if(init) {
  rNA2 = ResDataClass(2547072505,94610); # no tout dep,but troubling kink in 2 cp model
  rNA  = ResDataClass(2846844910,94610); # no tout dep
  rCO  = ResDataClass(553991005,93304);  # cooling only
  #rU   = ResDataClass(1882681258,93304); # U shape
  rU   = ResDataClass(6481381805,93304); # U shape
  #rV   = ResDataClass(6481399605,93304); # V shape
  rV   = ResDataClass(6502182810,93304); # V shape
  # todo: find more without temp dep or heating only...
}
r = rTest
tic()
runOut = testModelRun(cfg,r)
toc(prefix='run time')
summaries = runOut$summaries
others = runOut$others
d_summaries = runOut$d_summaries
d_others = runOut$d_others
rm(runOut)
#print(others$toutPiecesL[,'data']) # all 'other' data is accessed this way
#plot(r,estimates=toutDoubleChangePoint(df=rDFA(r)))

