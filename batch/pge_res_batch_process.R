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
library(cvTools)

# run 'source' on all includes to load them 
source(file.path(getwd(),'localConf.R'))         # Local computer specific configuration, especially db account info 
source(file.path(getwd(),'dbUtil.R'))            # generic database support functions for things like connection management
source(file.path(getwd(),'DataClasses.R'))       # Object code for getting meter and weather data 
source(file.path(getwd(),'batchSupport.R'))      # Object code for getting meter and weather data 
source(file.path(getwd(),'ksc.R'))               # k-Spectral Clustering (via Jungsuk)
source(file.path(getwd(),'basicFeatures.R'))     # typical max, min, mean, range
source(file.path(getwd(),'regressionSupport.R')) # mostly regressor manipulation
source(file.path(getwd(),'solaRUtil.R'))         # solar geometry
source(file.path(getwd(),'timer.R'))             # adds tic() and toc() functions

# run this if you need it, but everything should be installed by a setup script
# if (!require("RColorBrewer")) { install.packages("RColorBrewer") }

#library('DAAG') # Data Analysis and Graphics package has k-fold cross validation
# cv.lm(df=mydata, model, m=5) # 5 fold cross-validation


cfg = list()
cfg$outDir = 'results_daily_standard'

cfg$SKIP_EXISTING_RDATA = F # don't run models if the RData file for their zip is present
cfg$PLOT_INVALID = F # create png plots for residences that fail validaiton
cfg$PLOT_VALID   = F # create png plots for residences that pass validaiton

cfg$CACHE_QUERY_DATA   = T  # use file cache to store energy data and other query results


cfg$RUN_FEATURES       = F  # run basic features
cfg$RUN_HOURLY_MODELS  = F  # run hourly models
cfg$RUN_DAILY_MODELS   = T  # run daily summary data models (moderate time consuming)
cfg$RUN_MONTHLY_MODELS = F  # run monthly summary data models (moderate time consuming)
cfg$RUN_STEP_SELECTION = F  # run nested model selection algorithm (time consuming)

cfg$INVALID_IDS = NULL

invalidIdsFile = file.path(getwd(),'invalidIds.RData')
if(file.exists(invalidIdsFile)) {
  load(invalidIdsFile)
  cfg$INVALID_IDS = invalids
  rm('invalids')
  print('Using master list of invalid ids. No further validation will occur.')
} else {
  print('No list of invalid sp_ids. They will be checked one at a time.')
}

#cfg$INVALID_IDS = NULL

# generate the string values that will identify the desired subset of a data.frame
# using the command subset(df,subset=str,...)
cfg$subset = list(
  summer= paste("MOY %in% c(",paste("'M",5:10,"'",sep='',collapse=','),")"),
  winter= paste("MOY %in% c(",paste("'M",c(11:12,1:4),"'",sep='',collapse=','),")"),
  day   = paste("HOD %in% c(",paste("'H",8:19,"'",sep='',collapse=','),")")
  
)
cfg$subset$summer_day=paste(cfg$subset$summer,'&',cfg$subset$day)

# passing in forula objects directly creates lots of problems
# formulas are context specific and cant be stored in dataframes
# to get formula from string, call as.formula(str)
# to get a string repr of a formula object, call deparse(fmla)
cfg$models.hourly = list(
  #MOY         = ModelDescriptor(name='MOY',"kw ~ tout + MOY"),             
  #DOW         = ModelDescriptor(name='DOW',"kw ~ tout + DOW"),
  #DOW_HOD65    = ModelDescriptor(name='DOW_HOD65',formula="kw ~ tout65 + DOW + HOD",subset=list(all="TRUE",summer=cfg$subset$summer)),
  #HOW65        = ModelDescriptor(name='HOW65',formula="kw ~ tout65 + HOW",subset=list(all="TRUE",summer=cfg$subset$summer)),
  #  lagPieces    = DescriptorGenerator(name='toutPiecesL',genImpl=toutPieces24LagGenerator,subset=list(all="TRUE")),
  #  maPieces     = DescriptorGenerator(name='toutPiecesMA',genImpl=toutPieces24MAGenerator,subset=list(all="TRUE")),
  #  pieces       = DescriptorGenerator(name='toutPieces',genImpl=toutPieces24Generator,subset=list(all="TRUE")),
  #  lag          = DescriptorGenerator(name='lag',genImpl=lagGenerator,subset=list(all="TRUE"))
  #wea          = ModelDescriptor(name='wea',formula="kw ~ tout   + pout + rh + HOW + MOY",subset=list(all="TRUE",summer=cfg$subset$summer)), 
  #wea65        = ModelDescriptor(name='wea65',formula="kw ~ tout65 + pout + rh + HOW + MOY",subset=list(all="TRUE")),
  #HOW          = ModelDescriptor(name='HOW',"kw ~ tout + HOW")
  #toutTOD_WKND = ModelDescriptor(name='toutTOD_WKND',formula="kw ~ 0 + tout65:HODWK + HODWK",subset=list(all="TRUE"))
  #toutTOD      = ModelDescriptor(name='toutTOD',formula="kw ~ 0 + tout65:HOD + HOD",subset=list(all="TRUE")),
  #toutTOD_d1   = ModelDescriptor(name='toutTOD_d1',formula="kw ~ 0 + tout65:HOD + pout + rh + tout_d1 + HOD",subset=list(summer=cfg$subset$summer)),
  #toutTOD_65d1 = ModelDescriptor(name='toutTOD_65d1',formula="kw ~ 0 + tout65:HOD + pout + rh + tout65_d1 + HOD",subset=list(summer=cfg$subset$summer)),
  #toutTOD_l1   = ModelDescriptor(name='toutTOD_l1',formula="kw ~ 0 + tout65:HOD + pout + rh + tout65_l1 + HOD",subset=list(summer=cfg$subset$summer)),
  #toutTOD_l3   = ModelDescriptor(name='toutTOD_l3',formula="kw ~ 0 + tout65:HOD + pout + rh + tout65_l3 + HOD",subset=list(all="TRUE"))
  #toutTOD_min = ModelDescriptor(name='toutTOD_min',"kw_min ~ 0 + tout:HOD + HOW" # no intercept
  #parts       = DescriptorGenerator(name='parts',genImpl=partsGenerator,subset=list(all="TRUE"))
)

# todo: integration vacation days into regression
cfg$models.daily = list(
#   #tout           = "kwh ~ tout.mean",
#   DOW            = ModelDescriptor(name='DOW',formula="kwh ~ DOW",subset=list(all="TRUE")),
  tout             = "kwh ~ tout.mean",
  WKND             = "kwh ~ WKND",
  DOW              = "kwh ~ DOW",
  DOW_tout         = "kwh ~ DOW + tout.mean",
  DOW_tout_DL      = "kwh ~ DOW + tout.mean + day.length",
  DOW_tout_DL_l1   = "kwh ~ DOW + tout.mean + day.length + tout.mean.65.l1",
  DOW_tout.min_DL  = "kwh ~ DOW + tout.min  + day.length",
  DOW_tout.max_DL  = "kwh ~ DOW + tout.max  + day.length",
  DOW_DD_DL        = "kwh ~ DOW + CDH + HDH + day.length",
  DOW_tout_DL_vac  = "kwh ~ DOW + tout.mean + day.length + vac",
  DOW_toutCP_DL    = DescriptorGenerator(name='DOW_toutCP_DL',  genImpl=toutDailyCPGenerator,    subset=list(all="TRUE"), terms='+ DOW + day.length'), # 1 CP
  DOW_toutCP_DL_l1 = DescriptorGenerator(name='DOW_toutCP_DL_l1',  genImpl=toutDailyCPGenerator,    subset=list(all="TRUE"), terms='+ DOW + day.length + tout.mean.65.l1') # 1 CP
  
#   wea_mean       = "kwh ~ tout.mean + pout.mean + rh.mean + WKND + vac",
#   dailyCPFixed   = DescriptorGenerator(name='toutFixed',genImpl=toutDailyFixedCPGenerator,subset=list(all="TRUE")),
#   dailyCP        = DescriptorGenerator(name='tout',genImpl=toutDailyCPGenerator,subset=list(all="TRUE")),
#  tout_mean_WKND = "kwh ~ tout.mean + WKND + day.length -1"
  #tout_DL        = "kwh ~ day.length"
  #dailyWKND      = ModelDescriptor(    name='dailyWKDN',formula="kwh ~ WKND - 1",subset=list(all="TRUE"),cvReps=8), # no CP
  #dailyDOW       = ModelDescriptor(    name='dailyDOW', formula="kwh ~ DOW - 1",subset=list(all="TRUE"),cvReps=8), # no CP
  #dailyTout      = ModelDescriptor(    name='dailyTout',formula="kwh ~ tout.mean + DOW - 1",subset=list(all="TRUE"),cvReps=8), # no CP
  #dailyDL        = ModelDescriptor(    name='dailyDL',  formula="kwh ~ tout.mean + DOW + day.length -1",subset=list(all="TRUE"),cvReps=8), # no CP
  #dailyCP        = DescriptorGenerator(name='tout1CP',  genImpl=toutDailyCPGenerator,    subset=list(all="TRUE"),cvReps=8), # 1 CP
  #dailyFlexCP    = DescriptorGenerator(name='tout2CP',  genImpl=toutDailyFlexCPGenerator,subset=list(all="TRUE"),cvReps=8)  # 2 CPs
)

cfg$models.monthly = list(
  #tout       = "kwh ~ tout.mean",
  #CDD        = "kwh ~ CDD"
  #CDD_HDD    = "kwh ~ HDD + CDD"
)

cfg$triggerZip=NULL # default to NULL
cfg$truncateAt=-1 # default to -1

tic('batchRun')

cfg$allZips = c()
args = commandArgs(TRUE)
if (length(args) > 0) {
  print('Initializing batch run with command line zips:')
  cfg$allZips = args
  print(cfg$allZips)
} else { # no command line args so do a full run
  print('Initializing batch run with list of all zips in the database')
  cfg$allZips  <- db.getZips(useCache=cfg$CACHE_QUERY_DATA)
}
# bakersfield, oakland
#cfg$allZips = c(93304,94610)

#cfg$allZips = c(94923,94503,94574,94559,94028,94539,94564,94702,94704,94085,
#               95035,94041,95112,95113,95765,95648,95901,94531,94585,95205,
#               95202,93619,93614,93304,93701,95631,95726,95223,95666)

TEST_SINGLE = F
if(TEST_SINGLE) {
  source(file.path(getwd(),'testModelRun.R')) # test harness for regression code
  testModelRun(cfg)
} else {
  print('Beginning batch run')
  runDetails = runModelsByZip(cfg)
  runResultFile <- file.path(getwd(),cfg$outDir,paste('runDetails.RData',sep=''))
  save(runDetails,file=runResultFile)
  summarizeRun(runDetails,listFailures=F)
}

plot.modelResults = function(modelResults,vars=NULL,...) {
  op <- par(no.readonly = TRUE)
  if(is.null(vars)) {
    vars = c('kw.mean','kw.var','min.3%','max.97%','kw.tout.cor','kw.pout.cor','daily.kw.var','daily.kw.min.var','daily.kw.max.var')
  }
  a = ceiling(sqrt(length(vars)))
  b = floor(sqrt(length(vars)))
  par(mfrow=c(b,a), oma=c(2,0,3,0),mar=c(2,2,2,2))# Room for the title
  for(var in vars){
    hist(modelResults$features.basic[,var],breaks=100,main=var,xlab=var,mgp=c(1,0,0),tcl=0.5,...)
  }
  zip = 'unknown'
  if(length(modelResults$inputs) > 0) {
    zip = modelResults$inputs$zip
  }
  mtext(paste('Summary stats from',zip), line=0, font=2, cex=1.2,outer=TRUE)
  par(op)
}

if(F) {
  zip = cfg$allZips[2]
  load(file.path(getwd(),cfg$outDir,paste(zip,'_modelResults.RData',sep='')))
  print(names(modelResults))
  plot(modelResults)
  print(modelResults$summaries[1,]$coefficients)
  plot(modelResults$others$toutPiecesMA[[1,'data']]$maCorrs) # correlations for moving averages
  plot(modelResults$others$toutPiecesL[[1,'data']]$lagCorrs) # correlations for lags
  # matrix of values
  a = t(apply(modelResults$others$toutPiecesMA,1,function(x) c(x[['id']][[1]],x[['data']]$maCorrs)))
  a = cbind(a[,1],a[,-1]/apply(abs(a[,-1]),1,max)) # divide through by the max corr value
  colnames(a) <- c('id',paste('maCorr_H',1:(dim(a)[2]-1),sep=''))
  adf = data.frame(a)
  am = melt(adf,id.vars=c('id'))
  ggplot(am,aes(x=variable,y=value)) + geom_line(aes(group=id),alpha=0.1)
  
  r = ResDataClass(820735863,94610);  plot(r,type='hourly',estimates=hourlyChangePoint(regressorDF(r),as.list(1:24),reweight=F))  # heat, no cooling
  r = ResDataClass(553991005,93304);  plot(r,type='hourly',estimates=hourlyChangePoint(regressorDF(r),as.list(1:24),reweight=F))  # very clear cooling 24x7
  r = ResDataClass(554622151,93304);  plot(r,type='hourly',estimates=hourlyChangePoint(regressorDF(r),as.list(1:24),reweight=F))  # very clear cooling possible timed setback
  r = ResDataClass(637321210,93304);  plot(r,type='hourly',estimates=hourlyChangePoint(regressorDF(r),as.list(1:24),reweight=F))  # very clear cooling some bad data in july?
  r = ResDataClass(1064423310,93304); plot(r,type='hourly',estimates=hourlyChangePoint(regressorDF(r),as.list(1:24),reweight=F))  # cooling with high outliers. Unclear setpoint
  r = ResDataClass(1366549405,93304); plot(r,type='hourly',estimates=hourlyChangePoint(regressorDF(r),as.list(1:24),reweight=F))  # heat and cooling. slight tout slopes, bimodal
  r = ResDataClass(1064429605,93304); plot(r,type='hourly',estimates=hourlyChangePoint(regressorDF(r),as.list(1:24),reweight=F))  # VERY high change points. Change with HOD
  
  estimates=hourlyChangePoint(regressorDF(r),as.list(1:24))
  plot(estimates['cp',],type='l',col='blue')
  points(colMeans(r$toutMat),type='l',xlab='HOD',ylab='Mean temperature')
  
  legend(3,75,c('mean Tout','CP fit'),fill=c('black','blue'))
}
