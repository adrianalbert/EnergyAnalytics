#!/usr/bin/Rscript

# todo: is there a better way to detect the current directory?
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
source(file.path(getwd(),'batchSupport.R'))      # Object code for getting meter and weather data 
source(file.path(getwd(),'ksc.R'))               # k-Spectral Clustering (via Jungsuk)
source(file.path(getwd(),'basicFeatures.R'))     # typical max, min, mean, range
source(file.path(getwd(),'regressionSupport.R')) # mostly regressor manipulation
source(file.path(getwd(),'timer.R'))             # adds tic() and toc() functions

library(reshape)
library(timeDate)
library(RColorBrewer)
# run this if you need it, but everything should be installed by a setup script
# if (!require("RColorBrewer")) { install.packages("RColorBrewer") }

#library('DAAG') # Data Analysis and Graphics package has k-fold cross validation
# cv.lm(df=mydata, model, m=5) # 5 fold cross-validation
#library('cvTools') # cross validation tools

cfg = list()
cfg$outDir = 'results_basics'

cfg$PLOT_INVALID = FALSE # create png plots for residences that fail validaiton
cfg$PLOT_VALID   = TRUE  # create png plots for residences that pass validaiton

cfg$RUN_HOURLY_MODELS     = FALSE # run hourly models
cfg$RUN_AGGREGATED_MODELS = FALSE # run daily and monthly summary data models (moderate time consuming)
cfg$RUN_STEP_SELECTION    = FALSE # run nested model selection algorithm (time consuming)
cfg$RUN_PIECES_24         = FALSE # run piecewise model for every HOD (within hourly)
cfg$RUN_CP_24             = FALSE # fit changepoint for every HOD (within hourly)

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
  #MOY         = "kw ~ tout + MOY",             
  #DOW         = "kw ~ tout + DOW",
  #DOW_HOD65    = list(formula="kw ~ tout65 + DOW + HOD",subset=list(all="TRUE",summer=cfg$subset$summer)),
  #HOW65        = list(formula="kw ~ tout65 + HOW",subset=list(all="TRUE",summer=cfg$subset$summer)),
  #wea          = list(formula="kw ~ tout   + pout + rh + HOW + MOY",subset=list(all="TRUE",summer=cfg$subset$summer)), 
  wea65        = list(formula="kw ~ tout65 + pout + rh + HOW + MOY",subset=list(all="TRUE")),
  #HOW         = "kw ~ tout + HOW",
  toutTOD_WKND = list(formula="kw ~ 0 + tout65:HODWK + HODWK",subset=list(all="TRUE")),
  #toutTOD      = list(formula="kw ~ 0 + tout65:HOD + HOD",subset=list(all="TRUE")),
  #toutTOD_d1   = list(formula="kw ~ 0 + tout65:HOD + pout + rh + tout_d1 + HOD",subset=list(summer=cfg$subset$summer)),
  #toutTOD_65d1 = list(formula="kw ~ 0 + tout65:HOD + pout + rh + tout65_d1 + HOD",subset=list(summer=cfg$subset$summer)),
  #toutTOD_l1   = list(formula="kw ~ 0 + tout65:HOD + pout + rh + tout65_l1 + HOD",subset=list(summer=cfg$subset$summer)),
  toutTOD_l3   = list(formula="kw ~ 0 + tout65:HOD + pout + rh + tout65_l3 + HOD",subset=list(all="TRUE"))
  #toutTOD_min = "kw_min ~ 0 + tout:HOD + HOW" # no intercept
)

# todo: integration vacation days into regression
cfg$models.daily = list(
  #tout           = "kwh ~ tout.mean",
  DOW            = "kwh ~ DOW",
  tout_mean      = "kwh ~ tout.mean + DOW",
  tout_mean_WKND = "kwh ~ tout.mean + WKND",
  #tout_mean_vac  = "kwh ~ tout.mean + WKND + vac",
  tout_max       = "kwh ~ tout.max  + DOW",
  #tout_CDD       = "kwh ~ CDD + HDD + DOW",
  tout_CDD_WKND  = "kwh ~ CDD + HDD + WKND",
  wea_mean       = "kwh ~ tout.mean + pout.mean + rh.mean + WKND + vac"
)

cfg$models.monthly = list(
  #tout       = "kwh ~ tout.mean",
  #CDD        = "kwh ~ CDD"
  #CDD_HDD    = "kwh ~ HDD + CDD"
)

cfg$triggerZip=93304 #NULL # default to NULL
cfg$truncateAt=-1 #-1 # default to -1

tic('batchRun')

cfg$allZips = c()
args = commandArgs(TRUE)
if (length(args) > 0) {
  print('Initializing batch run with command line zips:')
  cfg$allZips = args
  print(cfg$allZips)
} else { # no command line args so do a full run
  print('Initializing batch run with list of all zips')
  cfg$allZips  <- db.getZips()
}
# bakersfield, oakland
cfg$allZips = c(94610,93304)

#cfg$allZips = c(94923,94503,94574,94559,94028,94539,94564,94702,94704,94085,
#               95035,94041,95112,95113,95765,95648,95901,94531,94585,95205,
#               95202,93619,93614,93304,93701,95631,95726,95223,95666)

print('Beginning batch run')
runResult = runModelsByZip(cfg)
summarizeRun(runResult,listFailures=F)



plot.modelResults = function(modelResults,vars=NULL) {
  op <- par(no.readonly = TRUE)
  if(is.null(vars)) {
    vars = c('kw.mean','kw.var','min.3%','max.97%','kw.tout.cor','kw.pout.cor','daily.kw.var','daily.kw.min.var','daily.kw.max.var')
  }
  a = ceiling(sqrt(length(vars)))
  b = floor(sqrt(length(vars)))
  par(mfrow=c(b,a), oma=c(2,0,3,0),mar=c(2,2,2,2))# Room for the title
  for(var in vars){
    hist(modelResults$features.basic[,var],breaks=100,main=var,xlab=var,mgp=c(1,0,0),tcl=0.5)
  }
  zip = 'unknown'
  if(length(modelResults$inputs) > 0) {
    zip = modelResults$inputs$zip
  }
  mtext(paste('Summary stats from',zip), line=0, font=2, cex=1.2,outer=TRUE)
  par(op)
}


if(F) {
  zip = cfg$allZips[1]
  load(file.path(getwd(),cfg$outDir,paste(zip,'_modelResults.RData',sep='')))
  print(names(modelResults))
  plot(modelResults)
  print(modelResults$summaries[1,]$coefficients)


  r = ResDataClass(820735863,94610);  plot(r,type='hourly',estimates=hourlyChangePoint(regressorDF(r),as.list(1:24)))  # heat, no cooling
  r = ResDataClass(553991005,93304);  plot(r,type='hourly',estimates=hourlyChangePoint(regressorDF(r),as.list(1:24)))  # very clear cooling 24x7
  r = ResDataClass(554622151,93304);  plot(r,type='hourly',estimates=hourlyChangePoint(regressorDF(r),as.list(1:24)))  # very clear cooling possible timed setback
  r = ResDataClass(637321210,93304);  plot(r,type='hourly',estimates=hourlyChangePoint(regressorDF(r),as.list(1:24)))  # very clear cooling some bad data in july?
  r = ResDataClass(1064423310,93304); plot(r,type='hourly',estimates=hourlyChangePoint(regressorDF(r),as.list(1:24)))  # cooling with high outliers. Unclear setpoint
  r = ResDataClass(1366549405,93304); plot(r,type='hourly',estimates=hourlyChangePoint(regressorDF(r),as.list(1:24)))  # heat and cooling. slight tout slopes, bimodal
  r = ResDataClass(1064429605,93304); plot(r,type='hourly',estimates=hourlyChangePoint(regressorDF(r),as.list(1:24)))  # VERY high change points. Change with HOD
  
}

estimates=hourlyChangePoint(regressorDF(r),as.list(1:24))
plot(estimates['cp',],type='l',col='blue')
points(colMeans(r$toutMat),type='l',xlab='HOD',ylab='Mean temperature')

legend(3,75,c('mean Tout','CP fit'),fill=c('black','blue'))

