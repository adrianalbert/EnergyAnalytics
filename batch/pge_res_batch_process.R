#!/usr/bin/Rscript

.libPaths('~/R/library') # use my local R library even from the comand line

# todo: is there a better way to detect the current directory?
conf.basePath = file.path('~/EnergyAnalytics/batch')
if(Sys.info()['sysname'] == 'Windows') {
  conf.basePath = file.path('c:/dev/pge_collab/EnergyAnalytics/batch')
}
# run 'source' on all includes to load them 
#source(file.path(conf.basePath,'localConf.R'))         # Sam's local computer specific configuration
source(file.path(conf.basePath,'stanfordConf.R'))     # Stanford on site specific configuration
source(file.path(conf.basePath,'DataClasses.R'))       # Object code for getting meter and weather data 
source(file.path(conf.basePath,'ksc.R'))               # k-Spectral Clustering (via Jungsuk)
source(file.path(conf.basePath,'basicFeatures.R'))     # typical max, min, mean, range
source(file.path(conf.basePath,'regressionSupport.R')) # mostly regressor manipulation
source(file.path(conf.basePath,'timer.R'))             # adds tic() and toc() functions

#library('DAAG') # Data Analysis and Graphics package has k-fold cross validation
# cv.lm(df=mydata, model, m=5) # 5 fold cross-validation
#library('cvTools') # cross validation tools

runModelsByZip = function(zipArray,triggerZip=NULL) {
  nZip     <- length(zipArray)
  zipCount <- 0
  outDir = 'results'
  
  print(paste('This batch process will run models for',nZip,'zip codes'))
  dir.create(file.path(conf.basePath,outDir),showWarnings=FALSE)
  print(paste('<zipcode>_modelResults.RData files for each can be found in',file.path(conf.basePath,outDir)))
  res = list(attemptedZip=zipArray,
             completedZip=c(),
             attemptedSP=c(),
             completedSP=c())
  triggered = FALSE
  for (zip in zipArray) { # i.e. 94610
    zipCount <- zipCount + 1
    if (is.null(triggerZip) || zip == triggerZip) triggered = TRUE
    if (! triggered) {
      next
    }
    
    print(paste('Running models for ',zip,' (',zipCount,'/',nZip,')', sep=''))
    resultsFile <- file.path(conf.basePath,outDir,paste(zip,'_modelResults.RData',sep=''))
    tryCatch( {
      tic('allDataForZip')
      zipData <- db.getAllData(zip)
      toc('allDataForZip')
      sp_ids <- db.getSPs(zip)
      sp_ids <- unique(zipData[,'sp_id']) # db.getSPs(zip)
      res$attemptedSP = c(res$attemptedSP,sp_ids)
      # we know that all sp's in the same zipcode share the same weather
      # so we speed execution by looking it up once and passing it in
      weather <- WeatherClass(zip)
      tic('modelsBySP')
      modelResults <- runModelsBySP(sp_ids,zip=zip,data=zipData,weather=weather)
      rm(zipData,weather,sp_ids)
      save(modelResults,file=resultsFile)
      res$completedZip <- rbind(res$completedZip,zip)
      res$completedSP  <- rbind(res$completedSP,modelResults$ids)
      rm(modelResults)
      toc('modelsBySP')
    }, 
    error = function(e) {
      print(e)
      modelResults = NA
      # put an empty file in place to differentiate between un-run zips and unsuccessful runs
      save(modelResults,e,file=resultsFile)   }, 
    finally={})
    
    # just for testing
    #break
  } # zip loop
  return(res)
}

base.models = list(
  MOY        = formula(kw ~ tout + MOY),
  DOW        = formula(kw ~ tout + DOW),
  DOW_HOD    = formula(kw ~ tout + DOW + HOD),
  standard   = formula(kw ~ tout + DOW + HOD + MOY),
  HOW        = formula(kw ~ tout + HOW),
  toutTOD    = formula(kw ~ tout:HOD + HOW)
)

daily.models = list(
  #MOY        = formula(kw.mean ~ tout.mean),
  tout_mean     = formula(kw.mean ~ tout.mean + DOW),
  tout_max      = formula(kw.mean ~ tout.max  + DOW),
  tout_CDD      = formula(kw.mean ~ CDD + HDD + DOW),
  tout_CDD_WKND = formula(kw.mean ~ CDD + HDD + WKND)
)

monthly.models = list(
  MOY        = formula(kw.mean ~ tout.mean + MOY),
  CDD        = formula(kw.mean ~ tout.mean + CDD),
  CDD_HDD    = formula(kw.mean ~ tout.mean + HDD + CDD)
)


regressorDF = function(residence,norm=TRUE) {
  hStr       = paste('H',sprintf('%02i',residence$dates$hour),sep='')
  dStr       = paste('D',residence$dates$wday,sep='')
  mStr       = paste('M',residence$dates$mon,sep='')
  howStrs    = paste(dStr,hStr,sep='')
  MOY = factor(mStr,levels=sort(unique(mStr)))       # month of year
  DOW = factor(dStr,levels=sort(unique(dStr)))       # day of week
  HOD = factor(hStr,levels=sort(unique(hStr)))       # hour of day
  HOW = factor(howStrs,levels=sort(unique(howStrs))) # hour of week
  tout = residence$w('tout')
  # todo: we need the names of the toutPIECES to build the model
  # but those names aren't returned form here
  # add special data to the data frame: piecewise tout data
  #toutPIECES = regressor.piecewise(r$tout,c(40,50,60,70,80,90))
  
  kw = residence$kw
  if(norm) kw = residence$norm(kw)
  df = data.frame(
      kw=kw,
      tout=tout,
      pout=residence$w('pout'),
      rain=residence$w('rain'),
      dates=residence$dates,
      wday=residence$dates$wday,
      MOY,DOW,HOD,HOW   )
  #df = cbind(df,toutPIECES) # add the columns with names from the matrix to the df
  return df
}

regressorDFAggregated = function(residence,norm=TRUE,bp=65) {
  # uses melt and cast to reshape and aggregate data
  df = residence$df() # kw, tout, dates
  if(norm) df$kw = residence$norm(df$kw)
  df$day   = format(df$dates,'%y-%m-%d') # melt has a problem with dates
  df$wday  = as.POSIXlt(df$dates)$wday   # raw for subsetting Su=0 ... Sa=6
  df$DOW   = paste('D',df$wday,sep='')  # Su=0 ... Sa=6
  df$WKND  = (df$wday == 0 | df$wday == 6) * 1 # weekend indicator
  df$DOW   = factor(df$DOW, levels=sort(unique(df$DOW)))
  month    = format(df$dates,'%y-%m')   # as.POSIXlt(df$dates)$mon
  df$mon   = as.POSIXlt(df$dates)$mon   # raw month data for subset functions Jan=0 ... Dec=11
  df$MOY   = factor(month, levels=sort(unique(month)))
  df <- subset(df, select = -c(dates) )  # melt has a problem with dates but  we don't need anymore
  
  # melt and cast to reshape data into monthly and daily time averages
  dfm = melt(df,id=c("day",'DOW','MOY','mon','wday','WKND'),na.rm=TRUE)
  monthly = cast(dfm,MOY + mon ~ variable,fun.aggregate=c(sum,mean,function(ar1) sum(ar1 > bp),function(ar2) sum(ar2 < bp)))
  colnames(monthly) <- c('month','mon','kwh','kw.mean','junk1','junk2','junk3','tout.mean','CDD','HDD')
  monthly <- subset(monthly, select = -c(junk1, junk2, junk3) )
  daily = cast(dfm, MOY + day + DOW + mon + wday + WKND ~ variable,fun.aggregate=c(sum,mean,max,function(ar1) sum(ar1 > bp),function(ar2) sum(ar2 < bp)))
  colnames(daily) <- c('MOY','day','DOW','mon','wday','WKND','kwh','kw.mean','kw.max','junk1','junk2','junk3','tout.mean','tout.max','CDD','HDD')
  daily <- subset(daily, select = -c(junk1, junk2, junk3) )
  return(list(daily,monthly))
}

runModelsBySP = function(sp_ids,zip=NULL,data=NULL,weather=NULL) {
  features.basic  <- c()
  ids             <- c()
  betas           <- c()
  summaries       <- c()
  
  # TODO:
  # use step for forward selection
  # find a faster method of cross validation than cvTools or DAAG provide
  # add other weather covariates.
  splen  <- length(sp_ids)
  i <- 0
  for (sp_id in sp_ids) { # i.e. "820735863"
    i <- i+1
    print(paste('  ',sp_id,' (',i,'/',splen,') in ',zip,sep=''))
    resData = NULL
    if (!is.null(data)) { 
      resData = data[data[,'sp_id']== sp_id,]
    }
    r <- tryCatch(ResDataClass(sp_id,zip=zip,weather=weather,data=resData,db=conf.meterDB()), 
                  error = function(e) {print(e)}, 
                  finally={} )
    if ( ! "ResDataClass" %in% class(r) ) {
      # do nothing: likely some sort of error, but not including it is sufficient
      print('    not found (see error)')
    }
    else { # it worked!
      r <- tryCatch( {
        tic('model run')
        df = regressorDF(r) # see also regressorDFAggregated
        models = base.models
        
        # add special data to the data frame: piecewise tout data
        toutPIECES = regressor.piecewise(r$tout,c(40,50,60,70,80,90))
        df = cbind(df,toutPIECES) # add the columns with names from the matrix to the df
        # define regression formula that uses the piecewise pieces
        models$toutPIECES = as.formula(paste('kw ~',paste(colnames(toutPIECES), collapse= "+"),'+ HOW'))
        
        WKND = (df$wday == 0 | df$wday == 6) # 0 is Sun, 6 is Sat
        WKDY = ! WKND
        # TODO: define other filter vectors, like seasonal
 
        features.basic  <- rbind(features.basic,basicFeatures(r$kwMat))
        ids             <- rbind(ids,sp_id)
        
        results   = lapply(models,lm, data=df)
        betas     = rbind(betas,lapply(results,coef))
        summaries = rbind(summaries,lapply(results,function(x) { summary(x)$sigma }))
        # TODO: write our own cross validation code - these are slow and picky! 
        #summaries = rbind(summaries,lapply(results,function(res) {cvFit(res,data=model.frame(res),y=model.frame(res)$kw)$cv}))
        #summaries = rbind(summaries,lapply(results,function(fit) {cv.lm(df=model.frame(fit),form.lm=fit,m=5,plotit=FALSE,printit=TRUE)$ss}))
        toc('model run',prefixStr='    model runs') 
      }, 
       error = function(e){
         print(e)
       },
       finally={
         # pass
       })
    }
    #break
  } # sp_id loop
  out <- list(
      features.basic  = features.basic,
      betas           = betas,
      ids             = ids,
      summaries       = summaries   )
  print(names(out))
  return(out)
}

# helper to format a summary of a batchRun
summarizeRun = function(runResult,listFailures=FALSE) {
  print('')
  print('----------- Batch run summary ------------')
  print('')
  toc('batchRun',prefixStr='Batch execution took ')
  print(paste('Zip code completion: ',length(runResult$completedZip),'/',length(runResult$attemptedZip)))
  print(paste('SP completion: ',length(runResult$completedSP),'/',length(runResult$attemptedSP)))
  if(listFailures) {
    print('Failed zips:')
    print(setdiff(runResult$attemptedZip,runResult$completedZip))
    print('Failed sp_ids:')
    print(setdiff(runResult$attemptedSP,runResult$completedSP))
  }
}

tic('batchRun')
allZips = c()
args = commandArgs(TRUE)
if (length(args) > 0) {
  print('Initializing batch run with command line zips:')
  allZips = args
  print(allZips)
} else { # no command line args so do a full run
  print('Initializing batch run with list of all zips')
  allZips  <- db.getZips()
}
# bakersfield, fresno, oakland
allZips = c('93304','93727','94610')
print('Beginning batch run')
runResult = runModelsByZip(allZips,triggerZip=94610)
summarizeRun(runResult,listFailures=FALSE)



