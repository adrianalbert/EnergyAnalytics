#!/usr/bin/Rscript

.libPaths('~/R/library') # use my local R library even from the comand line

# todo: is there a better way to detect the current directory?
conf.basePath = file.path('~/EnergyAnalytics/batch')
if(Sys.info()['sysname'] == 'Windows') {
  conf.basePath = file.path('c:/dev/pge_collab/pge_R')
}

# run 'source' on all includes to load them 
#source(file.path(conf.basePath,'localConf.R'))         # Sam's local computer specific configuration
source(file.path(conf.basePath,'stanfordConf.R'))     # Stanford on site specific configuration
source(file.path(conf.basePath,'DataClasses.R'))       # Object code for getting meter and weather data 
source(file.path(conf.basePath,'ksc.R'))               # k-Spectral Clustering (via Jungsuk)
source(file.path(conf.basePath,'basicFeatures.R'))     # typical max, min, mean, range
source(file.path(conf.basePath,'regressionSupport.R')) # mostly regressor manipulation
source(file.path(conf.basePath,'timer.R'))             # adds tic() and toc() functions

runModelsByZip = function(zipArray) {
  nZip     <- length(zipArray)
  zipCount <- 0
  outDir = 'test_results'
  
  print(paste('This batch process will run models for',nZip,'zip codes'))
  dir.create(file.path(conf.basePath,outDir),showWarnings=FALSE)
  print(paste('<zipcode>_modelResults.RData files for each can be found in',file.path(conf.basePath,outDir)))
  res = list(attemptedZip=zipArray,
             completedZip=c(),
             attemptedSP=c(),
             completedSP=c())
  triggerZip = '94610';
  triggered = TRUE;
  # oakland, sacramento, bakersfield, fresno
  #94610, 95823, 93304, 93727
  for (zip in zipArray) { # i.e. 94610
    if (zip == triggerZip) triggered = TRUE
    if (! triggered) {
      next
    }
    zipCount <- zipCount + 1
    print(paste('Running models for ',zip,' (',zipCount,'/',nZip,')', sep=''))
    resultsFile <- file.path(conf.basePath,outDir,paste(zip,'_modelResults.RData',sep=''))
    tryCatch( {
      tic('allDataForZip')
      zipData <- db.getAllData(zip)
      toc('allDataForZip')
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
    break
    
    
  } # zip loop
  return(res)
}

# allows for passing in weather and meter table data
# by looking up onece and running many times, this can save a lot of time
# if these data are not passed in, regular (slower) queries are performed
runModelsBySP = function(sp_ids,zip=NULL,data=NULL,weather=NULL) {
  coef.MOY        <- c()
  coef.DOW        <- c()
  coef.DOW_HOD    <- c()
  coef.standard   <- c()
  coef.HOW        <- c()
  coef.toutTOD    <- c()
  coef.toutPIECES <- c()
  features.basic  <- c()
  ids             <- c()
  summary         <- c()
  
  
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
        toutTOD    = NULL #regressor.split(r$tout,r$dates$hour)
        toutPIECES = regressor.piecewise(r$tout,c(40,50,60,70,80,90))
        hStr       = paste('H',sprintf('%02i',r$dates$hour),sep='')
        dStr       = paste('D',r$dates$wday,sep='')
        mStr       = paste('M',r$dates$mon,sep='')
        howStrs    = paste(dStr,hStr,sep='')
        MOY = factor(mStr,levels=sort(unique(mStr)))       # month of year
        DOW = factor(dStr,levels=sort(unique(dStr)))       # day of week
        HOD = factor(hStr,levels=sort(unique(hStr)))       # hour of day
        HOW = factor(howStrs,levels=sort(unique(howStrs))) # hour of week
        
        WKND = (r$dates$wday == 0 | r$dates$wday == 6) # 0 is Sun, 6 is Sat
        WKDY = ! WKND
        
        df = data.frame(MOY,DOW,HOD,HOW,toutTOD,toutPIECES,kw=r$norm(r$kw),tout=r$tout)
    
        tic('regression')
        model.MOY        <- NULL #rxLinMod(kw ~ tout       + MOY, data=df, verbose=0)
        model.DOW        <- NULL #rxLinMod(kw ~ tout       + DOW, data=df, verbose=0)
        model.DOW_HOD    <- rxLinMod(kw ~ tout       + DOW + HOD, data=df, verbose=0)
        model.standard   <- rxLinMod(kw ~ tout       + DOW + HOD + MOY, data=df, verbose=0)
        model.HOW        <- rxLinMod(kw ~ tout       + HOW, data=df, verbose=0)
        model.toutTOD    <- NULL #rxLinMod(kw ~ toutTOD    + HOW, data=df, verbose=0)
        model.toutPIECES <- rxLinMod(kw ~ toutPIECES + HOW, data=df, verbose=0)
        
        model.toutPIECES.WKND <- NULL #rxLinMod(kw ~ toutPIECES + HOW,subset=WKND, data=df, verbose=0)
        model.toutPIECES.WKDY <- NULL #rxLinMod(kw ~ toutPIECES + HOW,subset=WKDY, data=df, verbose=0)
        
        toc('regression',prefixStr='    regression')
        
        #coef.MOY        <- rbind(coef.MOY,coef(model.MOY))
        #coef.DOW        <- rbind(coef.DOW,coef(model.DOW))
        coef.DOW_HOD    <- rbind(coef.DOW_HOD,coef(model.DOW_HOD))
        
        coef.standard   <- rbind(coef.standard,coef(model.standard))
        coef.HOW        <- rbind(coef.HOW,coef(model.HOW))
        #coef.toutTOD    <- rbind(coef.toutTOD,coef(model.toutTOD))
        coef.toutPIECES <- rbind(coef.toutPIECES,coef(model.toutPIECES))
        features.basic  <- rbind(features.basic,basicFeatures(r$kwMat))
        ids             <- rbind(ids,sp_id)
        summary         <- rbind( summary,list(#MOY=summary(model.MOY)$kw$sigma,
                                               #DOW=summary(model.DOW)$kw$sigma,
                                               DOW_HOD=summary(model.DOW_HOD)$kw$sigma,
                                               standard=summary(model.standard)$kw$sigma,
                                               HOW=summary(model.HOW)$kw$sigma,
                                               #toutTOD=summary(model.toutTOD)$kw$sigma,
                                               toutPIECES=summary(model.toutPIECES)$kw$sigma,
                                               #toutPIECES_WKND=(summary(model.toutPIECES.WKND)$kw$sigma +
                                               #                 summary(model.toutPIECES.WKDY)$kw$sigma)/2    
                                               )
                                  )
        # dump the memory intensive parts
        rm(list = c('r','DOW','HOD','toutTOD','toutPIECES',
                    'model.MOY','model.DOW','model.DOW_HOD',
                    'model.standard','model.HOW','model.toutTOD',
                    'model.toutPIECES','model.toutPIECES.WKND','model.toutPIECES.WKDY'))
        #gc()
        toc('model run',prefixStr='    model run') 
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
  out <- list(coef.MOY        = coef.MOY,
              coef.DOW        = coef.DOW,
              coef.DOW_HOD    = coef.DOW_HOD,
              coef.standard   = coef.standard, 
              coef.HOW        = coef.HOW, 
              coef.toutTOD    = coef.toutTOD, 
              coef.toutPIECES = coef.toutPIECES, 
              features.basic  = features.basic,
              ids             = ids,
              summary         = summary)
  print(names(out))
  return(out)
}

# helper to format a summary of a batchRun
summarizeRun = function(runResult,listFailures=FALSE) {
  print('')
  print('----------- Batch run summary ------------')
  print('')
  toc('batchRun',prefixStr='Batch execution took (secs)')
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


print('Beginning batch run')
runResult = runModelsByZip(allZips)
summarizeRun(runResult,listFailures=FALSE)
toc('batchRun')



