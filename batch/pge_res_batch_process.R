#!/usr/bin/Rscript

.libPaths('~/R/library') # use my local R library even from the comand line

conf.basePath = file.path('~/EnergyAnalysis/batch')

# run 'source' on all includes to load them 
#source(file.path(conf.basePath,'localConf.R'))         # Sam's local computer specific configuration
source(file.path(conf.basePath,'stanfordConf.R'))     # Stanford on site specific configuration
source(file.path(conf.basePath,'DataClasses.R'))       # Object code for getting meter and weather data 
source(file.path(conf.basePath,'ksc.R'))               # k-Spectral Clustering (via Jungsuk)
source(file.path(conf.basePath,'basicFeatures.R'))     # typical max, min, mean, range
source(file.path(conf.basePath,'regressionSupport.R')) # mostly regressor manipulation
source(file.path(conf.basePath,'timer.R'))             # adds tic() and toc() functions

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
      sp_ids <- db.getSPs(zip)
      res$attemptedSP = c(res$attemptedSP,sp_ids)
      # we know that all sp's in the same zipcode share the same weather
      # so we speed execution by looking it up once and passing it in
      weather <- WeatherClass(zip)
      tic('modelsBySP')
      modelResults <- runModelsBySP(sp_ids,zip=zip,weather=weather)
      save(modelResults,file=resultsFile)
      res$completedZip <- rbind(res$completedZip,zip)
      res$completedSP  <- rbind(res$completedSP,modelResults$ids)
      rm(weather)
      rm(modelResults)
      toc('modelsBySP')
    }, 
    error = function(e) {
      print('runModelByZip')
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

runModelsBySP = function(sp_ids,zip=NULL,weather=NULL) {
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
    r <- tryCatch(ResDataClass(sp_id,zip,weather,conf.meterDB()), 
                  error = function(e) {
                    print('runModelBySP')
                    print(e)
                    next
                    }, 
                  finally={})
    if ( ! "ResDataClass" %in% class(r) ) {
      # do nothing: could be some sort of error, but not including it is sufficient
      print('    not found (see error)')
    }
    else { # it worked!
      r <- tryCatch( {
        
        tic()
        #tic('model run')
        hour = r$dates$hour
        #toutTOD = regressor.split(r$tout,hour)
        toutPIECES = regressor.piecewise(r$tout,c(40,50,60,70,80,90))
        hStr = paste('H',sprintf('%02i',r$dates$hour),sep='')
        dStr = paste('D',r$dates$wday,sep='')
        mStr = paste('M',r$dates$mon,sep='')
        howStrs = paste(dStr,hStr,sep='')
        MOY = factor(mStr,levels=sort(unique(mStr)))       # month of year
        DOW = factor(dStr,levels=sort(unique(dStr)))       # day of week
        HOD = factor(hStr,levels=sort(unique(hStr)))       # hour of day
        #HOW = factor(howStrs,levels=sort(unique(howStrs))) # hour of week
        
        #WKND = (r$dates$wday == 0 | r$dates$wday == 6) # 0 is Sun, 6 is Sat
        #WKDY = ! WKND
        df = data.frame(hour,MOY,DOW,HOD,toutPIECES,kw=r$norm(r$kw),tout=r$tout)

#         model.MOY        <- lm(r$norm(r$kw) ~ r$tout     + MOY)
#         model.DOW        <- lm(r$norm(r$kw) ~ r$tout     + DOW)
#         model.DOW_HOD    <- lm(r$norm(r$kw) ~ r$tout     + DOW + HOD)
#         model.standard   <- lm(r$norm(r$kw) ~ r$tout     + DOW + HOD + MOY)
#         model.HOW        <- lm(r$norm(r$kw) ~ r$tout     + HOW)
#         model.toutTOD    <- lm(r$norm(r$kw) ~ toutTOD    + HOW)
#         model.toutPIECES <- lm(r$norm(r$kw) ~ toutPIECES + HOW)
#         
#         model.toutPIECES.WKND <- lm(r$norm(r$kw) ~ toutPIECES + HOW,subset=WKND)
#         model.toutPIECES.WKDY <- lm(r$norm(r$kw) ~ toutPIECES + HOW,subset=WKDY)
        
        model.MOY        <- NULL #rxLinMod(kw ~ tout       + MOY, data=df, verbose=0, reportProgress=0)
        model.DOW        <- NULL #rxLinMod(kw ~ tout       + DOW, data=df, verbose=0, reportProgress=0)
        model.DOW_HOD    <- rxLinMod(kw ~ tout       + DOW + HOD, data=df, verbose=0, reportProgress=0)
        model.standard   <- rxLinMod(kw ~ tout       + DOW + HOD + MOY, data=df, verbose=0, reportProgress=0)
        model.HOW        <- rxLinMod(kw ~ tout       + DOW:HOD, data=df, verbose=0, reportProgress=0)
        model.toutTOD    <- NULL #rxLinMod(kw ~ tout:hour  + DOW:HOD, data=df, verbose=0, reportProgress=0)
        model.toutPIECES <- rxLinMod(kw ~ tout0_40 + tout40_50 + tout50_60 + tout60_70 + tout70_80 + tout80_90 + tout90_Inf + DOW:HOD, data=df, verbose=0, reportProgress=0)
        
        model.toutPIECES.WKND <- NULL #rxLinMod(kw ~ toutPIECES + DOW:HOD,subset=WKND, data=df, verbose=0, reportProgress=0)
        model.toutPIECES.WKDY <- NULL #rxLinMod(kw ~ toutPIECES + DOW:HOD,subset=WKDY, data=df, verbose=0, reportProgress=0)
        
        coef.MOY        <- rbind(coef.standard,coef(model.MOY))
        coef.DOW        <- rbind(coef.standard,coef(model.DOW))
        coef.DOW_HOD    <- rbind(coef.standard,coef(model.DOW_HOD))
        coef.standard   <- rbind(coef.standard,coef(model.standard))
        coef.HOW        <- rbind(coef.HOW,coef(model.HOW))
        coef.toutTOD    <- rbind(coef.toutTOD,coef(model.toutTOD))
        coef.toutPIECES <- rbind(coef.toutPIECES,coef(model.toutPIECES))
        features.basic  <- NULL #rbind(features.basic,basicFeatures(r$kwMat))
        ids             <- rbind(ids,sp_id)
        summary         <- rbind( summary,list(#MOY=summary(model.MOY)$kw$sigma,
                                               #DOW=summary(model.DOW)$kw$sigma,
                                               DOW_HOD=summary(model.DOW_HOD)$kw$sigma,
                                               standard=summary(model.standard)$sigma,
                                               HOW=summary(model.HOW)$kw$sigma,
                                               #toutTOD=summary(model.toutTOD)$kw$sigma
                                               toutPIECES=summary(model.toutPIECES)$kw$sigma
                                               #toutPIECES_WKND=(summary(model.toutPIECES.WKND)$kw$sigma +
                                               #                summary(model.toutPIECES.WKDY)$kw$sigma)/2    
                                               )
                                  )
        # dump the memory intensive parts
        #rm(list = c('r','DOW','HOD','toutTOD','toutPIECES',
        #            'model.MOY','model.DOW','model.DOW_HOD',
        #            'model.standard','model.HOW','model.toutTOD',
        #            'model.toutPIECES','model.toutPIECES.WKND','model.toutPIECES_WKDY'))
        #gc()
        toc()
        #toc('model run',prefixStr='    model run') 
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
#allZips = c('93304','93727','94610')
print('Beginning batch run')
runResult = runModelsByZip(allZips,triggerZip=NULL)
summarizeRun(runResult,listFailures=FALSE)



