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
source(file.path(getwd(),'ksc.R'))               # k-Spectral Clustering (via Jungsuk)
source(file.path(getwd(),'basicFeatures.R'))     # typical max, min, mean, range
source(file.path(getwd(),'regressionSupport.R')) # mostly regressor manipulation 
source(file.path(getwd(),'solaRUtil.R'))         # solar geometry
source(file.path(getwd(),'timer.R'))             # adds tic() and toc() functions

library(reshape)
library(timeDate)
library(RColorBrewer)
# run this if you need it, but everything should be installed by a setup script
# if (!require("RColorBrewer")) { install.packages("RColorBrewer") }

#library('DAAG') # Data Analysis and Graphics package has k-fold cross validation
# cv.lm(df=mydata, model, m=5) # 5 fold cross-validation
#library('cvTools') # cross validation tools


runModelsByZip = function(cfg) {
  zipArray <- cfg$allZips
  nZip     <- length(zipArray)
  zipCount <- 0
  
  print(paste('This batch process will run models for',nZip,'zip codes'))
  dir.create(file.path(getwd(),cfg$outDir),showWarnings=FALSE)
  print(paste('<zipcode>_modelResults.RData files for each can be found in',file.path(getwd(),cfg$outDir)))
  res = list(attemptedZip = zipArray,
             completedZip = c(),
             attemptedSP  = c(),
             completedSP  = c(),
             invalid.ids  = c(),
             cfg          = cfg )
  triggered = FALSE
  for (zip in zipArray) { # i.e. 94610
    zipCount <- zipCount + 1
    if (is.null(cfg$triggerZip) || zip == cfg$triggerZip) triggered = TRUE
    if (! triggered) {
      next
    }
    
    print(paste('Running models for ',zip,' (',zipCount,'/',nZip,')', sep=''))
    resultsFile <- file.path(getwd(),cfg$outDir,paste(zip,'_modelResults.RData',sep=''))
    if(file.exists(resultsFile) & cfg$SKIP_EXISTING_RDATA) {
      print(paste(resultsFile,'exists. Skipping to next zip.'))
      next
    }
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
      modelResults <- runModelsBySP(sp_ids,cfg,zip=zip,data=zipData,weather=weather)
      rm(zipData,weather,sp_ids)
      save(modelResults,file=resultsFile)
      res$completedZip <- rbind(res$completedZip,zip)
      res$completedSP  <- c(res$completedSP,modelResults$features.basic$id)
      res$invalid.ids  <- rbind.fill(res$invalid.ids,modelResults$invalid.ids)
      rm(modelResults)
      toc('modelsBySP')
    }, 
    error = function(e) {
      print(e)
      modelResults = NA
      # put an empty file in place to differentiate between un-run zips and unsuccessful runs
      save(modelResults,e,file=resultsFile)   }, 
    finally={})
  } # zip loop
  return(res)
}

# TODO: test rh, pout regressors
# TODO: define junk shot model for use with "step"
# TODO: aks Ram about his k-fold idea and look at serial folds as well as random
# TODO: consolidate after each model run, before serializing to disk

runModelsBySP = function(sp_ids,cfg,zip=NULL,data=NULL,weather=NULL) {
  naList = as.list(rep(NA,length(sp_ids)))
  features.basic  <- naList
  summaries       <- naList
  others          <- list() # special model specific data not found in the summaries
  d_summaries     <- naList
  d_others        <- list() # special model specific data not found in the summaries
  m_summaries     <- naList
  m_others        <- list() # special model specific data not found in the summaries
  steps           <- c()  
  invalid_ids     <- data.frame()
  inputs          <- list(zip=zip,ids=sp_ids,config=cfg)
  
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
  # TODO:
  # find a faster method of cross validation than cvTools or DAAG provide
  splen  <- length(sp_ids)
  i <- 0
  #skip = FALSE
  for (sp_id in sp_ids) { # i.e. "820735863" or "6502353905"
    i <- i+1
    print(paste('  ',sp_id,' (',i,'/',splen,') in ',zip,sep=''))
    #if(sp_id == 6502353905) { skip = FALSE } # debug shortcut
    #if(skip) { next }
    resData = NULL
    if (!is.null(data)) {  resData = data[data[,'sp_id']== sp_id,] }
    r <- tryCatch(ResDataClass(sp_id,zip=zip,weather=weather,data=resData,db=conf.meterDB()), 
                  error = function(e) {print(e)}, finally={} )
    if ( ! "ResDataClass" %in% class(r) ) { # constructor returns the error class if it has a problem
      print('    not found (see error)')    # ignore residences that produce errors & continue
    }
    else { # viable residence!
      r <- tryCatch( {
        issues = validateRes(r)
        if(length(issues)>1) # all issues will return with an id, so > 1 indicates a problem
        { 
          invalid_ids <- rbind.fill(invalid_ids,issues)
          print(paste('Bad or insufficient data:',paste(colnames(issues),collapse=', ')))
          if(cfg$PLOT_INVALID) {
            save.png.plot(r,file.path(getwd(),cfg$outDir,paste(r$zip,r$id,'invalid.png',sep='_')))
          }
          next # no further processing of invalid sp_ids
        }
        
        tic('model run')
        if(cfg$PLOT_VALID) {
          save.png.plot(r,file.path(getwd(),cfg$outDir,paste(r$zip,'_',r$id,'.png',sep='')))
        }
        basics = basicFeatures(r)
        
        tic()
        features.basic[[i]]  <- basics # max, min, etc.
        toc(prefix='features.basic')
        # hourly regressions
        if(cfg$RUN_HOURLY_MODELS) {
          df = regressorDF(r,norm=FALSE) # see also regressorDFAggregated
          # careful. lapply doesn't pass the right data to re-use the lm model
          # see https://stat.ethz.ch/pipermail/r-help/2007-August/138724.html
          
          for(mdName in names(cfg$models.hourly)) {
            #print(mdName)
            md = cfg$models.hourly[[mdName]]
            runOut = md$run(r,df)
            summaries[[i]] = runOut$summaries
            if(! empty(runOut$other)) {
              others[[md$name]] = rbind(others[[md$name]],list(id=r$id,data=runOut$other))
            }
          }
          if(cfg$RUN_STEP_SELECTION) { 
            df = regressorDF(r,norm=FALSE,rm.na=TRUE)
            stepped = step(lm(kw ~ 0,df),direction='forward',k=2,trace=0,
                           scope=kw ~ 0 + tout65:HOD + pout + rh + tout65_d1 + HOD)
            ano = stepped$anova
            matchOrder = ano[,'Step']
            matchOrder = matchOrder[matchOrder != ""]
            row = data.frame(t(matchOrder != ""))
            row[1,] = 1:length(matchOrder)
            colnames(row) <- matchOrder # b is now ready to be rbound with regressors as col heads.
            row$id = r$id
            steps = rbind.fill(steps,row)
          }
        }
        if (cfg$RUN_DAILY_MODELS) {
          #dfl = regressorDFAggregated(r,norm=FALSE)
          #dfd=dfl$daily
          dfd = rDFA(r)
          # daily regressions
          for(mdName in names(cfg$models.daily)) {
            md = cfg$models.daily[[mdName]]
            #print(mdName)
            runOut = md$run(r,dfd)
            tic()
            d_summaries[[i]] = runOut$summaries
            toc(prefix='rbind d_summaries')
            if(! empty(runOut$other)) {
              d_others[[md$name]] = rbind(d_others[[md$name]],list(id=r$id,data=runOut$other))
            }
          }
        }
        if (cfg$RUN_MONTHLY_MODELS) {
          # monthly regressions
          dfl = regressorDFAggregated(r,norm=FALSE)
          models = cfg$models.monthly     
          for(nm in names(models)) {
            #print(nm)
            fmla = models[[nm]]
            lm.result = lm(fmla,dfl$monthly, x=TRUE)
            m_summaries[[i]] = summarizeModel(lm.result,dfl$monthly,models,nm,
                                                           id=sp_id,zip=zip,subnm="all",
                                                           fold=FALSE,formula=fmla)
          }
        }
        # TODO: write our own cross validation code - these are slow and picky! 
        #summaries = rbind(summaries,lapply(results,function(lm.result) {cvFit(lm.result,data=model.frame(lm.result),y=model.frame(lm.result)$kw)$cv}))
        #summaries = rbind(summaries,lapply(results,function(fit) {cv.lm(df=model.frame(fit),form.lm=fit,m=5,plotit=FALSE,printit=TRUE)$ss}))
        toc('model run',prefixStr='    model runs') 
      }, 
       error = function(e){
         print(e)
         traceback()
       },
       finally={
         # pass
       })
    }
    # this allows for a config directive to halt execution after a cetain number
    # of model fits. If cfg$truncateAt is -1, all fits run. If it is a positive 
    # integer, only that many run. Note that i is only incremented for residences
    # that pass validation, so this is the truncatAt directive is the count of 
    # successful fits
    if (cfg$truncateAt > 0 & i >= cfg$truncateAt) break
  } # sp_id loop
  # strip out NAs from the lists
  fbArray        <- features.basic[! is.na(features.basic)] 
  summaryArray   <- summaries[! is.na(summaries)]
  d_summaryArray <- d_summaries[! is.na(d_summaries)]
  m_summaryArray <- m_summaries[! is.na(m_summaries)]
  out <- list(
      inputs          = inputs,
      features.basic  = as.data.frame(do.call(rbind,fbArray)),
      others          = others,
      summaries       = as.data.frame(do.call(rbind,summaryArray)),
      d_others        = d_others,
      d_summaries     = as.data.frame(do.call(rbind,d_summaryArray)),
      m_others        = m_others,
      m_summaries     = as.data.frame(do.call(rbind,m_summaryArray)),
      steps           = as.data.frame(steps),
      invalid.ids     = invalid_ids
      )
  class(out) <- 'modelResults'
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
  print(paste('Validation errors: ',dim(runResult$invalid.ids)[1],'/',length(runResult$attemptedSP)))
  if(listFailures) {
    print('Failed zips:')
    print(setdiff(runResult$attemptedZip,runResult$completedZip))
    print('Failed sp_ids:')
    print(setdiff(runResult$attemptedSP,runResult$completedSP))
  }
}
