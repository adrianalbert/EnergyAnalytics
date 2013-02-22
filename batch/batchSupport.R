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
             invalid.ids  = c() )
  triggered = FALSE
  for (zip in zipArray) { # i.e. 94610
    zipCount <- zipCount + 1
    if (is.null(cfg$triggerZip) || zip == cfg$triggerZip) triggered = TRUE
    if (! triggered) {
      next
    }
    
    print(paste('Running models for ',zip,' (',zipCount,'/',nZip,')', sep=''))
    resultsFile <- file.path(getwd(),cfg$outDir,paste(zip,'_modelResults.RData',sep=''))
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


# summarize regression model for future use
# designed to be friendly to storing a lot of results in sequence
# so it separates out the simple scalar metrics from the more complicated
# coefficients
summarizeModel = function(m,df,models,nm,id,zip,subnm=NULL,fold=FALSE,formula='',subset='') {
  #lm(m,subset=m$y > 1)
  s <- as.list(summary(m, correlation=FALSE)) # generic list is more friendly for adding to a data.frame
  class(s) <- 'list'         # make sure the class is no longer summary.lm
  s$call          <- c()     # lm model call (depends on variable scope and can be junk)
  s$terms         <- c()     # lm model terms (depends on variable scope and can be junk)
  s$residuals     <- c()     # residuals scaled by weights assigned to the model
  s$cov.unscaled  <- c()     # p x p matrix of (unscaled) covariances
  #s$aliased      <- c()     # named logical vector showing if the original coefficients are aliased
  s$na.action     <- c()     # get rid of extra meta info from the model
  #s$sigma                   # the square root of the estimated variance of the random error
  #s$adj.r.squared           # penalizing for higher p
  #s$r.squared               # 'fraction of variance explained by the model'
  #s$fstatistic              # 3-vector with the value of the F-statistic with its numerator and denominator degrees of freedom
  #s$df                      # degs of freedom, vector (p,n-p,p*), the last being the number of non-aliased coefficients.
  #s$coefficients            # a p x 4 matrix: cols = coefficient, standard error, t-stat, (two-sided) p-value
  s$id          <- id
  s$zip         <- zip
  s$model.name  <- nm        # name of model run
  s$subset.name <- subnm     # name of sample subset criteria (i.e. "summer" or "afternoon" or "weekdays")
  s$formula     <- formula   # string of lm model call
  s$subset      <- subset    # string of lm subset argument
  s$logLik      <- logLik(m) # log liklihood for the model
  s$AIC         <- AIC(m)    # Akaike information criterion
  
  # net contribution of each coefficient
  # todo: there are a lot of negative numbers in this, so the contributions aren't directly interpretable
  s$contribution <- colSums(t(m$coefficients * t(m$x)))
  s$total        <- sum(predict(m))
  
  # k-fold prediction error
  if (fold) {
    s$fold.rmse <- kFold(df,models,nm,nfolds=5)
  }
  return(s)
}

runModelsBySP = function(sp_ids,cfg,zip=NULL,data=NULL,weather=NULL) {
  features.basic  <- c()
  summaries       <- c()
  d_summaries     <- c()
  m_summaries     <- c()
  steps           <- c()
  
  invalid_ids     <- data.frame()
  # TODO:
  # find a faster method of cross validation than cvTools or DAAG provide
  splen  <- length(sp_ids)
  i <- 0
  #skip = FALSE
  for (sp_id in sp_ids) { # i.e. "820735863"
    i <- i+1
    print(paste('  ',sp_id,' (',i,'/',splen,') in ',zip,sep=''))
    #if(sp_id == 6502353905) { skip = FALSE }
    #if(skip) { next }
    resData = NULL
    if (!is.null(data)) {  resData = data[data[,'sp_id']== sp_id,] }
    r <- tryCatch(ResDataClass(sp_id,zip=zip,weather=weather,data=resData,db=conf.meterDB()), 
                  error = function(e) {print(e)}, finally={} )
    if ( ! "ResDataClass" %in% class(r) ) { # constructor returns the error class if it has a problem
      print('    not found (see error)') # ignore residences that produce errors
    }
    else { # viable residence!
      r <- tryCatch( {
        issues = validateRes(r)
        savePlot = FALSE
        if(length(issues)>1) # all issues will return with an id set, so > 1 is a problem
        { 
          invalid_ids <- rbind.fill(invalid_ids,issues)
          print(paste('Bad or insufficient data:',paste(colnames(issues),collapse=', ')))
          if(cfg$PLOT_INVALID) {
            png(file.path(getwd(),cfg$outDir,paste(r$zip,r$id,'invalid.png',sep='_')))
            colorMap = rev(colorRampPalette(brewer.pal(11,"RdBu"))(100))
            plot(  r, colorMap=colorMap, 
                   main=paste(r$zip, r$id, paste(colnames(issues),collapse=', ')) )
            dev.off()
          }
          next # no further processing
        }
        if(cfg$PLOT_VALID) {
          png(file.path(getwd(),cfg$outDir,paste(r$zip,'_',r$id,'.png',sep='')))
          redblue = rev(colorRampPalette(brewer.pal(11,"RdBu"))(100))
          plot( r,colorMap=redblue,main=paste(r$zip, r$id) )
          dev.off()
        }
        tic('model run')
        features.basic  <- rbind(features.basic,basicFeatures(r$kwMat,sp_id))
        # hourly regressions
        if(TRUE) {
          df = regressorDF(r,norm=FALSE) # see also regressorDFAggregated
          models = cfg$models.hourly
          if (FALSE) {
            # add special data to the data frame: piecewise tout data
            #toutPIECES = regressor.piecewise(r$tout,c(40,50,60,70,80,90))
            toutPIECES = regressor.piecewise(r$tout,c(60,75))
            df = cbind(df,toutPIECES) # add the columns with names from the matrix to the df
            # define regression formula that uses the piecewise pieces
            piecesf = paste('kw ~',paste(colnames(toutPIECES), collapse= "+"),'+ HOW')
            models$toutPIECES = list(formula=piecesf,subset=list(all="TRUE",summer=cfg$subset$summer))
          }
          if (TRUE) {
            # add special data to the data frame: piecewise tout data
            #toutPIECES = regressor.piecewise(r$tout,c(40,50,60,70,80,90))
            toutPIECES = regressor.piecewise(r$tout,c(55,65,75))
            df = cbind(df,toutPIECES) # add the columns with names from the matrix to the df
            # define regression formula that uses the piecewise pieces
            piecesf24 = paste('kw ~',paste(colnames(toutPIECES),':HOD',collapse=" + ",sep=''),'+ HOD - 1')
            models$toutPIECES24 = list(formula=piecesf24,subset=list(all="TRUE"))          
          }
          # careful. lapply doesn't pass the right data to re-use the lm model
          # that's why summarizeModel gets passed everything it needs to repeat the lm
          # with subsets...
          # see https://stat.ethz.ch/pipermail/r-help/2007-August/138724.html
          for(nm in names(models)) {
            fld = FALSE
            #print(nm) 
            fmla = models[[nm]]
            subs = list(all="TRUE") # this will include all obs
            if (class(fmla) == 'list') { # a list object will contain a formula and list of named subset commands
              subs = fmla$subset  # list of named subset code snippets
              fmla = fmla$formula
            }
            #if(class(fmla) == 'function') {
            #  # call the function assuming that it will pass back exactly what we need
            #}
            for(snm in names(subs)) {
              df$sub = eval(parse(text=subs[[snm]]),envir=df) # load the subset flags into the data.frame
              lm.result = lm(fmla,df, x=TRUE, subset=sub==TRUE) # run the lm on the subset indicated in the data frame
              summaries = rbind(summaries,
                                summarizeModel(lm.result,df,models,nm,
                                               id=sp_id,zip=zip,subnm=snm,fold=fld,
                                               formula=fmla,subset=subs[[snm]]) )
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

        if (cfg$RUN_AGGREGATED_MODELS) {
          dfl = regressorDFAggregated(r,norm=FALSE)
          
          # daily regressions
          models = cfg$models.daily
          for(nm in names(models)) {
            #print(nm)
            fmla = models[[nm]]
            lm.result = lm(fmla,dfl$daily, x=TRUE)
            d_summaries = rbind(d_summaries,
                                summarizeModel(lm.result,
                                               dfl$daily,models,nm,
                                               id=sp_id,zip=zip,subnm="all",
                                               fold=FALSE,formula=fmla)  )
          }
          #results = lapply(models,FUN=lm,dfl$daily,x=TRUE)
          #d_summaries = lapply(results,FUN=summarizeModel,dfl$daily,models,nm,id=sp_id,zip=zip,fold=FALSE)
          
          # monthly regressions
          models = cfg$models.monthly     
          for(nm in names(models)) {
            #print(nm)
            fmla = models[[nm]]
            lm.result = lm(fmla,dfl$monthly, x=TRUE)
            m_summaries = rbind(m_summaries,summarizeModel(lm.result,dfl$monthly,models,nm,
                                                           id=sp_id,zip=zip,subnm="all",
                                                           fold=FALSE,formula=fmla)  )
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
  out <- list(
      features.basic  = as.data.frame(features.basic),
      summaries       = as.data.frame(summaries),
      d_summaries     = as.data.frame(d_summaries),
      m_summaries     = as.data.frame(m_summaries),
      steps           = as.data.frame(steps),
      invalid.ids     = invalid_ids
      )
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
