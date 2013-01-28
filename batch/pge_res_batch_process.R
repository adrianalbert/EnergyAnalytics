#!/usr/bin/Rscript

# todo: is there a better way to detect the current directory?
conf.basePath = file.path('~/EnergyAnalytics/batch')
if(Sys.info()['sysname'] == 'Windows') {
  conf.basePath = file.path('f:/dev/pge_collab/EnergyAnalytics/batch')
} else {
  .libPaths('~/R/library') # use my local R library even from the comand line
}
# run 'source' on all includes to load them 
source(file.path(conf.basePath,'localConf.R'))         # Sam's local computer specific configuration
#source(file.path(conf.basePath,'stanfordConf.R'))     # Stanford on site specific configuration
source(file.path(conf.basePath,'DataClasses.R'))       # Object code for getting meter and weather data 
source(file.path(conf.basePath,'ksc.R'))               # k-Spectral Clustering (via Jungsuk)
source(file.path(conf.basePath,'basicFeatures.R'))     # typical max, min, mean, range
source(file.path(conf.basePath,'regressionSupport.R')) # mostly regressor manipulation
source(file.path(conf.basePath,'timer.R'))             # adds tic() and toc() functions

library(reshape)
#library('DAAG') # Data Analysis and Graphics package has k-fold cross validation
# cv.lm(df=mydata, model, m=5) # 5 fold cross-validation
#library('cvTools') # cross validation tools

runModelsByZip = function(zipArray,triggerZip=NULL,truncateAt=-1) {
  nZip     <- length(zipArray)
  zipCount <- 0
  outDir = 'results_summer'
  
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
      modelResults <- runModelsBySP(sp_ids,zip=zip,data=zipData,weather=weather,truncateAt=truncateAt)
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

# TODO: test rh, pout regressors
# TODO: define junk shot model for use with "step"
# TODO: define summer only runs
# TODO: aks Ram about his k-fold idea and look at serial folds as well as random
# TODO: normalize the output format for model runs so that consolidate works again
# TODO: develop consolidate across multiple zips
# TODO: consolidate after each model run, before serializing to disk

# passing in forula objects directly creates lots of problems
# formulas are context specific and cant be stored in dataframes
# to get formula from string, call as.formula(str)
# to get a string repr of a formula object, call deparse(fmla)
subset = list(
  summer=paste("MOY %in% c(",paste("'M",6:9,"'",sep='',collapse=','),")"),
  day=paste("HOD %in% c(",paste("'H",8:19,"'",sep='',collapse=','),")"),
  summer_day=paste(summer,'&',day)
)

models.hourly = list(
  MOY        = "kw ~ tout + MOY",             
  DOW        = "kw ~ tout + DOW",             
  DOW_prh    = "kw ~ tout + pout + rh + DOW", 
  DOW_HOD    = "kw ~ tout + DOW + HOD",
  standard   = "kw ~ tout + DOW + HOD + MOY",
  HOW        = "kw ~ tout + HOW",
  toutTOD    = list(formula="kw ~ tout:HOD + HOW",subset=list(all="TRUE",summer=subset$summer,day=subset$day)),
  toutTOD_min = "kw_min ~ 0 + tout:HOD + HOW" # no intercept
)

models.daily = list(
  tout           = "kwh ~ tout.mean",
  DOW            = "kwh ~ DOW",
  tout_mean      = "kwh ~ tout.mean + DOW",
  tout_mean_WKND = "kwh ~ tout.mean + WKND",
  tout_max       = "kwh ~ tout.max  + DOW",
  tout_CDD       = "kwh ~ CDD + HDD + DOW",
  tout_CDD_WKND  = "kwh ~ CDD + HDD + WKND"
)

models.monthly = list(
  tout       = "kwh ~ tout.mean",
  CDD        = "kwh ~ CDD",
  CDD_HDD    = "kwh ~ HDD + CDD"
)


regressorDF = function(residence,norm=TRUE,folds=1) {
  WKND       = (residence$dates$wday == 0 | residence$dates$wday == 6) * 1 # weekend indicator
  hStr       = paste('H',sprintf('%02i',residence$dates$hour),sep='')
  dStr       = paste('D',residence$dates$wday,sep='')
  mStr       = paste('M',residence$dates$mon,sep='')
  howStrs    = paste(dStr,hStr,sep='')
  MOY = factor(mStr,levels=sort(unique(mStr)))       # month of year
  DOW = factor(dStr,levels=sort(unique(dStr)))       # day of week
  HOD = factor(hStr,levels=sort(unique(hStr)))       # hour of day
  HOW = factor(howStrs,levels=sort(unique(howStrs))) # hour of week
  tout = residence$w('tout')
  pout = residence$w('pout')
  rain = residence$w('rain')
  dp   = residence$w('dp')
  rh   = residence$weather$rh(tout,dp)
  # todo: we need the names of the toutPIECES to build the model
  # but those names aren't returned form here
  # add special data to the data frame: piecewise tout data
  #toutPIECES = regressor.piecewise(r$tout,c(40,50,60,70,80,90))
  
  kw = residence$kw
  if(norm) kw = residence$norm(kw)
  kw_min = kw - quantile(kw,na.rm=TRUE,c(0.02)) # remove the min for regression w/o const term
  
  df = data.frame(
      kw=kw,
      kw_min=kw_min,
      tout=tout,
      pout=pout,
      rain=rain,
      rh=rh,
      dates=residence$dates,
      wday=residence$dates$wday,
      MOY,DOW,HOD,HOW,WKND   )
  #df = cbind(df,toutPIECES) # add the columns with names from the matrix to the df
  return(df)
}

regressorDFAggregated = function(residence,norm=TRUE,bp=65) {
  # uses melt and cast to reshape and aggregate data
  df = residence$df() # kw, tout, dates
  if(norm) df$kw_norm = residence$norm(df$kw)
  df$kw_min = df$kw - quantile(df$kw,na.rm=TRUE,c(0.02)) # remove the min for regression w/o const term
  
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
  monthly = cast(dfm,MOY + mon ~ variable,fun.aggregate=c(sum,mean,function(ar1,bp=65) sum(ar1 > bp),function(ar2,bp=65) sum(ar2 < bp)))
  colnames(monthly) <- c('MOY','mon','kwh','kw.mean','junk1','junk2','junk3','tout.mean','CDD','HDD')
  monthly <- subset(monthly, select = -c(junk1, junk2, junk3) )
  daily = cast(dfm, MOY + day + DOW + mon + wday + WKND ~ variable,fun.aggregate=c(sum,mean,max,function(ar1,bp=65) sum(ar1 > bp),function(ar2,bp=65) sum(ar2 < bp)))
  colnames(daily) <- c('MOY','day','DOW','mon','wday','WKND','kwh','kw.mean','kw.max','junk1','junk2','junk3','tout.mean','tout.max','CDD','HDD')
  daily <- subset(daily, select = -c(junk1, junk2, junk3) )
  return(list(daily=daily,monthly=monthly))
}

kFold = function(df,models,nm,nfolds=5) {
  folds <- sample(1:nfolds, dim(df)[1], replace=T)
  residuals = c()
  for (i in 1:nfolds) {
    fld = folds == i
    df$fold = fld
    subm = lm(models[[nm]], data=df, subset=(!fold),na.action=na.omit)
    yhat = predict(subm,newdata=df[fld,])
    #print(length(yhat))
    ynm = as.character(models[[nm]])[[2]]
    residuals = c(residuals,df[,ynm][fld] - yhat) # accumulate the errors from all predictions
  }
  #plot(residuals)
  rmnsqerr = sqrt(mean(residuals^2,na.rm=TRUE)) # RMSE
  return(rmnsqerr)
}

# summarize regression model for future use
# designed to be friendly to storing a lot of results in sequence
# so it separates out the simple scalar metrics from the more complicated
# coefficients
summarizeModel = function(m,df,models,nm,id,zip,subnm=NULL,fold=FALSE,formula='',subset='') {
  #lm(m,subset=m$y > 1)
  s <- as.list(summary(m, correlation=FALSE)) # generic list is more friendly for adding to a data.frame
  s$call          <- c()     # lm model call (depends on variable scope and can be junk)
  s$terms         <- c()     # lm model terms (depends on variable scope and can be junk)
  s$residuals     <- c()     # residuals scaled by weights assigned to the model
  s$cov.unscaled  <- c()     # p x p matrix of (unscaled) covariances
  s$aliased       <- c()     # named logical vector showing if the original coefficients are aliased
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

summerlm = function(fmla,df,x) {
  return(lm(fmla,df,x=x,subset=MOY %in% c('M6','M7','M8','M9')))
}

runModelsBySP = function(sp_ids,zip=NULL,data=NULL,weather=NULL,truncateAt=-1) {
  features.basic  <- c()
  summaries       <- c()
  d_summaries     <- c()
  m_summaries     <- c()
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
    if (!is.null(data)) {  resData = data[data[,'sp_id']== sp_id,] }
    r <- tryCatch(ResDataClass(sp_id,zip=zip,weather=weather,data=resData,db=conf.meterDB()), 
                  error = function(e) {print(e)}, finally={} )
    if ( ! "ResDataClass" %in% class(r) ) { # constructor returns the error class if it has a problem
      print('    not found (see error)') # ignore residences that produce errors
    }
    else { # viable residence!
      r <- tryCatch( {
        tic('model run')
        features.basic  <- rbind(features.basic,basicFeatures(r$kwMat,sp_id))
        
        # hourly regressions
        if(TRUE) {
          df = regressorDF(r,norm=FALSE) # see also regressorDFAggregated
          models = models.hourly
          
          # add special data to the data frame: piecewise tout data
          toutPIECES = regressor.piecewise(r$tout,c(40,50,60,70,80,90))
          df = cbind(df,toutPIECES) # add the columns with names from the matrix to the df
          # define regression formula that uses the piecewise pieces
          piecesf = paste('kw ~',paste(colnames(toutPIECES), collapse= "+"),'+ HOW')
          models$toutPIECES = list(formula=piecesf,subset=list(all="TRUE",summer=subset$summer,day=subset$day)),
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
            for(snm in names(subs)) {
              df$sub = eval(parse(text=subs[[snm]]),envir=df)
              res = lm(fmla,df, x=TRUE, subset=sub==TRUE)
              summaries = rbind(summaries,summarizeModel(res,df,models,nm,id=sp_id,zip=zip,subnm=snm,fold=fld,formula=fmla,subset=subs[[snm]]))
            }
          }
        }

        if (TRUE) {
          dfl = regressorDFAggregated(r,norm=FALSE)
          
          # daily regressions
          models = models.daily
          for(nm in names(models)) {
            #print(nm)
            fmla = models[[nm]]
            res = lm(fmla,dfl$daily, x=TRUE)
            d_summaries = rbind(d_summaries,summarizeModel(res,dfl$daily,models,nm,id=sp_id,zip=zip,subnm="all",fold=FALSE,formula=fmla))
          }
          #results = lapply(models,FUN=lm,dfl$daily,x=TRUE)
          #d_summaries = lapply(results,FUN=summarizeModel,dfl$daily,models,nm,id=sp_id,zip=zip,fold=FALSE)
          
          # monthly regressions
          models = models.monthly     
          for(nm in names(models)) {
            #print(nm)
            fmla = models[[nm]]
            res = lm(fmla,dfl$monthly, x=TRUE)
            m_summaries = rbind(m_summaries,summarizeModel(res,dfl$monthly,models,nm,id=sp_id,zip=zip,subnm="all",fold=FALSE,formula=fmla))
          }
        }
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
    if (i == truncateAt) break
  } # sp_id loop
  out <- list(
      features.basic  = as.data.frame(features.basic),
      summaries       = as.data.frame(summaries),
      d_summaries     = as.data.frame(d_summaries),
      m_summaries     = as.data.frame(m_summaries)    
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
allZips = c('94610','93304')
print('Beginning batch run')
runResult = runModelsByZip(allZips,triggerZip=NULL,truncateAt=-1)
summarizeRun(runResult,listFailures=FALSE)

load(file.path(conf.basePath,'results_summer',paste(zip,'_modelResults.RData',sep='')))
print(names(modelResults))
print(modelResults$summaries[1,]$coefficients)


