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

outDir = 'results_climate_test'

PLOT_INVALID = FALSE # create png plots for residences that fail validaiton
PLOT_VALID   = TRUE  # create png plots for residences that pass validaiton

RUN_AGGREGATED_MODELS = FALSE # run daily and monthly summary data models (moderate time consuming)
RUN_STEP_SELECTION    = FALSE # run nested model selection algorithm (time consuming)

# generate the string values that will identify the desired subset of a data.frame
# using the command subset(df,subset=str,...)
subset = list(
  summer=paste("MOY %in% c(",paste("'M",5:10,"'",sep='',collapse=','),")"),
  winter=paste("MOY %in% c(",paste("'M",c(11:12,1:4),"'",sep='',collapse=','),")"),
  day=paste("HOD %in% c(",paste("'H",8:19,"'",sep='',collapse=','),")")
  
)
subset$summer_day=paste(subset$summer,'&',subset$day)

# passing in forula objects directly creates lots of problems
# formulas are context specific and cant be stored in dataframes
# to get formula from string, call as.formula(str)
# to get a string repr of a formula object, call deparse(fmla)
models.hourly = list(
  #MOY         = "kw ~ tout + MOY",             
  #DOW         = "kw ~ tout + DOW",
  #DOW_HOD65    = list(formula="kw ~ tout65 + DOW + HOD",subset=list(all="TRUE",summer=subset$summer)),
  #HOW65        = list(formula="kw ~ tout65 + HOW",subset=list(all="TRUE",summer=subset$summer)),
  #wea          = list(formula="kw ~ tout   + pout + rh + HOW + MOY",subset=list(all="TRUE",summer=subset$summer)), 
  #wea65        = list(formula="kw ~ tout65 + pout + rh + HOW + MOY",subset=list(all="TRUE",summer=subset$summer)),
  #HOW         = "kw ~ tout + HOW",
  toutTOD_WKND = list(formula="kw ~ 0 + tout65:HODWK + HODWK",subset=list(summer=subset$summer))
  #toutTOD      = list(formula="kw ~ 0 + tout65:HOD + HOD",subset=list(summer=subset$summer)),
  #toutTOD_d1   = list(formula="kw ~ 0 + tout65:HOD + pout + rh + tout_d1 + HOD",subset=list(summer=subset$summer)),
  #toutTOD_65d1 = list(formula="kw ~ 0 + tout65:HOD + pout + rh + tout65_d1 + HOD",subset=list(summer=subset$summer)),
  #toutTOD_l1   = list(formula="kw ~ 0 + tout65:HOD + pout + rh + tout65_l1 + HOD",subset=list(summer=subset$summer)),
  #toutTOD_l3   = list(formula="kw ~ 0 + tout65:HOD + pout + rh + tout65_l3 + HOD",subset=list(summer=subset$summer))
  #toutTOD_min = "kw_min ~ 0 + tout:HOD + HOW" # no intercept
)

# todo: integration vacation days into regression
models.daily = list(
  #tout           = "kwh ~ tout.mean",
  #DOW            = "kwh ~ DOW",
  #tout_mean      = "kwh ~ tout.mean + DOW",
  #tout_mean_WKND = "kwh ~ tout.mean + WKND",
  #tout_mean_vac  = "kwh ~ tout.mean + WKND + vac",
  #tout_max       = "kwh ~ tout.max  + DOW",
  #tout_CDD       = "kwh ~ CDD + HDD + DOW",
  #tout_CDD_WKND  = "kwh ~ CDD + HDD + WKND",
  #wea_mean       = "kwh ~ tout.mean + pout.mean + rh.mean + WKND + vac"
)

models.monthly = list(
  #tout       = "kwh ~ tout.mean",
  #CDD        = "kwh ~ CDD"
  #CDD_HDD    = "kwh ~ HDD + CDD"
)

lag   = function(v,n=1) { return(c(rep(NA,n),head(v,-n))) } # prepend NAs and truncate to preserve length
diff2 = function(v,n=1) { return(c(rep(NA,n),diff(v, n))) } # prepend NAs to preserve length of standard diff

regressorDF = function(residence,norm=TRUE,folds=1,rm.na=FALSE) {
  WKND       = c('WK','ND')[(residence$dates$wday == 0 | residence$dates$wday == 6) * 1 + 1] # weekend indicator
  dateDays   = as.Date(residence$dates)
  # holidaysNYSE is a function from the dateTime package
  hdays      = as.Date(holidayNYSE(min(residence$dates$year + 1900):max(residence$dates$year + 1900)))
  vac        = factor(dateDays %in% hdays)
  hStr       = paste('H',sprintf('%02i',residence$dates$hour),sep='')
  hwkndStr   = paste(WKND,hStr,sep='')
  dStr       = paste('D',residence$dates$wday,sep='')
  mStr       = paste('M',residence$dates$mon,sep='')
  howStrs    = paste(dStr,hStr,sep='')
  MOY        = factor(mStr,levels=sort(unique(mStr)))         # month of year
  DOW        = factor(dStr,levels=sort(unique(dStr)))         # day of week
  HOD        = factor(hStr,levels=sort(unique(hStr)))         # hour of day
  HODWK      = factor(hwkndStr,levels=sort(unique(hwkndStr))) # hour of day for weekdays and weekends
  HOW        = factor(howStrs,levels=sort(unique(howStrs))) # hour of week
  tout       = residence$w('tout')
  tout65     = pmax(0,tout-65)
  
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
    tout65=tout65,
    tout65_l1 = lag(tout65,1),
    tout65_l3 = lag(tout65,3),
    tout_d1 = diff2(tout,1),
    tout65_d1 = diff2(tout,1)*(tout65 > 0),
    #tout_d3 = diff2(tout,3),
    pout=pout,
    rain=rain,
    rh=rh,
    dates=residence$dates,
    vac=vac,
    wday=residence$dates$wday,
    MOY,DOW,HOD,HODWK,HOW,WKND   )
  if(rm.na) { df = df[!rowSums(is.na(df)),] }
  #df = cbind(df,toutPIECES) # add the columns with names from the matrix to the df
  return(df)
}

regressorDFAggregated = function(residence,norm=TRUE,bp=65,rm.na=FALSE) {
  # uses melt and cast to reshape and aggregate data
  df = residence$df() # kw, tout, dates
  if(norm) df$kw_norm = residence$norm(df$kw)
  df$kw_min = df$kw - quantile(df$kw,na.rm=TRUE,c(0.02)) # remove the min for regression w/o const term
  
  df$pout  = residence$w('pout')
  dp       = residence$w('dp')
  df$rh    = residence$weather$rh(df$tout,dp)
  df$day   = format(df$dates,'%Y-%m-%d') # melt has a problem with dates
  df$wday  = as.POSIXlt(df$dates)$wday   # raw for subsetting Su=0 ... Sa=6
  df$DOW   = paste('D',df$wday,sep='')  # Su=0 ... Sa=6
  df$WKND  = (df$wday == 0 | df$wday == 6) * 1 # weekend indicator
  df$DOW   = factor(df$DOW, levels=sort(unique(df$DOW)))
  month    = format(df$dates,'%y-%m')   # as.POSIXlt(df$dates)$mon
  df$mon   = as.POSIXlt(df$dates)$mon   # raw month data for subset functions Jan=0 ... Dec=11
  df$MOY   = factor(month, levels=sort(unique(month)))
  df <- subset(df, select = -c(dates) )  # melt has a problem with dates but we don't need anymore
  # melt and cast to reshape data into monthly and daily time averages
  dfm = melt(df,id.vars=c("day",'DOW','MOY','mon','wday','WKND'),measure.vars=c('kw','tout','pout','rh'),na.rm=TRUE)
  
  monthly = cast(dfm,MOY + mon ~ variable,fun.aggregate=c(sum,mean,function(ar1,bp=65) sum(ar1 > bp),function(ar2,bp=65) sum(ar2 < bp)),subset= variable %in% c('kw','tout'))
  colnames(monthly) <- c('MOY','mon','kwh','kw.mean','junk1','junk2','junk3','tout.mean','CDD','HDD')
  monthly <- subset(monthly, select = -c(junk1, junk2, junk3) )
  
  daily = cast(dfm, MOY + day + DOW + mon + wday + WKND ~ variable,fun.aggregate=c(sum,mean,max,function(ar1,bp=65) sum(ar1 > bp),function(ar2,bp=65) sum(ar2 < bp)),subset= variable %in% c('kw','tout','pout','rh'))
  colnames(daily) <- c('MOY','day','DOW','mon','wday','WKND','kwh','kw.mean','kw.max','junk1','junk2','junk3','tout.mean','tout.max','CDD','HDD','junk4','pout.mean','pout.max','junk5','junk6','junk7','rh.mean','rh.max','junk8','junk9')
  daily <- subset(daily, select = grep("^junk", colnames(daily), invert=TRUE) )
  
  # add vacation days flags
  dateDays   = as.Date(daily$day)
  # holidaysNYSE is a function from the dateTime package
  hdays      = as.Date(holidayNYSE(min(as.POSIXlt(dateDays)$year+1900):max(as.POSIXlt(dateDays)$year+1900)))
  daily$vac  = factor(dateDays %in% hdays)
  
  if(FALSE) {
    M <- rbind(c(1, 2), c(3, 4), c(5, 6))
    layout(M)
    pacf(daily$kwh)
    acf(daily$kwh)
    plot(daily$pout.mean, daily$kwh)
    plot(daily$rh.mean, daily$kwh)
    plot(daily$tout.max, daily$kwh)
    plot(daily$kwh,type='l',main=paste('',residence$id))
    Sys.sleep(1)
  }
  if(rm.na) { 
    daily   = daily[!rowSums(is.na(daily)),] 
    monthly = monthly[!rowSums(is.na(monthly)),] 
  }
  return(list(daily=daily,monthly=monthly))
}


runModelsByZip = function(zipArray,triggerZip=NULL,truncateAt=-1) {
  nZip     <- length(zipArray)
  zipCount <- 0
  
  
  print(paste('This batch process will run models for',nZip,'zip codes'))
  dir.create(file.path(getwd(),outDir),showWarnings=FALSE)
  print(paste('<zipcode>_modelResults.RData files for each can be found in',file.path(getwd(),outDir)))
  res = list(attemptedZip =zipArray,
             completedZip =c(),
             attemptedSP  =c(),
             completedSP  =c(),
             invalid.ids  =c() )
  triggered = FALSE
  for (zip in zipArray) { # i.e. 94610
    zipCount <- zipCount + 1
    if (is.null(triggerZip) || zip == triggerZip) triggered = TRUE
    if (! triggered) {
      next
    }
    
    print(paste('Running models for ',zip,' (',zipCount,'/',nZip,')', sep=''))
    resultsFile <- file.path(getwd(),outDir,paste(zip,'_modelResults.RData',sep=''))
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
    
    # just for testing
    #break
  } # zip loop
  return(res)
}

# TODO: test rh, pout regressors
# TODO: define junk shot model for use with "step"
# TODO: aks Ram about his k-fold idea and look at serial folds as well as random
# TODO: normalize the output format for model runs so that consolidate works again
# TODO: develop consolidate across multiple zips
# TODO: consolidate after each model run, before serializing to disk


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

summerlm = function(fmla,df,x) {
  return(lm(fmla,df,x=x,subset=MOY %in% c('M6','M7','M8','M9')))
}

validate = function(r) {
  issues = data.frame(id=r$id)
  timeDiffs = diff(r$dates)
  units(timeDiffs) <- "hours"
  maxtd = max(timeDiffs) / 24
  span = difftime(tail(r$dates, n=1),r$dates[1],units='days')
  zerospct = sum((r$kw == 0)*1,na.rm=TRUE) / length(r$kw)
  kwmean = mean(r$kw,na.rm=TRUE)
  daylen = length(r$days)
  if( daylen < 180)     issues$days180    = daylen # less than 180 days (could be non-consecutive)
  #if( span < 270 )      issues$span270    = span # spanning less than 270 days total
  #if( maxtd > 60 )      issues$bigdiff    = maxtd # more than 2 months of missing data
  if( kwmean < 0.110 )  issues$lowmean    = kwmean # mean less than 150W is almost always empty or bad readings
  if( zerospct > 0.15 ) issues$zerospct15 = zerospct # over 15% of readings are zero
  return(issues)
}

runModelsBySP = function(sp_ids,zip=NULL,data=NULL,weather=NULL,truncateAt=-1) {
  features.basic  <- c()
  summaries       <- c()
  d_summaries     <- c()
  m_summaries     <- c()
  steps           <- c()
  
  invalid_ids     <- data.frame()
  # TODO:
  # use step for forward selection
  # find a faster method of cross validation than cvTools or DAAG provide
  # add other weather covariates.
  splen  <- length(sp_ids)
  i <- 0
  skip = FALSE
  for (sp_id in sp_ids) { # i.e. "820735863"
    i <- i+1
    print(paste('  ',sp_id,' (',i,'/',splen,') in ',zip,sep=''))
    if(sp_id == 6502353905) { skip = FALSE }
    if(skip) { next }
    resData = NULL
    if (!is.null(data)) {  resData = data[data[,'sp_id']== sp_id,] }
    r <- tryCatch(ResDataClass(sp_id,zip=zip,weather=weather,data=resData,db=conf.meterDB()), 
                  error = function(e) {print(e)}, finally={} )
    if ( ! "ResDataClass" %in% class(r) ) { # constructor returns the error class if it has a problem
      print('    not found (see error)') # ignore residences that produce errors
    }
    else { # viable residence!
      r <- tryCatch( {
        issues = validate(r)
        savePlot = FALSE
        if(length(issues)>1) # all issues will return with an id set, so > 1 is a problem
        { 
          invalid_ids <- rbind.fill(invalid_ids,issues)
          print(paste('Bad or insufficient data:',paste(colnames(issues),collapse=', ')))
          if(PLOT_INVALID) {
            png(file.path(getwd(),outDir,paste(r$zip,r$id,'invalid.png',sep='_')))
            colorMap = rev(colorRampPalette(brewer.pal(11,"RdBu"))(100))
            plot(  r, colorMap=colorMap, 
                   main=paste(r$zip, r$id, paste(colnames(issues),collapse=', ')) )
            dev.off()
          }
          next # no further processing
        }
        if(PLOT_VALID) {
          png(file.path(getwd(),outDir,paste(r$zip,'_',r$id,'.png',sep='')))
          redblue = rev(colorRampPalette(brewer.pal(11,"RdBu"))(100))
          plot( r,colorMap=redblue,main=paste(r$zip, r$id) )
          dev.off()
        }
        tic('model run')
        features.basic  <- rbind(features.basic,basicFeatures(r$kwMat,sp_id))
        
        # hourly regressions
        if(TRUE) {
          df = regressorDF(r,norm=FALSE) # see also regressorDFAggregated
          models = models.hourly
          if (TRUE) {
            # add special data to the data frame: piecewise tout data
            #toutPIECES = regressor.piecewise(r$tout,c(40,50,60,70,80,90))
            toutPIECES = regressor.piecewise(r$tout,c(60,75))
            df = cbind(df,toutPIECES) # add the columns with names from the matrix to the df
            # define regression formula that uses the piecewise pieces
            piecesf = paste('kw ~',paste(colnames(toutPIECES), collapse= "+"),'+ HOW')
            models$toutPIECES = list(formula=piecesf,subset=list(all="TRUE",summer=subset$summer))
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
            for(snm in names(subs)) {
              df$sub = eval(parse(text=subs[[snm]]),envir=df)
              res = lm(fmla,df, x=TRUE, subset=sub==TRUE)
              summaries = rbind(summaries,summarizeModel(res,df,models,nm,id=sp_id,zip=zip,subnm=snm,fold=fld,formula=fmla,subset=subs[[snm]]))
            }
          }
          if(RUN_STEP_SELECTION) { 
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

        if (RUN_AGGREGATED_MODELS) {
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
    if (truncateAt > 0 & i >= truncateAt) break
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
# bakersfield, oakland
testZips = c(94610,93304)

allZips = c(94923,94503,94574,94559,94028,94539,94564,94702,94704,94085,
            95035,94041,95112,95113,95765,95648,95901,94531,94585,95205,
            95202,93619,93614,93304,93701,95631,95726,95223,95666)

print('Beginning batch run')
runResult = runModelsByZip(testZips,triggerZip=93304,truncateAt=200)
summarizeRun(runResult,listFailures=FALSE)

zip = allZips[1]
load(file.path(getwd(),outDir,paste(zip,'_modelResults.RData',sep='')))
print(names(modelResults))
print(modelResults$summaries[1,]$coefficients)

r = ResDataClass(820735863,94610)  # heat, no cooling
r = ResDataClass(553991005,93304)  # very clear cooling 24x7
r = ResDataClass(554622151,93304)  # very clear cooling possible timed setback
r = ResDataClass(637321210,93304)  # very clear cooling some bad data in july?
r = ResDataClass(1064423310,93304) # cooling with high outliers. Unclear setpoint
r = ResDataClass(1366549405,93304) # heat and cooling. slight tout slopes


