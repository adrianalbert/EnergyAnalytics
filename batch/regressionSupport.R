# returns a matrix of length(sort(unique(membership))) columns, that is all 
# zeros except for the regressor values that match the membership values 
# corresponding to the column. This supports regression with separate 
# coefficints for each group defined by the membership
# for example, splitRegressor(Tout,dates$hour) will return a matrix with 24 
# columns, where the only non-zero entry per row contains the Tout value in
# the column corresponding to the hour of day it was recorded 
regressor.split = function(regressor,membership) {
  mat <- c()
  nm  <- c()
  for (i in sort(unique(membership))) {
    mat <- cbind(mat,ifelse(membership==i,1,0)) # add a colunm of 1's and 0's
    nm <- c(nm,i)
  }
  colnames(mat) <- nm
  return(mat * regressor)
}

# break a vector out for continuous piecewise regression (i.e. the fitted
# segments will join eachother at each end) into a matrix whose row
# totals are the original values, but whose columns divide the value across
# a set of bins, so 82 across bins with boundaries c(50,60,65,70,80,90) 
# becomes the row 50,10,5,5,10,2,0, which sums to 82...
# This is very useful for finding rough change points in thermal response
# as is expected for buildings with clear setpoints

# TODO: This can create a column of zeros, which should break the regression
# so we might need to prune the columns when we're done and keep track of
# which bins are in play when comparing across regressions
regressor.piecewise = function(regressor,bins) {
  if(any(is.na(bins))) return(as.matrix(regressor)) # if bins itself is NA or any of its values are NA, return the original data
  binLower = 0
  mat <- c()
  nm = c()
  for(binUpper in c(bins,Inf)) {
    col = regressor - binLower
    col[col < 0] = 0
    col[(col > binUpper-binLower)] = binUpper-binLower
    mat = cbind(mat,col)
    nm = c(nm,paste('tout',binLower,'_',binUpper,sep=''))
    binLower = binUpper
  }
  colnames(mat) <- nm
  return(mat)
}

# convienience function putting the bins first for apply style calls...
piecewise.regressor = function(bins,regressor) return(regressor.piecewise(regressor, bins))

lag   = function(v,n=1) { 
  if(n==0) return(v)
  return(c(rep(NA,n),head(v,-n))) 
} # prepend NAs and truncate to preserve length
diff2 = function(v,n=1) { return(c(rep(NA,n),diff(v, n))) } # prepend NAs to preserve length of standard diff
ma    = function(v,n=5) { filter(v,rep(1/n,n), sides=1)   } # calculate moving averages note adds n-1 NAs to beginning

regressorDF = function(residence,norm=FALSE,folds=1,rm.na=FALSE) {
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

hourlyChangePoint = function(df,hourBins=list(1:24),trange=c(50:80),fast=T,reweight=F) {
  # hourBins should be a list of n numeric arrays such that each member of the list 
  # is 1 or more hours of the day to use with the subset command to get 
  # n change point estimates. So the default list(1:24) corresponds to using
  # all the data. Whereas as.list(1:24) would fit 24 subsets...
  if(class(hourBins) != 'list') hourBins = as.list(hourBins)
  if(fast) {
    cps = sapply(hourBins,toutChangePointFast,subset(df),trange,reweight)
  } 
  else {
    cps = sapply(hourBins,toutChangePoint,subset(df),trange,reweight)
  }
  return(cps)
}

# this runs all the models and chooss the minimum one.
# likely waste of CPU on models past the min. See faster impl below
toutChangePoint = function(hrs,df,trange=c(50:85),reweight=F) {
  sub = df$HOD %in% paste('H',sprintf('%02i',(hrs-1)),sep='') # pull out hrs subset (#0-23 in the df)
  df = df[sub,]
  df = df[!is.na(df$kw),]
  steps = sapply(trange,FUN=evalCP,df,reweight)  # run all the models in the range
  #print(steps['SSR',])
  bestFit = steps[,which.min(steps['SSR',])] # find and return the min SSR in the range
  return(bestFit)
}

# this takes advantage of the fact that for one change point, 
# the SSR will be convex so it stops when the change in SSR is positive
# note that each cp in trange could be a list c(40,50,60,70) or just a number
toutChangePointFast = function(hrs,df,trange=c(50:85),reweight) {
  sub = df$HOD %in% paste('H',sprintf('%02i',(hrs-1)),sep='')  # define hrs subset (#0-23 in the df)
  df = df[sub,]                                                # pull out the subset from the df
  df = df[!is.na(df$kw),]                                      # no NA's. Breaks the algorithm
  rng = floor(quantile(df$tout,c(0.1,0.90)))
  trange = c( rng[1]:rng[2]  )
  prev = c(cp=-1,SSR=Inf)                                      # init the compare options
  warnMulti = F
  for(cp in trange) {
    if(length(cp)>1) warnMulti = T # there is no guarantee against global minima
    out = evalCP(cp,df,reweight)                               # run piecewise regression
    if(out['SSR'] > prev['SSR']) { # the previous value was the min
      #plot(df$tout,df$kw,main=paste('Hr',paste(hrs,collapse=',')))
      #points(quickEst(cp,prev['(Intercept)'],prev['lower'],prev['upper']),type='l')
      if(warnMulti) {
        print(paste('warning: toutChangePointFast returning multiple change points: ',
                    paste(prev[grep('^cp[0-9]+',names(prev))],collapse=','),'. ',
                    'May be a local minima. Consider toutChangePoint instead.',sep=''))
      }
      return(prev)
    }
    prev = out
  }
  # failed search
  print(paste('Warning. SSR min not found for hr ',paste(hrs,collapse=','),
                  '. Increase your temperature range? Returning higest value: ',
                  paste(prev[grep('^cp',names(prev),value=T)],collapse=','),sep=''))
  return(prev)
}

evalCP = function(cp,df,reweight=F) {
  out = list()
  tPieces = regressor.piecewise(df$tout,c(cp))
  middle = c()
  nSegs = dim(tPieces)[2]
  if(nSegs > 2) middle = paste('middle_',1:(nSegs-2),sep='')
  colnames(tPieces) <- c('lower',middle,'upper')
  # calculate weights such that each partition of Tout data makes the same 
  # potential contribution to the SSR. So if, for example the data is partitioned
  # into 1/4 and 3/4 of obs in each, the minority data would be weighted at 3x
  n = dim(tPieces)[1]
  m = dim(tPieces)[2]
  df$w = rep(1,n) # default weights
  if(reweight) { 
    # find the index of the last column with a non-zero value
    # for a piecewise regressor, this happens to be the number of non-zero values per row
    highestCol = rowSums(tPieces != 0) 
    colCounts  = table(highestCol)
    # only alter the weights when all columns are participating
    if(length(colCounts) == m) { 
      names(colCounts) <- colnames(tPieces)
      nobs       = sum(colCounts)
      ncols      = length(colCounts)
      colWeights = (nobs/ncols) / colCounts
      #print(colWeights)
      #print(cp)
      #print(colCounts)
      if(colWeights[1] < 1) { # only re-weight if it improves cooling estimate
        df$w = colWeights[highestCol]
      }
    }
  }
  colSums(tPieces != 0)
  df = cbind(df,tPieces) # add the columns from the matrix to the df
  # define regression formula that uses the piecewise pieces
  f_0  = 'kw ~ tout'
  f_cp = paste('kw ~',paste(colnames(tPieces),collapse=" + ",sep=''))
  fit_0  = lm(f_0, df,weights=w)
  fit_cp = lm(f_cp,df,weights=w)
  s_cp = summary(fit_cp)
  s_0  = summary(fit_0)
  coefficients = s_cp$coefficients[,'Estimate'] # regression model coefficients
  pvals        = s_cp$coefficients[,'Pr(>|t|)'] # coefficient pvalues
  # if one of our partitions has no data, we need to fill in the missing columns 
  # in the returned data to ensure that sapply comes back as a matrix, not a list
  if(any(colSums(tPieces != 0)==0)) { 
    missingCols = setdiff(colnames(tPieces),names(coefficients))
    coefficients[missingCols] = NA
    pvals[missingCols]        = NA
  }
  names(pvals) = paste('pval',names(pvals))
  
  # weights to enforce a bayesian idea that higher temps should matter more
  # t > 75 gets a double weighting in the SSR
  w = (df$tout > 75)*3 + 1 
  SSR_cp = (w * fit_cp$residuals) %*% (w * fit_cp$residuals)
  SSR_0  = (w * fit_0$residuals ) %*% (w * fit_0$residuals )
  k = 1 # we have one regressor, so k = 1
  n = length(fit_cp$residuals)
  out$SSR = SSR_cp
  out$RMSE = sqrt(out$SSR/(n - 1))
  out$AIC_cp <- AIC(fit_cp)
  out$AIC_0  <- AIC(fit_0)
  k_cp = s_cp$df[1] #- 1 # degs of freedom, assuming intercept doesn't count
  k_0  = s_0$df[1]  #- 1
  # for comparison of models, see also f-test http://en.wikipedia.org/wiki/F-test#Regression_problems
  fstat = ( (SSR_0 - SSR_cp) / (k_cp - k_0) ) / ( SSR_cp / (n - k_cp) )
  #print(paste((SSR_0 - SSR_cp),(k_cp - k_0),SSR_cp,(n - k_cp)))
  plower = pf(fstat,k_cp - k_0, n - k_cp) # single sided
  nullModelTest = 2 * min(plower, 1 - plower)      # double sided p test for whether the cp model improves on a non-cp model
  return(c(cp=cp,SSR=out$SSR,AIC_cp=out$AIC_cp,AIC_0=out$AIC_0,nullModelTest=nullModelTest,coefficients,pvals))
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
