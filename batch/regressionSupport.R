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
  binLower = 0
  mat <- c()
  nm = c()
  for(binUpper in c(bins,Inf)) {
    col = regressor - binLower
    col[col < 0] = 0
    col[(col> binUpper-binLower)] = binUpper-binLower
    mat = cbind(mat,col)
    nm = c(nm,paste('tout',binLower,'_',binUpper,sep=''))
    binLower = binUpper
  }
  colnames(mat) <- nm
  return(mat)
}

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
