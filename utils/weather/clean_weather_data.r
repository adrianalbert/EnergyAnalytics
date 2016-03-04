# clean_weather_data.r
#
# Perform imputation on weather data
#
# Adrian Albert
#
# Last modified: May 2014.

library('dummies')
library('timeDate')
library('lubridate')
library('Amelia')
#library('imputation')
library('zoo')
library('VIM')

# _______________________________________________________
# Function to imputate missing observations using EM/SVD

imputateMissingValues = function(X.imp, verbose=T, method = 'AMELIA'){
  sd_col   = apply(X.imp, 2, function(x) sd(x,na.rm=T))
  rm.vars  = which(sd_col == 0 | is.na(sd_col))
  
  if (length(rm.vars) == ncol(X.imp)) {
    
    print('All non-NA covariate values are constant - cannot impute!')
    X.ok.imp = X.imp
    return(X.ok.imp)
  }
  
  if (length(rm.vars)>0) X.imp = X.imp[,-rm.vars]
  
  if (method == 'AMELIA') {
    bds        = matrix(nrow=ncol(X.imp), ncol=3)
    bds[,1]    = 1:ncol(X.imp)
    bds[,2:ncol(bds)] = t(sapply(1:ncol(X.imp), function(i) range(X.imp[,i], na.rm=T)))    
    res   = try(amelia(m=1, x = X.imp, bounds = bds)$imputations[[1]])
    if (class(res) == 'try-error') X.ok.imp = X.imp else X.ok.imp = res
  }
  
  if (method == 'IRMI') {
    X.ok.imp = irmi(X.imp)#, robust = T, mi = 3)
  }
  
  if (method == 'NONE') {
    X.ok.imp = X.imp
  }
  
  return(X.ok.imp)
}

# ________________________________________________________
# Function to perform analysis and do plots of imputation

perform_analysis = function(orig, clean, timestamps) {
    
  # plots of densities for each covariate (before/after)
  covars = names(clean)
  orig   = orig[,covars]
  orig$Type  = 'orig'
  clean$Type = 'clean'
  df = rbind(melt(orig, id.vars = c('Type')), melt(clean, id.vars = c('Type')))
  plt1 = ggplot(df, aes(x = value, color = Type)) + 
    facet_wrap(~variable, nrow=2, scales = 'free') + 
    geom_density(size=2) + 
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.text.x      = element_text(size=18),                  
          axis.text.y      = element_text(size=18),  
          strip.text.x     = element_text(size=18),                                   
          legend.text      = element_text(size=16),                         
          plot.title       = element_text(size=18),          
          axis.title.x = element_blank(),
          panel.background = element_blank(),
          axis.ticks = element_blank()) + 
    ggtitle('Density: Raw and Cleaned')
  
  # table of # of points imputed per covariate
  nr.NA = sapply(covars, function(v) length(which(is.na(orig[,v]))))
  levels(df$variable) = paste(levels(df$variable), nr.NA[levels(df$variable)])
  
  # time series plot: original and estimate
  
  orig[,dateCol] = timestamps
  clean[,dateCol]= timestamps
  df = rbind(melt(orig, id.vars = c('Type', 'date')), melt(clean, id.vars = c('Type', 'date')))
  plt2 = ggplot(df, aes(x = date, y = value, color = Type)) + 
    facet_wrap(~variable, ncol = 1, scales = 'free') + 
    geom_line() + 
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.text.x      = element_text(size=18),                  
          axis.text.y      = element_text(size=18),                         
          legend.text      = element_text(size=16),                         
          strip.text.x     = element_text(size=18),                         
          plot.title       = element_text(size=18),
          axis.title.x     = element_blank(),
          panel.background = element_blank(),
          axis.ticks = element_blank()) + 
   ggtitle(paste('Weather: Actual and Imputed'))
  
  return(list(nr.NA/nrow(orig), plt1, plt2))
}

# __________________________
# Function to clean up data 

clean_weather_data = function(wthr_data, dateCol = 'date', addToD = F) {
    
  # select covariates of interest
  wthr_times  = as.character(wthr_data[,dateCol])
  wthr_posix  = as.POSIXct(wthr_times)
  wthr_names  = c('TemperatureF', 'DewpointF', 'Pressure', 'WindSpeed', 'Humidity', 
                  'HourlyPrecip', 'SolarRadiation')
  wthr_data   = wthr_data[,intersect(names(wthr_data), wthr_names)]            
    
  # remove duplicate timestamps in weather data, if any
  idx_dup     = which(duplicated(wthr_times))
  if (length(idx_dup) > 0) {
    wthr_times = wthr_times[-idx_dup]
    wthr_data  = wthr_data[-idx_dup,]
  }  
  
  # clean up weather covariates (remove unreasonable values)
  ranges = list(TemperatureF = c(0, 125), 
                Pressure     = c(25, 35),
                Humidity     = c(5, 100),
                DewpointF    = c(0, 125))
  for (r in intersect(names(ranges),names(wthr_data))){
    idx = which(wthr_data[,r]<ranges[[r]][1] | wthr_data[,r]>ranges[[r]][2]) 
    if (length(idx)>0) wthr_data[idx,r] = NA
  }

  # remove covariates that have less than 1/2 of the values
  na.number = sapply(names(wthr_data), function(v) {
    return(length(which(is.na(wthr_data[,v]))))
  })
  idx = which(na.number > nrow(wthr_data) / 2)                   
  if (length(idx) > 0) wthr_data = wthr_data[,-idx]              
    
  # imputate missing values (at least one covariate value in an observation tuple needs to be defined)
  wthr_clean = imputateMissingValues(wthr_data)
  wthr_covar = names(wthr_clean)
  
  # add temporal indicators?
  if (addToD) {
    Day.Of.Week = as.factor(dayOfWeek(timeDate(wthr_times)))
    Month       = as.factor(month(timeDate(wthr_times)))
    Hour.Of.Day = as.factor(as.numeric(sapply( sapply( sapply(strsplit(as.character(wthr_times),' '), '[', 2), strsplit, ':' ), '[',1)))
    Is.Holiday  = as.factor(as.numeric(isHoliday( timeDate(wthr_times) )))
    wthr_clean = cbind(data.frame(Month = Month, HourOfDay = Hour.Of.Day), wthr_clean)
  }
  
  # use gbm to estimate gaps (no covariate value available for a given observation)
  X.dum = dummy.data.frame(wthr_clean)
  X.imp = gbmImpute(X.dum, max.iters = 2, cv.fold = 4, verbose = T)$x 
  wthr_clean = X.imp[,wthr_covar]  
  wthr_clean = cbind(date = wthr_times, wthr_clean)
  
  # add in dummies
  if (addToD) {
    wthr_clean = cbind(data.frame(Month = Month, HourOfDay = Hour.Of.Day, 
                                  IsHoliday = Is.Holiday, DayOfWeek = Day.Of.Week), wthr_clean)
  } 
  
  return(wthr_clean)
}

# # plots 
# res = perform_analysis(wthr_data, wthr_clean, wthr_posix)
# png(paste(plots_path, paste(zip, 'density.png', sep='_'), sep=''), width=1200, height = 800)
# print(res[[2]])
# dev.off()
# png(paste(plots_path, paste('zip', 'ts.png', sep='_'), sep=''), width=2000, height = 1200)
# print(res[[3]])
# dev.off()
# 
# 
