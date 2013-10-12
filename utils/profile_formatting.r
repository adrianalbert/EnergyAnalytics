# ----------------------------------------------------
# Formatting utilities
# ----------------------------------------------------

# ____________________________________________________________________
# Function to convert hourly data to day profile (multiple variables)

convert_long_to_profile = function(data) {
  # assumes that time stamps are in the first column
  days     = as.character(as.Date(data[,1]))
  data.mat = matrix(as.matrix(data[,-1]), ncol = 24*ncol(data[,-1]))
  colnames(data.mat) = c(t(outer(names(data[,-1]), 1:24, FUN = function(x,y) paste(x, y, sep='.'))))
  data.mat = as.data.frame(data.mat)
  data.mat = cbind(DateTime = unique(days), data.mat)
  data.mat$DateTime = as.character(data.mat$DateTime)
  
  return(data.mat)
}

# ____________________________________________________________________
# Function to convert daily data to hourly profile (univariate)

convert_day_to_long = function(kwh.profile) {
  
  # convert to long format
  cons_names  = paste('hkw', 1:24, sep='')
  days        = kwh.profile$date              
  hours       = paste(rep(formatC(0:23, flag=0, width=2), length(days)), ':00:00', sep='')
  days_times  = rep(days, each = 24)
  days_times  = paste(days_times, hours, sep=' ')
  # days_times  = c(days_times, days_times[length(days_times)] + 3600)
  idx.dup     = which(duplicated(days_times))
  if (length(idx.dup) > 0) days_times = days_times[-idx.dup]
  kwh.long    = as.vector(t(kwh.profile[,cons_names]))  
  
  result      = data.frame(DateTime = days_times, kWh = kwh.long)
  result$DateTime = as.character(result$DateTime)
  
  return(result)
}



# ______________________________________________________________
# Function to convert stacked day profile to multiple profiles

convert_day_to_long = function(kwh.profile) {
  
  # convert to long format
  cons_names  = paste('hkw', 1:24, sep='')
  days        = kwh.profile$date              
  hours       = paste(rep(formatC(0:23, flag=0, width=2), length(days)), ':00:00', sep='')
  days_times  = rep(days, each = 24)
  days_times  = paste(days_times, hours, sep=' ')
  # days_times  = c(days_times, days_times[length(days_times)] + 3600)
  idx.dup     = which(duplicated(days_times))
  if (length(idx.dup) > 0) days_times = days_times[-idx.dup]
  kwh.long    = as.vector(t(kwh.profile[,cons_names]))  
  
  result      = data.frame(DateTime = days_times, kWh = kwh.long)
  result$DateTime = as.character(result$DateTime)
  
  return(result)
}
