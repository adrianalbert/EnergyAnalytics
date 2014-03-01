toHourly = function(data,dateCol='date'){
  data$dayhr = strftime(data[,dateCol],'%Y-%m-%d %H',tz='CST6CDT')
  hourly = aggregate(. ~ dayhr,data=data, mean, na.action = na.pass )
  hourly$date = strptime(hourly$dayhr,'%Y-%m-%d %H')
  hourly$dayhr = NULL
  hourly = hourly[,-which(names(hourly) == dateCol)]
  return(hourly)
}

toXmin = function(data,dateCol='date', min = 15){
  dayhr    = strftime(data[,dateCol],'%Y-%m-%d %H',tz='CST6CDT')
  minute   = as.numeric(strftime(data[,dateCol],'%M',tz='CST6CDT'))
  minute   = (minute %/% 15) * 15
  data$dayhrXm = as.character(paste(dayhr, minute))
  data$date       = NULL
  toXmin_cur      = aggregate(. ~ dayhrXm, data=data, mean,  na.action = na.pass )
  toXmin_cur$date = strptime(toXmin_cur$dayhrXm, '%Y-%m-%d %H%M')
  toXmin_cur$dayhrXm = NULL
  toXmin_cur = toXmin_cur[,-which(names(toXmin_cur) == dateCol)]
  return(toXmin_cur)
}

