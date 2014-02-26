toHourly = function(data,dateCol='date'){
  data$dayhr = strftime(data[,dateCol],'%Y-%m-%d.%H',tz='CST6CDT')
  hourly = aggregate(. ~ dayhr,data=data, mean, na.rm=T )
  hourly$date = strptime(hourly$dayhr,'%Y-%m-%d.%H')
  return(hourly)
}

toXmin = function(data,dateCol='date', min = 15){
  dayhr    = strftime(data[,dateCol],'%Y-%m-%d.%H',tz='CST6CDT')
  minute   = as.numeric(strftime(data[,dateCol],'%M',tz='CST6CDT'))
  minute   = (minute %/% 15) * 15
  data$dayhrXm = as.character(paste(dayhr, minute))
  data$date       = NULL
  toXmin_cur      = aggregate(. ~ dayhrXm, data=data, mean, na.rm=T )
  toXmin_cur$date = strptime(toXmin_cur$dayhrXm, '%Y-%m-%d.%H%M')
  return(toXmin_cur)
}

