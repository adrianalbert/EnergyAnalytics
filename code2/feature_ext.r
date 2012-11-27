ext.daily.feature.96=function(ddata){ 
  ## ddata: sp_id+date+96points usage data
  days = c('Monday','Tuesday','Wednesday','Thursday','Friday','Saturday','Sunday')
  day.ind = which(days==weekdays(as.Date(as.character(ddata[2]))))
  udata = as.numeric(ddata[3:98])
  peak.ind = which.max(apply(matrix(udata,nrow=4),2,sum))
  
  ## make one row to return: it should be 'list' type to cover uncertain length of ramp info
  arow = list(sp_id = ddata[1],date = ddata[2], ## sp_id , date
           dmin = min(udata), dmax = max(udata), dmean = mean(udata), dsd = sd(udata), ## min, max, mean, sd
           dmaxmean = max(udata)/mean(udata), drange = max(udata)-min(udata), ## max/mean, max-min
           dq25 = quantile(udata,0.25), dq50 = quantile(udata,0.50), dq75 = quantile(udata,0.75),  ## quantile 25,50,75%   
           ## didn't put 1,5,95,99% quantile, total is 96points=> it would be so close to max, min
           dpeak = peak.ind, dind = day.ind, dmaxdiff = max((diff(udata))), dmindiff = min((diff(udata))), ## peak hour indicator, day indicator(if <6, it's weekday), max and min of usage change
           d3du = sum(udata>max(udata)*0.97)/4, d6du = sum(udata>max(udata)*0.94)/4 ## 3% duration, 6% duration: hour measure
           )
  arow
}

ext.daily.feature.24=function(ddata){ 
  ## ddata: sp_id+date+24points usage data
  days = c('Monday','Tuesday','Wednesday','Thursday','Friday','Saturday','Sunday')
  day.ind = which(days==weekdays(as.Date(as.character(ddata[2]))))
  udata = as.numeric(ddata[3:26])
  peak.ind = which.max(udata)
  ## one row has min, max, mean, max/mean, sd, peak hour ind, day ind, label ind, simple detected ramp, ramp info, max use change
  arow = list(sp_id = ddata[1],date = ddata[2], ## sp_id , date
         dmin = min(udata), dmax = max(udata), dmean = mean(udata), dsd = sd(udata), ## min, max, mean, sd
         dmaxmean = max(udata)/mean(udata), drange = max(udata)-min(udata), ## max/mean, max-min
         dq25 = quantile(udata,0.25), dq50 = quantile(udata,0.50), dq75 = quantile(udata,0.75),  ## quantile 25,50,75%   
         ## didn't put 1,5,95,99% quantile, total is 96points=> it would be so close to max, min
         dpeak = peak.ind, dind = day.ind, dmaxdiff = max((diff(udata))), dmindiff = min((diff(udata))), ## peak hour indicator, day indicator(if <6, it's weekday), max and min of usage change
         d3du = sum(udata>max(udata)*0.97), d6du = sum(udata>max(udata)*0.94) ## 3% duration, 6% duration: hour measure
         )
  arow
}
