## Put date info and get the day of week info
get.day = function(date){
  ## date format is assumed like '2011-12-31'
  ## 1:Mon 2:Tue 3:Wed 4:Thu 5:Fri 6:Sat 7:Sun
  ## if day.ind is less than 6 => weekdays, if not => weekends
  days = c('Monday','Tuesday','Wednesday','Thursday','Friday','Saturday','Sunday')
  day.ind = which(days==weekdays(as.Date(as.character(date),format='%Y-%m-%d')))
  day.ind
}

## will give you qkw1,qkw2,...,qkw96
q.cols = paste(paste('qkw',1:96,sep=''),collapse=',')

## will give you qkw1+qkw2+...+qkw96
q.cols = paste(paste('qkw',1:96,sep=''),collapse='+')

## will give you hkw1,hkw2,...,hkw24
h.cols = paste(paste('hkw',1:24,sep=''),collapse=',')

## will give you hkw1+hkw2+...+hkw24
h.cols = paste(paste('hkw',1:24,sep=''),collapse='+')
