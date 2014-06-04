library(lubridate)

select_data = function(df, dateCol = 'date',
                       include.weekends = TRUE,
                       seasons = c('Summer', 'Winter')) {
  
  # add month and season
  df$month  = month(df[,dateCol], label = TRUE, abbr = FALSE)
  df$season = 'None'; 
  summer.months = c('May', 'June', 'July', 'August', 'September')
  idx.summer= which(df$month %in% summer.months)
  idx.winter= which(!(df$month %in% summer.months))
  if (length(idx.summer)>0) df$season[idx.summer] = 'Summer'
  if (length(idx.winter)>0) df$season[idx.winter] = 'Winter'
  
  # add weekdays
  df$weekday = wday(df$date, label = TRUE, abbr = FALSE)  
  
  # select just weekdays?
  df.sel = df
  if (!include.weekends) df.sel = subset(df.sel, !(weekday %in% c('Saturday', 'Sunday')))
  
  if (nrow(df.sel) == 0) return(NULL)
  
  # which seasons to analyze?
  df.sel = subset(df.sel, season %in% seasons)

  if (nrow(df.sel) == 0) return(NULL)
  
  return(df.sel)
}