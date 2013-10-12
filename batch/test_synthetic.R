source(file.path(conf.basePath,'regressionSupport.R')) # mostly regressor manipulation
source(file.path(conf.basePath,'timer.R'))             # adds tic() and toc() functions

base = 'C:/dev/pge_collab/'
file.simOut     = paste(base,'simulation/CZ016-Sched7to6Wkdys-4CSB-GainsMeter.csv',sep='')
file.simWeather = paste(base,'simulation/CZ01RV2.epw',sep='')

rawOut = read.csv(file.simOut)
p = rawOut[,4]
# R assumes the current year, since none is give. 
# However, this means that days of the week will move around
# so we should assign a year with the same days of the week as the 
# simulation. For example, 2010.
dates = as.POSIXlt(paste('2010',rawOut[,1]),format='%Y %m/%d  %H:%M:%S')
rm(rawOut)

rawWeather = read.csv(file.simWeather,skip=8,header=FALSE)
wDates = as.POSIXlt(paste('2010',rawWeather[,2],rawWeather[,3],rawWeather[,4]-1,rawWeather[,5]),format='%Y %m %d %H %M')
tOut = approx(wDates, rawWeather[,7], dates, method="linear")[[2]]
rm(wDates,rawWeather,base,file.simOut,file.simWeather)

df = data.frame(p,dates,tOut)
rm(p,dates,tOut)

# todo:
# define linear model
model.MOY <- lm(r$norm(r$kw) ~ r$tout     + MOY)
# assign each row of data to 1 of N pools of data
# fit N models with their assigned data
# for every row, and each of the N models, calculate error for model fit
# assign row to minimal error model
# repeat until no points change models.