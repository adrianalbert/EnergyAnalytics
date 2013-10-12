# test_segmented.r
#
# Test segmented model with HMM transitions.
# 
# Adrian Albert
# Last modified: April 2013.
# -----------------------------------------------------------------------

rm(list=ls())

setwd('~/Dropbox/OccupancyStates/')

options(error = recover)
library('timeDate')
library('lubridate')

# ______________
# Get test data

load('data/selection_consumption.RData')
raw_data  = subset(consumption.ok, UID == unique(consumption.ok$UID)[1800])
wthr_data = weather.ok[[as.character(raw_data$ZIP5[1])]]
names(wthr_data)[1] = 'date'
TemperatureF = aggregate(data = wthr_data, TemperatureF ~ as.Date(date), FUN = mean, na.rm = T)
df.day    = data.frame(date = raw_data$date, 
                       DoW  = as.factor(weekdays(as.Date(raw_data$date))), 
                       kWh  = rowSums(raw_data[,-(1:5)]),                        
                       TemperatureF = TemperatureF$TemperatureF)
df.hour   = data.frame(date = wthr_data$date, 
                       DoW  = as.factor(weekdays(as.Date(wthr_data$date))), 
                       ToD  = as.factor(hour(wthr_data$date)),
                       kWh  =  as.vector(t(raw_data[,paste('hkw',1:24,sep='')])),                        
                       TemperatureF = wthr_data$TemperatureF)
# library('dummies')
# df.day    = dummy.data.frame(df.day, names = 'DoW', sep = '.')

# ___________________________________________
# Test single breakpoint model on daily data

library('segmented')
fit.lm <- lm(kWh ~ TemperatureF, data = df.day)
fit.seg<-segmented(fit.lm, seg.Z = ~TemperatureF, psi=c(65, 75))

test = davies.test(fit.lm, ~TemperatureF, k=5)

fit.lm<-update(fit.lm,.~.-TemperatureF)
fit.seg1<-update(fit.seg)

plot(df.hour$TemperatureF, df.hour$kWh)

# ________________________________________
# Construct depmixS4 model for breakpoint

library('depmixS4')
source('code/responseSelect.r')

mod = makeModelSelect(df.hour, kWh ~ TemperatureF + ToD + DoW, ~TemperatureF, K = 5)
fm3 = fit(mod,emc=em.control(rand=FALSE, maxit = 50, tol = 1e-3))

summary(fm3)

library('ggplot2')
df = data.frame(TemperatureF = df.hour$TemperatureF, kWh = df.hour$kWh, State = as.factor(posterior(fm3)$state))
ggplot(df, aes(TemperatureF, kWh, color = State)) + geom_point()  + facet_wrap(~State)



