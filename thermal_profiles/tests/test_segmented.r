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
raw_data  = subset(consumption.ok, UID == unique(consumption.ok$UID)[1])
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
plot(df.hour$TemperatureF, df.hour$kWh)
# library('dummies')
# df.day    = dummy.data.frame(df.day, names = 'DoW', sep = '.')

# ___________________________________________
# Test single breakpoint model on daily data

library('segmented')
fit.lm <- lm(kWh ~ TemperatureF, data = df.hour)
fit.seg<-segmented(fit.lm, seg.Z = ~TemperatureF, psi=65)

test = davies.test(fit.lm,~TemperatureF,k=5)

fit.lm<-update(fit.lm,.~.-TemperatureF)
fit.seg1<-update(fit.seg)

# ________________________________________
# Construct depmixS4 model for breakpoint

library('depmixS4')
source('code/responseBreakpoint.r')

K = 2
bp = list()
pstart = runif(5)
bp[[1]] = breakpoint(formula = kWh~TemperatureF, data = df.hour, pstart=c(0.1,   0, 0.9, 75, 0.1))
bp[[2]] = breakpoint(formula = kWh~TemperatureF, data = df.hour, pstart=c(0.1,   0, 0.9, 70, 0.1))
bp[[3]] = breakpoint(formula = kWh~TemperatureF, data = df.hour, pstart=c(0.1,   0, 0.9, 80, 0.1))

rModels <- list(
  list(
    bp[[1]]
  ),
  list(
    bp[[2]]
  ),
  list(
    bp[[3]]
  )
)

trstart  = rbind(c(0.9,0.1),
                 c(0.1,0.9))
trstart3 = rbind(c(0.8,0.1,0.1),
                 c(0.1,0.8,0.1),
                 c(0.1,0.1,0.8))
transition <- list()
transition[[1]] <- transInit(~1,nstates=3,data=df.hour,pstart=trstart3[1,])
transition[[2]] <- transInit(~1,nstates=3,data=df.hour,pstart=trstart3[2,])
transition[[3]] <- transInit(~1,nstates=3,data=df.hour,pstart=trstart3[3,])

instart=c(1,1,1)/3
inMod <- transInit(~1,ns=3,ps=instart,data=data.frame(rep(1,1,1)))

mod <- makeDepmix(response=rModels,transition=transition,prior=inMod,stat=FALSE)

fm3 <- fit(mod,emc=em.control(rand=FALSE, maxit = 20, tol = 1e-3))
summary(fm3)

library('ggplot2')
df = data.frame(TemperatureF = df.hour$TemperatureF, kWh = df.hour$kWh, State = as.factor(posterior(fm3)$state))
ggplot(df, aes(TemperatureF, kWh, color = State)) + geom_point() + facet_wrap(~State)



