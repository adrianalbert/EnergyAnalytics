

allZips = c(93304)
summary = combineSummaries(allZips)
estimates.s = cf(summary,'toutPIECES','summer')
estimates.std = cf(summary,'toutPIECES','summer',col='Std. Error')[,c('id','tout75_Inf')]
colnames(estimates.std) <- c('id','std')
estimates.s = merge(estimates.s,estimates.std,by='id')
estimates.HOD.s = cf(summary,'toutTOD_WKND','summer')
HODts = estimates.HOD.s[,c('id','tout65.HODWKWKH15')]
piecests = estimates.s[c('id','tout75_Inf','std')]
combined = merge(piecests,HODts,by='id')
ord = order(combined[,'tout75_Inf'])
combined = combined[ord,]
piece = combined[,'tout75_Inf']
lim = c(min(piece),max(piece))
plot(piece,xlab='Index',ylab='kW/deg F',mgp=c(1,0,0),tcl=0.5,main='Sorted Tout slopes from 93304')
par(new=T)
plot(combined[,'tout65.HODWKWKH15'],col='grey',axes=F,ylim=lim,xlab='',ylab='')
par(new=T)
plot(piece+ 2*combined[,'std'],type='l',col='blue',axes=F,ylim=lim,xlab='',ylab='')
par(new=T)
plot(piece- 2*combined[,'std'],type='l',col='blue',axes=F,ylim=lim,xlab='',ylab='')
legend(0,max(piece)*0.9,c('T > 75 piece','3pm HOD','5% - 95%'),col=c('black','grey','blue'),pch=c(1,1,NA),lwd=c(NA,NA,1))


r = ResDataClass(553991005,93304)  # very clear cooling 24x7
plot(r)
r = ResDataClass(554622151,93304)  # very clear cooling possible timed setback
plot(r)
r = ResDataClass(637321210,93304)  # very clear cooling some bad data in july?
plot(r)
r = ResDataClass(1064423310,93304) # cooling with high outliers. Unclear setpoint
plot(r)
r = ResDataClass(1366549405,93304) # heat and cooling. slight tout slopes
plot(r)

r_max = ResDataClass(combined[dim(combined)[1]-2,'id'],93304)
plot(r_max)
r_max = ResDataClass(combined[dim(combined)[1]-1,'id'],93304)
plot(r_max)