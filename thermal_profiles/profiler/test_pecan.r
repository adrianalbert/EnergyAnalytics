# test_pecan.r
#
# Applies HMM decoding on Pecan Street data. 
#
# Adrian Albert
# Last modified: December 2013.
# ---------------------------------------------------------

rm(list = ls())
options(error = recover)
library('segmented')

# __________________________________________________
# Initializations...

setwd('~/Dropbox/EnergyAnalytics/thermal_profiles/profiler/')
source('stateProcessorWrapper.r')

load('~/energy-data/pecan_street/pecan_plus_weather.RData')

PLOTS_PATH = '~/Dropbox/OccupancyStates/plots/pecan-street'
dir.create(PLOTS_PATH)

# __________________________________________________
# Helper functions

usage = function(data,labels) {
  sub = match(labels,names(data))
  sub = sub[!is.na(sub)]
  return(apply(as.matrix(data[,sub]),1,sum))
}

toHourly = function(data,dateCol='date'){
  data$dayhr = strftime(data[,dateCol],'%Y-%m-%d.%H',tz='CST6CDT')
  hourly = aggregate(. ~ dayhr,data=data, mean, na.rm=T )
  hourly$date = strptime(hourly$dayhr,'%Y-%m-%d.%H')
  return(hourly)
}

to10min = function(data,dateCol='date'){
  dayhr    = strftime(data[,dateCol],'%Y-%m-%d.%H',tz='CST6CDT')
  minute   = as.numeric(strftime(data[,dateCol],'%M',tz='CST6CDT'))
  minute   = (minute %/% 10) * 10
  data$dayhr10m = as.character(paste(dayhr, minute))
  data$date     = NULL
  to10min       = aggregate(. ~ dayhr10m,data=data, mean, na.rm=T )
  to10min$date  = strptime(to10min$dayhr10m,'%Y-%m-%d.%H%M')
  return(to10min)
}

# __________________________________________________
# Select appliances of interest

HVAC   = c('AIR1','AIR2','FURNACE1','FURNACE2')
HV     = c('FURNACE1','FURNACE2')
AC     = c('AIR1','AIR2')
fridge = c('REFRIGERATOR1')
dhw    = c('WATERHEATER1')
aux    = c('SPRINKLER1')
user   = c('DININGROOM1','DISHWASHER1','DRYG1','KITCHEN1','KITCHEN2','KSAC1','KSAC2','SECURITY1',
           'LIGHTING_PLUGS1','LIGHTING_PLUGS2','LIGHTING_PLUGS3','LIGHTING_PLUGS4',
           'LIVINGROOM1','MASTERBATH1','MICROWAVE1','OVEN1','RANGE1','WASHER1','WASHINGMACHINE1',
           'BATHROOM1','BEDROOM1','BEDROOM2','DISPOSAL1','DRYE1','FAMROOM1','BATH1','DRYER2','GARAGE','GARAGE1',
           'GENLIGHT1','SMALLAPPLIANCE1','SMALLAPPLIANCE2','SMALLAPPLIANCE3','COOKTOP1',
           'BACKYARD1','OFFICE1','TVROOM1','THEATER1','MASTERBED1')
total   = c('use')
solar   = c('gen')
ev      = c('CAR1')
net     = c('Grid')

# __________________________________________________
# Plot usage by minute & by hour

ppw = pecanCombinedData
if (1 == 1) {
for(homeName in names(ppw)){
  print(homeName)
  homeData = ppw[[homeName]]
  #print(names(homeData))
  print(dim(homeData))
	homeName = trim(homeName)

  # define plot parameters
  op <- par(no.readonly = TRUE)
  cols = c('black','#FE2E2E','#0040FF','#ff9999', '#bbbbbb')
  m <- matrix(c(1,2,3,4),nrow=4,ncol=1,byrow=T)
  layout(mat = m,heights = c(0.3,0.3,0.3,0.1))
  par(oma=c(2,2,2,0),mar=c(2,4,2,1)) # Room for the title
  
  # print minute-by-minute data
  maxUsage = max(usage(homeData,total))
  plot(homeData$date,usage(homeData,total),type='l',col=cols[1],main='minute',ylab='kW',ylim=c(0,1.1*maxUsage),xaxt='n')
  points(homeData$date,usage(homeData,AC),type='l',col=cols[3])
  points(homeData$date,usage(homeData,HV),type='l',col=cols[2])
  points(homeData$date,usage(homeData,user),type='l',col=cols[4])
  points(homeData$date,homeData$TemperatureF/100*maxUsage,type='l',col=cols[5])
  d = homeData$date
  dts = seq(d[1],d[length(d)],by=3600*24) # one per day from one per minute
  axis(1, dts, format(dts, "%a, %m/%d"), cex.axis=1)
  grid(nx=NA,ny=NULL)
  abline(v=dts,col="black",lty=3)
  
  # aggregate data to 10 minutes  
  min10Data = to10min(homeData)
  maxUsage = max(usage(min10Data,total))
  plot(min10Data$date,usage(min10Data,total),type='l',col=cols[1],main='10-minute',ylab='kW',ylim=c(0,1.1*maxUsage),xaxt='n')
  points(min10Data$date,usage(min10Data,AC),type='l',col=cols[3])
  points(min10Data$date,usage(min10Data,HV),type='l',col=cols[2])
  points(min10Data$date,usage(min10Data,user),type='l',col=cols[4])
  points(min10Data$date,min10Data$TemperatureF/100*maxUsage,type='l',col=cols[5])
  d = min10Data$date
  dts = seq(d[1],d[length(d)],by=3600*24) # one per day from one per minute
  axis(1, dts, format(dts, "%a, %m/%d"), cex.axis=1)
  grid(nx=NA,ny=NULL)
  abline(v=dts,col="black",lty=3)
  
  # print hourly data
  hrData = toHourly(homeData)
  maxUsage = max(usage(hrData,total))
  plot(hrData$date,usage(hrData,total),type='l',col=cols[1],main='hourly',ylim=c(0,1.1*maxUsage),ylab='kW',xaxt='n')
  points(hrData$date,usage(hrData,HV),type='l',col=cols[2])
  points(hrData$date,usage(hrData,AC),type='l',col=cols[3])
  points(hrData$date,usage(hrData,user),type='l',col=cols[4])
  points(hrData$date,hrData$TemperatureF/100*maxUsage,type='l',col=cols[5])
  d = as.POSIXct(hrData$date)
  dts = seq(d[1],d[length(d)],by=3600*24) # one per day from one per hour
  axis(1, dts, format(dts, "%a, %m/%d"), cex.axis=1)
  grid(nx=NA,ny=NULL)
  abline(v=dts,col="black",lty=3)
  mtext(paste('Pecan Street example','kW demand'), line=0, font=2, cex=1.2, outer=TRUE)
  par(mar=c(1,4,1,1))
  plot.new()
  legend("center", lty=1,cex=1,lwd=2,
         legend=c('total','HV', 'AC', 'user','temperature'), 
         col=cols,horiz=T)
  # par(op)
  #break
	dir.create(file.path(PLOTS_PATH, homeName))
  dev.copy2pdf(file=paste(paste(PLOTS_PATH, homeName, sep = '/'), paste(homeName,'.pdf',sep=''), sep='/'),width=10,height=6)
}
}

#plot(homeData$date,usage(homeData,'THEATER1'),type='l')
#points(homeData$date,usage(homeData,fridge),type='l',col='green')

# __________________________________________________
# Apply Thermal States model to Pecan data

train.frac = 0.95
psi    = matrix(ncol = 2, nrow = length(ppw))
colnames(psi) = c('mins10', 'minute')
rownames(psi) = as.character(sapply(names(ppw), trim))
df.tot = list()
if (1 == 1) {
for(homeName in names(ppw)){
  print(homeName)
  homeData = ppw[[homeName]]
  print(dim(homeData))
  homeName = trim(homeName)
  homeData$TemperatureD = c(0, diff(homeData$TemperatureF))
    
  datasets = list(mins10 = to10min(homeData),
                  minute = homeData)
  results  = list()
  
  k = 0
  for (d in names(datasets)) {
    data     = datasets[[d]]
    nTrain   = trunc(nrow(data) * train.frac)
    k        = k + 1
   	cur_data = subset(data, select = c('date', 'use'))
   	names(cur_data)[2] = 'obs'
   	cur_data$date = as.character(cur_data$date)
   	cur_covar = subset(data, select = c('date', 'TemperatureF', 'TemperatureD'))
   	cur_covar$date = as.character(cur_covar$date)   
    
    # estimate a "breakpoint" indoors temperature
    fit  = lm('use ~ TemperatureF', data = data)    
    fmla = as.formula(paste('~', 'TemperatureF'))
    fit.seg <- try(segmented(fit, seg.Z = fmla, psi = 75))
    psi[homeName, d]  = fit.seg$psi[1,2]        
    cur_covar$TemperatureD = cur_covar$TemperatureF - psi[homeName, d]
    data$TemperatureD = cur_covar$TemperatureD 
    
    # compute breakpoint model temperature contributions
    temp.diff   = data$TemperatureF[1:nTrain] - fit.seg$psi[1,2] 
    temp.diff.p = sapply(temp.diff, function(x) max(c(x,0)))
    temp.diff.n = sapply(temp.diff, function(x) min(c(x,0)))
    bpm.hvac    = fit.seg$coefficients['TemperatureF'] * temp.diff.n + 
      (fit.seg$coefficients['TemperatureF'] + fit.seg$coefficients['U1.TemperatureF']) * temp.diff.p
    bpm.totl = predict(fit.seg)[1:nTrain]
    
    # define model learning controls
    controls = list(
      Kmin = 4, Kmax = 4, 
      maxit = 50, nRestarts = 5, tol = 1e-6,
      thresh.R2 = 0.85, thresh.MAPE = 0.10,
      test.periods = 12,
      vis.interval = 10^(k-1) * 20 * 6)
  
    dir.create(file.path(PLOTS_PATH, paste(homeName, d, sep='/')))
    # learn model
    results[[d]] = stateProcessorWrapper(cur_data, cur_covar, homeName, 
   				  		  	  					         controls = controls,
                                         train.frac = train.frac, 
   					    	  						         verbose = F, plots_path = paste(PLOTS_PATH, '/', homeName, '/', d, '/', sep=''))
    # print(bla)
    
    # format decoded data 
    obs.hvac = usage(data, HVAC)[1:nTrain]
    obs.ac   = usage(data,AC)[1:nTrain]
    obs.hv   = usage(data,HV)[1:nTrain]
    obs.totl = data$use[1:nTrain]
    date.time= data$date[1:nTrain]
    params   = results[[d]]$decoded_data$response
    states   = results[[d]]$decoded_data$states
    st.types = results[[d]]$interp_data$regime.types
    fit.hvac = as.numeric(params$means[2,states]) * data$TemperatureD[1:nTrain]
    fit.totl = colSums(as.matrix(params$means)[,states] * rbind(rep(1,nTrain), data$TemperatureD[1:nTrain]))
    df       = data.frame(date.time, obs.hvac, obs.ac, obs.hv, obs.totl,
                          state = as.factor(states), fit.hvac, fit.totl, 
                          bpm.totl, bpm.hvac)
    
    # compute performance metrics
    df$obs.pr= abs(df$obs.hvac) / df$obs.totl
    df$fit.pr= abs(df$fit.hvac) / df$fit.totl
    df$bpm.pr= abs(df$bpm.hvac) / df$bpm.totl
    df$fit.re= abs(df$obs.hvac - df$fit.hvac) / df$obs.hvac
    df$bpm.re= abs(df$obs.hvac - df$bpm.hvac) / df$obs.hvac
#     is.neutr = sapply(st.types, function(s) length(grep('N', s)>0))
#     df$match.fit  = 1*(df$obs.pr > 0.25 & !(df$state %in% is.neutr))
    df$premise = homeName
    df$resolution = d
    df.tot[[length(df.tot)+1]] = df
  }  
}
}

# plots & aggregation
df.fin   = do.call('rbind', df.tot)
df.fin$hour  = hour(df.fin$date)    
df.hour  = aggregate(data = df.fin, FUN = mean,  
                     cbind(bpm.re, fit.re) ~ hour + premise + resolution)
# df.hour  = aggregate(data = df.fin, FUN = mean,  
#                      cbind(bpm.pr, fit.pr, obs.pr) ~ hour + premise + resolution)
df.hour.m= melt(subset(df.hour, resolution == 'mins10'), id.vars = c('hour', 'premise', 'resolution'))

# plots

plt = ggplot(df.hour.m, aes(hour, value, color = variable, shape = resolution))
plt = plt + geom_point(size = 3) + geom_line(size=1.5)
plt = plt + facet_wrap(~premise, ncol = 2, scales = 'free')
plt

