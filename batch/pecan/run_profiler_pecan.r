# run_profiler_pecan.r
#
# Applies HMM decoding on Pecan Street data. 
#
# Adrian Albert
# Last modified: May 2014.
# ---------------------------------------------------------

rm(list = ls())
options(error = recover)
library('segmented')

# __________________________________________________
# Initializations...

setwd('~/EnergyAnalytics/thermal_profiles/profiler/')
source('stateProcessorWrapper.r')
source('stateVisualizerWrapper.r')

DATA_PATH  = '~/energy-data/pecan_street/usage-processed/2012/'
PLOTS_PATH = '~/Dropbox/OccupancyStates/plots/pecan-street'
dir.create(PLOTS_PATH)

# load weather data
weather.hourly = read.csv('~/energy-data/pecan_street/weather/weather_hourly.csv')
weather.15mins = read.csv('~/energy-data/pecan_street/weather/weather_15mins.csv')
weather.hourly$date = as.POSIXct(as.character(weather.hourly$date))
weather.15mins$date = as.POSIXct(as.character(weather.15mins$date))

# load baby names
baby_names  = read.csv('~/Dropbox/OccupancyStates/data/baby-names.csv')

# __________________________________________________
# Select appliances of interest

HVAC   = c('AIR1','AIR2','FURNACE1','FURNACE2')
HV     = c('FURNACE1','FURNACE2')
AC     = c('AIR1','AIR2')
lights = c('LIGHTS_PLUGS1','LIGHTS_PLUGS2','LIGHTS_PLUGS3','LIGHTS_PLUGS4')
total   = c('USE')

# __________________________________________________
# Load up user data

# list all data files for 2013
files    = list.files(path=DATA_PATH, full.names = T)
files_01 = files[grep('minute',files)]
files_15 = files[grep('15mins',files)]
files_60 = files[grep('hourly', files)]

# __________________________________________________
# Plot usage by minute & by hour for each user

# aggregate data from columns defined by labels
usage = function(data,labels) {
  sub = match(labels,names(data))
  sub = sub[!is.na(sub)]
  return(apply(as.matrix(data[,sub]),1,sum))
}

# plot ground truth components
plot_user = function(homeData, main = 'minute') {
  
  names(homeData) = toupper(names(homeData))
  
  # aggregate components
  AC_kwh        = usage(homeData,AC)
  HV_kwh        = usage(homeData,HV)
  total_kwh     = usage(homeData,total)
  occupancy_kwh = total_kwh - AC_kwh - HV_kwh
  
  maxUsage = max(total_kwh)
  plot(homeData$DATE,total_kwh,type='l',col=cols[1],main=main,ylab='kW',ylim=c(0,1.1*maxUsage),xaxt='n')
  points(homeData$DATE,AC_kwh,type='l',col=cols[3])
  points(homeData$DATE,HV_kwh,type='l',col=cols[2])
  points(homeData$DATE,occupancy_kwh,type='l',col=cols[4])
  if ('TEMPERATUREF' %in% names(homeData))
    points(homeData$DATE,homeData$TEMPERATUREF/100*maxUsage,type='l',col=cols[5])
  d = homeData$DATE
  dts = seq(d[1],d[length(d)],by=3600*24) # one per day from one per minute
  axis(1, dts, format(dts, "%a, %m/%d"), cex.axis=1)
  grid(nx=NA,ny=NULL)
  abline(v=dts,col="black",lty=3)
  
}

all_data = list()
for(i in 1:length(files_01)){

  # some gymnastics to get user ID
  userID = strsplit(rev(strsplit(files_01[i], '/')[[1]])[1], '_')[[1]][1]
  
  # read in data at different resolutions...
  cat(paste(files_01[i], files_15[i], files_60[i], '...\n'))
  data_01 = read.csv(files_01[i]); data_01$localminute = as.POSIXct(data_01$localminute); names(data_01)[1] = 'date';
  data_15 = read.csv(files_15[i]); data_15$date = as.POSIXct(data_15$date);
  data_60 = read.csv(files_60[i]); data_60$date = as.POSIXct(data_60$date);
  
  # add in weather (temperature) data
  data_15 = merge(data_15, subset(weather.15mins, select = c('date', 'TemperatureF')), by = 'date')
  data_60 = merge(data_60, subset(weather.hourly, select = c('date', 'TemperatureF')), by = 'date')  
  
  # define plot parameters
  op <- par(no.readonly = TRUE)
  m <- matrix(c(1,2,3,4),nrow=4,ncol=1,byrow=T)
  layout(mat = m,heights = c(0.3,0.3,0.3,0.1))
  par(oma=c(2,2,2,0),mar=c(2,4,2,1)) # Room for the title
  
  cols = c('black','#FE2E2E','#0040FF', '#088A08', '#424242')
  # print minute-by-minute data
  print(plot_user(subset(data_01, date > as.POSIXct('2012-08-01') & date < as.POSIXct('2012-08-14')), main = 'minute'))
  # print 15 min data
  print(plot_user(subset(data_15, date > as.POSIXct('2012-08-01') & date < as.POSIXct('2012-08-14')), main = '15 minute'))
  # print hourly data
  print(plot_user(subset(data_60, date > as.POSIXct('2012-08-01') & date < as.POSIXct('2012-08-14')), main = 'hourly'))
  
  mtext(paste('Pecan Street Experiment User', userID), line=0, font=2, cex=1.2, outer=TRUE)
  par(mar=c(1,4,1,1))
  plot.new()
  legend("center", lty=1,cex=1,lwd=2,
         legend=c('total','HV', 'AC', 'occupant', 'temperature (F)'), 
         col=cols,horiz=T)

  # save plot to file
  dir.create(file.path(PLOTS_PATH, userID))
  dev.copy2pdf(file=paste(paste(PLOTS_PATH, userID, sep = '/'), paste(userID,'.pdf',sep=''), sep='/'),width=10,height=6)
  
  # store data for later processsing
  all_data[[userID]] = list(min_15 = data_15, min_60 = data_60)
}

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

