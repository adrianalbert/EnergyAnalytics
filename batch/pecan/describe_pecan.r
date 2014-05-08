# describe_pecan.r
#
# Performs a basic statistical description of the Pecan St dataset.
#
# Adrian Albert
# Last modified: May 2014.
# ---------------------------------------------------------

rm(list = ls())
options(error = recover)

# __________________________________________________
# Initializations...

setwd('~/EnergyAnalytics/batch/pecan/')
source('../../utils/aggregate_data.r')
source('../../batch/pecan/define_categories_pecan.r')

# use data from 2013
DATA_PATH_SEL = '~/energy-data/pecan_street/usage-select/'
METADATA_PATH = '~/energy-data/pecan_street/metadata/'
PLOTS_PATH    = '~/Dropbox/OccupancyStates/plots/pecan-street'
dir.create(PLOTS_PATH)

# __________________________________________________
# Access usage data

# list all data files for the selected subset of data
files_prc    = list.files(path=DATA_PATH_PRC, full.names = T, recursive = T, pattern = '*hourly*')
files_prc_60 = list.files(path=DATA_PATH_PRC, full.names = T, pattern = '*hourly*')
files_15 = files[grep('15mins',files)]
files_60 = files[grep('hourly', files)]

# extract user IDs 
users_01 = sapply(files_01, function(s) strsplit(tail(strsplit(s, '/')[[1]], 1), '\\.')[[1]][1])
users_15 = sapply(files_15, function(s) strsplit(tail(strsplit(s, '/')[[1]], 1), '_')[[1]][1])
users_sel= intersect(users_01, users_15)


# load IDs and names
usr_name = read.csv(paste(METADATA_PATH, 'user_names_ids.csv', sep = '/'))


# __________________________________________________
# Plot usage by minute & by hour for each user

# plot ground truth components
plot_user = function(homeData, main = 'minute') {
  
  names(homeData) = tolower(names(homeData))
  
  # aggregate components
  AC_kwh        = add.columns(homeData,AC)
  HV_kwh        = add.columns(homeData,HV)
  total_kwh     = add.columns(homeData,total)
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
for(i in 1:nrow(user_names)){

  # current user info
  user   = user_names[i,]
  
  # some gymnastics to get appropriate data files for this user  
  files_01_i = files_01[grep(paste('/', user$ID, '.csv', sep=''), files_01)]
  files_15_i = files_15[grep(paste('/', user$ID, '_15mins.csv', sep=''), files_15)]
  files_60_i = files_60[grep(paste('/', user$ID, '_hourly.csv', sep=''), files_60)]

  # read in data at different resolutions...
  cat(paste(user$ID, '...\n'))
  data_01 = read.csv(files_01_i); data_01$localminute = as.POSIXct(data_01$localminute); names(data_01)[1] = 'date';
  data_15 = read.csv(files_15_i); data_15$date = as.POSIXct(data_15$date);
  data_60 = read.csv(files_60_i); data_60$date = as.POSIXct(data_60$date);
  
  dir.create(file.path(PLOTS_PATH, userID))
  pdf(file=paste(paste(PLOTS_PATH, userID, sep = '/'), paste(userID,'.pdf',sep=''), sep='/'),width=10,height=6)
  
  # define plot parameters
  op <- par(no.readonly = TRUE)
  m <- matrix(c(1,2,3,4),nrow=4,ncol=1,byrow=T)
  layout(mat = m,heights = c(0.3,0.3,0.3,0.1))
  par(oma=c(2,2,2,0),mar=c(2,4,2,1)) # Room for the title
  
  cols = c('black','#FE2E2E','#0040FF', '#088A08', '#424242')
  # print minute-by-minute data
  print(plot_user(data_01[1:(7*24*60),], main = 'minute'))
  # print 15 min data
  print(plot_user(data_15[1:(7*24*4),], main = '15 minute'))
  # print hourly data
  print(plot_user(data_60[1:(7*24),], main = 'hourly'))
  
  mtext(paste('Pecan Street Experiment User', userID), line=0, font=2, cex=1.2, outer=TRUE)
  par(mar=c(1,4,1,1))
  plot.new()
  legend("center", lty=1,cex=1,lwd=2,
         legend=c('total','HV', 'AC', 'occupant', 'temperature (F)'), 
         col=cols,horiz=T)

  # save plot to file
  dev.off()
  
  # store data for later processsing
  all_data[[userID]] = list(min_15 = data_15, min_60 = data_60)
}

# save data to RData file
save(list = c('all_data'), file = '~/energy-data/pecan_street/')