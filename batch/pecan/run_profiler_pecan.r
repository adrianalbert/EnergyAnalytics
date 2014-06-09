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
source('../../batch/pecan/define_categories_pecan.r')

DATA_PATH = '~/energy-data/pecan_street/usage-select/'
DUMP_PATH = '~/energy-data/pecan_street/models/'
PLOT_PATH = '~/Dropbox/OccupancyStates/plots/pecan-street'

# load user names
user_names = read.csv('~/energy-data/pecan_street/metadata/user_names_ids.csv')

# __________________________________________________
# Load up user data

# list all data files by year/uid
files.input = list.files(path=DUMP_PATH, pattern = '*_decoded*', full.names = T, recursive = T)
already_done  = lapply(files.input, function(x) {
  tmp = strsplit(x, '/')[[1]]
  ret = data.frame(res = as.character(tmp[length(tmp)-1]),                   
                   yr  = as.character(tmp[length(tmp)-3]),
                   uid = as.character(tmp[length(tmp)-2]))
  rownames(ret) = NULL
  return(ret)
})
already_done = do.call('rbind', already_done)
already_done$yr = as.character(already_done$yr)
already_done$uid= as.character(already_done$uid)
already_done$res= as.character(already_done$res)

# list all data files
files    = list.files(path=DATA_PATH, full.names = T, recursive = T)
files_01 = files[grep('01min',files)]
files_15 = files[grep('15min',files)]
files_60 = files[grep('60min', files)]

# extract ID and year
usersVec = data.frame(UID = as.character(sapply(files_60, function(s) strsplit(tail(strsplit(s, '/')[[1]], 1), '\\.')[[1]][1])),
                      year= as.character(sapply(files_60, function(s) tail(strsplit(s, '/')[[1]], 2)[1])))
rownames(usersVec) = NULL

# __________________________________________________
# Apply Thermal States model to Pecan data

# format data in the way it's expected by the HMM package
format_data = function(homeData) {
  
  # temperature above reference
  homeData$TemperatureD = homeData$TemperatureF - 65
  
  # format data as expected by the HMM package
  cur_data = subset(homeData, select = c('date', 'use'))
  names(cur_data)[2] = 'obs'
  cur_data$date = as.character(cur_data$date)
  cur_covar = subset(homeData, select = c('date', 'TemperatureF', 'TemperatureD'))
  cur_covar$date = as.character(cur_covar$date)
  
  return(list(cur_data, cur_covar))
}

apply_thermal_model = function(cur_data, cur_covar, userName, 
                               dump_path = NULL, 
                               plot_path = NULL,
                               train.frac = 0.9) {
  
  nTrain   = trunc(nrow(cur_data) * train.frac)
  
  # define model learning controls
  controls = list(
    Kmin = 4, Kmax = 4, 
    maxit = 50, nRestarts = 5, tol = 1e-6,
    thresh.R2 = 0.85, thresh.MAPE = 0.10,
    test.periods = 12,
    vis.interval = 3 * 24
  )
  
  # generate visualization interval; make sure there's data in there
  # TODO: there was an error generated here (indices for subsetting were messed up)
  ok = FALSE
  no.secs = controls$vis.interval * 3600
  while (!ok) {
    start_date = sample(cur_data$date[-((nrow(cur_data)-controls$vis.interval):nrow(cur_data))], 1)
    stop_date  = as.character(as.POSIXct(start_date) + no.secs)
    dat        = subset(cur_data, date >= start_date & date < stop_date)
    if (nrow(na.omit(dat)) > 0) 
      ok = TRUE
  }      
  
  
  # learn model
  res = try(stateProcessorWrapper(cur_data, cur_covar, userName, 
                              controls = controls,
                              train.frac = train.frac, 
                              verbose = F, 
                              dump_path = dump_path))
  if (class(res) == 'try-error') {
    cat('Error in learning model for current user!\n')
  }
  # produce visualizations
  res = try(stateVisualizerWrapper(res$decoder, 
                               res$interpreter, 
                               plots_path = plot_path, 
                               interval = c(start_date, stop_date)))
  if (class(res) == 'try-error') {
    cat('Error in visualizing current user!\n')
  }
  
  return(NULL)
}

res = mclapply(1:nrow(usersVec), 
               mc.cores = 5,
               function(i) {
  # load data             
  user     = usersVec[i,] 
  userName = as.character(user_names[which(user_names$ID == user$UID),'name'])
  cat(paste('Processing user', user$UID, '/', user$year, ':', i, '/', nrow(usersVec)))  
  
  idx = which(user$UID == already_done$uid & as.character(user$year) == already_done$yr)
  if (length(idx)>0) {
    cat('Already processed!\n')
    return(NULL)
  }
  
  homeData15 = read.csv(files_15[i])     
  homeData60 = read.csv(files_60[i])     
  
  # create directory to store models
  dump_path_15 = file.path(DUMP_PATH, paste(user$year, user$UID, '15min/', sep='/')); 
  dir.create(dump_path_15, recursive = T)
  dump_path_60 = file.path(DUMP_PATH, paste(user$year, user$UID, '60min/', sep='/')); 
  dir.create(dump_path_60, recursive = T)
  
  # create directory to store plots
  plot_path_15 = file.path(PLOT_PATH, paste(user$year, user$UID, '15min/', sep='/')); 
  dir.create(plot_path_15, recursive = T)
  plot_path_60 = file.path(PLOT_PATH, paste(user$year, user$UID, '60min/', sep='/')); 
  dir.create(plot_path_60, recursive = T)
  
  # format datasets
  res = format_data(homeData15); cur_data15 = res[[1]]; cur_covar15 = res[[2]];
  res = format_data(homeData60); cur_data60 = res[[1]]; cur_covar60 = res[[2]];

  # apply model to data
  res = apply_thermal_model(cur_data15, cur_covar15, userName, 
                            dump_path = dump_path_15, 
                            plot_path = plot_path_15)
  res = apply_thermal_model(cur_data60, cur_covar60, userName, 
                            dump_path = dump_path_60, 
                            plot_path = plot_path_60)  
  return(NULL)
})

