# ---------------------------
# Format REDD data (CA/MA)
# ---------------------------

rm(list=ls())

# CA data
data_path = '/home/adrian/REDD/CA_Data/redd_pre'
dump_path = '/home/adrian/REDD/CA_Data/redd_pre_processed'

# MA data
# data_path = '/home/adrian/REDD/MA_Data/low_freq/'
# dump_path = '/home/adrian/REDD/MA_Data/redd_pre_processed'

if (!file.exists(dump_path)) dir.create(dump_path, recursive = TRUE)

library('lubridate')

# read in directory names
dirs.list = list.files(path = data_path, full.names = T, include.dirs = TRUE)

# choose some devices of interest
sel.devices = c('ac_plugs', 'furnace', 'space_heater', 'floor_heater', 'electric_heat', 'furance', 
                'air_conditioning' )

for (cur_dir in dirs.list) {
  # read in device labels
  labels  = read.csv(paste(cur_dir, '/labels.dat', sep=''), header = F, sep = ' ')

  # does current house have interesting data?
  ok.vars = subset(labels, V2 %in% sel.devices)
  if (nrow(ok.vars) == 0) {
    cat(paste(cur_dir, 'has no intersting data!\n'))
    next
  }
  ok.vars = subset(labels, V2 %in% c('mains',sel.devices))
  
  # read in data for current house
  files.ok = paste('channel_', ok.vars$V1, '.dat', sep='')
  files.ok = paste(cur_dir, files.ok, sep='/')  
  current_house = list()
  i = 0
  for (f in files.ok) {
    i = i + 1
    cat(paste('Processing device', f, '\n'))
    
    # read in data for current device
    cur_data           = try(read.csv(f, header = F, sep = ' '))
    if (class(cur_data) == 'try-error') {
      cat(paste('no data for', f, '\n'))
      next
    }
    processed          = cur_data
    names(processed)   = c('UTC', 'Watt')
    processed$timestamp= as.POSIXct(processed$UTC, origin='1970-01-01')
    processed$DateTime = paste(paste( year(processed$timestamp), sprintf('%02d',month(processed$timestamp)), 
                                      sprintf('%02d',mday(processed$timestamp)), sep = '-'), 
                               paste( sprintf('%02d',hour(processed$timestamp)), ':00:00', sep=''), sep=" ")
      
    # aggregate watts to kWh
    processed$kWh      = processed$Watt / 3600 / 1000 # * c(1,diff(processed$UTC))
    noObs              = table(processed$DateTime)
    final              = aggregate(data = processed, kWh ~ DateTime, FUN = mean)
    final$kWh          = final$kWh * as.numeric(noObs)
    
    # save data to list
    cur_name           = paste(ok.vars$V2[i], ok.vars$V1[i], sep='.')
    current_house[[cur_name]] = final
  }

  # align all devices to same time vector
  cur_names     = names(current_house)
  timestamps    = lapply(current_house, function(d) d$DateTime)
  timestamps.ok = timestamps[[1]]
  for (i in 2:length(timestamps))
    timestamps.ok = intersect(timestamps[[i]], timestamps.ok)
  current_house = lapply(current_house, function(h) subset(h, DateTime %in% timestamps.ok, select = 'kWh'))
  current_house = do.call('cbind', current_house)   
  names(current_house) = cur_names
  current_house$DateTime = timestamps.ok  
  
  # save data to file
  if (nrow(current_house) == 0) {
    cat(paste(cur_dir, 'has no intersting data!\n'))
    next    
  }
  house_id = strsplit(cur_dir, '/')[[1]]
  house_id = house_id[[length(house_id)]]
  write.csv(current_house, file = paste(dump_path, '/', house_id, '.csv', sep = ''), row.names = F)
}

# ______________________________________
# Format weather data for REDD zipcodes

wthr_path = '/home/adrian/REDD/CA_Data/weather/'
wdmp_path = '/home/adrian/REDD/CA_Data/weather_processed'
if (!file.exists(wdmp_path)) dir.create(wdmp_path, recursive = TRUE)
wthr.list = list.files(path = wthr_path, full.names = T)

for (f in wthr.list) {
  
  cat(paste('Processing:', f, '\n'))
  
  cur_data = read.csv(f)  
  cur_data = subset(cur_data, select = c('Time', 'TemperatureF', 'PressureIn'))
  cur_data$Time = as.character(cur_data$Time)
  cur_data$date = paste(paste( year(cur_data$Time), 
                               sprintf('%02d',month(cur_data$Time)), 
                               sprintf('%02d',mday(cur_data$Time)), sep = '-'), 
                             paste( sprintf('%02d',hour(cur_data$Time)), ':00:00', sep=''), sep=" ")
  cur_data      = aggregate(data = cur_data, cbind(TemperatureF, PressureIn)~ date, FUN = mean)  
  f_save        = strsplit(f, '\\.')[[1]][1]
  write.csv(cur_data , file = paste(f_save, '_processed.csv', sep=''), row.names = F)
}

