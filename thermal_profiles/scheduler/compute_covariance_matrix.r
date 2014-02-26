# compute_covariance_matrix.r
#
# Estimate empirical covariance matrix.
# 
# Adrian Albert
# Last modified: February 2014.
# -----------------------------------------------------------------------

# ____________________________________________________
# Initializations....

rm(list = ls())
# options(warn=2)
options(error = recover)
setwd('~/EnergyAnalytics/thermal_profiles/scheduler/')

library('lubridate')
library('parallel')

PLOTS_PATH = '~/Dropbox/OccupancyStates/plots/scheduling/'
DATA_PATH  = '~/Dropbox/OccupancyStates/fits/bakersfield/'

D   = 20 # use at least D days
rho = 5  # total temperature discrepancy has to be smaller than rho for two profiles to be "similar"

# ____________________________________________________
# Load response data

# data files for users
dec_files  = list.files(path = DATA_PATH, pattern = '*decoded*', full.names = T, recursive = T)

# load & format data
data = mclapply(1:length(dec_files), mc.cores = 5, 
                function(i) {                        
                  if (i %% 10 == 0) cat(paste(i, '/', length(dec_files), '...'))      
                  if (i %% 100 == 0) cat('\n')                  
                  load(dec_files[i])
                  a_vec = as.numeric(data$response$means[2,data$states])
                  no.days = length(a_vec) %/% 24
                  a_vec = as.data.frame(matrix(a_vec[1:(no.days * 24)], ncol = 24, byrow=T))
                  names(a_vec) = 1:24
                  UID = strsplit(strsplit(dec_files[i], '_')[[1]][1], '//')[[1]][2]
                  a_vec$UID = UID
                  return(a_vec)
                })
uids = unlist(sapply(data, function(d) unique(d$UID)))
names(data) = uids

# ____________________________________________________
# Load weather data

# load historic weather data
load('~/Dropbox/OccupancyStates/data/selection_consumption.RData')

# load in temperature forecast
weather  = read.csv('~/Dropbox/OccupancyStates/data/weather_stanford_10_29_2013.csv')
weather  = subset(weather, select = c('Time', 'TemperatureF'))
weather$Hour = hour(as.POSIXct(weather$Time))
weather  = aggregate(TemperatureF ~ Hour, data = weather, FUN = mean)
weather$Hour = NULL
weather$TemperatureF = weather$TemperatureF + 25
weather = rbind(weather, weather$TemperatureF[nrow(weather)]*0.95)

# map users to zipcodes
raw_data  = subset(consumption.ok, UID %in% uids)
uid.zips  = subset(raw_data, UID %in% uids)
uid.zips  = paste(uid.zips$ZIP5, uid.zips$UID)
uid.zips  = unique(uid.zips)
uid.zips  = sapply(uid.zips, function(s) t(strsplit(s, ' ')[[1]]))
uid2zip   = t(uid.zips)[,1]
names(uid2zip) = t(uid.zips)[,2]

# for each user, find D days with temperature profiles "similar" to the input profile
days.resp = mclapply(uids, mc.cores = 5, 
             function(u) {
               
               # get weather data for current user
               cur_zip  = uid2zip[u]
               no.days  = nrow(data[[u]])
               cur_wthr = weather.ok[[cur_zip]]$TemperatureF[1:(no.days * 24)]
               cur_wthr = as.data.frame(matrix(cur_wthr, ncol = 24, byrow=T))
               
               # find "similar weather" days 
               DTemp = cur_wthr - matrix(rep(weather$TemperatureF, nrow(cur_wthr)), nrow = nrow(cur_wthr), byrow=T)               
               DTemp = sqrt(rowSums(DTemp^2))
               DTemp = order(DTemp)
               
               ret = t(as.matrix(DTemp[1:D]))
               rownames(ret) = u
               return(ret)
             })
days.resp = do.call('rbind', days.resp)

# ____________________________________________________
# Compute covariance matrix

# compute E[aiaj]
idx = expand.grid(i = uids, j = uids)
covmat = mclapply(1:nrow(idx), 
                  mc.cores = 5,
                  function(k) {
                    
                    ui = as.character(idx[k,'i']); uj = as.character(idx[k,'j']); #ui = uids[i]; uj = uids[j]
                    
                    if (k %% length(uids) == 0) cat(paste(k/length(uids), '/', length(uids), '\n'))      

                    a_i = data[[ui]][days.resp[ui,],-25]
                    a_j = data[[uj]][days.resp[uj,],-25]
                    wij = lapply(1:D, function(d) t(as.matrix(a_i[d,])) %*% as.matrix(a_j[d,]))
                    mat = wij[[1]]; for(p in 2:D) mat = mat + wij[[p]]
                    mat = (1 / D) * mat
                    dm  = t(as.matrix(diag(mat)))
                    rownames(dm) = k
                    
                    return(dm)
                  })

save.image(file = '~/Dropbox/OccupancyStates/data/bakersfield_covmat.RData')

