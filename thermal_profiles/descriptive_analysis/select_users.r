# select_users.r
#
# Select consumption and weather data for users in the PGE residential dataset. 
# 
# Adrian Albert
#
# Last modified: June 2013.

# -------------------------------
# Initializations ... 
# -------------------------------

rm(list=ls())

options(error = recover)

library('lubridate')
library('zoo')
library('ggplot2')
library('reshape')
library('dummies')
library('RMySQL')
library('imputation')
library('parallel')

# ------------------------------------------
# Connect to databases, get data
# ------------------------------------------

# disconnect all MySQL connections
library('RMySQL')
all_cons <- dbListConnections(MySQL())
for(con in all_cons) dbDisconnect(con)

# get user info
db.con <- dbConnect(dbDriver("MySQL"), 
                    host = "sssl-cyclops.stanford.edu",
                    user = 'adrian', password = 'xmenxmen', 
                    dbname = 'pge_res')
userInfo   = dbGetQuery(db.con, 'SELECT ZIP5, UID, SP_ID, PER_ID, PREM_ID, AREA, CECCLMZN FROM pge_res_electric')
dbDisconnect(db.con)

# only select users in zones 3,4,12,13 
userInfo = subset(userInfo, CECCLMZN %in% paste('Z', c('03', '04', '12', '13'), sep = ''))

# only consider zipcodes with at least 50 users 
zip_tab = table(userInfo$ZIP5)
zip_tab = names(zip_tab[which(zip_tab > 100)])
userInfo= subset(userInfo, ZIP5 %in% zip_tab)

# select data for a subset of users
db.con <- dbConnect(dbDriver("MySQL"), 
                    host = "sssl-cyclops.stanford.edu",
                    user = 'adrian', password = 'xmenxmen', 
                    dbname = 'pge_res')
qspid  = paste('SELECT DISTINCT sp_id FROM pge_res_1hour_aligned')
spids  = dbGetQuery(db.con, qspid)

# sample from valid sp_ids 
userInfo = subset(userInfo, SP_ID %in% spids$sp_id)
sel.usr  = sample(unique(userInfo$SP_ID), 10000)
userInfo = subset(userInfo, SP_ID %in% sel.usr)
query = paste('SELECT * FROM pge_res_1hour_aligned WHERE sp_id IN (', paste(sel.usr, collapse = ','), ')')
consumption  = dbGetQuery(db.con, query)

dbDisconnect(db.con)

spid_zip         = userInfo$ZIP5
names(spid_zip)  = userInfo$SP_ID
consumption      = cbind(ZIP5 = spid_zip[as.character(consumption$sp_id)], consumption)

  # ------------------------------------------
# Pull & process corresponding weather data
# ------------------------------------------

zips = unique(userInfo$ZIP5)

db.con <- dbConnect(dbDriver("MySQL"), 
                    host = "sssl-cyclops.stanford.edu",
                    user = 'adrian', password = 'xmenxmen', 
                    dbname = 'PGE_WEATHER')
query = paste("SELECT zip5, date, TemperatureF, Humidity from weather_60 where zip5 IN (", paste(zips, collapse=','),
              ") and date between '2010-08-01 00:00:00' and '2011-07-31 23:00:00' order by zip5, date;")
weather.data = dbGetQuery(db.con, query)
dbDisconnect(db.con)

## impute NA values in temp.info by linear regression
process_weather = function(zip) {  
  cat(paste('ZIP=', zip, '\n'))  
  temp.info = subset(weather.data, zip5 == zip)  
  col.na    = colSums(is.na(temp.info))
  if (length(which(col.na>0))>0) {
    temp.info = temp.info[,colSums(is.na(temp.info))<0.8 * nrow(temp.info)]
    varnames  = names(temp.info)
    if (ncol(temp.info) == 3) vars = c(3,3) else vars = -c(1,2)
    res       = gbmImpute(temp.info[,vars], verbose = F)
    if (ncol(temp.info) == 3) imp = res$x[,1] else imp = res$x
    temp.info = cbind(temp.info[,c(1,2)], imp)
    names(temp.info) = varnames
  }
  return(temp.info[,c(1,2,3)])
}
weather = mclapply(zips, process_weather, 
                   mc.cores = 8, mc.silent = F, mc.preschedule = FALSE)
names(weather) = zips
idx.rm = which(sapply(weather, class) == 'try-error')

# concatenate weather data
weather = do.call('rbind', weather)

# save data to disc
save(file = '~/Dropbox/OccupancyStates/data/selection_10k_1yr.RData', list = c('consumption', 'weather'))

