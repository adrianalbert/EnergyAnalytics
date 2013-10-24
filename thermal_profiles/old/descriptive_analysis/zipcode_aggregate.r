# zipcode_aggregate.r
#
# Analyze aggregation/users at the block group level in the PGE data.
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

plots_path = '~/Dropbox/OccupancyStates/plots/aggregate/'

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
userInfo   = dbGetQuery(db.con, 'SELECT ZIP5, UID, SP_ID, PER_ID, PREM_ID, AREA FROM pge_res_electric')
dbDisconnect(db.con)

# get data for each zipcode
zips  = unique(userInfo$ZIP5)
aggregate_zip = function(zip) {  
  cat(paste('Processing Zip', zip, '\n'))

  # connect to database
  db.con <- dbConnect(dbDriver("MySQL"), 
                      host = "sssl-cyclops.stanford.edu",
                      user = 'adrian', password = 'xmenxmen', 
                      dbname = 'pge_res')
  # pull data
  spids = subset(userInfo, ZIP5 == zip, select = 'SP_ID')
  spids = unique(spids$SP_ID)
  query = paste('SELECT * FROM pge_res_1hour_aligned WHERE sp_id IN (', paste(spids, collapse = ','), ')')
  data  = dbGetQuery(db.con, query)
  dbDisconnect(db.con)
  
  # aggregate
  data$Count = 1
  aggr  = aggregate(data = data, 
                    cbind(hkw1, hkw2, hkw3, hkw4, hkw5, hkw6, hkw7, hkw8, hkw9, hkw10, hkw11, hkw12, 
                          hkw13, hkw14, hkw15, hkw16, hkw17, hkw18, hkw19, hkw20, hkw21, hkw22, hkw23, hkw24, Count) ~ date, FUN = sum)
  aggr[-c(1,26)] = aggr[-c(1,26)] / aggr$Count
  return(aggr)
}

# compute in parallel
aggrZip = mclapply(zips, FUN = aggregate_zip, 
                   mc.cores = 8, mc.silent = F, mc.preschedule = FALSE)
names(aggrZip) = as.character(zips)

# remove zipcodes with no data in the processed table
idx.rm = which(sapply(aggrZip, class) == 'try-error')
aggrZip = aggrZip[-idx.rm]

# save data to disc
save(file = '~/Dropbox/OccupancyStates/data/zip_aggregate.RData', list = c('aggrZip'))

# ------------------------------------------
# Pull & process corresponding weather data
# ------------------------------------------

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
save(file = '~/Dropbox/OccupancyStates/data/zip_aggregate_120k_.RData', list = c('aggrZip', 'weather'))

# ------------------------------------------
# Descriptive plots
# ------------------------------------------

noUsers = sapply(aggrZip, function(l) unique(l$Count))
df = data.frame(noUsers = noUsers)

# plot distribution of #users / zipcode
plt = ggplot(df, aes(x = noUsers, y = ecdf(noUsers)(noUsers)))
plt = plt + geom_step(color = 'black', size = 1.5) + geom_point(size = 4, color = 'blue')
plt = plt + 
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.text.x     = element_text(size=18),
        axis.text.y      = element_text(size=15),
        axis.text.x      = element_text(size=15),
        axis.title.y     = element_text(size=18),
        axis.title.x     = element_text(size=18),
        plot.title       = element_text(size=20),
        legend.text      = element_text(size=15),
        legend.title     = element_text(size=15),
        axis.ticks       = element_blank() ) + 
  ylab('CDF') + xlab('Users per Zip Code') + 
  theme(plot.title=element_text(family="Times", face="bold", size=20)) + 
  ggtitle( 'Distribution of Users per Zip Code' )

png('~/Dropbox/OccupancyStates/plots/aggregate/spid_distribution_aggregate.png', width = 1000, height = 800, res = 200)
print(plt)
dev.off()

